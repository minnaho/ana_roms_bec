!   $Log: wrt_phys_flux_avg.F,v $
!   Revision 1.2  2002/04/19 00:53:19  frenzel
!   Added code for the analysis of mixing and nudging fluxes.
!
!   Revision 1.1  2002/01/16 19:58:53  frenzel
!   New routines related to flux analysis.
!
!
#include "cppdefs.h"
#if defined PHYS_FLUX_ANALYSIS && defined AVERAGES
! Write time-averaged physical flux fields into averages netCDF file
      subroutine wrt_phys_flux_avg

#ifdef PFA_BY_TRACER
      use phys_flux
#endif
      implicit none
#include "netcdf.inc"
      integer ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), ibuff(4), nf_fwrite, iflux, iflux_out            

# if defined MPI && !defined PARALLEL_FILES
#  ifdef HOFFMAN2
      include "mpif.h"
#  else
#  include "mpif.h"
#  endif
      integer status(MPI_STATUS_SIZE), blank
# endif
# include "param.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "ncvars.h"
# include "physflux.h"
# include "ncroms_physflux.h"

#if defined MPI && !defined PARALLEL_FILES
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
#endif
!
! Create/open averages file; write grid arrays, if so needed.
!
      call def_phys_flux_avg(ncid_phys_flux_avg, nrec_phys_flux_avg, 
     &     ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(phys_flux_avg_name)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index
      nrec_phys_flux_avg=max(nrec_phys_flux_avg,1)!!! nrec_phys_flux_avg 
      if (nrpf_phys_flux_avg.eq.0) then      !!! is already advanced by main.
        record=nrec_phys_flux_avg
      else
        record=1+mod(nrec_phys_flux_avg-1, nrpf_phys_flux_avg)
      endif
!
! Write out time-averaged variables:
! ----- --- ------------- ----------
!
! Time step and record indices.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrec_phys_flux_avg

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int(ncid_phys_flux_avg, phys_flux_avgTstep,  
     &     start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,1) 'time_step', record,ierr MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged time
!
      ierr=nf_put_var1_double (ncid_phys_flux_avg, phys_flux_avgTime, 
     &     record, time_phys_flux_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
! Barotropic mode variable: free-surface
! Always include this in the output so that depths of sigma levels
! can be computed
      ierr=nf_fwrite (zeta_phys_flux_avg(START_2D_ARRAY), 
     &     ncid_phys_flux_avg, phys_flux_avgZ,record, r2dvar)
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxZ))
         write(*,1) vname(1,indxZ)(1:lvar), record, ierr
     &        MYID
         goto 99                !--> ERROR
      endif
! horizontal advective fluxes (xi-direction)
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(HorXAdvFlux_avg(START_2D_ARRAY,1,iflux_out), 
     &        ncid_phys_flux_avg, avgHorXAdvFlux(iflux_out), record, 
     &        r3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_phys(1,indxHorXAdvFlux+iflux-1))
            write(*,1) 
     &           vname_phys(1,indxHorXAdvFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
! horizontal advective fluxes (eta-direction)
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(HorYAdvFlux_avg(START_2D_ARRAY,1,iflux_out), 
     &        ncid_phys_flux_avg, avgHorYAdvFlux(iflux_out), record, 
     &        r3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_phys(1,indxHorYAdvFlux+iflux-1))
            write(*,1) 
     &           vname_phys(1,indxHorYAdvFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
! vertical advective fluxes
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(VertAdvFlux_avg(START_2D_ARRAY,
     &        0,iflux_out), 
     &        ncid_phys_flux_avg, avgVertAdvFlux(iflux_out), record, 
     &        w3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_phys(1,indxVertAdvFlux+iflux-1))
            write(*,1) 
     &           vname_phys(1,indxVertAdvFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
#ifdef HOR_DIFF_ANALYSIS
! horizontal mixing fluxes (xi-direction)
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(HorXMixFlux_avg(START_2D_ARRAY,1,iflux_out), 
     &        ncid_phys_flux_avg, avgHorXMixFlux(iflux_out), record, 
     &        r3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxHorXMixFlux+iflux-1))
            write(*,1) vname(1,indxHorXMixFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
! horizontal mixing fluxes (eta-direction)
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(HorYMixFlux_avg(START_2D_ARRAY,1,iflux_out), 
     &        ncid_phys_flux_avg, avgHorYMixFlux(iflux_out), record, 
     &        r3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxHorYMixFlux+iflux-1))
            write(*,1) vname(1,indxHorYMixFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
! vertical mixing fluxes
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(VertMixFlux_avg(START_2D_ARRAY,
     &        0,iflux_out), 
     &        ncid_phys_flux_avg, avgVertMixFlux(iflux_out), record, 
     &        w3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxVertMixFlux+iflux-1))
            write(*,1) vname(1,indxVertMixFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
#endif /* HOR_DIFF_ANALYSIS */
#ifdef VERT_DIFF_ANALYSIS
! vertical diffusive fluxes
      do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
        if (wrt_pfa(iflux)) then
         iflux_out = pfa_out(iflux)
# else 
         iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
         ierr = nf_fwrite(VertDiffFlux_avg(START_2D_ARRAY,0,iflux_out), 
     &        ncid_phys_flux_avg, avgVertDiffFlux(iflux_out), 
     &        record, w3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_phys(1,indxVertDiffFlux+iflux-1))
            write(*,1) 
     &           vname_phys(1,indxVertDiffFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
# ifdef PFA_BY_TRACER
        endif
# endif /* PFA_BY_TRACER */ 
      end do
#endif /* VERT_DIFF_ANALYSIS */
# ifdef FULL_PHYS_FLUX_ANALYSIS
! top fluxes
      do iflux = 1, NT_PFA
         ierr = nf_fwrite(TopFlux_avg(START_2D_ARRAY,iflux), 
     &        ncid_phys_flux_avg, avgTopFlux(iflux), record, r2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTopFlux+iflux-1))
            write(*,1) vname(1,indxTopFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
      end do
# ifdef BOTTOM_FLUX
! bottom fluxes
      do iflux = 1, NT_PFA
         ierr = nf_fwrite(BottomFlux_avg(START_2D_ARRAY,iflux), 
     &        ncid_phys_flux_avg, avgBottomFlux(iflux), record, r2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxBottomFlux+iflux-1))
            write(*,1) vname(1,indxBottomFlux+iflux-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
      end do
# endif /* BOTTOM_FLUX */
! nudging divergences
      do iflux = 1, NT_PFA
         ierr = nf_fwrite(NudgingFlux_avg(START_2D_ARRAY,
     &        1,iflux), 
     &        ncid_phys_flux_avg, avgNudgingFlux(iflux), record, r3dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxNudgingFlux+iflux-1))
            write(*,1) vname(1,indxNudgingFlux+iflux-1)(1:lvar), 
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
         end if
      end do
# endif /* FULL_PHYS_FLUX_ANALYSIS */
  1   format(/' WRT_PHYS_FLUX_AVG - ERROR while writing variable(',
     &     1x,a,1x,
     &               ')into averages file.',/,11x,'Time record:',
     &                      i6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
  99  may_day_flag=3
 100  continue

!
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!
# if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close(ncid_phys_flux_avg)
      if (nrpf_phys_flux_avg.gt.0 .and. record.ge.nrpf_phys_flux_avg) 
     &     ncid_phys_flux_avg=-1
# else
      if (nrpf_phys_flux_avg.gt.0 .and. 
     &     record.ge.nrpf_phys_flux_avg) then
        ierr=nf_close(ncid_phys_flux_avg)
        ncid_phys_flux_avg=-1
      else
        ierr=nf_sync(ncid_phys_flux_avg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
        mpi_master_only write(*,'(6x,A,2(A,I4,1x),A,I3)') 
     &        'WRT_PHYS_FLUX_AVG -- wrote ',
     &        'averaged fields into time record =', record, '/',
     &        nrec_phys_flux_avg  MYID
      else
        write(*,'(/1x,2A/)') 
     &        'WRT_PHYS_FLUX_AVG ERROR: Cannot ',
     &        'synchronize/close averages netCDF file.'
        may_day_flag=3
      endif

# if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) then
         blank=mynode
         call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_WORLD,  ierr)
      endif
# endif
      return
      end
#else /* PHYS_FLUX_ANALYSIS && AVERAGES */
      subroutine wrt_phys_flux_avg_empty
      end
#endif /* PHYS_FLUX_ANALYSIS && AVERAGES */
