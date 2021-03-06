#include "cppdefs.h"
                               ! Writes requested model
      subroutine wrt_his       ! fields at requested levels
                               ! into history netCDF file.
      implicit none
#include "param.h"
#include "ncvars.h"
      integer ierr,  record, i, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(iaux), nf_fwrite
      character(len=18) tstring
#include "scalars.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "eos_vars.h"
#include "mixing.h"
#include "netcdf.inc"
#ifdef SEDIMENT_BIOLOGY
# include "bgcflux.h"
#endif
#ifdef BIOLOGY_BEC
# include "ecosys_bec.h"
#endif
#ifdef SOLVE3D
      integer tile, itrc, j,k
# include "work.h"
#endif
#if defined MPI && !defined PARALLEL_FILES
# ifdef HOFFMAN2
      include "mpif.h"
# else
# include "mpif.h"
# endif
      integer status(MPI_STATUS_SIZE), blank
      if (mynode.gt.0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif

! Create/open history file; write grid arrays, if so needed.

      call def_his (nrechis, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(hisname)
                                        !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
                                        !!! started, it is assumed
      nrechis=max(nrechis,1)            !!! that the global history
      if (nrpfhis.eq.0) then            !!! record index "nrechis"
        record=nrechis                  !!! is advanced by main.
      else
        record=1+mod(nrechis-1, nrpfhis)
      endif
 
 
c-#define CR
CR      write(*,*) 'wrt_his: Entry  ' MYID

! Write out evolving model variables:
! ----- --- -------- ----- ----------

      ibuff(1)=iic-1               ! Time step and
      ibuff(2)=nrecrst             ! record numbers.
      ibuff(3)=nrechis
      ibuff(4:iaux)=0
#ifdef AVERAGES
      ibuff(4)=nrecavg
#endif
#ifdef STATIONS
      ibuff(5)=nrecstn
#endif
#ifdef FLOATS
      ibuff(6)=nrecflt
#endif
      start(1)=1
      start(2)=record
      count(1)=iaux
      count(2)=1
      ierr=nf_put_vara_int (ncidhis, hisTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_his :: ',
     &       'Cannot write variable ''time_step'' into history file, ',
     &       'rec =', record, nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif

! Time

      ierr=nf_put_var1_FTYPE (ncidhis, hisTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_his :: ',
     &        'Cannot write variable ''',    vname(1,indxTime)(1:lvar),
     &        ''' into history file, rec =', record, nf_strerror(ierr)
     &            MYID
        goto 99                                           !--> ERROR
      endif
CR      write(*,*) 'wrt_his:  time  ' MYID

! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.

      if (wrthis(indxZ)) then
        ierr=nf_fwrite (zeta(START_2D_ARRAY,knew), ncidhis, hisZ,
     &                                              record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          write(*,1) vname(1,indxZ)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtZ  ' MYID
      endif
      if (wrthis(indxUb)) then
        ierr=nf_fwrite (ubar(START_2D_ARRAY,knew), ncidhis, hisUb,
     &                                              record, u2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          write(*,1) vname(1,indxUb)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtUBAR' MYID
      endif
      if (wrthis(indxVb)) then
        ierr=nf_fwrite (vbar(START_2D_ARRAY,knew), ncidhis, hisVb,
     &                                            record, v2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          write(*,1) vname(1,indxVb)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtVBAR' MYID
      endif
 
#ifdef SOLVE3D

! 3D momentum components in XI- and ETA-directions.

      if (wrthis(indxU)) then
        ierr=nf_fwrite (u(START_2D_ARRAY,1,nstp), ncidhis,
     &                             hisU,  record,  u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          write(*,1) vname(1,indxU)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtU  ' MYID
      endif
      if (wrthis(indxV)) then
        ierr=nf_fwrite (v(START_2D_ARRAY,1,nstp), ncidhis,
     &                             hisV,  record,  v3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          write(*,1) vname(1,indxV)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtV  ' MYID
      endif

! Tracer variables.

# if defined OUTPUT_AVG_PHYS_BIO
          do itrc=1,1+ntrc_salt
# else
      do itrc=1,NT
# endif
        if (wrthis(indxT+itrc-1)) then
          ierr=nf_fwrite (t(START_2D_ARRAY,1,nstp,itrc), ncidhis,
     &                                hisT(itrc), record, r3dvar)


c          write(*,*) 'Replacing T,S with rho1,qp1.'
c          if (itrc.eq.1) then
c            ierr=nf_fwrite (rho1, ncidhis, hisT(itrc), record,r3dvar)
c          elseif (itrc.eq.2) then
c            ierr=nf_fwrite  (qp1, ncidhis, hisT(itrc), record,r3dvar)
c          endif

          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            write(*,1) vname(1,indxT+itrc-1)(1:lvar),
     &                 record MYID
            goto 99                                       !--> ERROR
          endif
CR      write(*,*) 'wrt_his:  wrtT  ' MYID
        endif
      enddo
#if defined SEDIMENT_BIOLOGY && !defined OUTPUT_AVG_PHYS_BIO
!
! Sediment concentrations
!
      do itrc=1,NT_sed
        if (wrthis(indxSedFirst+itrc-1)) then
          ierr=nf_fwrite(t_sed(START_2D_ARRAY,itrc), ncidhis,
     &                                hisTsed(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            write(*,1) vname(1,indxSedFirst+itrc-1)(1:lvar),
     &                      record, ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
CR      write(*,*) 'wrt_his:  wrtTsed  ' MYID
        endif
      enddo
#endif /* SEDIMENT_BIOLOGY */
!
! Density anomaly.

      if (wrthis(indxR)) then
# ifdef SPLIT_EOS
cxxxxxxxxxx
cxxxxxxxxxx
cxxxxxxxxxx
        ierr=nf_fwrite (rho1, ncidhis, hisR, record, r3dvar)
# else
        ierr=nf_fwrite (rho, ncidhis, hisR, record, r3dvar)
# endif
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          write(*,1) vname(1,indxR)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtRHO' MYID
      endif

! S-coordinate omega vertical velocity (m/s).

      if (wrthis(indxO)) then
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
            enddo
          enddo
        enddo
        ierr=nf_fwrite (work, ncidhis, hisO, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxO))
          write(*,1) vname(1,indxO)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtO  ' MYID
      endif

! True vertical velocity (m/s).

      if (wrthis(indxW)) then
        do tile=0,NSUB_X*NSUB_E-1
          call wvlcty (tile, work)
        enddo
        ierr=nf_fwrite (work, ncidhis, hisW, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          write(*,1) vname(1,indxW)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtW  ' MYID
      endif

! Vertical viscosity coefficient.

      if (wrthis(indxAkv)) then
        ierr=nf_fwrite (Akv, ncidhis, hisAkv, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          write(*,1) vname(1,indxAkv)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAkv' MYID
      endif

! Vertical diffusion coefficient for potential temperature.

      if (wrthis(indxAkt)) then
        ierr=nf_fwrite (Akt(START_2D_ARRAY,0,itemp), ncidhis,
     &                                hisAkt, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          write(*,1) vname(1,indxAkt)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAkt' MYID
      endif
# ifdef SALINITY

! Vertical diffusion coefficient for salinity.

      if (wrthis(indxAks)) then
        ierr=nf_fwrite (Akt(START_2D_ARRAY,0,isalt), ncidhis,
     &                                hisAks, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          write(*,1) vname(1,indxAks)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAks' MYID
      endif
# endif
# ifdef LMD_KPP
! Depth of planetary boundary layer [m].

      if (wrthis(indxHbls)) then
        ierr=nf_fwrite (hbls(START_2D_ARRAY,nstp), ncidhis,
     &                             hisHbls, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbls))
          write(*,1) vname(1,indxHbls)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtHBL' MYID
      endif
# endif
# ifdef LMD_BKPP
! Thickness of KPP bottom boundary layer [m].

      if (wrthis(indxHbbl)) then
        ierr=nf_fwrite (hbbls(START_2D_ARRAY,nstp), ncidhis,
     &                              hisHbbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbbl))
          write(*,1) vname(1,indxHbbl)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtHbBL' MYID
      endif
# endif


# if (defined BIOLOGY_BEC || (defined BIOLOGY_NPZDOC && \
     defined OXYGEN && defined CARBON)) && ! defined OUTPUT_AVG_PHYS_BIO
!
! Write out averaged pH:
!
      ierr=nf_fwrite (pH_hist(START_2D_ARRAY), ncidhis,
     &                              hisPH, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPH_rst))
        write(*,1) vname(1,indxPH_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! Write out averaged pCO2:
!
      ierr=nf_fwrite (pCO2sw(START_2D_ARRAY), ncidhis,
     &                              hisPCO2, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPCO2_rst))
        write(*,1) vname(1,indxPCO2_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! Write out averaged pCO2air:
!
      ierr=nf_fwrite (pCO2air(START_2D_ARRAY), ncidhis,
     &                              hisPCO2air, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPCO2air_rst))
        write(*,1) vname(1,indxPCO2air_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
# endif /* BIOLOGY_BEC || (BIOLOGY_NPZDOC && OXYGEN && CARBON) */
#  if (defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC) && ! defined OUTPUT_AVG_PHYS_BIO
!
! PARinc:
!
      ierr=nf_fwrite (PARinc(START_2D_ARRAY), ncidhis,
     &     hisPARinc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPARinc_rst))
        write(*,1) vname(1,indxPARinc_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! PAR:
!
      ierr=nf_fwrite (PAR(START_2D_ARRAY,1), ncidhis,
     &                              hisPAR, record, r3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPAR_rst))
        write(*,1) vname(1,indxPAR_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
# endif /* BIOLOGY_BEC || BIOLOGY_NPZDOC */
#endif /* SOLVE3D */
  1   format(/1x, 'WRT_HIS ERROR while writing variable ''', A,
     &         ''' into history file.'  /11x, 'Time record:', I6,
     &         3x,'netCDF error code',i4 /11x,'Cause of error: ',
     &                                            A, 3x, A, i4)
      goto 100
  99  if (may_day_flag.eq.0) may_day_flag=3
 100  continue 
 
! Close or synchronize netCDF file to the disk to allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close (ncidhis)
      if (nrpfhis.gt.0 .and. record.ge.nrpfhis) ncidhis=-1
#else
      if (nrpfhis.gt.0 .and. record.ge.nrpfhis) then
        ierr=nf_close (ncidhis)
        ncidhis=-1
      else
        ierr=nf_sync(ncidhis)
      endif
#endif
      if (ierr .eq. nf_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode.eq.0) then
#endif
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i.lt.18 .and. tstring(i:i).eq.' ')
            i=i+1
          enddo
          write(*,'(7x,A,1x,A,2x,A,I7,1x,A,I4,A,I4,1x,A,I4)')
     &      'wrt_his :: wrote history, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', record, '/',  nrechis
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')      '### ERROR: wrt_his :: Cannot ',
     &                       'synchronize or close history file.'
        if (may_day_flag.eq.0) may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                           mynode+1, 2, ocean_grid_comm, ierr)
#endif
      return
      end
