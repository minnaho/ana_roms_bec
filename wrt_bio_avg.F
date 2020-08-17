!
! $Id: wrt_avg.F,v 1.9 2005/09/23 15:56:31 pmarches Exp $
!
#include "cppdefs.h"
#if defined AVERAGES && defined OUTPUT_AVG_PHYS_BIO
                                     ! Write time-averaged 
      subroutine wrt_bio_avg        ! fields into averages
                                     ! netCDF file.
      implicit none
      integer ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), ibuff(4), nf_fwrite
# ifdef SOLVE3D
     &            , itrc            
# endif
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
# include "averages.h"
# include "ncvars.h"
# include "netcdf.inc"
#ifdef SEDIMENT_BIOLOGY
# include "bgcflux.h"
# include "bec_avg.h"
#endif

#if defined MPI && !defined PARALLEL_FILES
      if (mynode.gt.0) then
#if defined UCLA_COUPLED
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_MODEL, status, ierr)
#else
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, ocean_grid_comm, status, ierr)
#endif
      endif
#endif
!
! Create/open averages file; write grid arrays, if so needed.
!
      call def_bio_avg (ncid_bio_avg, nrec_bio_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(bio_avgname)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index
      nrec_bio_avg=max(nrec_bio_avg,1)       !!! nrecrst is already
      if (nrpf_bio_avg.eq.0) then            !!! advanced by main.
        record=nrec_bio_avg
      else
        record=1+mod(nrec_bio_avg-1, nrpf_bio_avg)
      endif
!
! Write out time-averaged variables:
! ----- --- ------------- ----------
!
! Time step and record indices.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=0 ! nrechis
      ibuff(4)=nrec_bio_avg

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (ncid_bio_avg, avg_bio_Tstep, start, count, 
     &     ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,1) 'time_step', record,ierr MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged time
!
      ierr=nf_put_var1_double (ncid_bio_avg, avg_bio_time, record, 
     &     time_bio_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.
!
      if (wrtavg(indxZ)) then
        ierr=nf_fwrite (zeta_bio_avg(START_2D_ARRAY), ncid_bio_avg,
     &        avg_bio_Z, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          write(*,1) vname(1,indxTime)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif

# ifdef SOLVE3D
!
! Tracer variables.
!
      do itrc=2+ntrc_salt,NT
        if (wrtavg(indxT+itrc-1)) then
          ierr=nf_fwrite (t_avg(START_2D_ARRAY,1,itrc), ncid_bio_avg,
     &                              avgT(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            write(*,1) vname(1,indxT+itrc-1)(1:lvar),
     &                      record, ierr MYID
            goto 99                                       !--> ERROR
          endif
        endif
      enddo
#ifdef SEDIMENT_BIOLOGY
!
! Sediment concentrations
!
      do itrc=1,NT_sed
        if (wrtavg(indxSedFirst+itrc-1)) then
          ierr=nf_fwrite (t_sed_avg(START_2D_ARRAY,itrc), ncid_bio_avg,
     &                                avgTsed(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            write(*,1) vname(1,indxSedFirst+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
CR      write(*,*) 'wrt_avg:  wrtT_sed  ' MYID
        endif
      enddo
#endif /* SEDIMENT_BIOLOGY */

!
! Write out averaged pH:
!
      ierr=nf_fwrite (pH_avg(START_2D_ARRAY), ncid_bio_avg,
     &                              avgPH, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPH_rst))
        write(*,1) vname(1,indxPH_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! Write out averaged pCO2:
!
      ierr=nf_fwrite (pCO2_avg(START_2D_ARRAY), ncid_bio_avg,
     &                              avgPCO2, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPCO2_rst))
        write(*,1) vname(1,indxPCO2_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! Write out averaged pCO2air:
!
      ierr=nf_fwrite (pCO2air_avg(START_2D_ARRAY), ncid_bio_avg,
     &                              avgPCO2air, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPCO2air_rst))
        write(*,1) vname(1,indxPCO2air_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! PARinc:
!
      ierr=nf_fwrite (PARinc_avg(START_2D_ARRAY), ncid_bio_avg,
     &                              avgPARinc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPARinc_rst))
        write(*,1) vname(1,indxPARinc_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
!
! PAR:
!
      ierr=nf_fwrite (PAR_avg(START_2D_ARRAY,1), ncid_bio_avg,
     &                              avgPAR, record, r3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPAR_rst))
        write(*,1) vname(1,indxPAR_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
# ifdef OUTPUT_SPECDOM
! frequency of small phytoplankton dominance at the surface
      ierr=nf_fwrite (dom_sp_sfc_avg(START_2D_ARRAY), ncid_bio_avg,
     &     avgFreqDomSP_sfc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        write(*,1) vname(1,indxFreqDomSP_sfc)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diatom dominance at the surface
      ierr=nf_fwrite (dom_diat_sfc_avg(START_2D_ARRAY), ncid_bio_avg,
     &     avgFreqDomDIAT_sfc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        write(*,1) vname(1,indxFreqDomDIAT_sfc)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diazotroph dominance at the surface
      ierr=nf_fwrite (dom_diaz_sfc_avg(START_2D_ARRAY), ncid_bio_avg,
     &     avgFreqDomDIAZ_sfc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        write(*,1) vname(1,indxFreqDomDIAZ_sfc)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of small phytoplankton dominance, integrated
      ierr=nf_fwrite (dom_sp_int_avg(START_2D_ARRAY), ncid_bio_avg,
     &     avgFreqDomSP_int, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        write(*,1) vname(1,indxFreqDomSP_int)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diatom dominance, integrated
      ierr=nf_fwrite (dom_diat_int_avg(START_2D_ARRAY), ncid_bio_avg,
     &     avgFreqDomDIAT_int, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        write(*,1) vname(1,indxFreqDomDIAT_int)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diazotroph dominance, integrated
      ierr=nf_fwrite (dom_diaz_int_avg(START_2D_ARRAY), ncid_bio_avg,
     &     avgFreqDomDIAZ_int, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        write(*,1) vname(1,indxFreqDomDIAZ_int)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
# endif /* OUTPUT_SPECDOM */


#  ifdef BIOLOGY
!
!  Write out depth of the euphotic layer (m).
!
      if (wrtavg(indxHel)) then
        ierr=nf_fwrite(hel_avg, ncid_bio_avg, avgHel, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHel))
          write(*,1) vname(1,indxHel)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
!
!  Write out Chlorophyll to Carbon ratio (m).
!
#   ifdef BIO_NPZD
      if (wrtavg(indxChC)) then
        ierr=nf_fwrite(theta_avg, ncid_bio_avg, avgChC, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxChC))
          write(*,1) vname(1,indxChC)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
#     ifdef OXYGEN
      if (wrtavg(indxU10)) then
        ierr=nf_fwrite (u10_avg, ncid_bio_avg, avgU10, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxU10))
          write(*,1) vname(1,indxU10)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxKvO2)) then
        ierr=nf_fwrite (Kv_O2_avg, ncid_bio_avg, avgKvO2, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKvO2))
          write(*,1) vname(1,indxKvO2)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxO2sat)) then
        ierr=nf_fwrite (O2satu_avg, ncid_bio_avg, avgO2sat, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxO2sat))
          write(*,1) vname(1,indxO2sat)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
      endif
#     endif /* OXYGEN */
#   elif defined BIO_N2P2Z2D2 && defined VAR_CHL_C
      if (wrtavg(indxChC1)) then
        ierr=nf_fwrite(theta1_avg, ncid_bio_avg, avgChC1, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxChC1))
          write(*,1) vname(1,indxChC1)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif

      if (wrtavg(indxChC2)) then
        ierr=nf_fwrite(theta2_avg, ncid_bio_avg, avgChC2, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxChC2))
          write(*,1) vname(1,indxChC2)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
#   endif
#  endif /* BIOLOGY */
# endif /* SOLVE3D */
  1   format(/' WRT_BIO_AVG - ERROR while writing variable(',1x,a,1x,
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
      ierr=nf_close(ncid_bio_avg)
      if (nrpf_bio_avg.gt.0 .and. record.ge.nrpf_bio_avg) 
     &     ncid_bio_avg=-1
# else
      if (nrpf_bio_avg.gt.0 .and. record.ge.nrpf_bio_avg) then
        ierr=nf_close(ncid_bio_avg)
        ncid_bio_avg=-1
      else
        ierr=nf_sync(ncid_bio_avg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
        write(*,'(6x,A,2(A,I4,1x),A,I3)') 'WRT_BIO_AVG -- wrote ',
     &            'averaged fields into time record =', record, '/',
     &             nrec_bio_avg  MYID
      else
        write(*,'(/1x,2A/)') 'WRT_BIO_AVG ERROR: Cannot ',
     &             'synchronize/close averages netCDF file.'
        may_day_flag=3
      endif

# if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) then
#if defined UCLA_COUPLED
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_MODEL,  ierr)
#else
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, ocean_grid_comm,  ierr)
#endif
      endif
# endif
      return
      end
#else
      subroutine wrt_bio_avg_empty
      return
      end
#endif /* AVERAGES && OUTPUT_AVG_PHYS_BIO */
