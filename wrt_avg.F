#include "cppdefs.h"
#ifdef AVERAGES

      subroutine wrt_avg             ! Write time-averaged fields
                                     ! into averages netCDF file.
      implicit none
      integer record, i, ierr, lstr, lvar, lenstr, nf_fwrite
# ifdef SOLVE3D
     &             , itrc
# endif
      character(len=18) tstring
# include "param.h"
# include "scalars.h"
# include "averages.h"
# include "ncvars.h"
# include "netcdf.inc"
#ifdef SEDIMENT_BIOLOGY
!# include "bgcflux.h"
#endif
#ifdef BIOLOGY_BEC
!# include "ecosys_bec.h"
# include "bec_avg.h"
#endif
# if defined MPI && !defined PARALLEL_FILES
# ifdef HOFFMAN2
      include "mpif.h"
# else
#  include "mpif.h"
# endif
      integer status(MPI_STATUS_SIZE), blank
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,  3,
     &                       ocean_grid_comm, status, ierr)
      endif
# endif

! Create/open averages file; write grid arrays, if so needed,
! check for errors, then set record within the file "record".

      call def_avg (nrecavg, ierr)
      lstr=lenstr(avgname)
      if (ierr .ne. nf_noerr) goto 99       !!! WARNING: Here it is
      nrecavg=max(nrecavg,1)                !!! assumed that global
      if (nrpfavg.eq.0) then                !!! averages record index
        record=nrecavg                      !!! "nrecavg" is already
      else                                  !!! advanced by main.
        record=1+mod(nrecavg-1, nrpfavg)
      endif

! Write out time-averaged variables:
! ----- --- ------------- ----------
! Averaged time

      ierr=nf_put_var1_FTYPE (ncidavg, avgTime, record, time_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) vname(1,indxTime)(1:lvar), record MYID
        goto 99                                           !--> ERROR
      endif

! Barotropic mode variables: free-surface and 2D momentum components
! in XI-,ETA-directions.

      if (wrtavg(indxZ)) then
        ierr=nf_fwrite (zeta_avg(START_2D_ARRAY), ncidavg, avgZ,
     &                                            record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          write(*,1) vname(1,indxTime)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
 
      if (wrtavg(indxUb)) then
        ierr=nf_fwrite (ubar_avg(START_2D_ARRAY), ncidavg, avgUb,
     &                                            record, u2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          write(*,1) vname(1,indxUb)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
 
      if (wrtavg(indxVb)) then
        ierr=nf_fwrite (vbar_avg(START_2D_ARRAY), ncidavg, avgVb,
     &                                            record, v2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          write(*,1) vname(1,indxVb)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
# ifdef SOLVE3D

! 3D momentum components in XI- and ETA-directions.

      if (wrtavg(indxU)) then
        ierr=nf_fwrite (u_avg(START_2D_ARRAY,1), ncidavg, avgU,
     &                                           record, u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          write(*,1) vname(1,indxU)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
 
      if (wrtavg(indxV)) then
        ierr=nf_fwrite (v_avg(START_2D_ARRAY,1), ncidavg, avgV,
     &                                          record, v3dvar)
        if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxV))
          write(*,1) vname(1,indxV)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif

! Tracer variables.

# if defined OUTPUT_AVG_PHYS_BIO
          do itrc=1,1+ntrc_salt
# else
      do itrc=1,NT
# endif
        if (wrtavg(indxT+itrc-1)) then
          ierr=nf_fwrite (t_avg(START_2D_ARRAY,1,itrc), ncidavg,
     &                              avgT(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            write(*,1) vname(1,indxT+itrc-1)(1:lvar),record MYID
            goto 99                                       !--> ERROR
          endif
        endif
      enddo
#if defined SEDIMENT_BIOLOGY && !defined OUTPUT_AVG_PHYS_BIO
!
! Sediment concentrations
!
      do itrc=1,NT_sed
        if (wrtavg(indxSedFirst+itrc-1)) then
          ierr=nf_fwrite (t_sed_avg(START_2D_ARRAY,itrc), ncidavg,
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
#endif /* SEDIMENT_BIOLOGY... */

! Density anomaly.

      if (wrtavg(indxR)) then
        ierr=nf_fwrite (rho_avg(START_2D_ARRAY,1), ncidavg, avgR,
     &                                             record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          write(*,1) vname(1,indxR)(1:lvar), record  MYID
          goto 99                                         !--> ERROR
        endif
      endif

! S-coordinate omega vertical velocity (m/s).

      if (wrtavg(indxO)) then
        ierr=nf_fwrite (w_avg, ncidavg, avgO, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxO))
          write(*,1) vname(1,indxO)(1:lvar), record  MYID
          goto 99                                         !--> ERROR
        endif
      endif

! True vertical velocity (m/s).

      if (wrtavg(indxW)) then
        ierr=nf_fwrite (wvl_avg, ncidavg, avgW, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          write(*,1) vname(1,indxW)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif

! Vertical eddy viscosity coefficient for momentum.

      if (wrtavg(indxAkv)) then
        ierr=nf_fwrite (akv_avg(START_2D_ARRAY,0), ncidavg,
     &                                avgAkv, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          write(*,1) vname(1,indxAkv)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif

! Vertical heat conductivity

      if (wrtavg(indxAkt)) then
        ierr=nf_fwrite (akt_avg(START_2D_ARRAY,0), ncidavg,
     &                                avgAkt, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          write(*,1) vname(1,indxAkt)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
#  ifdef SALINITY

! Vertical diffusion coefficient for salinity.

      if (wrtavg(indxAks)) then
        ierr=nf_fwrite (aks_avg(START_2D_ARRAY,0), ncidavg,
     &                                avgAks, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          write(*,1) vname(1,indxAks)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
#  endif
#  ifdef LMD_KPP

! Depth of planetary boundary layer [m].

      if (wrtavg(indxHbls)) then
        ierr=nf_fwrite (hbl_avg(START_2D_ARRAY), ncidavg,
     &                            avgHbls, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbls))
          write(*,1) vname(1,indxHbls)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
#  endif
#  ifdef LMD_BKPP

! Thickness of KPP bottom boundary layer [m].

      if (wrtavg(indxHbbl)) then
        ierr=nf_fwrite (hbbl_avg(START_2D_ARRAY), ncidavg,
     &                              avgHbbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbbl))
          write(*,1) vname(1,indxHbbl)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
#  endif

# if ((defined BIOLOGY_NPZDOC && defined CARBON) || defined BIOLOGY_BEC) && ! defined OUTPUT_AVG_PHYS_BIO
!
! Write out averaged pH:
!
      ierr=nf_fwrite (pH_avg(START_2D_ARRAY), ncidavg,
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
      ierr=nf_fwrite (pCO2_avg(START_2D_ARRAY), ncidavg,
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
      ierr=nf_fwrite (pCO2air_avg(START_2D_ARRAY), ncidavg,
     &                              avgPCO2air, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPCO2air_rst))
        write(*,1) vname(1,indxPCO2air_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
#  endif /* BIOLOGY_NPZDOC ... */
#  if (defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC) && ! defined OUTPUT_AVG_PHYS_BIO
!
! PARinc:
!
      ierr=nf_fwrite (PARinc_avg(START_2D_ARRAY), ncidavg,
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
      ierr=nf_fwrite (PAR_avg(START_2D_ARRAY,1), ncidavg,
     &                              avgPAR, record, r3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxPAR_rst))
        write(*,1) vname(1,indxPAR_rst)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
# ifdef OUTPUT_SPECDOM
! frequency of small phytoplankton dominance at the surface
      ierr=nf_fwrite (dom_sp_sfc_avg(START_2D_ARRAY), ncidavg,
     &     avgFreqDomSP_sfc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        write(*,1) vname(1,indxFreqDomSP_sfc)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diatom dominance at the surface
      ierr=nf_fwrite (dom_diat_sfc_avg(START_2D_ARRAY), ncidavg,
     &     avgFreqDomDIAT_sfc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        write(*,1) vname(1,indxFreqDomDIAT_sfc)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diazotroph dominance at the surface
      ierr=nf_fwrite (dom_diaz_sfc_avg(START_2D_ARRAY), ncidavg,
     &     avgFreqDomDIAZ_sfc, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        write(*,1) vname(1,indxFreqDomDIAZ_sfc)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of small phytoplankton dominance, integrated
      ierr=nf_fwrite (dom_sp_int_avg(START_2D_ARRAY), ncidavg,
     &     avgFreqDomSP_int, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        write(*,1) vname(1,indxFreqDomSP_int)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diatom dominance, integrated
      ierr=nf_fwrite (dom_diat_int_avg(START_2D_ARRAY), ncidavg,
     &     avgFreqDomDIAT_int, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        write(*,1) vname(1,indxFreqDomDIAT_int)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
! frequency of diazotroph dominance, integrated
      ierr=nf_fwrite (dom_diaz_int_avg(START_2D_ARRAY), ncidavg,
     &     avgFreqDomDIAZ_int, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        write(*,1) vname(1,indxFreqDomDIAZ_int)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                         !--> ERROR
      endif
# endif /* OUTPUT_SPECDOM */
# endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC */
# endif /* SOLVE3D */

  1   format(/1x, '### ERROR: wrt_avg :: Cannot write variable ''',
     &             A, ''' into averages file, rec =', i6, 3x,A,i4)
      goto 100
  99  if (may_day_flag.eq.0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to disk to allow other 
! processes to access data immediately after it is written.

# if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close(ncidavg)
      if (nrpfavg.gt.0 .and. record.ge.nrpfavg) ncidavg=-1
# else
      if (nrpfavg.gt.0 .and. record.ge.nrpfavg) then
        ierr=nf_close(ncidavg)
        ncidavg=-1
      else
        ierr=nf_sync(ncidavg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
# ifdef MPI_SILENT_MODE
        if (mynode.eq.0) then
# endif
          write(tstring,'(F18.8)') time_avg*sec2day
          i=1
          do while (i.lt.18 .and. tstring(i:i).eq.' ')
            i=i+1
          enddo
          write(*,'(7x,2A,1x,A,1x,A,I4,A,I4,1x,A,I4)') 'wrt_avg :: ',
     &      'wrote averages for time_avg =', tstring(i:i+8), 'rec =',
     &       record, '/', nrecavg  MYID
# ifdef MPI_SILENT_MODE
        endif
# endif
      else
        write(*,'(/1x,2A/)')      '### ERROR: wrt_avg :: Cannot ',
     &                         'synchronize/close averages file.'
        if (may_day_flag.eq.0) may_day_flag=3
      endif
# if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER, 
     &                           mynode+1, 3, ocean_grid_comm, ierr)
# endif
      return
      end
#else
      subroutine wrt_avg_empty
      end
#endif /* AVERAGES */
