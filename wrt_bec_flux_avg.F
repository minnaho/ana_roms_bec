#include "cppdefs.h"
#if defined BIOLOGY_BEC && defined BGC_FLUX_ANALYSIS && defined AVERAGES

       subroutine wrt_bec_flux_avg
       implicit none
#include "param.h"
#include "ncroms_bgcflux.h"
#include "bec_avg.h"
#include "netcdf.inc"
#include "scalars.h"
!for iic
#include "ncvars.h"

      integer ierr,  record, lstr, lvar, lenstr,itrc
     &   , start(2), count(2), ibuff(4), nf_fwrite
     &   ,varid
      REAL,DIMENSION(GLOBAL_2D_ARRAY,1)::WORK
   
#if defined MPI && !defined PARALLEL_FILES
# ifdef HOFFMAN2
      include "mpif.h"
# else
# include "mpif.h"
# endif
      integer status(MPI_STATUS_SIZE), blank
      if (mynode.gt.0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 4, MPI_COMM_WORLD, status, ierr)
#ifdef ASCII_OUTPUT
      write(*,*)"bec_flux_avg recv,mynode,blank",mynode,blank
#endif
#endif

      call def_bec_flux_avg (ncid_bgc_flux_avg, nrec_bgc_flux_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(bgc_flux_avg_name)
!                                                !!! WARNING: Here it is
! Set record within the file.                    !!! assumed that global
!                                                !!! restart record index
      nrec_bgc_flux_avg=max(nrec_bgc_flux_avg,1) !!! that global avg
      if (nrpf_bgc_flux_avg.eq.0) then           !!! record index is
        record=nrec_bgc_flux_avg                 !!! advanced by main.  
      else
        record=1+mod(nrec_bgc_flux_avg-1, nrpf_bgc_flux_avg)
      endif
!
! Time step and record indices.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrec_bgc_flux_avg

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int(ncid_bgc_flux_avg, bgc_flux_avgTstep,
     &     start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,1) 'time_step', record,ierr MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged time
!
      ierr=nf_put_var1_double (ncid_bgc_flux_avg, bgc_flux_avgTime,
     &     record, time_bgc_flux_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

! Barotropic mode variable: free-surface
! Always include this in the output so that depths of sigma levels
! can be computed
      ierr=nf_fwrite (zeta_bgc_flux_avg(START_2D_ARRAY), 
     &     ncid_bgc_flux_avg, bgc_flux_avgZ,record, r2dvar)
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxZ))
         write(*,1) vname(1,indxZ)(1:lvar), record, ierr
     &        MYID
         goto 99                !--> ERROR
      endif

       itrc=1
          ierr=nf_fwrite(ws_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                 vid_bec_flux_avg(itrc), record, r2dvar)
            if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif

       itrc=itrc+1
          ierr=nf_fwrite(xkw_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif

       itrc=itrc+1
          ierr=nf_fwrite(ap_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(schmidt_o2_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       

        itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(o2sat_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_o2_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# ifdef DMS_CYCLE
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_dms_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif

# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_n2o_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_n2_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(denitr_sed_avg(START_2D_ARRAY), 
     &         ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
#  ifdef N15_CYCLE
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_n2o_n15_avg(START_2D_ARRAY), 
     &         ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record,r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_n2_n15_avg(START_2D_ARRAY), 
     &         ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record,r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(denitr_sed_n15_avg(START_2D_ARRAY), 
     &         ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif

#  endif /* N15_CYCLE */
# endif /* NCYCLE_ANOXIC */

# ifdef NOT_USED___NEEDS_RECODING_IF_USED
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(ph_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(pco2surf_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dpco2_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# endif /* NOT_USED___NEEDS_RECODING_IF_USED */

       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(schmidt_co2_avg(START_2D_ARRAY),ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(co2star_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(dco2star_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fg_co2_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(iron_flux_avg(START_2D_ARRAY), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(PARinc_flux_avg(START_2D_ARRAY),ncid_bgc_flux_avg,
     &      vid_bec_flux_avg(itrc), record, r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname_bgcflux(1,itrc))
          write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &         ierr MYID
          goto 99               !--> ERROR
       endif
#endif
#ifdef SEDIMENT_BIOLOGY 
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(bot_flux_poc_avg(START_2D_ARRAY), 
     &      ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record, r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname_bgcflux(1,itrc))
          write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &         ierr MYID
          goto 99               !--> ERROR
       endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(bot_flux_caco3_avg(START_2D_ARRAY), 
     &      ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record, r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname_bgcflux(1,itrc))
          write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &         ierr MYID
          goto 99               !--> ERROR
       endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(bot_flux_si_avg(START_2D_ARRAY), 
     &      ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record, r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname_bgcflux(1,itrc))
          write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &         ierr MYID
          goto 99               !--> ERROR
       endif
#endif
# ifdef N15_CYCLE
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(bot_flux_pon15_avg(START_2D_ARRAY), 
     &      ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record, r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname_bgcflux(1,itrc))
          write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &         ierr MYID
          goto 99               !--> ERROR
       endif
#endif
# endif
#endif /* SEDIMENT_BIOLOGY */

!write 3d fluxes

       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(po4_restore_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg, vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(no3_restore_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(sio3_restore_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
        itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
        ierr=nf_fwrite(par_flux_avg(START_2D_ARRAY,1),ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(po4star_avg(START_2D_ARRAY,1), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(poc_flux_in_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
        itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(poc_prod_avg(START_2D_ARRAY,1), 
     &       ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(poc_remin_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(caco3_flux_in_avg(START_2D_ARRAY,1),
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(caco3_prod_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(caco3_remin_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sio2_flux_in_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
        itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sio2_prod_avg(START_2D_ARRAY,1), 
     &       ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(sio2_remin_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(dust_flux_in_avg(START_2D_ARRAY,1), 
     &      ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
         itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
         ierr=nf_fwrite(dust_remin_avg(START_2D_ARRAY,1), 
     &        ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(p_iron_flux_in_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
        itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
        ierr=nf_fwrite(p_iron_prod_avg(START_2D_ARRAY,1), 
     &       ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(p_iron_remin_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_sp_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_diat_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_tot_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_loss_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_loss_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(zoo_loss_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
       itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
       ierr=nf_fwrite(sp_agg_avg(START_2D_ARRAY,1), ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_agg_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(photoc_sp_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(f_ratio_sp_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(photoc_diat_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(f_ratio_diat_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(tot_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(no3_v_sp_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nh4_v_sp_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(no3_v_diat_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nh4_v_diat_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# ifdef PHAEOCYSTIS
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(no3_v_phaeo_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nh4_v_phaeo_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# endif /* PHAEOCYSTIS */
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(doc_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(doc_remin_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fe_scavenge_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_n_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_fe_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_po4_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_light_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_n_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_fe_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_po4_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_sio3_lim_avg(START_2D_ARRAY,1),
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_light_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(caco3_form_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
#ifdef PHAEOCYSTIS
! the "Nfix" variable is used for N limitation instead
#endif
          ierr=nf_fwrite(diaz_nfix_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_diaz_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diaz_loss_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(photoc_diaz_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diaz_p_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diaz_fe_lim_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diaz_light_lim_avg(START_2D_ARRAY,1),
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(fe_scavenge_rate_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(don_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
           ierr=nf_fwrite(don_remin_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dofe_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dofe_remin_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dop_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
           ierr=nf_fwrite(dop_remin_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(bsi_form_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(photofe_diaz_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(photofe_diat_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(photofe_sp_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nitrif_avg(START_2D_ARRAY,1),ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(j_o2_avg(START_2D_ARRAY,1),ncid_bgc_flux_avg,
     &                          vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# ifdef DMS_CYCLE     
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dmsp_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif       
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dmsp_p_uptake_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dmsp_lysis_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dms_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(dms_photlys_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(ammox_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(denitr_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(n2o_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(n2_prod_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
# endif /* NCYCLE_ANOXIC */

# ifdef N15_CYCLE
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(no3_v_sp_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nh4_v_sp_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(no3_v_diat_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nh4_v_diat_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_sp_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_diat_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(graze_diaz_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_loss_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_loss_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diaz_loss_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(zoo_loss_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(sp_agg_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diat_agg_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(diaz_Nfix_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(pon_remin_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(don_prod_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(don_remin_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(nitrif_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif
#  ifdef NCYCLE_ANOXIC
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(ammox_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(denitr_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(n2o_prod_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
          itrc=itrc+1
#ifndef REDUCED_BGC_OUTPUT
          ierr=nf_fwrite(n2_prod_n15_avg(START_2D_ARRAY,1), 
     &         ncid_bgc_flux_avg,vid_bec_flux_avg(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_bgcflux(1,itrc))
            write(*,1) vname_bgcflux(1,itrc)(1:lvar), record,
     &                      ierr MYID
            goto 99                                       !--> ERROR
          endif
#endif          
#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */


  1   format(/1x,'WRT_bec_flux_avg ERROR while writing variable ''', A,
     &          ''' into avg file.', /11x, 'Time record:',
     &                   I6,3x,'netCDF error code',i4,3x,a,i4)
       goto 100 
  99  may_day_flag=3
 100  continue
!
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!
#if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close (ncid_bgc_flux_avg)
      if (nrpf_bgc_flux_avg.gt.0 .and. record.ge.nrpf_bgc_flux_avg) 
     &     ncid_bgc_flux_avg=-1
#else
      if (nrpf_bgc_flux_avg.gt.0 .and. record.ge.nrpf_bgc_flux_avg) then
        ierr=nf_close (ncid_bgc_flux_avg)
        ncid_bgc_flux_avg=-1
      else
        ierr=nf_sync(ncid_bgc_flux_avg)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        mpi_master_only write(*,'(6x,A,2(A,I4,1x),A,I3)') 
     &             'WRT_BEC_FLUX_AVG -- wrote ',
     &            'avg fields into time record =', record, '/',
     &             nrec_bgc_flux_avg  MYID
      else
        write(*,'(/1x,2A/)') 'WRT_BEC_FLUX_AVG ERROR: Cannot ',
     &             'synchronize/close avg netCDF file.'
        may_day_flag=3
      endif

#if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) then
         blank=40+mynode
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        4, MPI_COMM_WORLD,  ierr)
#ifdef ASCII_OUTPUT
        write(*,*)"wrt_bec_flux_avg send,mynode,blank",mynode,blank
#endif
      endif
#endif

      return
       end subroutine wrt_bec_flux_avg
#else
       subroutine ecosys_wrt_bec_flux_avg_empty ()
       
         end
#endif  /* BIOLOGY_BEC && BGC_FLUX_ANALYSIS && AVERAGES */
