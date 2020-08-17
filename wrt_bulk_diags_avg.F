#include "cppdefs.h"
#if defined BULK_DIAGS && defined AVERAGES
! Write time-averaged physical flux fields into averages netCDF file
      subroutine wrt_bulk_diags_avg

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
# include "forces.h"
# include "averages.h"

#if defined MPI && !defined PARALLEL_FILES
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
#endif
!
! Create/open averages file; write grid arrays, if so needed.
!
      call def_bulk_diags_avg(ncid_bulk_diags_avg, nrec_bulk_diags_avg, 
     &     ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(bulk_diags_avg_name)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index
      nrec_bulk_diags_avg=max(nrec_bulk_diags_avg,1)!!! nrec_bulk_diags_avg 
      if (nrpf_bulk_diags_avg.eq.0) then      !!! is already advanced by main.
        record=nrec_bulk_diags_avg
      else
        record=1+mod(nrec_bulk_diags_avg-1, nrpf_bulk_diags_avg)
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
      ibuff(4)=nrec_bulk_diags_avg

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int(ncid_bulk_diags_avg, bulk_diags_avgTstep,  
     &     start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,1) 'time_step', record,ierr MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged time
!
      ierr=nf_put_var1_double (ncid_bulk_diags_avg, bulk_diags_avgTime, 
     &     record, time_bulk_diags_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
! Barotropic mode variable: free-surface
! Always include this in the output so that depths of sigma levels
! can be computed
      ierr=nf_fwrite (zeta_bulk_diags_avg(START_2D_ARRAY), 
     &     ncid_bulk_diags_avg, bulk_diags_avgZ,record, r2dvar)
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxZ))
         write(*,1) vname(1,indxZ)(1:lvar), record, ierr
     &        MYID
         goto 99                !--> ERROR
      endif

       ierr=nf_fwrite(sustr_blk_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgSustr_blk, record,
# ifdef WND_AT_RHO_POINTS
     &                                                     r2dvar)
# else
     &                                                     u2dvar)
# endif
       ierr=nf_fwrite(svstr_blk_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                        avgSvstr_blk, record,
# ifdef WND_AT_RHO_POINTS
     &                                                     r2dvar)
# else
     &                                                     v2dvar)
# endif

       ierr=nf_fwrite(shflx_net_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgShflx_net, record,
     &                                                     r2dvar)

       ierr=nf_fwrite(shflx_lat_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgShflx_lat, record,
     &                                                     r2dvar)

       ierr=nf_fwrite(shflx_sen_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgShflx_sen, record,
     &                                                     r2dvar)

        ierr=nf_fwrite(shflx_rad_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                        avgShflx_rad, record,
     &                                                     r2dvar)

        ierr=nf_fwrite(swflx_emp_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                        avgSwflx_emp, record,
     &                                                     r2dvar)

       ierr=nf_fwrite(shflx_wwk_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgShflx_wwk, record,
     &                                                     r2dvar)
       ierr=nf_fwrite(surf_u_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgSurf_u, record,
     &                                                     r2dvar)
       ierr=nf_fwrite(surf_v_avg(START_2D_ARRAY),ncid_bulk_diags_avg,
     &                                       avgSurf_v, record,
     &                                                     r2dvar)


  1   format(/' WRT_BULK_DIAGS_AVG - ERROR while writing variable(',
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
      ierr=nf_close(ncid_bulk_diags_avg)
      if (nrpf_bulk_diags_avg.gt.0 .and. record.ge.nrpf_bulk_diags_avg) 
     &     ncid_bulk_diags_avg=-1
# else
      if (nrpf_bulk_diags_avg.gt.0 .and. 
     &     record.ge.nrpf_bulk_diags_avg) then
        ierr=nf_close(ncid_bulk_diags_avg)
        ncid_bulk_diags_avg=-1
      else
        ierr=nf_sync(ncid_bulk_diags_avg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
        mpi_master_only write(*,'(6x,A,2(A,I4,1x),A,I3)') 
     &        'WRT_BULK_DIAGS_AVG -- wrote ',
     &        'averaged fields into time record =', record, '/',
     &        nrec_bulk_diags_avg  MYID
      else
        write(*,'(/1x,2A/)') 
     &        'WRT_BULK_DIAGS_AVG ERROR: Cannot ',
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
#else /* BULK_DIAGS && AVERAGES */
      subroutine wrt_bulk_diags_avg_empty
      end
#endif /* BULK_DIAGS && AVERAGES */
