!
!
#include "cppdefs.h"
#if defined SOLVE3D && defined OUTPUT_SST_AVG
                                     ! Write time-averaged SST fields  
      subroutine wrt_sst_avg         ! into separate SST averages
                                     ! netCDF file.
      implicit none
      integer ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), ibuff(4), nf_fwrite
     &            , itrc            
# if defined MPI && !defined PARALLEL_FILES
#  ifdef HOFFMAN2
      include "mpif.h"
#  else
# include "mpif.h"
#  endif
      integer status(MPI_STATUS_SIZE), blank
# endif
# include "param.h"
# include "scalars.h"
# include "averages.h"
# include "ncvars.h"
# include "netcdf.inc"

      character*70 text

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
      call def_sst_avg (ncid_sst_avg, nrec_sst_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(avgname)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index
      nrec_sst_avg=max(nrec_sst_avg,1)       !!! nrecrst is already
      if (nrpf_sst_avg.eq.0) then            !!! advanced by main.
        record=nrec_sst_avg
      else
        record=1+mod(nrec_sst_avg-1, nrpf_sst_avg)
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
      ibuff(4)=nrec_sst_avg

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (ncid_sst_avg, avg_sst_Tstep, start, count, 
     &     ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,1) 'time_step', record,ierr MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged time
!
      ierr=nf_put_var1_double (ncid_sst_avg, avg_sst_time, record, 
     &     time_sst_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr()
        write(*,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! SST
!
      ierr=nf_fwrite (sst_avg(START_2D_ARRAY), ncid_sst_avg,
     &     avgSST, record, r2dvar)
      if (ierr .ne. nf_noerr) then
         text='SST'
         lvar=lenstr(text)
         write(*,1) text(1:lvar),
     &        record, ierr MYID
         goto 99                !--> ERROR
      endif
  1   format(/' WRT_SST_AVG - ERROR while writing variable(',1x,a,1x,
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
      ierr=nf_close(ncid_sst_avg)
      if (nrpf_sst_avg.gt.0 .and. record.ge.nrpf_sst_avg) 
     &     ncid_sst_avg=-1
# else
      if (nrpf_sst_avg.gt.0 .and. record.ge.nrpf_sst_avg) then
        ierr=nf_close(ncid_sst_avg)
        ncid_sst_avg=-1
      else
        ierr=nf_sync(ncid_sst_avg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
        write(*,'(6x,A,2(A,I4,1x),A,I3)') 'WRT_SST_AVG -- wrote ',
     &            'averaged fields into time record =', record, '/',
     &             nrec_sst_avg  MYID
      else
        write(*,'(/1x,2A/)') 'WRT_SST_AVG ERROR: Cannot ',
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
#else  /* SOLVE3D && OUTPUT_SST_AVG */
      subroutine wrt_sst_avg_empty
      end
#endif /* SOLVE3D && OUTPUT_SST_AVG */
