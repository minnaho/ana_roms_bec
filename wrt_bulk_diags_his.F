!   $Log: wrt_bulk_diags_his.F,v $
!   Revision 1.1  2002/04/19 00:59:56  frenzel
!   New routine: Output of physical flux histories.
!
!
!
#include "cppdefs.h"
#ifdef BULK_DIAGS
! Write physical flux fields into histories netCDF file
      subroutine wrt_bulk_diags_his
      implicit none
#include "netcdf.inc"
      integer ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), ibuff(4), nf_fwrite, iflux, iflux_out

# include "param.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "ncvars.h"
# include "forces.h"

#if defined MPI && !defined PARALLEL_FILES
# ifdef HOFFMAN2
      include "mpif.h"
# else
# include "mpif.h"
# endif
      integer status(MPI_STATUS_SIZE), blank

      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
#endif
!
! Create/open histories file; write grid arrays, if so needed.
!
      call def_bulk_diags_his(ncid_bulk_diags_his, nrec_bulk_diags_his, 
     &     ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(bulk_diags_his_name)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index
      nrec_bulk_diags_his=max(nrec_bulk_diags_his,1)!!! nrec_bulk_diags_his 
      if (nrpf_bulk_diags_his.eq.0) then      !!! is already advanced by main.
        record=nrec_bulk_diags_his
      else
        record=1+mod(nrec_bulk_diags_his-1, nrpf_bulk_diags_his)
      endif
!
! Write out variables:
! ----- --- ------------- ----------
!
! Time step and record indices.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrec_bulk_diags_his

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int(ncid_bulk_diags_his, bulk_diags_hisTstep,  
     &     start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(*,1) 'time_step', record,ierr, nf_strerror(ierr)
     &        MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_double (ncid_bulk_diags_his, bulk_diags_hisTime, 
     &     record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,1) vname(1,indxTime)(1:lvar), record, ierr, 
     &       nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif
! Barotropic mode variable: free-surface
! Always include this in the output so that depths of sigma levels
! can be computed
      ierr=nf_fwrite (zeta(START_2D_ARRAY,knew), 
     &     ncid_bulk_diags_his, bulk_diags_hisZ,record, r2dvar)
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxZ))
         write(*,1) vname(1,indxZ)(1:lvar), record, ierr,
     &        nf_strerror(ierr) MYID
         goto 99                !--> ERROR
      endif


# if defined BULK_DIAGS
        ierr=nf_fwrite(sustr_blk, ncid_bulk_diags_his, hisSustr_blk 
     &                 ,record,
# ifdef WND_AT_RHO_POINTS        
     &                                                     r2dvar)
# else
     &                                                     u2dvar)
#endif
        ierr=nf_fwrite(svstr_blk, ncid_bulk_diags_his, hisSvstr_blk
     &                 ,record,
# ifdef WND_AT_RHO_POINTS
     &                                                     r2dvar)
# else
     &                                                     v2dvar)
#endif     
!
!  Write out surface heat flux component (degC/s)
!
        ierr=nf_fwrite(shflx_net, ncid_bulk_diags_his, hisShflx_net, 
     &                                               record, r2dvar)

        ierr=nf_fwrite(shflx_lat, ncid_bulk_diags_his, hisShflx_lat,
     &                                               record, r2dvar)

        ierr=nf_fwrite(shflx_sen, ncid_bulk_diags_his, hisShflx_sen,
     &                                               record, r2dvar)

        ierr=nf_fwrite(shflx_rad, ncid_bulk_diags_his, hisShflx_rad,
     &                                               record, r2dvar)

         ierr=nf_fwrite(swflx_emp, ncid_bulk_diags_his, hisSwflx_emp,
     &                                               record, r2dvar)

        ierr=nf_fwrite(shflx_wwk, ncid_bulk_diags_his, hisShflx_wwk,
     &                                               record, r2dvar)
        ierr=nf_fwrite(surf_u, ncid_bulk_diags_his, hisSurf_u,
     &                                               record, r2dvar)
        ierr=nf_fwrite(surf_v, ncid_bulk_diags_his, hisSurf_v,
     &                                               record, r2dvar)


# endif /* BULK_DIAGS */


  1   format(/' WRT_PHYS_FLUX_HIS - ERROR while writing variable(',
     &     1x,a,1x,
     &     ')into phys flux histories file.',/,11x,
     &     'Time record:',
     &     i6,3x,'netCDF error code',i4,3x,a,1x,a,i4)
      goto 100
  99  may_day_flag=3
 100  continue


! to access data immediately after it is written.
!
# if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close(ncid_bulk_diags_his)
      if (nrpf_bulk_diags_his.gt.0 .and. record.ge.nrpf_bulk_diags_his) 
     &     ncid_bulk_diags_his=-1
# else
      if (nrpf_bulk_diags_his.gt.0 .and. 
     &     record.ge.nrpf_bulk_diags_his) then
        ierr=nf_close(ncid_bulk_diags_his)
        ncid_bulk_diags_his=-1
      else
        ierr=nf_sync(ncid_bulk_diags_his)
      endif
# endif
      if (ierr .eq. nf_noerr) then
        mpi_master_only write(*,'(6x,A,2(A,I4,1x),A,I3)') 
     &        'WRT_BULK_DIAGS_HIS -- wrote ',
     &        'bulk_diags fields into time record =', record, '/',
     &        nrec_bulk_diags_his  MYID
      else
        write(*,'(/1x,2A/)') 
     &        'WRT_BULK_DIAGS_HIS ERROR: Cannot ',
     &        'synchronize/close phys fluxes netCDF file.'
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
#else /* BULK_DIAGS */
      subroutine wrt_bulk_diags_his_empty
      end
#endif /* BULK_DIAGS */
