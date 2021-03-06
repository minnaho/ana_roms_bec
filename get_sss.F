#include "cppdefs.h"
#if defined SOLVE3D && (defined QCORRECTION || defined SALINITY_RESTORING) && !defined ANA_SST
 
      subroutine get_sss (ierr)     ! Read sea surface salinity
      implicit none
# define SST_DATA
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncvars.h"
      integer ierr, ncid, ifrc, lfrc, lenstr, nf_fread
#ifdef BULK_FLUX
     $,       iblk
      character(len=max_name_size) blkname
#else
      character(len=max_name_size) frcname
#endif
      character(len=3), parameter :: sss_var='SSS'
      character(lEN=8), parameter :: sss_var_time='sss_time'
      real cff

! Initialization:  Check, whether forcing netCDF file is already
!----------------  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables, and determine whether sss data exists
! as 2D field or scalar. After that determine whether there is
! cycling to reuse the input data and what is cycling period
! "sss_cycle", set initial cycling index "sss_ncycle", record index
! "sss_rec", and time index "itsss" and both time record bounds to
! large negative artificial values, so that it will trigger logic
! in the reading part below.
# ifdef VERBOSE
      write(*,*) 'enter get_sss' MYID
# endif

#ifdef BULK_FLUX
      iblk=max(sss_file_id,1)
 100  ncid=ncidbulk(iblk)
      blkname=blkfile(iblk)
      ierr = nf_noerr
      lfrc=lenstr(blkname)
#else
      ifrc=max(sss_file_id,1)
 100  ncid=ncidfrc(ifrc)
      frcname=frcfile(ifrc)
      lfrc=lenstr(frcname)
      ierr=nf_noerr
#endif


      if (iic.eq.ntstart) then
        if (ncid.eq.-1) then
#ifdef BULK_FLUX
          ierr=nf_open(blkname(1:lfrc), nf_nowrite, ncid)
#else
          ierr=nf_open(frcname(1:lfrc), nf_nowrite, ncid)
#endif
          if (ierr.eq.nf_noerr) then
#ifdef BULK_FLUX
            ncidbulk(iblk)=ncid
#else
            ncidfrc(ifrc)=ncid
#endif
          else
#ifdef BULK_FLUX
            write(*,'(/1x,4A/)')  '### ERROR: get_sss :: Cannot ',
     &              'open netCDF file ''', blkname(1:lfrc), '''.'
#else
            write(*,'(/1x,4A/)')  '### ERROR: get_sss :: Cannot ',
     &              'open netCDF file ''', frcname(1:lfrc), '''.'
#endif
          endif
        endif
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varid (ncid, sss_var, sss_id)
          if (ierr .eq. nf_noerr) then
#ifdef BULK_FLUX
            sss_file_id=iblk
            mpi_nonexit_warn write(*,'(1x,5A,I3)') 'Found ''',sss_var,
     &                ''' in ''', blkname(1:lfrc), ''', ifrc =', iblk
#else
            sss_file_id=ifrc
            mpi_nonexit_warn write(*,'(1x,5A,I3)') 'Found ''',sss_var,
     &                ''' in ''', frcname(1:lfrc), ''', ifrc =', ifrc
#endif


            ierr=nf_inq_varid (ncid, sss_var_time, sss_tid)
            if (ierr .eq. nf_noerr) then
              call set_cycle (ncid, sss_tid, ntsss, sss_cycle,
     &                              sss_ncycle, sss_rec, ierr)
              itsss=1
              sss_time(1)=-1.E+20
              sss_time(2)=-1.E+20
            else
#ifdef BULK_FLUX
              mpi_nonexit_warn write(*,2) sss_var_time, blkname(1:lfrc)
#else
              mpi_nonexit_warn write(*,2) sss_var_time, frcname(1:lfrc)
#endif
            endif

#ifdef BULK_FLUX
          elseif (iblk .lt. max_blk) then
            mpi_nonexit_warn write(*,1) sss_var, blkname(1:lfrc)
            iblk=iblk+1
            goto 100
          else
            write(*,2) sss_var, blkname(1:lfrc)
          endif
#else
# if defined BULK_FLUX
          elseif (iblk .lt. max_blk) then
            mpi_nonexit_warn write(*,1) sss_var, blkname(1:lfrc)
            iblk=iblk+1
# else
          elseif (ifrc .lt. max_frc) then
            mpi_nonexit_warn write(*,1) sss_var, frcname(1:lfrc)
            ifrc=ifrc+1
# endif
            goto 100
          else
# if defined BULK_FLUX
            write(*,2) sss_var, blkname(1:lfrc)
# else
            write(*,2) sss_var, frcname(1:lfrc)
# endif
          endif
#endif
        endif
      endif
  1   format(1x,  'WARNING: get_sss :: Cannot find variable ''',
     &                             A, ''' in file ''', A, '''.')
  2   format(1x, '### ERROR: get_sss :: Cannot find variable ''',
     &                             A, ''' in file ''', A, '''.')

! Read data from the file: Check if model time is bounded by the past
!----- ---- ---- --- ----- and the future data times: if not, increment
! record and cycling indices, flip time index and read a new portion of
! the data. Repeat until model time falls between the two data times.

      do while (sss_time(itsss).lt.time+0.5*dt .and.ierr.eq.nf_noerr)
        call advance_cycle (sss_cycle, ntsss, sss_ncycle,
     &                                 sss_rec,     ierr)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, sss_tid, sss_rec, cff)
          if (ierr .eq. nf_noerr) then
            itsss=3-itsss
            sss_time(itsss)=cff*day2sec + sss_cycle*sss_ncycle
            ierr=nf_fread (sssg(START_2D_ARRAY,itsss), ncid,
     &                                 sss_id, sss_rec, r2dvar)
            if (ierr .eq. nf_noerr) then
# ifdef MPI_SILENT_MODE
              if (mynode.eq.0) then
# endif
                write(*,'(7x,A,10x,A,F12.4,1x,A,I4,1x,A,I4)')
     &                  'get_sss :: read sea surface salinity',
     &                  'time =', cff, 'rec =', sss_rec  MYID
# ifdef MPI_SILENT_MODE
              endif
# endif
            else
              write(*,3) sss_var, sss_rec, nf_strerror(ierr)
            endif
          else
            write(*,3) sss_var_time, sss_rec, nf_strerror(ierr)
          endif
        else
# if defined BULK_FLUX
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &     '### ERROR: get_sss :: requested time record ',  sss_rec,
     &     'exeeds the last record', ntsss, 'available in forcing ',
     &     'netCDF file ''',   blkname(1:lfrc),   '''',  'tdays = ',
     &      tdays,        '  but the last available  sss_time =',
     &                                      sss_time(itsss)*sec2day
# else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &     '### ERROR: get_sss :: requested time record ',  sss_rec,
     &     'exeeds the last record', ntsss, 'available in forcing ',
     &     'netCDF file ''',   frcname(1:lfrc),   '''',  'tdays = ',
     &      tdays,        '  but the last available  sss_time =',
     &                                      sss_time(itsss)*sec2day
# endif

        endif
      enddo
  3   format(/1x,'### ERROR: get_sss: cannot read variable ''',
     &                         A, ''' for record ', I4/12x,A/)
# ifdef VERBOSE
      write(*,*) 'return from get_sss' MYID
# endif
      return
      end
#else
      subroutine get_sss_empty
      end
#endif /* QCORRECTION && !ANA_SST */
