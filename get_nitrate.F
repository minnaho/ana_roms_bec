#include "cppdefs.h"
#ifdef BIOLOGY_BEC

      subroutine get_nitrate(ierr)
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncvars.h"
      real cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, ifrc, ncid
      character*(max_name_size) frcname
      character(len=7), parameter :: nitrate_var='nitrate'

!
! Initialization: Check, whether forcing netCDF file is already
!===============  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables, and determine whether sst data exists
! as 2D field or scalar. After that determine whether there is
! cycling to reuse the input data and what is cycling period
! "nitrate_cycle", set initial cycling index "nitrate_ncycle", record index
! "nitrate_rec", and time index "itnitrate" and both time record bounds to
! large negative artificial values, so that it will trigger logic in
! reading part below.
!
# if defined BULK_FLUX
      ifrc=max(nitrate_file_id,1)
 100  ncid=ncidbulk(ifrc)
      frcname=blkfile(ifrc)
# else
      ifrc=max(nitrate_file_id,1)
 100  ncid=ncidfrc(ifrc)
      frcname=frcfile(ifrc)
# endif
      ierr=nf_noerr

      if (iic.eq.ntstart) then
        lstr=lenstr(frcname)
        if (ncid.eq.-1) then
# if defined BULK_FLUX
          ierr=nf_open(frcname(1:lstr), nf_nowrite, ncid)
          if (ierr.eq.nf_noerr) then
            ncidbulk(ifrc)=ncid
          else
            write(*,'(/1x,4A/)') '### ERROR ',
     &                   'in get_nitrate: can not open netCDF file ''',
     &                                       frcname(1:lstr), '''.'
          endif
# else
          ierr=nf_open(frcname(1:lstr), nf_nowrite, ncid)
          if (ierr.eq.nf_noerr) then
            ncidfrc(ifrc)=ncid
          else
            write(*,'(/1x,4A/)') '### ERROR ',
     &                   'in get_nitrate: can not open netCDF file ''',
     &                                       frcname(1:lstr), '''.'
          endif
# endif
        endif
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varid (ncid, 'nitrate_time', nitrate_tid)
          if (ierr .eq. nf_noerr) then
            nitrate_file_id = ifrc
            lvar=lenstr(nitrate_var)
            ierr=nf_inq_varid (ncid, nitrate_var(1:lvar),
     &                                                   nitrate_id)
            if (ierr .eq. nf_noerr) then
              call set_cycle (ncid, nitrate_tid, ntnitrate, nitrate_cycle,
     &                                 nitrate_ncycle, nitrate_rec, ierr)
              itnitrate=1
              nitrate_time(1)=-1.E+20
              nitrate_time(2)=-1.E+20
            else
              write(*,1) nitrate_var(1:lvar),
     &                                  frcname(1:lstr)
            endif
          else
            mpi_nonexit_warn write(*,1) 'nitrate_time', frcname(1:lstr)
# if defined BULK_FLUX
            if (ifrc .lt. max_blk_file) then
# else
            if (ifrc .lt. max_frc_file) then
# endif
              ifrc=ifrc+1
              goto 100
            endif
          endif
        endif
      endif
  1   format(/1x,'WARNING: get_nitrate cannot find variable ''',
     &     A, ''''/7x, 'in forcing netCDF file ''', A, '''.')
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
      do while (nitrate_time(itnitrate).lt.time+0.5*dt .and. ierr.eq.nf_noerr)
        call advance_cycle (nitrate_cycle, ntnitrate, nitrate_ncycle,
     &                                 nitrate_rec,     ierr)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, nitrate_tid, nitrate_rec, cff)
          if (ierr .eq. nf_noerr) then
            itnitrate=3-itnitrate
            nitrate_time(itnitrate)=cff*day2sec + nitrate_cycle*nitrate_ncycle
            ierr=nf_fread (nitrateg(START_2D_ARRAY,itnitrate), ncid,
     &                                 nitrate_id, nitrate_rec, r2dvar)
            if (ierr .eq. nf_noerr) then
              mpi_master_only write(*,
     &             '(6x,2A,G12.4,1x,A,I4)') 'get_nitrate    --',
     &                 ' read nitrate                       for time =',
     &                   cff MYID
            else
              write(*,2) 'nitrate', nitrate_rec
            endif
          else
            write(*,2) 'nitrate_time', nitrate_rec
          endif
        else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &   'ERROR in get_nitrate: requested time record ',      nitrate_rec,
     &   'exeeds the last record', ntnitrate, 'available in forcing ',
     &   'netCDF file ''',   frcname(1:lstr),   '''',  'tdays = ',
     &    tdays,        '  but the last available  nitrate_time =',
     &                                    nitrate_time(itnitrate)*sec2day
        endif
      enddo
  2   format(/1x,'ERROR in get_nitrate: cannot read variable ''', A,
     &                                    ''' for record ', I4/)
      return
      end
 
      subroutine set_nitrate_tile (Istr,Iend,Jstr,Jend,ierr)
!
! Set-up sea surface nitrate data for current tile.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2,ierr
      real cff, cff1,cff2, cff3,cff4, val1

# include "param.h"
# include "forces.h"
# include "scalars.h"
!
# include "compute_extended_bounds.h"
 

! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.
!
 
      it1=3-itnitrate
      it2=itnitrate
      cff1=nitrate_time(it2)-time
      cff2=time-nitrate_time(it1)


      if (ZEROTH_TILE  .and. ( cff1.lt.dt )) synchro_flag=.true.
 
      if (cff1.ge.0. .and. cff2.ge.0. ) then

        cff=1./(cff1+cff2)
        cff1=cff1*cff                ! Interpolate SST and dQdSST
        cff2=cff2*cff                ! data to model time.  Convert
                                     ! units of dQdSST from

        do j=jstrR,jendR             !
          do i=istrR,iendR
             ! convert units from CCSM to ROMS (nmol/cm2 to mmol/m2)
              nitrate(i,j)=(cff1*nitrateg(i,j,it1)+cff2*nitrateg(i,j,it2))
     &            * nmol_cm2_to_mmol_m2
          enddo
        enddo
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,2A/3(1x,A,F16.10)/)')
     &       'SET_NITRATE_TILE - current model time is outside bounds ',
     &       'of ''nitrate_time''.', 'NITRATE_TSTART=',nitrate_time(it1)*sec2day,
     &                           'TDAYS=',     tdays,
     &                           'NITRATE_TEND=',  nitrate_time(it2)*sec2day
        ierr=ierr+1
      endif
      return
      end
#else /* BIOLOGY_BEC */
      subroutine get_nitrate_empty
      end
#endif /* BIOLOGY_BEC */
 
