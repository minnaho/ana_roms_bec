#include "cppdefs.h"
#ifdef BIOLOGY_BEC

      subroutine get_iron(ierr)
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncvars.h"
      real cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, ifrc, ncid
      character*(max_name_size) frcname
      character(len=4), parameter :: iron_var='iron'

!
! Initialization: Check, whether forcing netCDF file is already
!===============  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables, and determine whether sst data exists
! as 2D field or scalar. After that determine whether there is
! cycling to reuse the input data and what is cycling period
! "iron_cycle", set initial cycling index "iron_ncycle", record index
! "iron_rec", and time index "itiron" and both time record bounds to
! large negative artificial values, so that it will trigger logic in
! reading part below.
!
# if defined BULK_FLUX
      ifrc=max(iron_file_id,1)
 100  ncid=ncidbulk(ifrc)
      frcname=blkfile(ifrc)
# else
      ifrc=max(iron_file_id,1)
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
     &                   'in get_iron: can not open netCDF file ''',
     &                                       frcname(1:lstr), '''.'
          endif
# else
          ierr=nf_open(frcname(1:lstr), nf_nowrite, ncid)
          if (ierr.eq.nf_noerr) then
            ncidfrc(ifrc)=ncid
          else
            write(*,'(/1x,4A/)') '### ERROR ',
     &                   'in get_iron: can not open netCDF file ''',
     &                                       frcname(1:lstr), '''.'
          endif
# endif
        endif
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varid (ncid, 'iron_time', iron_tid)
          if (ierr .eq. nf_noerr) then
            iron_file_id = ifrc
            lvar=lenstr(iron_var)
            ierr=nf_inq_varid (ncid, iron_var(1:lvar),
     &                                                   iron_id)
            if (ierr .eq. nf_noerr) then
              call set_cycle (ncid, iron_tid, ntiron, iron_cycle,
     &                                 iron_ncycle, iron_rec, ierr)
              itiron=1
              iron_time(1)=-1.E+20
              iron_time(2)=-1.E+20
            else
              write(*,1) iron_var(1:lvar),
     &                                  frcname(1:lstr)
            endif
          else
            mpi_nonexit_warn write(*,1) 'iron_time', frcname(1:lstr)
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
  1   format(/1x,'WARNING: get_iron cannot find variable ''',
     &     A, ''''/7x, 'in forcing netCDF file ''', A, '''.')
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
      do while (iron_time(itiron).lt.time+0.5*dt .and. ierr.eq.nf_noerr)
        call advance_cycle (iron_cycle, ntiron, iron_ncycle,
     &                                 iron_rec,     ierr)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, iron_tid, iron_rec, cff)
          if (ierr .eq. nf_noerr) then
            itiron=3-itiron
            iron_time(itiron)=cff*day2sec + iron_cycle*iron_ncycle
            ierr=nf_fread (irong(START_2D_ARRAY,itiron), ncid,
     &                                 iron_id, iron_rec, r2dvar)
            if (ierr .eq. nf_noerr) then
              mpi_master_only write(*,
     &             '(6x,2A,G12.4,1x,A,I4)') 'get_iron    --',
     &                 ' read iron                       for time =',
     &                   cff MYID
            else
              write(*,2) 'iron', iron_rec
            endif
          else
            write(*,2) 'iron_time', iron_rec
          endif
        else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &   'ERROR in get_iron: requested time record ',      iron_rec,
     &   'exeeds the last record', ntiron, 'available in forcing ',
     &   'netCDF file ''',   frcname(1:lstr),   '''',  'tdays = ',
     &    tdays,        '  but the last available  iron_time =',
     &                                    iron_time(itiron)*sec2day
        endif
      enddo
  2   format(/1x,'ERROR in get_iron: cannot read variable ''', A,
     &                                    ''' for record ', I4/)
      return
      end
 
      subroutine set_iron_tile (Istr,Iend,Jstr,Jend,ierr)
!
! Set-up sea surface iron data for current tile.
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
 
      it1=3-itiron
      it2=itiron
      cff1=iron_time(it2)-time
      cff2=time-iron_time(it1)


      if (ZEROTH_TILE  .and. ( cff1.lt.dt )) synchro_flag=.true.
 
      if (cff1.ge.0. .and. cff2.ge.0. ) then

        cff=1./(cff1+cff2)
        cff1=cff1*cff                ! Interpolate SST and dQdSST
        cff2=cff2*cff                ! data to model time.  Convert
                                     ! units of dQdSST from

        do j=jstrR,jendR             !
          do i=istrR,iendR
             ! convert units from CCSM to ROMS (nmol/cm2 to mmol/m2)
              iron(i,j)=(cff1*irong(i,j,it1)+cff2*irong(i,j,it2))
     &            * nmol_cm2_to_mmol_m2
          enddo
        enddo
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,2A/3(1x,A,F16.10)/)')
     &       'SET_IRON_TILE - current model time is outside bounds ',
     &       'of ''iron_time''.', 'IRON_TSTART=',iron_time(it1)*sec2day,
     &                           'TDAYS=',     tdays,
     &                           'IRON_TEND=',  iron_time(it2)*sec2day
        ierr=ierr+1
      endif
      return
      end
#else /* BIOLOGY_BEC */
      subroutine get_iron_empty
      end
#endif /* BIOLOGY_BEC */
 
