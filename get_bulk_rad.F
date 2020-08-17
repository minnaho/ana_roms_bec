#include "cppdefs.h"
#if defined BULK_FLUX

#ifdef MPI
#  define mpi_nonexit_warn if (mynode.eq.0)
#else
# define mpi_nonexit_warn
#endif 

      subroutine get_bulk_rad(ierr) 
!
!  Read in LW down and SW down heat fluxes
!
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncvars.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread,iblk
      integer ncid
      character*(max_name_size) blkname
!
! Initialization: Check, whether bulk netCDF file is already
!===============  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables. Determine whether there is cycling to
! reuse the input data and what is cycling period "sms_cycle", find
! initial cycling index "sms_ncycle" and record index "sms_rec".
! Set initial value for time index "itsms" and set both time record
! bounds to large negative artificial values to trigger the logic
! in reading part below.
!      
      iblk=max(rad_file_id,1)
 100  ncid=ncidbulk(iblk)
      blkname=blkfile(iblk)
       ierr = nf_noerr
       lstr=lenstr(blkname)
       
       if (iic.eq.ntstart) then
        if (ncid.eq.-1) then
          ierr=nf_open (blkname(1:lstr), nf_nowrite, ncid)
          if (ierr.eq.nf_noerr) then
            ncidbulk(iblk)=ncid
          else
            write(*,'(/1x,4A/)') '### ERROR in get_bulk: ',
     &      'can not open netCDF file ''', blkname(1:lstr), '''.'
          endif            
        endif
        if (ierr .eq. nf_noerr) then
          
          ierr=
     &    nf_inq_varid (ncid,'radsw', radsw_id)

          ierr=ierr + 
     &    nf_inq_varid (ncid,'radlw', radlw_id)        
          
          if (ierr .eq. nf_noerr) then
            rad_file_id=iblk
            mpi_nonexit_warn write(*,*) 'Read all fields in ',blkname(1:lstr)
            ierr=nf_inq_varid (ncid, 'rad_time', rad_tid)
              if (ierr .eq. nf_noerr) then
                call set_cycle (ncid, rad_tid, ntrad, rad_cycle,
     &                                   rad_ncycle, rad_rec, ierr)
                itrad=1
                rad_time(1)=-1.E+20
                rad_time(2)=-1.E+20
              else
                write(*,1) 'rad_time', blkname(1:lstr)
              endif
          else
            mpi_nonexit_warn write(*,1) 'the heat flux variables ',blkname(1:lstr)            
            if (iblk .lt. max_blk) then
              iblk=iblk+1
              goto 100
            endif            
          endif
        endif
      endif  

  1   format(/1x,'WARNING: get_bulk cannot find ''',
     &              A, ''''/7x, 'in netCDF file ''', A, '''.')

!
! Read data from the file:  Check if model time is bounded by past
!===== ==== ==== === =====  and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion
! of data. Repeat until model time falls between the two data times.
!
      do while (rad_time(itrad).lt.time+0.5*dt
     &  .and.ierr.eq.nf_noerr)
        call advance_cycle (rad_cycle, ntrad, rad_ncycle,
     &                                    rad_rec, ierr)
        if (ierr.eq.nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, rad_tid, rad_rec, cff)
          if (ierr .eq. nf_noerr) then
            itrad=min(3-itrad,ntrad)
            rad_time(itrad)=rad_cycle*rad_ncycle + cff*day2sec                       
            ierr=
     &      nf_fread (radswg(START_2D_ARRAY,itrad), ncid, 
     &      radsw_id,rad_rec, r2dvar)
            ierr=ierr + 
     &      nf_fread (radlwg(START_2D_ARRAY,itrad), ncid, 
     &      radlw_id,rad_rec,  r2dvar)        
              if (ierr .eq. nf_noerr) then
       mpi_nonexit_warn write(*,'(6x,A,3x,A,F12.4,2(1x,A,I6))')
     &          'get_bulk -- read bulk for all fields',
     &          'rad_time =', cff, 'rec =',rad_rec MYID
                 if (ntrad.eq.1) return 
              else
                write(*,2) ' one field ', rad_rec
              endif
          else
            write(*,2) 'rad_time', rad_rec
          endif
        else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &      'ERROR in get_bulk: requested time record ', rad_rec,
     &      'exeeds the last record', ntrad,  'available in netCDF ',
     &      'file ''',   blkname(1:lstr),   '''',  'tdays = ', tdays,
     &      '  but the last available  rad_time =',
     &                                       rad_time(itrad)*sec2day
        endif
      enddo
  2   format(/1x,'ERROR in get_bulk: cannot read variable ''',
     &                                 A, ''' for record ', I4/)

      return
      end



      subroutine set_bulk_rad_tile (Istr,Iend,Jstr,Jend,ierr)
!
! Set-up bulk data for current tile.
!
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "grid.h"
      integer Istr,Iend,Jstr,Jend, ierr,i,j, it1,it2
      real cff,cff1,cff2, cff3,cff4
      real cff5,cff6    

      real val1,val2,val3,val4,val5,val6,val7,val8
!
# include "compute_extended_bounds.h"
!

      it1=3-itrad
      it2=itrad
      cff=time+0.5*dt
      cff1=rad_time(it2)-cff
      cff2=cff-rad_time(it1)   
      
      srf_scale = 1.        

      if (cff1.lt.dt) synchro_flag=.true.
 
      if (ntrad.eq.1) then                ! Load time-invariant
        if (iic.eq.ntstart) then          ! wind stress components
          do j=jstrR,jendR
            do i=IstR,iendR
              radlw(i,j)= radlwg(i,j,itrad)              
              radsw(i,j)= radswg(i,j,itrad)
            enddo                                ! Interpolate wind
          enddo                                  ! stress in time,

        endif                                    
      elseif (cff1.ge.0. .and. cff2.ge.0.) then               
        cff=srf_scale/(cff1+cff2)
        cff3=cff1*cff
        cff4=cff2*cff      
        do j=jstrR,jendR                         ! flux.
          do i=istrR,iendR
            radlw(i,j)=cff3*radlwg(i,j,it1)+cff4*radlwg(i,j,it2)
            radsw(i,j)=cff3*radswg(i,j,it1)+cff4*radswg(i,j,it2)
          enddo
        enddo
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,2A/3(1x,A,F16.10)/)')
     &     'SET_BULK_RAD_TILE - current model time is outside bounds ',
     &     'of ''rad_time''.',  'BULK_TSTART=', rad_time(it1)*sec2day,
     &                          'TDAYS=',      tdays,
     &                          'BULK_TEND=',   rad_time(it2)*sec2day
        ierr=ierr+1
      endif       

      return
      end
#else
      subroutine get_bulk_rad_empty
      return
      end
#endif /* BULK_FLUX */
