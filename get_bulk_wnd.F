#include "cppdefs.h"
#if defined BULK_FLUX

      subroutine get_bulk_wnd(ierr) 
!
!  Read in Atmospheric wind
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
      iblk=max(wnd_file_id,1)
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
            write(*,'(/1x,4A,1x,2A,I4 /)') '### ERROR: get_smflux :: ',
     &            'Cannot open netCDF file ''', blkname(1:lstr), '''.',
     &             nf_strerror(ierr) MYID
          endif
        endif
        if (ierr .eq. nf_noerr) then
          
          ierr=nf_inq_varid (ncid, 'uwnd', uwnd_id)
          
          ierr=ierr + 
     &         nf_inq_varid (ncid, 'vwnd', vwnd_id)
          
          
          
          if (ierr .eq. nf_noerr) then

            wnd_file_id=iblk
            mpi_nonexit_warn write(*,*) 'Read all fields in ',blkname(1:lstr) 

!              ierr=nf_get_att_FTYPE (ncid, uwnd_id, 'level',zref)
!              if (ierr .gt. nf_noerr) then
!                 write(*,*) 'Error when reading the wnd level reference, 10m default'
                 zref=10
!              endif


            ierr=nf_inq_varid (ncid, 'wnd_time', wnd_tid)
              if (ierr .eq. nf_noerr) then
                
                call set_cycle (ncid, wnd_tid, ntwnd, wnd_cycle,
     &                                   wnd_ncycle, wnd_rec, ierr)
                itwnd=1
                wnd_time(1)=-1.E+20
                wnd_time(2)=-1.E+20
              else
                write(*,1) 'wind variables', blkname(1:lstr)
              endif
          elseif (iblk .lt. max_blk) then
            mpi_nonexit_warn write(*,1) 'uwnd ',blkname(1:lstr)   
            iblk=iblk+1
            goto 100
          else
                write(*,2) 'wnd_time', blkname(1:lstr)

          endif
        endif
      endif
  1   format(1x,   'WARNING: get_wnd    :: Cannot find variable ''',
     &                                A, ''' in file ''', A, '''.')
  2   format(1x, '###EMFMDS ERROR: get_wnd :: Cannot find variable ''',
     &                                A, ''' in file ''', A, '''.')

!
! Read data from the file:  Check if model time is bounded by past
!===== ==== ==== === =====  and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion
! of data. Repeat until model time falls between the two data times.
!
!      if (ntwnd.gt.1) then


        do while (wnd_time(itwnd).lt.time+0.5*dt  .and.
     &                                    ierr.eq.nf_noerr)
        call advance_cycle (wnd_cycle, ntwnd, wnd_ncycle,
     &                                    wnd_rec, ierr)
          if (ierr.eq.nf_noerr) then
            ierr=nf_get_var1_FTYPE (ncid, wnd_tid, wnd_rec, cff)
            if (ierr .eq. nf_noerr) then
              itwnd=min(3-itwnd,ntwnd)
              wnd_time(itwnd)=wnd_cycle*wnd_ncycle + cff*day2sec 
# ifdef WND_AT_RHO_POINTS                     
            ierr=
     &      nf_fread (uwndg(START_2D_ARRAY,itwnd), ncid, 
     &      uwnd_id,wnd_rec, r2dvar)
            ierr=ierr + 
     &      nf_fread (vwndg(START_2D_ARRAY,itwnd), ncid, 
     &      vwnd_id,wnd_rec,  r2dvar)
#else
            ierr=
     &      nf_fread (uwndg(START_2D_ARRAY,itwnd), ncid, 
     &      uwnd_id,wnd_rec, u2dvar)
            ierr=ierr + 
     &      nf_fread (vwndg(START_2D_ARRAY,itwnd), ncid, 
     &      vwnd_id,wnd_rec,  v2dvar)
#endif
              if (ierr .eq. nf_noerr) then
                 mpi_nonexit_warn write(*,'(6x,A,3x,A,F12.4,2(1x,A,I6))')
     &          'get_bulk -- read bulk for all fields',
     &          'wnd_time =', cff, 'rec =',rad_rec MYID

!#                 write(*,
!#     &        '(6x,A,3x,A,F12.4,2(1x,A,I4))')
!     &          'get_bulk -- read bulk for all fields',
!     &          'wnd_time =', cff, 'rec =',wnd_rec MYID
                 if (ntwnd.eq.1) return
              else
                write(*,2) ' one field ', wnd_rec
              endif

          else
            write(*,2) 'wnd_time', wnd_rec
          endif

          else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &     '### ERROR: get_wnd :: requested time record ',  wnd_rec,
     &     'exeeds the last record', ntwnd, 'available in forcing ',
     &     'netCDF file ''',   blkname(1:lstr),   '''',  'tdays = ',
     &      tdays,        '  but the last available  sss_time =',
     &                                      wnd_time(itwnd)*sec2day
          endif
        enddo
!      elseif (iic.eq.ntstart) then
!        wnd_rec=1
!              ierr=nf_fread (uwndg(START_2D_ARRAY,itwnd),ncid, 
!     &                                        uwnd_id,wnd_rec, 
!# ifdef WND_AT_RHO_POINTS       
!     &         r2dvar 
!# else
!     &         u2dvar 
!# endif   
!     &                     )
!        if (ierr .eq. nf_noerr) then
!              ierr=nf_fread (vwndg(START_2D_ARRAY,itwnd), ncid, 
!     &                                         vwnd_id,wnd_rec,  
!# ifdef WND_AT_RHO_POINTS       
!     &         r2dvar 
!# else
!     &         vdvar 
!# endif  
!     &                      )
!          if (ierr .eq. nf_noerr) then
!            mpi_nonexit_warn write(*,'(6x,2A,3x,1x,A,I4)')
!     &         'get_wnd    :: read time-invariant surface ',
!     &                             'momentum stress' MYID
!          else
!            write(*,3) 'vwnd', sms_rec, nf_strerror(ierr)
!          endif
!        else
!          write(*,3) 'uwnd',  sms_rec,  nf_strerror(ierr)
!        endif
!      endif
!  3   format(/1x,'### ERROR: get_wnd    :: Cannot read variable ''',
!     &                                A, ''' for rec =', I4/12x,A/)


      return
      end



      subroutine set_bulk_wnd_tile (Istr,Iend,Jstr,Jend,ierr)
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
      it1=3-itwnd
      it2=itwnd
      cff=time+0.5*dt
      cff1=wnd_time(it2)-cff
      cff2=cff-wnd_time(it1)   
      
      wnd_scale = 1.        

      if (cff1.lt.dt) synchro_flag=.true.
 
      if (ntwnd.eq.1) then                ! Load time-invariant
        if (iic.eq.ntstart) then          ! wind stress components
          do j=jstrR,jendR
            do i=IstR,iendR                             
              uwnd(i,j)= uwndg(i,j,itwnd)
              vwnd(i,j)= vwndg(i,j,itwnd)             
            enddo                                ! Interpolate wind
          enddo                                  ! stress in time,

        endif                                    
      elseif (cff1.ge.0. .and. cff2.ge.0.) then    
        cff=wnd_scale/(cff1+cff2)
        cff3=cff1*cff
        cff4=cff2*cff      
        do j=jstrR,jendR                         ! flux.
          do i=istrR,iendR 
            uwnd(i,j)=cff3*uwndg(i,j,it1)+cff4*uwndg(i,j,it2)
            vwnd(i,j)=cff3*vwndg(i,j,it1)+cff4*vwndg(i,j,it2)
          enddo
        enddo
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,2A/3(1x,A,F16.10)/)')
     &     'set_bulk_wnd_tile - current model time is outside bounds ',
     &     'of ''wnd_time''.',  'bulk_tstart=', wnd_time(it1)*sec2day,
     &                          'tdays=',      tdays,
     &                          'bulk_tend=',   wnd_time(it2)*sec2day
        ierr=ierr+1
      endif
      
      return
      end
#else
      subroutine get_bulk_wnd_empty
      return
      end
#endif /* BULK_FLUX */
