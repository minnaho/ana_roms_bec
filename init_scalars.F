#include "cppdefs.h"

       subroutine init_scalars (ierr)

! Set initial values for globally accessable (stored in common blocks)
! scalar variables.  Typically these are associated with time indices,
! global summation counters, model utility switches, i/o contlol, etc.
! Note that NONE of the model physical parameters is initialized here.

      implicit none
      integer ierr, i, itrc

#ifdef PASSIVE_TRACER
      integer itrc
      character*3 c_indx
#endif

#define ALL_DATA
#include "param.h"
#include "scalars.h"
#include "comm_vars.h"
#include "forces.h"
#include "sources.h"
#include "mixing.h"
#include "ncvars.h"
#include "ncvars_floats.h"
#include "boundary.h"
#ifdef DIAGNOSTICS_EDDY
# include "diags_eddy.h"
#endif

c      integer size, kmp_get_stacksize_s
c      size=kmp_get_stacksize_s()
c      write(*,*) 'kmp_stacksize =', size

      start_time=0.                  ! Set time and evolving time
      nfast=1                        ! stepping indices to their
      iic=0                          ! initial values needed for
      kstp=1                         ! proper logical initialization
      knew=1                         ! of model fields.
      ntstart=1
#ifdef SOLVE3D
      iif=1
      nstp=1
      nrhs=1
      nnew=1
#endif
#ifdef FLOATS
      nfp1=0                         ! set time step indices for
      nf=3                           ! 5th order AB4AM4 corrected
      nfm1=2                         ! scheme
      nfm2=1
      nfm3=0
#endif
                                     ! Reset intra-thread
      may_day_flag=0                 ! communication signal
      trd_count=0                    ! variables and counters.
      trd_count=0
      do i=1,4
        cpu_all(i)=0.
      enddo
      do i=1,16
        barr_count(i)=0
      enddo

      mcheck=0
      first_time=0
      tile_count=0                   ! Reset global sumation variables
      bc_count=0                     ! and counters needed for volume
      volume=0.                      ! and energy diagnostics, and for
c      avke=0.                        ! volume conservation constraint.
c      avpe=0.
c      prev_ke=0.
c      avke2b=0.
c      avke3bc=0.
c      v2d_max=0.
c      v3d_max=0.
c      v3bc_max=0.

      hmin=+1.E+20                   ! Set extremal values to large
      hmax=-1.E+20                   ! unrealistic values, positive or
      grdmin=+1.E+20                 ! negative, depending on context
      grdmax=-1.E+20
      Cg_min=+1.E+20
      Cg_max=-1.E+20
#ifdef SOLVE3D
      rx0=-1.E+20
      rx1=-1.E+20
#endif
      bc_crss=QuadZero
#ifdef OBC_VOLCONS
      bc_flux=QuadZero
      ubar_xs=QuadZero
#endif
#if defined BIOLOGY || defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
      do i=0,2*NT+1
        global_sum(i)=QuadZero
      enddo
#endif
                                     ! Set slipperness parameter to
      gamma2=1.                      ! "free-slip" along the straight
                                     ! outer boundaries: this may be
                                     ! over written by "read_inp".
      ncidrst=-1
      nrecrst=0                      ! Set netCDF file IDs to negative
      ncidhis=-1                     ! values to signal their initially
      nrechis=0                      ! closed status.
#ifdef AVERAGES
      ncidavg=-1                     ! Also reset record counters for
      nrecavg=0                      ! restart, history, average, and
# ifdef OUTPUT_AVG_PHYS_BIO
      ncid_bio_avg = -1
# endif 
#endif
#ifdef STATIONS
      ncidstn=-1                     ! stations netCDF output files.
      nrecstn=0
#endif
#ifdef DIAGNOSTICS_EDDY
      nciddiags_eddy=-1
      nrecdiags_eddy=0
# ifdef AVERAGES
      nciddiags_eddy_avg=-1
      nrecdiags_eddy_avg=0
# endif
# endif
#ifdef FLOATS
      ncidflt=-1
      nrecflt=0
#endif
#if defined PSOURCE || defined Q_PSOURCE && ! defined ANA_PSOURCE
      ncidfpsrc=-1
#endif
# ifdef BULK_FLUX     
      do i=1,max_blk_file
        ncidbulk(i)=-1
      enddo
# else      
      do i=1,max_frc_file
        ncidfrc(i)=-1
      enddo
# endif

#if !defined ANA_SMFLUX && !defined BULK_FLUX
      sms_file_id=0
#endif

#ifdef BULK_FLUX
      rad_file_id=0
      tra_file_id=0
      prec_file_id=0
      wnd_file_id=0
# if defined SALINITY_RESTORING
      sss_file_id=0
# endif
# ifdef BULK_DIAGS
      ncid_bulk_diags_his = -1
#ifdef AVERAGES
      ncid_bulk_diags_avg = -1
#endif
! Reset record counters for phys flux analysis averages/history output
!
      nrec_bulk_diags_his=0
#ifdef AVERAGES
      nrec_bulk_diags_avg=0
#endif

# endif

#endif




#ifdef SOLVE3D
# ifndef BULK_FLUX
#  if !defined ANA_STFLUX || !defined ANA_SSFLUX
      do i=1,NT
        stf_file_id(i)=0
      enddo
#  endif
#  ifndef ANA_SRFLUX
      srf_file_id=0
#  endif
# endif

# if defined QCORRECTION && !defined ANA_SST
      sst_file_id=0
      sss_file_id=0
      dqdt_file_id=0
# endif
# ifdef BIOLOGY_BEC
      dust_file_id=0
# endif
#endif
      ncidclm=-1
#ifndef ANA_BRY
       bry_id=-1
#endif

      call get_date (date_str)

! Names of variables in NetCDF output files. The first element
! is the name of the variable; the other two are are attributes.

c--#define CR
CR      write(*,'(1x,A,I3)') 'indxTime =',indxTime
      vname(1,indxTime)='ocean_time'
      vname(2,indxTime)='Time since initialization'
      vname(3,indxTime)='second'

CR      write(*,'(1x,A,I3)') 'indxZ =',indxZ
      vname(1,indxZ)='zeta'
      vname(2,indxZ)='free-surface elevation'
      vname(3,indxZ)='meter'

CR      write(*,'(1x,A,I3)') 'indxUb =',indxUb
      vname(1,indxUb)='ubar'
      vname(2,indxUb)='vertically averaged u-momentum component'
      vname(3,indxUb)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxVb = ',indxVb
      vname(1,indxVb)='vbar'
      vname(2,indxVb)='vertically averaged v-momentum component'
      vname(3,indxVb)='meter second-1'

#ifdef SOLVE3D
CR      write(*,'(1x,A,I3)') 'indxU = ', indxU
      vname(1,indxU)='u'
      vname(2,indxU)='u-momentum component'
      vname(3,indxU)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxV = ', indxV
      vname(1,indxV)='v'
      vname(2,indxV)='v-momentum component'
      vname(3,indxV)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxO =',indxO
      vname(1,indxO)='omega'
      vname(2,indxO)='S-coordinate vertical velocity'
      vname(3,indxO)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxW =', indxW
      vname(1,indxW)='w'
      vname(2,indxW)='vertical velocity'
      vname(3,indxW)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxR =', indxR
      vname(1,indxR)='rho'
      vname(2,indxR)='density anomaly'
      vname(3,indxR)='kilogram meter-3'

CR      write(*,'(1x,A,I3)') 'indxT = ', indxT
      vname(1,indxT)='temp'
      vname(2,indxT)='potential temperature'
      vname(3,indxT)='Celsius'

# ifdef SALINITY
CR      write(*,'(1x,A,I3)') 'indxS = ', indxS
      vname(1,indxS)='salt'
      vname(2,indxS)='salinity'
      vname(3,indxS)='PSU'
# endif

# ifdef PASSIVE_TRACER
      do itrc = 1, ntrc_pas
         itpas(itrc) = itemp + ntrc_salt + itrc
      end do
#  ifdef TRACER_BE7
! assuming that this is the only passive tracer
      itrc = indxT+ntrc_salt+1
      vname(1,itrc) = 'Be7'
      vname(2,itrc) = 'Beryllium-7'
      vname(3,itrc)='no unit'
      vname(4,itrc)='passive tracer Be-7, scalar, series'

#  else
! set up generic names      
      do itrc = indxT+ntrc_salt+1,indxT+ntrc_salt+ntrc_pas
        if (itrc-indxT-ntrc_salt .lt. 10) then
           write (c_indx, '(i1)') itrc-indxT-ntrc_salt
        else if (itrc-indxT-ntrc_salt .lt. 100) then
           write (c_indx, '(i2)') itrc-indxT-ntrc_salt
        else
           write (c_indx, '(i3)') itrc-indxT-ntrc_salt
        end if
        vname(1,itrc) = 'tpas'/ /trim(c_indx)
        vname(2,itrc) = 'passive tracer '/ /c_indx
        vname(3,itrc)='no unit'
      end do
#  endif /* TRACER_BE7 */
# endif /* PASSIVE_TRACER */

# ifdef BIOLOGY
CR      write(*,'(1x,A,I3)') 'indxNO3 =', indxNO3
      vname(1,indxNO3)='NO3'
      vname(2,indxNO3)='NO3 Nutrient'
      vname(3,indxNO3)='mMol N m-3'

CR      write(*,'(1x,A,I3)') 'indxNH4 =', indxNH4
      vname(1,indxNH4='NH4'
      vname(2,indxNH4)='NH4 Nutrient'
      vname(3,indxNH4)='mMol N m-3'

CR      write(*,'(1x,A,I3)') 'indxDet =', indxDet
      vname(1,indxDet)='DET'
      vname(2,indxDet)='Detritus'
      vname(3,indxDet)='mMol N m-3'

CR      write(*,'(1x,A,I3)') 'indxPhyt =', indxPhyt
      vname(1,indxPhyt)='PHYTO'
      vname(2,indxPhyt)='Phytoplankton'
      vname(3,indxPhyt)='mMol N m-3'

CR      write(*,'(1x,A,I3)') 'indxZoo =', indxZoo
      vname(1,indxZoo)='ZOO'
      vname(2,indxZoo)='Zooplankton'
      vname(3,indxZoo)='mMol N m-3'
# endif

CR      write(*,'(1x,A,I3)') 'indxAkv =', indxAkv
      vname(1,indxAkv)='AKv'
      vname(2,indxAkv)='vertical viscosity coefficient'
      vname(3,indxAkv)='meter2 second-1'

CR      write(*,'(1x,A,I3)') 'indxAkt =',indxAkt
      vname(1,indxAkt)='AKt'
      vname(2,indxAkt)='vertical thermal conductivity coefficient'
      vname(3,indxAkt)='meter2 second-1'
# ifdef SALINITY
CR      write(*,'(1x,A,I3)') 'indxAkt =',indxAks
      vname(1,indxAks)='AKs'
      vname(2,indxAks)='salinity vertical diffusion coefficient'
      vname(3,indxAks)='meter2 second-1'
# endif
# ifdef LMD_KPP
CR      write(*,'(1x,A,I3)') 'indxHbls =',indxHbls
      vname(1,indxHbls)='hbls'
      vname(2,indxHbls)='Thickness of KPP surface boundary layer'
      vname(3,indxHbls)='meter'
# endif
# ifdef LMD_BKPP
CR      write(*,'(1x,A,I3)') 'indxHbbl =',indxHbbl
      vname(1,indxHbbl)='hbbl'
      vname(2,indxHbbl)='Thickness of KPP bottom boundary layer'
      vname(3,indxHbbl)='meter'
# endif
#endif

#ifdef SOLVE3D
# ifdef SG_BBL96
#  ifndef ANA_WWAVE
      vname(1,indxWWA)='Awave'
      vname(2,indxWWA)='wind induced wave amplitude'
      vname(3,indxWWA)='meter'

      vname(1,indxWWD)='Dwave'
      vname(2,indxWWD)='wind induced wave direction'
      vname(3,indxWWD)='degrees'

      vname(1,indxWWP)='Pwave'
      vname(2,indxWWP)='wind induced wave Period'
      vname(3,indxWWP)='second'
#  endif
#  ifndef ANA_BSEDIM
      vname(1,indxBSS)='Ssize'
      vname(2,indxBSS)='bottom sediment grain diameter size'
      vname(3,indxBSS)='meter'

      vname(1,indxBSD)='Sdens'
      vname(2,indxBSD)='bottom sediment grain density'
      vname(3,indxBSD)='kilogram meter-3'
#  endif
# endif
#endif
#ifdef ICE
      vname(1,indxAi)='aice'
      vname(2,indxAi)='fraction of cell covered by ice'
      vname(3,indxAi)='nondimensional'

      vname(1,indxUi)='uice'
      vname(2,indxUi)='u-component of ice velocity'
      vname(3,indxUi)='meter sec-1'

      vname(1,indxVi)='vice'
      vname(2,indxVi)='v-component of ice velocity'
      vname(3,indxVi)='meter sec-1'

      vname(1,indxHi)='hice'
      vname(2,indxHi)='depth of ice cover'
      vname(3,indxHi)='meter'

      vname(1,indxHS)='snow_thick'
      vname(2,indxHS)='depth of snow cover'
      vname(3,indxHS)='meter'

      vname(1,indxTIsrf)='tisrf'
      vname(2,indxTIsrf)='temperature of ice surface'
      vname(3,indxTIsrf)='Kelvin'
#endif

      ierr=0
# ifdef BIOLOGY_BEC
      call init_scalars_bec(ierr)
      if (ierr.ne.0) goto 100
# endif
# ifdef PHYS_FLUX_ANALYSIS
      call init_scalars_physflux(ierr)
      if (ierr.ne.0) goto 100 
# endif /* PHYS_FLUX_ANALYSIS */



# ifdef DIAGNOSTICS_EDDY
      itrc=1

CR       write(stdout,'(1x,A,I3)') 'indxeddyuu = ',indxeddyuu+itrc-1
       vname(1,indxeddyuu+itrc-1)='uu'
CR       write(stdout,'(20x,A)') vname(1,indxeddyxadv+itrc-1)
       vname(2,indxeddyuu+itrc-1)='uu'
       vname(3,indxeddyuu+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddyvv = ',indxeddyvv+itrc-1
       vname(1,indxeddyvv+itrc-1)='vv'
CR       write(stdout,'(20x,A)') vname(1,indxeddyvv+itrc-1)
       vname(2,indxeddyvv+itrc-1)='vv'
       vname(3,indxeddyvv+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddyuv = ',indxeddyuv+itrc-1
       vname(1,indxeddyuv+itrc-1)='uv'
CR       write(stdout,'(20x,A)') vname(1,indxeddyuv+itrc-1)
       vname(2,indxeddyuv+itrc-1)='uv'
       vname(3,indxeddyuv+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddyub = ',indxeddyub+itrc-1
       vname(1,indxeddyub+itrc-1)='ub'
CR       write(stdout,'(20x,A)') vname(1,indxeddycor+itrc-1)
       vname(2,indxeddyub+itrc-1)='ub'
       vname(3,indxeddyub+itrc-1)='meter2 second-3'

CR       write(stdout,'(1x,A,I3)') 'indxeddyvb = ',indxeddyvb+itrc-1
       vname(1,indxeddyvb+itrc-1)='vb'
CR       write(stdout,'(20x,A)') vname(1,indxeddyvb+itrc-1)
       vname(2,indxeddyvb+itrc-1)='vb'
       vname(3,indxeddyvb+itrc-1)='meter2 second-3'

CR       write(stdout,'(1x,A,I3)') 'indxeddywb = ',indxeddywb+itrc-1
       vname(1,indxeddywb+itrc-1)='wb'
CR       write(stdout,'(20x,A)') vname(1,indxeddywb+itrc-1)
       vname(2,indxeddywb+itrc-1)='wb'
       vname(3,indxeddywb+itrc-1)='meter2 second-3'

CR       write(stdout,'(1x,A,I3)') 'indxeddynotused1 =
C',indxeddynotused1+itrc-1
        vname(1,indxeddynotused1+itrc-1)='uw'
CR       write(stdout,'(20x,A)') vname(1,indxeddynotused1+itrc-1)
       vname(2,indxeddynotused1+itrc-1)='uw'
       vname(3,indxeddynotused1+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddynotused2 =
C',indxeddynotused2+itrc-1
       vname(1,indxeddynotused2+itrc-1)='vw'
CR       write(stdout,'(20x,A)') vname(1,indxeddynotused2+itrc-1)
       vname(2,indxeddynotused2+itrc-1)='vw'
       vname(3,indxeddynotused2+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddynotused3 =
C',indxeddynotused3+itrc-1
       vname(1,indxeddynotused3+itrc-1)='b'
CR       write(stdout,'(20x,A)') vname(1,indxeddynotused3+itrc-1)
       vname(2,indxeddynotused3+itrc-1)='b'
       vname(3,indxeddynotused3+itrc-1)='None'


# endif /* DIAGNOSTICS_EDDY */


      call read_inp (ierr)
      if (ierr.ne.0) goto 100 
# ifdef BIOLOGY_BEC
      call init_scalars_bec_part2(ierr)
      if (ierr.ne.0) goto 100
# endif
#ifdef SOLVE3D
      call set_weights
#endif
#ifdef GRID_LEVEL
# if GRID_LEVEL < MAX_GRID_LEVEL
#  if GRID_LEVEL == 1
        call init_scalars_2 (ierr)
#  elif GRID_LEVEL == 2
        call init_scalars_3 (ierr)
#  elif GRID_LEVEL == 3
        call init_scalars_4 (ierr)
#  endif
# endif
#endif
      return
 100  stop
      end
