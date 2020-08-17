#include "cppdefs.h"

! ANALYTICAL PACKAGE:
!=========== ========
!
! This package is used to provide various analytical fields to the
! model when appropriate.
!
! Routines:
!
! ana_bsedim_tile   Analytical bottom sediment grain size
!                          and density.
! ana_meanRHO_tile  Analytical mean density anomaly.
! ana_smflux_tile   Analytical kinematic surface momentum flux
!                          (wind stress).
! ana_srflux_tile   Analytical kinematic surface shortwave
!                          radiation.
! ana_ssh_tile      Analytical sea surface height climatology.
! ana_sst_tile      Analytical sea surface temperature and dQdSST
!                         which are used during heat flux correction.
! ana_stflux_tile   Analytical kinematic surface flux of tracer type
!                          variables.
! ana_tclima_tile   Analytical tracer climatology fields.
! ana_uclima_tile   Analytical tracer climatology fields.
! ana_wwave_tile    Analytical wind induced wave amplitude,
!                         direction and period.
! modification: 25/04/2018
! PSOURCE AND QPSOURCE for rivers and anthropogenic pollution (AP)
!-------------------------------------------------------------------

#if defined ANA_BSEDIM && defined SG_BBL96
      subroutine ana_bsedim_tile (istr,iend,jstr,jend)

! Set initial bottom sediment grain diameter size [m]  and density
! used in the bottom boundary formulation [kg/m^3].

      implicit none
      integer istr,iend,jstr,jend, i,j
# include "param.h"
# include "bblm.h"
# include "grid.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          Ssize(i,j)=???
          Sdens(i,j)=???
        enddo
      enddo
      return
      end
#endif /* ANA_BSEDIM && SG_BBL96 */



#if defined SOLVE3D && defined ANA_MEANRHO
      subroutine ana_meanRHO_tile (istr,iend,jstr,jend)
!
!--------------------------------------------------------------------
!  This subroutine sets mean density anomaly [kg/m^3] using an
!  analytical expression.
!--------------------------------------------------------------------
!
      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff1, cff2 , cff3, cff4
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

!  Set mean density anomaly (kg/m^3) at horizontal and vertical
!  RHO-points.

# if defined BASIN
c     cff1=(44.690/39.382)**2
c     cff2=cff1*(rho0*800./g)*(5.0e-5/((42.689/44.690)**2))
c-song:
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c           rhobar(i,j,k)=R0-cff2*exp(z_r(i,j,k)/800.)*
            rhobar(i,j,k)=R0-
     &  (44.690/39.382)**2
     &  *(rho0*800./g)*(5.0e-5/((42.689/44.690)**2))
     &                           *exp(z_r(i,j,k)/800.)*
     &                    (0.6-0.4*tanh(z_r(i,j,k)/800.))
           enddo
        enddo
      enddo
# elif defined CANYON_A || defined GRAV_ADJ || defined OVERFLOW
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=R0
          enddo
        enddo
      enddo
# elif defined CANYON_B
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=R0-3.488*exp(z_w(i,j,k)/800.)*
     &                  (1.-(2./3.)*tanh(z_r(i,j,k)/800.))
          enddo
        enddo
      enddo
# elif defined NJ_BIGHT
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (z_r(i,j,k).ge.-15.0 ) then
              cff1=2.049264257728403e+01-z_r(i,j,k)*(
     &               2.640850848793918e-01+z_r(i,j,k)*(
     &               2.751125328535212e-01+z_r(i,j,k)*(
     &               9.207489761648872e-02+z_r(i,j,k)*(
     &               1.449075725742839e-02+z_r(i,j,k)*(
     &               1.078215685912076e-03+z_r(i,j,k)*(
     &               3.240318053903974e-05+
     &                1.262826857690271e-07*z_r(i,j,k)
     &                                           ))))))
              cff2=3.066489149193135e+01-z_r(i,j,k)*(
     &               1.476725262946735e-01+z_r(i,j,k)*(
     &               1.126455760313399e-01+z_r(i,j,k)*(
     &               3.900923281871022e-02+z_r(i,j,k)*(
     &               6.939014937447098e-03+z_r(i,j,k)*(
     &               6.604436696792939e-04+z_r(i,j,k)*(
     &               3.191792361954220e-05+
     &               6.177352634409320e-07*z_r(i,j,k)
     &                                           ))))))
            else
               cff1=14.6+6.7 *tanh(1.1*z_r(i,j,k)+15.9)
               cff2=31.3-0.55*tanh(1.1*z_r(i,j,k)+15.9)
            endif
            rhobar(i,j,k)=R0+Tcoef*cff1+Scoef*cff2
          enddo
        enddo
      enddo
# elif defined SEAMOUNT
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=28.-2.*exp(z_r(i,j,k)/1000.)
          enddo
        enddo
      enddo
# elif defined SHELFRONT
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=R0
          enddo
        enddo
      enddo
# elif defined UPWELLING
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=26.459472-2.24*exp(z_r(i,j,k)/50.)
          enddo
        enddo
      enddo
# elif defined JET
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
           cff1=TL1+
     &          TL2*z_r(i,j,k)+
     &          TL3*exp(z_r(i,j,k)/ZL1)
           cff2=TC1+
     &          TC2*z_r(i,j,k)+
     &          TC3*exp((z_r(i,j,k)-ZC1)/ZC2)
           cff3=YC1+YC2*
     &          tanh(YC3*xr(i,j)+YC4)
!           cff4=0.5*(
!     &          cff2*(1+tanh((cff3-yr(i,j))
!     &           /LCUR))+
!     &          cff1*(1+tanh((yr(i,j)-cff3)
!     &           /LCUR)))
           cff4=cff1
            rhobar(i,j,k)=R0+Tcoef*cff4
          enddo
        enddo
      enddo
# else
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=???
          enddo
        enddo
      enddo
# endif
# ifdef MASKING
!
!  Apply Land/Sea mask.
!
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=rhobar(i,j,k)*rmask(i,j)
          enddo
        enddo
      enddo
# endif
      return
      end
#endif /* ANA_MEANRHO && SOLVE3D */


#ifdef ANA_SMFLUX
      subroutine ana_smflux_tile (istr,iend,jstr,jend)

! Sets kinematic surface momentum flux (wind stress) "sustr"
! and "svstr"  [m^2/s^2] using an analytical expression.

      implicit none
      integer istr,iend,jstr,jend, i,j
# include "param.h"
# include "grid.h"
# include "forces.h"
# include "scalars.h"
      real Ewind, Nwind, dircoef, windamp
      data Ewind, Nwind, dircoef /0., 0., 0./
!        real wdir_mvco, cff, windday, wcos, wsin %% FAYCAL (In YU)

# include "compute_extended_bounds.h"

                          ! Set kinematic surface momentum flux
      do j=jstrR,jendR    ! (wind stress) components in XI- and ETA-
        do i=istrR,iendR  ! directions [m^2/s^2] at horizontal U-
          sustr(i,j)=0.   ! and V-points respectively.  For the sake
          svstr(i,j)=0.   ! of code compactness, both components are
        enddo             ! set to zero at first (default values)
      enddo               ! and then overwritten, when appropriate.

# ifdef BASIN
      cff1=0.0001 * 0.5*(1.+tanh((time-6.*86400.)/(3.*86400.)))
      cff2=2.*pi/el
      do j=jstrR,jendR
        do i=istr,iendR
          sustr(i,j)=-cff1*cos(cff2*yr(i,j))
        enddo
      enddo
# elif defined CANYON_A || defined CANYON_B
      do j=jstrR,jendR
        do i=istr,iendR
          sustr(i,j)=0.0001*0.5*sin(2.*pi*tdays/10.)*
     &               (1.-tanh((yr(i,j)-0.5*el)/10000.))
        enddo
      enddo
# elif defined EKMAN_SPIRAL
      Ewind=0.1/rho0
      do j=jstrR,jendR
        do i=istr,iendR
          sustr(i,j)=Ewind
        enddo
      enddo

# elif defined NJ_BIGHT
      windamp=0.086824313
      dircoef=0.5714286
      if (tdays.le.0.5) then
        Ewind=windamp*dircoef*sin(pi*tdays)/rho0
        Nwind=windamp*sin(pi*tdays)/rho0
      else
        Ewind=windamp*dircoef/rho0
        Nwind=windamp/rho0
      endif
      do j=jstrR,jendR
        do i=istr,iendR
          cff1=0.5*(angler(i-1,j)+angler(i,j))
          sustr(i,j)=Ewind*cos(cff1)+Nwind*sin(cff1)
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          cff1=0.5*(angler(i,j)+angler(i,j-1))
          svstr(i,j)=-Ewind*sin(cff1)+Nwind*cos(cff1)
        enddo
      enddo

# elif defined TASMAN_SEA
      windamp=0.5*(1.+tanh((time-6.*86400.)/(3.*86400.)))
      do j=jstrR,jendR
        do i=istr,iendR
          cff1=yr(i,j)/el
          sustr(i,j)=windamp*cff1*(1.-cff1)*( 1.7580e-3 +cff1*(
     &                       -2.7907e-3+cff1*( -4.0764e-3+cff1*(
     &                                  9.1903e-3-cff1*4.7915e-3
     &                                                      ))))
        enddo
      enddo

# elif defined UPWELLING
      if (tdays.le.2.) then
        windamp=-0.1*sin(pi*tdays/4.)/rho0
      else
        windamp=-0.1/rho0
      endif
c??
c??
c      windamp=0.
c??
c??
      do j=jstrR,jendR
        do i=istrR,iendR
          sustr(i,j)=windamp
          svstr(i,j)=0.2*windamp !<-- to make it ETALON_CHECK sensitive
        enddo
      enddo
# elif defined BALTIC

      do j=jstrR,jendR
        do i=istr,iendR
          sustr(i,j)=0.1/rho0  ! Wind blowing West --> East
        enddo
      enddo


c# elif defined USWEST
# elif defined USSWC1
      do j=jstrR,jendR
        do i=istr,iendR
c         cff1=(latr(i,j)-latr(Lm/2,Mm/2))/20.
c         sustr(i,j)=1.d-4.*cff1
c         sustr(i,j)=-1.d-4
          sustr(i,j)=0.
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
c          svstr(i,j)=-1.0d-4
          svstr(i,j)=0.
        enddo
      enddo

# elif defined FLAT_BOTTOM
      do j=jstrR,jendR
        do i=istr,iendR
c         cff1=(latr(i,j)-latr(Lm/2,Mm/2))/20.
c         sustr(i,j)=1.d-4.*cff1
c         sustr(i,j)=-1.d-4
          sustr(i,j)=0.
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
c          svstr(i,j)=-1.0d-4
          svstr(i,j)=0.
        enddo
      enddo
# endif
      return
      end
#endif /* ANA_SMFLUX */

#if defined ANA_SRFLUX && defined SOLVE3D

      subroutine ana_srflux_tile (istr,iend,jstr,jend)

! Set kinematic surface solar shortwave radiation flux "srflx"
! [degC m/s] using an analytical expression.

      implicit none
      integer istr,iend,jstr,jend, i,j
# include "param.h"
# include "grid.h"
# include "forces.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          srflx(i,j)=0.
        enddo
      enddo
      return
      end
#endif /* ANA_SRFLUX */


#if defined ANA_SSH
      subroutine ana_ssh (tile)           ! Set analytical sea-surface
      implicit none                       ! height climatology [m]
      integer tile
#include "param.h"
#include "compute_tile_bounds.h"
      call ana_ssh_tile   (istr,iend,jstr,jend)
      return
      end

      subroutine ana_ssh_tile (istr,iend,jstr,jend)
      implicit none
      integer istr,iend,jstr,jend, i,j
#ifdef JET
      real TL0, TC0, cff1, cff2, cff3, cff4
#endif
# include "param.h"
# include "grid.h"
# include "climat.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

# if defined JET
      TL0=-TL1 - 0.5*TL2*ZREF**2 - ZL1*TL3*exp(ZREF/ZL1)
      TC0=-TC1 - 0.5*TC2*ZREF**2 - ZC2*TC3*exp((ZREF-ZC1)/ZC2)
      do j=jstrR,jendR
        do i=istrR,iendR
         cff1=YC1+YC2*tanh(YC3*xr(i,j)+YC4)
         cff2=TL0+ZL1*TL3
         cff3=TC0+ZC2*TC3*exp(-ZC1/ZC2)
         cff4=-0.5*g*Tcoef*(cff2+cff3)
         ssh(i,j)=cff4/(rho0*g)
                +Tcoef*cff2*(1+tanh((cff1-yr(i,j))/LCUR))/(2*rho0)
                +Tcoef*cff3*(1+tanh((yr(i,j)-cff1)/LCUR))/(2*rho0)
        enddo
      enddo
# else
      do j=jstrR,jendR
        do i=istrR,iendR
          ssh(i,j)=???
        enddo
      enddo
#endif /* JET */
      return
      end
#endif
#ifdef SOLVE3D
# if defined ANA_STFLUX || defined ANA_SSFLUX
      subroutine ana_stflux_tile (istr,iend,jstr,jend, itrc)

! Set kinematic surface flux of tracer type variables "stflx" (tracer
! units m/s) using analytical expressions.

      implicit none
      integer istr,iend,jstr,jend, itrc, i,j
#  include "param.h"
#  include "grid.h"
#  include "forces.h"
#  include "scalars.h"

# include "compute_extended_bounds.h"

      if (itrc.eq.itemp) then

! Set kinematic surface heat flux [degC m/s] at horizontal
! RHO-points.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,itemp)=0.
          enddo
        enddo

#  ifdef SALINITY
      elseif (itrc.eq.isalt) then

! Set kinematic surface freshwater flux (m/s) at horizontal
! RHO-points, scaling by surface salinity is done in STEP3D.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,isalt)=0.
          enddo
        enddo
#  endif
      else

!  Set kinematic surface flux of additional tracers, if any.

      endif
      return
      end
# endif /* ANA_STFLUX || ANA_SSFLUX */



# if defined ANA_TCLIMA && defined TCLIMATOLOGY
      subroutine ana_tclima (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call ana_tclima_tile   (istr,iend,jstr,jend)
      return
      end

      subroutine ana_tclima_tile (istr,iend,jstr,jend)
!
! Set analytical tracer climatology fields.
!
      implicit none
      integer istr,iend,jstr,jend, i,j,k
#  ifdef JET
      real cff1, cff2, cff3
#  endif
#  include "param.h"
#  include "climat.h"
#  include "grid.h"
#  include "ocean3d.h"
#  include "scalars.h"

# include "compute_extended_bounds.h"

#  if defined NJ_BIGHT
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (z_r(i,j,k).ge.-15.0) then
              tclm(i,j,k,itemp)=2.049264257728403e+01-z_r(i,j,k)*(
     &                          2.640850848793918e-01+z_r(i,j,k)*(
     &                          2.751125328535212e-01+z_r(i,j,k)*(
     &                          9.207489761648872e-02+z_r(i,j,k)*(
     &                          1.449075725742839e-02+z_r(i,j,k)*(
     &                          1.078215685912076e-03+z_r(i,j,k)*(
     &                          3.240318053903974e-05+
     &                          1.262826857690271e-07*z_r(i,j,k)
     &                                                      ))))))
              tclm(i,j,k,isalt)=3.066489149193135e+01-z_r(i,j,k)*(
     &                          1.476725262946735e-01+z_r(i,j,k)*(
     &                          1.126455760313399e-01+z_r(i,j,k)*(
     &                          3.900923281871022e-02+z_r(i,j,k)*(
     &                          6.939014937447098e-03+z_r(i,j,k)*(
     &                          6.604436696792939e-04+z_r(i,j,k)*(
     &                          3.191792361954220e-05+
     &                          6.177352634409320e-07*z_r(i,j,k)
     &                                                       ))))))
            else
               tclm(i,j,k,itemp)=14.6+6.7 *tanh(1.1*z_r(i,j,k)+15.9)
               tclm(i,j,k,isalt)=31.3-0.55*tanh(1.1*z_r(i,j,k)+15.9)
            endif
          enddo
        enddo
      enddo
#  elif defined JET
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
           cff1=TL1+ TL2*z_r(i,j,k) +TL3*exp(z_r(i,j,k)/ZL1)
           cff2=TC1+ TC2*z_r(i,j,k) +TC3*exp((z_r(i,j,k)-ZC1)/ZC2)
           cff3=YC1+YC2*tanh(YC3*xr(i,j)+YC4)
           tclm(i,j,k,itemp)=0.5*(cff2*(1+tanh((cff3-yr(i,j))/LCUR))
     &                           +cff1*(1+tanh((yr(i,j)-cff3)/LCUR)))
          enddo
        enddo
      enddo

#  elif defined ISWAKE
!
! There is no code here because for this case climatological data
! is saved from the initial conditions and remained unchanged after
! that.
!
#  else
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            tclm(i,j,k,itemp)=???
#   ifdef SALINITY
            tclm(i,j,k,isalt)=???
#   endif
          enddo
        enddo
      enddo
#  endif
      return
      end
# endif /* ANA_TCLIMA && TCLIMATOLOGY */
#endif /* SOLVE3D */


#if defined ANA_UCLIMA && defined UCLIMATOLOGY
      subroutine ana_uclima (tile)
      implicit none
      integer tile
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      call ana_uclima_tile   (istr,iend,jstr,jend)!,A2d(1,1),A2d(1,2))
      return
      end

      subroutine ana_uclima_tile (istr,iend,jstr,jend )!, wrk1,wrk2)

!---------------------------------------------------------------------
!  This routine sets analytical momentum climatology fields.
!---------------------------------------------------------------------

      implicit none
      integer istr,iend,jstr,jend, i,j,k
c      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk1, wrk2
# ifdef JET
      real TL0, TC0, cff1, cff2, cff3, cff4, cff5
# endif
# include "param.h"
# include "climat.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

! Set climatology fields for tracer type variables.
!---- ----------- ------ --- ------ ---- ----------

# if defined USWEST
      do j=jstrR,jendR
        do i=istrR,iendR
          ubclm(i,j)=0.
          vbclm(i,j)=0.
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            uclm(i,j,k)=0.
            vclm(i,j,k)=0.
          enddo
        enddo
      enddo
#  endif
# elif defined JET
      TL0=-TL1-
     &    0.5*TL2*ZREF**2-
     &    ZL1*TL3*exp(ZREF/ZL1)
      TC0=-TC1-
     &    0.5*TC2*ZREF**2-
     &    ZC2*TC3*exp((ZREF-ZC1)/ZC2)
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
         cff1=TL0*h(i,j)-
     &        0.5*TL1*h(i,j)**2+
     &        TL2*(h(i,j)**3)/6-
     &        (ZL1**2)*TL3*
     &         (1-exp(-h(i,j)/ZL1))
         cff2=TC0*h(i,j)-
     &        0.5*TC1*h(i,j)**2+
     &        TC2*(h(i,j)**3)/6-
     &        (ZC2**2)*TC3*
     &         exp(-ZC1/ZC2)*
     &         (1-exp(-h(i,j)/ZC2))
         cff3=YC1+YC2*
     &        tanh(YC3*xr(i,j)+YC4)
         cff4=-g*Tcoef*
     &        (1-(tanh((cff3-yr(i,j))/LCUR))**2)/
     &        (2*rho0*f(i,j)*LCUR*h(i,j))
! ubar
         wrk1(i,j)=cff4*(cff1-cff2)
! vbar
         cff5=YC2*YC3*
     &        (1-(tanh(YC3*xr(i,j)
     &            +YC4))**2)
         wrk2(i,j)=cff4*cff5*(cff1-cff2)
        enddo
      enddo
      do j=jstrR,jendR
        do i=istr,iendR
          ubclm(i,j)=0.5*(wrk1(i,j)+wrk1(i-1,j))
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          vbclm(i,j)=0.5*(wrk2(i,j)+wrk2(i,j-1))
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
           cff1=TL0+
     &        TL1*z_r(i,j,k)+
     &        0.5*TL2*(z_r(i,j,k)**2)+
     &        ZL1*TL3*
     &        exp(z_r(i,j,k)/ZL1)
           cff2=TC0+
     &        TC1*z_r(i,j,k)+
     &        0.5*TC2*(z_r(i,j,k)**2)+
     &        ZC2*TC3*
     &        exp((z_r(i,j,k)-ZC1)/ZC2)
           cff3=YC1+YC2*
     &        tanh(YC3*xr(i,j)+YC4)
           cff4=g*Tcoef*
     &        (1-(tanh((cff3-yr(i,j))/LCUR))**2)/
     &        (2*rho0*f(i,j)*LCUR)
! u
           wrk1(i,j)=cff4*(cff1-cff2)
! v
           cff5=YC2*YC3*
     &          (1-(tanh(YC3*xr(i,j)
     &              +YC4))**2)
           wrk2(i,j)=cff4*cff5*(cff1-cff2)
          enddo
        enddo
        do j=jstrR,jendR
          do i=istr,iendR
           uclm(i,j,k)=0.5*(wrk1(i,j)+wrk1(i-1,j))
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
           vclm(i,j,k)=0.5*(wrk2(i,j)+wrk2(i,j-1))
          enddo
        enddo
      enddo
#  endif
# elif defined ISWAKE
!
! There is no code here because for this case climatological data
! is saved from the initial conditions and remained unchanged after
! that.
!
# else
      do j=jstrR,jendR
        do i=istrR,iendR
          ubclm(i,j)=???
          vbclm(i,j)=???
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            uclm(i,j,k)=???
            vclm(i,j,k)=???
          enddo
        enddo
      enddo
#  endif
# endif
      return
      end
#endif /* ANA_UCLIMA && UCLIMATOLOGY */
!
#if defined ANA_WWAVE && defined SG_BBL96
      subroutine ana_wwave_tile (istr,iend,jstr,jend)

! Set wind induced wave amplitude, direction and period for
! bottom boundary layer formulation.

      implicit none
      integer istr,iend,jstr,jend, i,j
# include "param.h"
# include "bblm.h"
# include "grid.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

!  Set wind induced wave amplitude (m), direction (radians) and
!  period (s) at RHO-points.

      do j=jstrR,jendR
        do i=istrR,iendR
          Awave(i,j)=???
          Dwave(i,j)=???
          Pwave(i,j)=???
        enddo
      enddo
      return
      end
#endif /* ANA_WWAVE && SG_BBL96 */


!#if (defined PSOURCE || defined Q_PSOURCE) && defined ANA_PSOURCE
#if defined ANA_PSOURCE
                                  ! Set analytical tracer and mass point
      subroutine ana_psource_tile (istr,iend,jstr,jend) ! sources and sinks. River runoff can
!      subroutine ana_psource
      implicit none               ! be consider as a point source.
!      integer istr,iend,jstr,jend
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "sources.h"
# include "ocean3d.h"
# include "grid.h"
# include "mixing.h"
#  ifdef MPI
#   include "hidden_mpi_vars.h"
#  endif

      integer is, k, itrc, psrcgrd
      real cff, fac, csum, flow

      integer i,j,itpas, i1,i2,j1,j2,Npas
      integer istr,iend,jstr,jend
      real Zcff, Scff, GFT(9)
!      real, dimension(:,:), allocatable :: D
!      if (.not.allocated(D)) allocate (D(Nsrc,N))
      if (iic.eq.ntstart) then
# if defined BIOLOGY_BEC

      if (iic.eq.ntstart) then
      do j=jstr,jend
        do i=istr,iend
!          ammonium(i,j)=1000
          anthpco2(i,j)=1000
        enddo
      enddo
      endif

# endif

# if defined FLAT_BOTTOM
      ! hardcoded values
      Nsrc = 4 ! number of point sources/sinks
      flow = 4.0 ! volume flux m3/s

      ! tracers all in mmol/m3
      Npas = 2 ! number of tracers
      temp = 10
      salt = 1
      NH4 = 1000

      ! location of input
      ! positive u/v need to be placed on water
      ! negative u/v need to be placed on land
      Isrc(1) = 6
      Jsrc(1) = 5
      Isrc(2) = 4
      Jsrc(2) = 5
      Isrc(3) = 5
      Jsrc(3) = 6
      Isrc(4) = 5
      Jsrc(4) = 4
!      Isrc(1) = 7
!      Jsrc(1) = 5
!      Isrc(3) = 5
!      Jsrc(3) = 7

      do is=1,Nsrc
        write(*,*)
     &  '--> PSOURCE No.',is,'is ','', 
     &  ': (i,j)=(',Isrc(is),',',Jsrc(is),')'
      enddo

!      do is=1,Nsrc
!        do itrc=1,Npas
!          Lsrc(is,itrc)=.true.
!        enddo
!        write(*,*)
!     &  '--> PSOURCE No.',is,'is ','', 
!     &  ': (i,j)=(',Isrc(is),',',Jsrc(is),')'
!      enddo

#  ifdef MPI
      do is=1,Nsrc
        if ((Isrc(is)-iSW_corn).lt.iwest  .or. (Isrc(is)-iSW_corn).gt.ieast  .or.
     &      (Jsrc(is)-jSW_corn).lt.jsouth .or. (Jsrc(is)-jSW_corn).gt.jnorth ) then
          do itrc=1,Npas
            Lsrc(is,itrc)=.false.
          enddo
        else
          write(*,'(2X,A,X,I4,X,A,X,I3,A,2(I4,A))')
     &       '--> PSOURCE No.',is,'is found in node',mynode,
     &       ' : (i,j)=(',Isrc(is),',',Jsrc(is),')'
        endif
      enddo
#  endif
      
      do is=1,Nsrc
        do k=1,N
          Qshape(is,k) = 1./N
        enddo
      enddo
 
      ! direction of input
      Dsrc(1) = 0
      Dsrc(2) = 0
      Dsrc(3) = 1
      Dsrc(4) = 1
      ! positive u/v need to have negative sign for Qbar
      Qbar(1) = -flow/Nsrc
      Qbar(2) = flow/Nsrc
      Qbar(3) = -flow/Nsrc
      Qbar(4) = flow/Nsrc
      do is=1,Nsrc
!        Dsrc(is) = 0
!        Qbar(is) = flow/Nsrc
        do k=1,N
          Qsrc(is,k) = Qbar(is)*Qshape(is,k)
          Qsrctemp(k)= temp    *Qshape(is,k)
          Qsrcsalt(k)= salt    *Qshape(is,k)
          QsrcNH4(k)= NH4    *Qshape(is,k)

                Tsrc(is,k,itemp)=Qsrctemp(k)
                Tsrc(is,k,isalt)=Qsrcsalt(k)
                Tsrc(is,k,iNH4)=QsrcNH4(k)
          write(*,*) 'Qshape(',is,',',k,') = ',Qshape(is,k)
          write(*,*) 'Qsrc(',is,',',k,') = ',Qsrc(is,k)
          write(*,*) 'Tsrc(',is,',',k,',itemp) = ',Tsrc(is,k,itemp)
          write(*,*) 'Tsrc(',is,',',k,',isalt) = ',Tsrc(is,k,isalt)
        enddo ! k

      enddo ! is
      endif ! iic.eq.ntstart

! end FLAT_BOTTOM

# endif ! ends psource scenarios
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! define the position of the rivers (i,j) in the global grid !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
! J-grid locations (Isrc,Jsrc).  Currently, the direction can be along
! XI- (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The mass sources
! are located at U- or V-points so the grid locations should range
! from 1 =< Isrc =< L  and  1 =< Jsrc =< M.

!
! Set-up vertically integrated mass transport [m3/s] of point
! sources (these may time-dependent; positive in the positive U-
! or V-direction and vice-versa) and vertically distribute them
! according to mass transport profile chosen above.
!

      return
      end        !/* SUBROUTINE ana_psource */
# endif
c# else ! defined ANA_PSOURCE
!#endif           !/* (defined PSOURCE || defined Q_PSOURCE) && defined ANA_PSOURCE */

c      subroutine ana_empty
c      return
c      end
c#endif
