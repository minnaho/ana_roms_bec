#include "cppdefs.h"
#ifdef ANA_INITIAL

      subroutine ana_init (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call ana_init_tile (istr,iend,jstr,jend, A3d(1,1),A3d(1,2))
      return
      end

      subroutine ana_init_tile (istr,iend,jstr,jend, AA)

! Set initial conditions for momentum and tracer variables using
! analytical expressions.

      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j,k, itrc, i0,j0
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: AA

      real Ampl, fac, x,y, x0,y0, cff, cff1,cff2,cff3
# ifdef COLD_FILAMENT

# elif defined ISWAKE
      integer jwake,mid_wake
      real AmpR,AmpU, cffS, Zthrm, dZthrm, xitm, alpha, csa,sna
# endif
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "scalars.h"
# include "eos_vars.h"
!# ifdef ISWAKE
# ifdef ANA_BRY
#  include "boundary.h"
# else
#  include "climat.h"
# endif
!# endif

# include "compute_auxiliary_bounds.h"

      do j=jstrR,jendR       ! Set everything (except temperature
        do i=istrR,iendR     ! and salinity) to all-zero state, then
          zeta(i,j,1)=0.     ! modify some of the variables, if a
          ubar(i,j,1)=0.     ! non-trivial initialization required.
          vbar(i,j,1)=0.     ! Note: A code to initialize T [and S]
        enddo                ! must always be supplied for 3D
      enddo                  ! applications.
# ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
          enddo
        enddo
      enddo
# endif

# if defined SOLITON

! References for Equatorial Rossby Soliton:
!
!   Boyd, J. P., 1980:  Equatorial solitary waves. part i: Rossby
!                solitons.  J. Phys. Oceanogr., 10:1699-1717.
!
!   Boyd, J. P., 1985:  Equatorial solitary waves. part 3: Westward-
!                traveling modons.  J. Phys. Oceanogr., 15:46-54.

! For the purpose of testing the EXACT symmetry of the code there
! is an option to reverse beta-effect in "ana_grid" and, consistently,
! reverse velocities here.  This results in eastward traveling
! soliton, however it has the same energy and theoretically this
! solution should pass ETALON CHECK with exactly the same check
! points. In either case the soliton is initialized 1/3 into the
! domain in the direction of its propagation (so that it travels
! the remaining 2/3 of the domain before exiting through the
! boundary.


      if (fomn(istr,jend)-fomn(istr,jstr) .lt. 0.) then
C$OMP MASTER
          mpi_master_only write(*,'(/1x,2A/)')  '>>> Eastward-',
     &      'traveling soliton: reversing initial velocities.'
C$OMP END MASTER
         x0=xl/3.
      else
         x0=2.*xl/3.
      endif

c**   x0=xl/2.

      y0=el/2.
      cff1=0.395
      cff2=0.771*(cff1*cff1)
      do j=jstrR,jendR
        do i=istrR,iendR
          x=xr(i,j)-x0
          y=yr(i,j)-y0
          cff3=exp(-cff1*x)
          fac=cff2*(2.*cff3/(1.+cff3*cff3))**2
          zeta(i,j,1)=0.25*fac*(6.*y*y+3.)*exp(-0.5*y*y)
        enddo
      enddo

      do j=jstrR,jendR
        do i=istr,iendR
          x=0.5*(xr(i-1,j)+xr(i,j))-x0
          y=0.5*(yr(i-1,j)+yr(i,j))-y0
          cff3=exp(-cff1*x)
          fac=cff2 * (2.*cff3/(1.+cff3*cff3))**2
          ubar(i,j,1)=0.25*fac*(6.*y*y-9.)*exp(-0.5*y*y)
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          x=0.5*(xr(i,j-1)+xr(i,j))-x0
          y=0.5*(yr(i,j-1)+yr(i,j))-y0
          cff3=exp(-cff1*x)
          fac=cff2 * (2.*cff3/(1.+cff3*cff3))**2
          vbar(i,j,1)=2.*fac*y*(-2.*cff1*tanh(cff1*x))
     &                                  *exp(-0.5*y*y)
        enddo
      enddo
      if (fomn(istr,jend)-fomn(istr,jstr) .lt. 0.) then
        do j=jstrR,jendR               ! Reverse velocities, if the
          do i=istrR,iendR             ! sign of gradient of Coriolis
            ubar(i,j,1)=-ubar(i,j,1)   ! parameter "beta" is negative
            vbar(i,j,1)=-vbar(i,j,1)   ! (see ana_grid). In this case
          enddo                        ! soliton should travel to the
        enddo                          ! east, while maintaining all
      endif                            ! energetics exactly the same.

# elif defined WAVE_RAD

      x0=xl/2.                     ! Set initial perturbation in
      y0=el/2.                     ! in free surface and velocities
      cff=64./xl                   ! for isotropically expanding
      do j=jstr,jend               ! circular wave.
        do i=istr,iend
          x=cff*(xr(i,j)-x0)
          y=cff*(yr(i,j)-y0)
          cff1=sqrt(x*x+y*y)
          if (cff1.gt.9. .and. cff1.lt.11.) then
            zeta(i,j,1)=1.+cos(pi*(cff1-10.))
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          x=cff*(0.5*(xr(i-1,j)+xr(i,j))-x0)
          y=cff*(               yr(i,j) -y0)
          cff1=sqrt(x*x+y*y)
          if (cff1.gt.9. .and. cff1.lt.11.) then
            ubar(i,j,1)=(1.+cos(pi*(cff1-10.)))
     &      *(x/cff1)*sqrt(2.*g/(h(i-1,j)+h(i,j)))
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          x=cff*(               xr(i,j) -x0)
          y=cff*(0.5*(yr(i,j-1)+yr(i,j))-y0)
          cff1=sqrt(x*x+y*y)
          if (cff1.gt.9. .and. cff1.lt.11.) then
            vbar(i,j,1)=(1.+cos(pi*(cff1-10.)))
     &      *(y/cff1)*sqrt(2.*g/(h(i,j-1)+h(i,j)))
          endif
        enddo
      enddo


# elif defined BASIN
#  ifdef SOLVE3D
      cff1=(44.690/39.382)**2
      cff2=cff1*(rho0*800./g)*(5.0e-5/((42.689/44.690)**2))
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=cff2*exp(z_r(i,j,k)/800.)
     &                 *(0.6-0.4*tanh(z_r(i,j,k)/800.))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined CANYON_A
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=T0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined CANYON_B
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=3.488*exp(z_r(i,j,k)/800.)
     &        *(1.-0.666666666666*tanh(z_r(i,j,k)/800.))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined DAMEE_B
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c            t(i,j,k,1,itemp)=0.2+6.*exp(z_r(i,j,k)/2500.)
c     &              +30./(1.+exp((-200.-z_r(i,j,k))/200.))

            t(i,j,k,1,itemp)=0.2+6.*exp(z_r(i,j,k)/2500.)
     &                    +20.*exp(z_r(i,j,k)/250.)

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=35.
#   endif
          enddo
        enddo
      enddo
#  endif
# elif defined EKMAN_SPIRAL

# elif defined GRAV_ADJ
#  ifdef SOLVE3D
      i0=LLm/2
#   ifdef MPI
      i0=i0 -iSW_corn
#   endif
      cff=Tcoef*(35.-5.)*20./(4.*rho0)  ! A provision is made to form
      do j=jstrR,jendR                  ! a 2-point transitional layer
        do i=istrR,iendR                ! between the left- and right-
          if (i.lt.i0) then             ! constant values in order to
            t(i,j,1,1,itemp)=5.         ! "soften" the step function.
            zeta(i,j,1)=-cff
          elseif (i.eq.i0) then         ! Setting free surface to
            t(i,j,1,1,itemp)=5. +3.75   !
            zeta(i,j,1)=-cff*0.75       !         Tcoef*DeltaT
          elseif (i.eq.i0+1) then       ! zeta = -------------- * h
            t(i,j,1,1,itemp)=35.-3.75   !           2.*rho0
            zeta(i,j,1)=+cff*0.75       !
          else                          ! balances the barotropic mode
            t(i,j,1,1,itemp)=35.        ! making the solution be almost
            zeta(i,j,1)=+cff            ! free of barotropic seiche
          endif                         ! mode.
          zeta(i,j,2)=zeta(i,j,1)
          t(i,j,1,2,itemp)=t(i,j,1,1,itemp)    ! Canonically, when
        enddo                                  ! initialized to zero,
        do k=2,N                               ! free surface
          do i=istrR,iendR                     ! oscillates with its
            t(i,j,k,1,itemp)=t(i,j,1,1,itemp)  ! own frequency,
            t(i,j,k,2,itemp)=t(i,j,1,1,itemp)  ! period = 7900 sec,
          enddo                                ! amplitude=-0.06...0.06
        enddo   !<-- k                         ! provided that
      enddo   !<-- j                           ! Tcoef=0.2.
#  endif
# elif defined COLD_FILAMENT

      do j=jstr,jend+1            ! Initialize density
        do i=istr,iend+1          ! integration used below
          AA(i,j,0)=0.            ! to achieve geostrophic
        enddo                     ! balance
      enddo

      x0=xl/2.
      do k=1,N,+1
        do j=jstr,jend+1
          do i=istr,iend+1
            cff1=0.25*( z_r(i,j,k)+z_r(i-1,j,k)+z_r(i,j-1,k)
     &                                         +z_r(i-1,j-1,k))
            cff=0.25*(xr(i,j)+xr(i-1,j)+xr(i,j-1)+xr(i-1,j-1))-x0

            cff=exp(-(cff/2.e+3)**2)

            AA(i,j,k)=0.3*tanh( (cff1+60.-60.*cff )/(20.*(1.+cff) ))
c     &               -0.1*(cff1/500.)**2
          enddo
        enddo

        do j=jstr,jend
          do i=istr,iend
            t(i,j,k,1,itemp)=0.25*( AA(i,j,k)+AA(i+1,j,k)
     &                          +AA(i,j+1,k)+AA(i+1,j+1,k)
     &                                                   )
          enddo
        enddo

        if (k.gt.1) then
          cff=0.125*Tcoef/rho0
          do j=jstr,jend+1
            do i=istr,iend+1
              AA(i,j,0)=AA(i,j,0) +cff*(AA(i,j,k)+AA(i,j,k-1))*(
     &            z_r(i,j,k)+z_r(i-1,j,k)+z_r(i,j-1,k)+z_r(i-1,j-1,k)
     &           -z_r(i,j,k-1)-z_r(i-1,j,k-1)
     &                               -z_r(i,j-1,k-1)-z_r(i-1,j-1,k-1)
     &                                                              )
            enddo
          enddo
        endif

        cff=0.5 *g/1.e-4
        do j=jstr,jend
          do i=istrU,iend
            u(i,j,k,1)=-cff*(pn(i,j)+pn(i-1,j))*(AA(i,j+1,0)-AA(i,j,0))
            ubar(i,j,1)=ubar(i,j,1) +(Hz(i-1,j,k)+Hz(i,j,k))*u(i,j,k,1)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            v(i,j,k,1)=+cff*(pm(i,j)+pm(i,j-1))*(AA(i+1,j,0)-AA(i,j,0))
            vbar(i,j,1)=vbar(i,j,1) +(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,1)
          enddo
        enddo
      enddo  !<-- k

      do j=jstr,jend
        do i=istr,iend
          zeta(i,j,1)=0.009 + 0.25*( AA(i,j,0)+AA(i+1,j,0)
     &                         +AA(i,j+1,0)+AA(i+1,j+1,0))
        enddo
      enddo
      do j=jstr,jend
        do i=istrU,iend
          ubar(i,j,1)=ubar(i,j,1)/( z_w(i,j,N)+z_w(i-1,j,N)
     &                             -z_w(i,j,0)-z_w(i-1,j,0))
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          vbar(i,j,1)=vbar(i,j,1)/( z_w(i,j,N)+z_w(i,j-1,N)
     &                             -z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo





#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          zeta(istr-1,j,1)=zeta(istr,j,1)
          ubar(istr-1,j,1)=0.
        enddo
        do k=1,N
          do j=jstr,jend
            t(istr-1,j,k,1,itemp)=t(istr,j,k,1,itemp)
            u(istrU-1,j,k,1)=0.
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          zeta(iend+1,j,1)=zeta(iend,j,1)
          ubar(iend+1,j,1)=0.
        enddo
        do k=1,N
          do j=jstr,jend
            t(iend+1,j,k,1,itemp)=t(iend,j,k,1,itemp)
            u(iend+1,j,k,1)=0.
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          zeta(i,jstr-1,1)=zeta(i,jstr,1)
          vbar(i,jstr-1,1)=0.
        enddo
        do i=istrR,iendR
          ubar(i,jstr-1,1)=ubar(i,jstr,1)
        enddo
        do k=1,N
          do i=istr,iend
            t(i,jstr-1,k,1,itemp)=t(i,jstr,k,1,itemp)
            v(i,jstrV-1,k,1)=0.
          enddo
          do i=istrR,iendR
            u(i,jstr-1,k,1)=u(i,jstr,k,1)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          zeta(i,jend+1,1)=zeta(i,jend,1)
          vbar(i,jend+1,1)=0.
        enddo
        do i=istrR,iendR
          ubar(i,jend+1,1)=ubar(i,jend,1)
        enddo
        do k=1,N
          do i=istr,iend
            t(i,jend+1,k,1,itemp)=t(i,jend,k,1,itemp)
            v(i,jend+1,k,1)=0.
          enddo
          do i=istrR,iendR
            u(i,jend+1,k,1)=u(i,jend,k,1)
          enddo
        enddo
      endif
#  endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstrR,jendR
          vbar(istr-1,j,1)=vbar(istr,j,1)
        enddo
        do k=1,N
          do j=jstrR,jendR
            v(istr-1,j,k,1)=v(istr,j,k,1)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstrR,jendR
          vbar(iend+1,j,1)=vbar(iend,j,1)
        enddo
        do k=1,N
          do j=jstrR,jendR
            v(iend+1,j,k,1)= v(iend,j,k,1)
          enddo
        enddo
      endif
#   ifndef NS_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        zeta(istr-1,jstr-1,1)=zeta(istr,jstr,1)
        do k=1,N
          t(istr-1,jstr-1,k,1,itemp)=t(istr,jstr,k,1,itemp)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        zeta(iend+1,jstr-1,1)=zeta(iend,jstr,1)
        do k=1,N
          t(iend+1,jstr-1,k,1,itemp)=t(iend,jstr,k,1,itemp)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        zeta(istr-1,jend+1,1)=zeta(istr,jend,1)
        do k=1,N
          t(istr-1,jend+1,k,1,itemp)=t(istr,jend,k,1,itemp)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        zeta(iend+1,jend+1,1)=zeta(iend,jend,1)
        do k=1,N
          t(iend+1,jend+1,k,1,itemp)=t(iend,jend,k,1,itemp)
        enddo
      endif
#   endif
#  endif


# elif defined OVERFLOW
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c           t(i,j,k,1,itemp)=T0-0.5*T0*(1.+tanh((yr(i,j)-25000.)
c     &                                                   /1000.))

            t(i,j,k,1,itemp)=T0*(0.5-0.5*tanh( yr(i,j)/1000.-25.))

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined NJ_BIGHT
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (z_r(i,j,k).ge.-15.0) then
              t(i,j,k,1,itemp)=2.049264257728403e+01-z_r(i,j,k)*(
     &                         2.640850848793918e-01+z_r(i,j,k)*(
     &                         2.751125328535212e-01+z_r(i,j,k)*(
     &                         9.207489761648872e-02+z_r(i,j,k)*(
     &                         1.449075725742839e-02+z_r(i,j,k)*(
     &                         1.078215685912076e-03+z_r(i,j,k)*(
     &                         3.240318053903974e-05+
     &                         1.262826857690271e-07*z_r(i,j,k)
     &                                                     ))))))
              t(i,j,k,1,isalt)=3.066489149193135e+01-z_r(i,j,k)*(
     &                         1.476725262946735e-01+z_r(i,j,k)*(
     &                         1.126455760313399e-01+z_r(i,j,k)*(
     &                         3.900923281871022e-02+z_r(i,j,k)*(
     &                         6.939014937447098e-03+z_r(i,j,k)*(
     &                         6.604436696792939e-04+z_r(i,j,k)*(
     &                         3.191792361954220e-05+
     &                         6.177352634409320e-07*z_r(i,j,k)
     &                                                     ))))))
            else
               t(i,j,k,1,itemp)=14.6+6.7 *tanh(1.1*z_r(i,j,k)+15.9)
               t(i,j,k,1,isalt)=31.3-0.55*tanh(1.1*z_r(i,j,k)+15.9)
            endif
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
# elif defined SEAMOUNT
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#   ifdef NONLIN_EOS
            t(i,j,k,1,itemp)=4.0  ! Just spatially uniform T,S to
            t(i,j,k,2,isalt)=33.0 ! illustrate compressibility effect.
#   else

! There are four possible choices here: exponential profile which can
! be initialized as either a set of values at locations z_r (to test
! prsgrd31, 32, 32A); or in finite volume sense (to test prsgrd40, 41,
! 42, 44). Also there are two similar versions of linear ! profiles.

#    define AMPL 3.
#    define STRAT 500.

           t(i,j,k,1,itemp)=AMPL*exp(z_r(i,j,k)/STRAT)

c**      t(i,j,k,1,itemp)=AMPL*exp((z_r(i,j,k)-0.0003*yr(i,j))/STRAT)

c            t(i,j,k,1,itemp)=AMPL*STRAT*( exp(z_w(i,j,k)/STRAT)
c     &                                   -exp(z_w(i,j,k-1)/STRAT)
c     &                                )/(z_w(i,j,k)-z_w(i,j,k-1))

! Linear profile to check exact error cancellation.


c            t(i,j,k,1,itemp)=AMPL*(1.+z_r(i,j,k)/5000.)

c           t(i,j,k,1,itemp)=0.5*(z_w(i,j,k)+z_w(i,j,k-1))

c            t(i,j,k,1,itemp)=AMPL*(1.+0.5*(z_w(i,j,k)
c     &                            +z_w(i,j,k-1))/5000.)

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)

#    undef STRAT
#    undef AMPL
#   endif
          enddo
        enddo
      enddo
/*
#   define AREA 7.e+9
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=(xr(i,j)-0.8e+5)**2 + (yr(i,j)-0.8e+5)**2
          if (cff .lt. AREA) then
            zeta(i,j,1)=1.0*( 1.-cff/AREA )**2
          endif

          cff=(xr(i,j)-2.4e+5)**2 + (yr(i,j)-2.4e+5)**2
          if (cff .lt. AREA) then
            zeta(i,j,1)=1.0*( 1.-cff/AREA )**2
          endif

          cff=(xr(i,j)-0.8e+5)**2 + (yr(i,j)-2.4e+5)**2
          if (cff .lt. AREA) then
            zeta(i,j,1)=-1.0*( 1.-cff/AREA )**2
          endif

          cff=(xr(i,j)-2.4e+5)**2 + (yr(i,j)-0.8e+5)**2
          if (cff .lt. AREA) then
            zeta(i,j,1)=-1.0*( 1.-cff/AREA )**2
          endif
        enddo
      enddo
*/

#  endif
# elif defined ISWAKE
#  ifdef SOLVE3D

! Analytical initial conditions for 3D Island Wake problem:
!-------------------- ------- ---------- --- -- ------ ---- --------
! free-surface, velocity and temperature (actually negative density --
! it is assumed that thermal expansion coefficient is equal to 1 for
! this problem) fields are set to satisfy thermal wind balance
! equation
!                      du      g       d rho
!                 f * ---- = ------ * -------     [here rho=-T]
!                      dz     rho0       dy
!
! subject to the conditions that surface velocity is in geostrophic
! balance with free surface, and velocity vanishes at infinitely large
! negative depth.  Both velocity and temperature field are set via
! the same hyperbolic shape function
!
!                  1     1         z-Zthrm(y)
!                 --- + --- * tanh ----------
!                  2     2           dZthrm
!
! where "Zthrm" and "dZthrm" set thermocline depth and sharpness
! respectively, and y-dependency of "Zthrm" creates thermocline slope
! needed for thermal balance; free surface is determined by vertical
! integral of temperature profile and the assumption that total
! pressure (baroclinic+barotropic) gradient vanishes at infinite
! negative depth. Integration is done analytically.

      AmpR=3.0           ! magnitude of density perturbation [kg/m^3]
      AmpU=0.2           ! velocity magnitude above thermocline [m/s]
      Zthrm=-120.        ! median depth of thermocline [m]
      dZthrm=60.         ! thickness of thermocline [m]

                               ! set thermocline slope parameter to
      cffS=AmpU*rho0/(AmpR*g)  ! achieve geostrophic balance between
      y0=el/2.                 ! XI-component of baroclinic velocity
                               ! and density field.

      cff=1./(1.+exp(-2.*Zthrm/dZthrm))
      do j=jstrR,jendR
        do i=istrR,iendR
          xitm=(-Zthrm  -cffS*f(i,j)*(yr(i,j)-y0))/dZthrm
          zeta(i,j,1)=(AmpR/rho0) *0.5*dZthrm
     &                       *log( cff*(1.+exp(2.*xitm)) )
#   ifdef MASKING
     &                                 *rmask(i,j)
#   endif
          ubar(i,j,1)=0.
          vbar(i,j,1)=0.
        enddo
      enddo
      do k=N,1,-1
        do j=jstrR,jendR
          do i=istrR,iendR
            xitm=(z_r(i,j,k)-Zthrm -cffS*f(i,j)*(yr(i,j)-y0))/dZthrm
            t(i,j,k,1,itemp)= 0.5*AmpR*tanh(xitm)
#   ifdef MASKING
     &                                 *rmask(i,j)
#   endif
          enddo
        enddo
        do j=jstrR,jendR
          do i=istr,IendR
            xitm=( 0.5*(z_r(i,j,k)+z_r(i-1,j,k)) -Zthrm
     &                    -cffS*0.5*(f(i,j)+f(i-1,j))
     &                  *(0.5*(yr(i,j)+yr(i-1,j))-y0) )/dZthrm

            u(i,j,k,1)=AmpU*(0.5 + 0.5*tanh(xitm))
#   ifdef MASKING
     &                                  *umask(i,j)
#   endif
            ubar(i,j,1)=ubar(i,j,1) +(Hz(i-1,j,k)+Hz(i,j,k))
     &                                           *u(i,j,k,1)
          enddo
        enddo
      enddo
      do j=jstrR,jendR
        do i=istr,IendR
          ubar(i,j,1)=ubar(i,j,1)/( z_w(i,j,N)+z_w(i-1,j,N)
     &                             -z_w(i,j,0)-z_w(i-1,j,0))
        enddo
      enddo

#   define SIGNELL
#   ifdef SIGNELL

      AmpU=0.15                       ! This setup is originally from
      cff=AmpU/g                      ! Rich Signell, who chose to pose
      x0=xl/2.                        ! a non-stratified problem to be
      y0=el/2.                        ! solved by 3D code.

      j0=MMm/2
#    ifdef MPI
      j0=j0 -jSW_corn
#    endif


      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=-cff*f(i,j)*(yr(i,j)-y0)
          ubar(i,j,1)=AmpU
          vbar(i,j,1)=0.
        enddo
      enddo
      do k=1,N
        do j=jstrR,jendR              ! Temperature field just below
          do i=istrR,iendR            ! is passive tracer initialized
            u(i,j,k,1)=ubar(i,j,1)    ! as a "smoothed" step function
            v(i,j,k,1)=vbar(i,j,1)    ! for the flow visualization
                                      ! purpose. It does not affect
            if (j.lt.j0) then         ! the dynamics of the flow.
              t(i,j,k,1,itemp)=-1.
            elseif (j.eq.j0) then
              t(i,j,k,1,itemp)=-0.333333333333333
            elseif (j.eq.j0+1) then
              t(i,j,k,1,itemp)=+0.333333333333333
            else
              t(i,j,k,1,itemp)=+1.
            endif

          enddo
        enddo
      enddo
#   endif
#  else
      AmpU=0.2            ! 2D setup for flow around cylinder
      cff=AmpU/g          ! problem: set u-velocity constant and
      alpha=0.3           ! free surface from geostrophic balance.
      csa=cos(alpha)
      sna=sin(alpha)
      x0=xl/2.
      y0=el/2.
      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=-cff*f(i,j)*((yr(i,j)-y0)*csa-(xr(i,j)-y0)*sna)
          ubar(i,j,1)=AmpU*csa
          vbar(i,j,1)=AmpU*sna
        enddo
      enddo
#  endif  /* SOLVE3D */

#  ifdef ANA_BRY
      if (WESTERN_EDGE) then
#   ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#   endif
#   ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#    ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#    endif
#    ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#    endif
          enddo
        enddo
#   endif
      endif

      if (EASTERN_EDGE) then
#   ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#   endif
#   ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#    ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#    endif
#    ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#    endif
          enddo
        enddo
#   endif
      endif

#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */


# elif defined SHELFRONT
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=T0+2.5*tanh((yr(i,j)-50000.0)/20000.0)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,1,isalt)=S0  +  tanh((yr(i,j)-50000.0)/20000.0)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
# elif defined TASMAN_SEA
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=17.4+z_r(i,j,k)*(1.727e-2+z_r(i,j,k)*
     &                           (5.707e-06+z_r(i,j,k)*(5.921e-10)))
            t(i,j,k,1,isalt)=35.08+z_r(i,j,k)*(7.56e-4+z_r(i,j,k)*
     &                           (3.185e-07+z_r(i,j,k)*(3.702e-11)))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
# elif defined UPWELLING
#  ifdef SOLVE3D
      j0=MMm/2
#   ifdef MPI
      j0=j0 -jSW_corn
#   endif

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c            t(i,j,k,1,itemp)=14.+8.*exp(z_r(i,j,k)/50.)


#   define Z0 (-35.)
#   define THKNSS 6.5
#   define Z1 (-75.)
#   define STRAT 150.

#   define THKNSS 25.0
#   define Z1 (-120.)
#   define STRAT 75.

c           t(i,j,k,1,itemp)=18.+4.*tanh((z_r(i,j,k)-Z0)/THKNSS)
c     &                                  +(z_r(i,j,k)-Z1)/STRAT

           t(i,j,k,1,itemp)=14.+4.*THKNSS*log(
     &                          cosh((z_w(i,j,k )-Z0)/THKNSS)
     &                         /cosh((z_w(i,j,k-1)-Z0)/THKNSS)
     &                             )/(z_w(i,j,k)-z_w(i,j,k-1))
     &               +((z_w(i,j,k)+z_w(i,j,k-1))/2.-Z1)/STRAT

#   undef Z0
#   undef THKNSS
#   undef Z1
#   undef STRAT





            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)

#   ifdef SALINITY
c*            t(i,j,k,1,isalt)=1.E-4*yr(i,j)-S0 +0.1
            t(i,j,k,1,isalt)=1.

            if (j.lt.j0) then          ! A provision is made to form
              t(i,j,k,1,isalt)=0.      ! a 2-point transitional layer
            elseif (j.eq.j0) then      ! between S=0 and S=+1, hence
              t(i,j,k,1,isalt)=0.2     ! "softening" the step function.
            elseif (j.eq.j0+1) then
              t(i,j,k,1,isalt)=0.8
            else
              t(i,j,k,1,isalt)=1.
            endif

            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          enddo
        enddo
      enddo
#  endif
# elif defined  RIVER
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=4.+10.*exp(z_r(i,j,k)/50.)
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
# elif defined PACIFIC_2D
      Ampl=1.0
      cff1=800./dble(LLm)

      if (LLm.eq.840 .and. MMm.eq.512) then
        i0=519
        j0=283
      elseif (LLm.eq.1680 .and. MMm.eq.1024) then
        i0=1036
        j0=565
      else
        i0=LLm/2
        j0=MMm/2
      endif
#  ifdef MPI
      i0=i0 -iSW_corn
      j0=j0 -jSW_corn
#  endif

      cff1=800./dble(LLm)
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=cff1*sqrt(dble(i-i0)**2+dble(j-j0)**2)
          if (cff.lt.pi) then
            zeta(i,j,1)=Ampl*(0.5+0.5*cos(cff))
          endif
        enddo
      enddo
# elif defined FLAT_BOTTOM
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            include "ideal_winter_temp.h"
#   ifdef SALINITY
            include "ideal_winter_salt.h"
#   endif
#   ifdef BIOLOGY
              do itrc=3,NT
                t(i,j,k,1,itrc)=0.
                t(i,j,k,2,itrc)=0.
              enddo
#   endif
          enddo
        enddo
      enddo

#  endif /* SOLVE3D */ 

#  ifdef ANA_BRY
      if (WESTERN_EDGE) then
#   ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#   endif
#   ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#    ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#    endif
#    ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#    endif
          enddo
        enddo
#   endif
      endif

      if (EASTERN_EDGE) then
#   ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#   endif
#   ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#    ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#    endif
#    ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#    endif
          enddo
        enddo
#   endif
      endif

#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */
!!!!!!!!!!!!!!!!!!!!!!!!! end FLAT_BOTTOM !!!!!!!!!!!!!!!!!!!!!!!!!!!!

# elif defined BALTIC
#  ifdef SOLVE3D
      i0=116
#   ifdef MPI
      i0=i0 -iSW_corn
#   endif

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (i.lt.i0) then

              t(i,j,k,1,itemp)=8.+8.*exp(z_r(i,j,k)/40.) !   Baltic sea code must go here.
              t(i,j,k,1,isalt)=16.
            else

            endif




            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif

# else
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=???
            t(i,j,k,2,itemp)=???
            t(i,j,k,1,isalt)=???
            t(i,j,k,2,isalt)=???
          enddo
        enddo
      enddo
#  endif
# endif  /* end of the switch over applications */

# ifdef EXCHANGE
      call exchange2d_3_tile (istr,iend,jstr,jend,
     &                        zeta(START_2D_ARRAY,1),
     &                        ubar(START_2D_ARRAY,1),
     &                        vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
      call exchange_2_tile (istr,iend,jstr,jend,
     &                      u(START_2D_ARRAY,1,1), N,
     &                      v(START_2D_ARRAY,1,1), N)
      do itrc=1,NT
        call exchange_tile (istr,iend,jstr,jend,
     &                 t(START_2D_ARRAY,1,1,itrc), N)
      enddo
#  endif
#  ifdef ISWAKE
#   ifdef UCLIMATOLOGY
      call exchange2d_2_tile (istr,iend,jstr,jend, ubclm,vbclm)
#   endif
#   ifdef SOLVE3D
#    ifdef UCLIMATOLOGY
      call exchange_2_tile (istr,iend,jstr,jend, uclm,N, vclm,N)
#    endif
#    ifdef TCLIMATOLOGY
      call exchange_tile (istr,iend,jstr,jend,
     &                    tclm(START_2D_ARRAY,1,itemp), N)
#    endif
#   endif
#  endif
# endif  /* EXCHANGE */
      return
      end
#else
      subroutine ana_init_empty
      end
#endif /* ANA_INITIAL */
