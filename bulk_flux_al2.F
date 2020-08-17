#include "cppdefs.h"
#ifdef BULK_FLUX
!
!==================================================================================
!
!  This routine computes the surface fluxes through a Bulk 
!  formulation
!
!==>  The Bulk formulation used here is adapted from the CORE/CCSM 
!  algorithm described in 
!
!  * Large, W.G. and S. G. Yeager, 2004: Diurnal to decadal global 
!    forcing for ocean and sea-ice models: The data sets and flux 
!    climatologies. Technical Report TN-460+STR, NCAR, 105pp.
!
!  * Large, W.G., 2006: Surface Fluxes for Practitioners of Global 
!    Ocean Data Assimilation", in Ocean Weather and Forecasting. 
!    Ed. E. Chassignet and J. Verron, Springer, 229-270.
!  
!    => T is given in Celsius
!    => Q is the specific humidity in kg/kg
!    => RAIN is in kg meter-2 s-1
!
!    The notations in this subroutine are meant to be consistent with
!    the CCSM subroutine shr_flux_mod.F90
!
!=================================================================================
!
      subroutine bulk_flux_tile (istr,iend,jstr,jend, ierr) 
      implicit none
      integer  istr,iend,jstr,jend, ierr,   imin,imax,jmin,jmax
# include "../param.h"
# include "../grid.h"
# include "../ocean3d.h"
# include "../forces.h"
# include "../scalars.h"
      real    tol, ustarold, q10, t10, cff, cff1,cff2,cff3
      real    ithv, iqhv, vmag, sst, ssq, al2, alz, ssu, ssv
      real    delt, delq, evap, lwup, hfsen, hflat
      real    rdn,ustar,qstar, tstar
      real    cph, hol, psixh, psimh, psixhu, psixh2
      real    rhn, ren, cdn, xqq, qsat, psimhu, cp_a
      
# define tind nstp    
      
       real   fwflx, cff7,cff8,al2
       integer it7,it8,niter,  kiter, i,j
       real, parameter ::
     &            umin   = 0.5,            
     &            zref   = 10.,                      !<== this variable should be read in the input file for winds   (attribute level)
     &            ztref  = 2.,                    !<== this variable should be read in the input file for tracers (attribute level)
     &            CtoK   = 273.16,
     &            cpvir  = 0.801640362327584,
     &            zvir   = 0.607792984078361,
     &           stebol  = 5.670400E-08,
     &            latvap = 2.5e6,
     &            irhoa  = 0.819672,                !<== 1/rhoA (rhoA ~ 1.22 kg.m^-3)
     &            rho_a  = 1.22,        
     &            a1     = 0.000142, 
     &            a2     = 0.0000764,
     &            a3     = 0.0027
          
#  ifdef EW_PERIODIC
      imin=istr-1
      imax=iend+1
#  else
      if (WESTERN_EDGE) then
        imin=istr
      else
        imin=istr-1
      endif
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=jstr-1
      jmax=jend+1
#  else
      if (SOUTHERN_EDGE) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=jend
      else
        jmax=jend+1
      endif
#  endif
#  define I_EXT_RANGE imin,imax
#  define J_EXT_RANGE jmin,jmax

#if defined SALINITY && defined SALINITY_RESTORING
      it7=3-itsss
      it8=itsss
      cff7=sss_time(it8)-time
      cff8=time-sss_time(it7)
      cff=1./(cff7+cff8)
      cff7=cff7*cff
      cff8=cff8*cff
#endif

!++      al2   = 1.60943791243410         !<== al2 = log(10./2.)

      al2   = log(zref/ztref)
      alz   = 0.                       !<== alz = log(zref/10.) = 0.

      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          if (rmask(i,j) .gt. 0.5 ) then                   
!
!++===========================================
! surface wind magnitude   (m/s)
!++===========================================               
!
# ifdef COUPLED_SURF_CURR
            ssu  = 0.5*(u(i,j,N,tind)+u(i+1,j,N,tind))
            ssv  = 0.5*(v(i,j,N,tind)+v(i,j+1,N,tind))
#else
            ssu  = 0.
            ssv  = 0.
#endif                  
      
            vmag = max( umin,sqrt( 
     &             (uwnd(i,j)-ssu)**2
     &            +(vwnd(i,j)-ssv)**2  ) )          !<-- assumes that winds are at rho points
!++===========================================
! Temperature / Humidity
!++===========================================
            sst  = t(i,j,N,tind,itemp) + CtoK              !<-- SST is needed in K for upward longwave flux and ssq 
            ssq  = irhoA*0.98*640380.*exp(-5107.4/sst)     !<== qsat = 640380. * exp(-5107.4/max(Tk,1.e-14)) 
            t10  = tair(i,j) + CtoK
            q10  = qair(i,j)                                   
            delt = t10 - sst
            delq = q10 - ssq
            ren  = 0.0346
!++===========================================           
!++ First guess for tstar,qstar,rdn and ustar  
!++===========================================                                   
            rdn   = sqrt( a1 + a2*vmag + a3 / vmag )
            ustar = rdn * vmag
            if (delt .gt. 0) then
              tstar=0.018*delt
            else
              tstar=0.0327*delt
            endif            
            qstar = ren*delq
!++===========================================
!++ Start iterations here  
!++===========================================                        

            do niter = 1,5   !--> start iterations here   
  
              cff  = 1./ ( 1.+ zvir*q10 )               !<== qair is strictly positive 
              ithv = cff / t10                          !<== tair can't be "zero" here
              iqhv = zvir*cff

              hol  = vonKar*g*zref*(tstar*ithv+qstar*iqhv)/ustar**2

              if (hol .gt. 0.) then   
                hol    = min(10., hol)
                psimh  = -5. * hol
                psixh  = -5. * hol 
                psixh2 = 0.2*psixh               
                rhn    = 0.018                       
              else
                hol    = max(-10., hol)               
                xqq    = sqrt( sqrt(1.-16.*hol) )
                psimh  = psimhu(xqq)
                psixh  = psixhu(xqq)
                xqq    = sqrt( sqrt(1.-3.2*hol) )       ! 2.*16 / 10.  = 3.2
                psixh2 = 2.*log( 0.5*(1.+xqq*xqq) )            
                rhn   = 0.0327            
              endif         
              cff = vonKar*vmag/sqrt( (vonKar + rdn*(alz-psimh))**2 ) !<== u10n                    
              rdn = sqrt( a1 + a2*cff + a3 / cff )                    !<== sqrt(CDN)

              cff  = ( psixh - psixh2 - al2 - alz ) / vonKar  
              t10  = tair(i,j) - cff*tstar + CtoK
              q10  = qair(i,j) - cff*qstar
              delt = t10 - sst
              delq = q10 - ssq
              
!++
!++ update ustar,qstar,tstar
!++
              ustar = vonKar * vmag * rdn / (vonKar + rdn*(alz-psimh)) 
              tstar = vonKar * delt * rhn / (vonKar + rhn*(alz-psixh)) 
              qstar = vonKar * delq * ren / (vonKar + ren*(alz-psixh)) 
!++
            enddo
!++===========================================
!++ End iteration  
!++=========================================== 






          
!++===========================================
!+++ Compute the air-sea  fluxes           
!++===========================================        
            cp_a       = 1004.64*(1. + cpvir*ssq)
            hfsen      =  cp_a*rho_a*ustar*tstar
            hflat      =  latvap*rho_a*ustar*qstar             
            evap       = -rho_a*ustar*qstar          
            lwup       = -stebol*sst**4
            fwflx      =  0.001*(evap-rain(i,j))                   !<== 0.001 is for the conversion from [kg.m-2.s-1]   to [m.s-1] 
!++===========================================        
              cff=1./(rho0*vmag)
              sustr(i,j) = cff*rho_a*ustar*ustar*(uwnd(i,j)-ssu)
              svstr(i,j) = cff*rho_a*ustar*ustar*(vwnd(i,j)-ssv)
              cff=1./(rho0*Cp)
              stflx(i,j,itemp) = cff*(
     &          radsw(i,j)+hflat+hfsen+radlw(i,j)+lwup
     &                                 )
              srflx(i,j)       = cff*radsw(i,j)          
#ifdef SALINITY
              stflx(i,j,isalt) = fwflx*t(i,j,N,tind,isalt)
#ifdef SALINITY_RESTORING
     &              -dSdt*dSdh*( t(i,j,N,tind,isalt)
     &              -(cff7*sssg(i,j,it7)+cff8*sssg(i,j,it8))
     &                                                          )
#endif
#endif              
!======================================================================           
#ifdef BULK_DIAGS
            shflx_net(i,j)= stflx(i,j,itemp)/cff            !<== diagnostics
            shflx_rad(i,j)= radsw(i,j)                      !<== diagnostics
# ifdef SALINITY
            swflx_emp(i,j)=(  stflx(i,j,isalt)
     &                    /   max(t(i,j,N,tind,isalt),1.e-14) )*8640000. !<== diagnostics
# endif
            sustr_blk(i,j)= sustr(i,j)*rho0  !<== diagnostics
            svstr_blk(i,j)= svstr(i,j)*rho0  !<== diagnostics
#endif
!======================================================================                     
!++ 8640 * rain[kg.m-2.s-1] = prec[cm.day-1] // prec[mm.s-1] = 0.01/86400 * prec[cm.day-1]                               
!====================================================================== 
          else
            stflx(i,j,isalt)=0.
            stflx(i,j,itemp)=0.
            srflx(i,j)=0.
            sustr(i,j) = 0.
            svstr(i,j) = 0.
#ifdef BULK_DIAGS
# ifdef SALINITY
            swflx_emp(i,j)=0.
# endif            
            shflx_rad(i,j)=0.
            shflx_net(i,j)=0.           
            sustr_blk(i,j)=0.
            svstr_blk(i,j)=0.      
#endif            
          endif
        enddo  !<-- i 
      enddo  !<-- j  
# undef I_EXT_RANGE
# undef J_EXT_RANGE     
      return
      end 

      function psimhu (xd)
      real psimhu, xd
      psimhu=1.570796326794896619D0 -2.D0*atan(xd)
     &         +log(0.125D0*(1.+xd*xd)*(1.+xd)**2)
      return
      end

      function psixhu (xd)
      real psixhu, xd
      psixhu = 2.*log( 0.5*(1.+xd*xd) )
      return
      end

#else
      subroutine bulk_flux_empty
      end      
#endif            
