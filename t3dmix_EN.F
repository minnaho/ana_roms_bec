#include "cppdefs.h"
 
      subroutine t3dmix (tile)
      implicit none
      integer tile, itrc
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      do itrc=1,NT
        call t3dmix_tile (istr,iend,jstr,jend, itrc, A3d(1,1),
     &                                    A2d(1, 1), A2d(1, 2),
     &                          A2d(1,3), A2d(1, 5), A2d(1, 7),
     &                          A2d(1,9), A2d(1,11), A2d(1,13))
      enddo
      return
      end
 
      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, LapT, FX,FE,
     &                                  FC,dTdr, dTdx,dTde,dRdx,dRde)
!
!  Compute isopycnic horizontal mixing terms for tracer equations.
!
      implicit none
#include "param.h"
      integer itrc, istr,iend,jstr,jend, i,j,k,k1,k2
      real LapT(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &       FX(PRIVATE_2D_SCRATCH_ARRAY),
     &       FE(PRIVATE_2D_SCRATCH_ARRAY),
     &       FC(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     dTdr(PRIVATE_2D_SCRATCH_ARRAY,2),   cff,
     &     dTdx(PRIVATE_2D_SCRATCH_ARRAY,2),   cff1,
     &     dTde(PRIVATE_2D_SCRATCH_ARRAY,2),   cff2,
     &     dRdx(PRIVATE_2D_SCRATCH_ARRAY,2),   cff3,
     &     dRde(PRIVATE_2D_SCRATCH_ARRAY,2),   cff4
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "climat.h"
#include "scalars.h"
!
#include "compute_auxiliary_bounds.h"
!
#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif
#ifdef SPLIT_EOS
# define rho rho1
#endif
 
 
c--#define rho z_r
c--#define MAX min
c--#define MIN max
 
#define EPS 0.5
c-#define MAX_SLOPE 1.E+4
c-#define MIN_STRAT  0.1
 
#ifdef TS_DIF2
!
! Compute and add in horizontal Laplacian     !!! WARNING: RECURSIVE
! diffusion along geopotential surfaces.            BLOCKING SEQUENCE
!--------------------------------------------------------------------
!                   !  Vertical    dTdx,dTde(:,:,k2) k+1   rho-points
      k2=1          !  placement     dTdr,FC(:,:,k2) k+1/2   W-points
      do k=0,N,+1   !  of fields:  dTdx,dTde(:,:,k1) k     pho-points
       k1=k2        !                dTdr,FC(:,:,k1) k-1/2   W-points
       k2=3-k1
 
        if (k.lt.N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dRdx(i,j,k2)=cff*(rho(i,j,k+1) -rho(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i-1,j,k+1,nrhs,itrc)
# ifdef CLIMAT_TS_MIXH
     &               -tclm(i,j,k+1,itrc)+tclm(i-1,j,k+1,itrc)
# endif
     &                                                      )
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dRde(i,j,k2)=cff*(rho(i,j,k+1) -rho(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i,j-1,k+1,nrhs,itrc)
# ifdef CLIMAT_TS_MIXH
     &               -tclm(i,j,k+1,itrc)+tclm(i,j-1,k+1,itrc)
# endif
     &                                                      )
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdr(i,j,k2)=0.0     ! Set no-flux bottom or
              FC(i,j,k2)=0.0       ! top boundary conditions.
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
#if defined MAX_SLOPE
              cff=MAX( rho(i,j,k)-rho(i,j,k+1),
     &                 1.E-14,
     &                 0.25*MAX_SLOPE*(z_r(i,j,k+1)-z_r(i,j,k))
     &                     *sqrt( dRdx(i,j,k2)**2+dRdx(i+1,j,k2)**2
     &                           +dRdx(i,j,k1)**2+dRdx(i+1,j,k1)**2
     &                           +dRde(i,j,k2)**2+dRde(i,j+1,k2)**2
     &                           +dRde(i,j,k1)**2+dRde(i,j+1,k1)**2)
     &                                                             )
              cff1=-1.0/cff
#elif defined MIN_STRAT
              cff1=-1.0/MAX( rho(i,j,k)-rho(i,j,k+1),
     &                       MIN_STRAT*(z_r(i,j,k+1)-z_r(i,j,k)))
#elif defined EPS
              cff1=-1.0/MAX( rho(i,j,k)-rho(i,j,k+1), EPS )
#endif
              dTdr(i,j,k2)=cff1*( t(i,j,k+1,nrhs,itrc)
     &                                            -t(i,j,k,nrhs,itrc)
# ifdef CLIMAT_TS_MIXH
     &                           -tclm(i,j,k+1,itrc)+tclm(i,j,k,itrc)
# endif
     &                                                              )
              FC(i,j,k2)=cff1*(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif
 
        if (k.gt.0) then
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             MAX(dRdx(i,j,k1),0.)*(dTdr(i-1,j,k1)+dTdr(i,j,k2))
     &            +MIN(dRdx(i,j,k1),0.)*(dTdr(i-1,j,k2)+dTdr(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=0.25*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
     &                        *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                                                  dTde(i,j,k1)
     &      -0.5*( MAX(dRde(i,j,k1),0.)*(dTdr(i,j-1,k1)+dTdr(i,j,k2))
     &            +MIN(dRde(i,j,k1),0.)*(dTdr(i,j-1,k2)+dTdr(i,j,k1))
     &                                                            ))
            enddo
          enddo
          if (k.lt.N) then
            do j=jstr,jend
              do i=istr,iend
c                cff1=MAX(dRdx(i,j,k1)+dRdx(i+1,j,k2),0.)
c                cff2=MIN(dRdx(i,j,k2)+dRdx(i+1,j,k1),0.)
c
c                cff3=MAX(dRde(i,j,k1)+dRde(i,j+1,k2),0.)
c                cff4=MIN(dRde(i,j,k2)+dRde(i,j+1,k1),0.)
c
c                FC(i,j,k2)=0.25*FC(i,j,k2)*diff2(i,j,itrc)*(
c     &                             dTdr(i,j,k2)*( (cff1+cff2)**2
c     &                                           +(cff3+cff4)**2)
c
c     &                           -cff1*(dTdx(i,j,k1)+dTdx(i+1,j,k2))
c     &                           -cff2*(dTdx(i,j,k2)+dTdx(i+1,j,k1))
c
c     &                           -cff3*(dTde(i,j,k1)+dTde(i,j+1,k2))
c     &                           -cff4*(dTde(i,j,k2)+dTde(i,j+1,k1))
c     &                                                             )
 
                cff1=MAX(dRdx(i  ,j,k1),0.)
                cff2=MAX(dRdx(i+1,j,k2),0.)
 
                cff3=MIN(dRdx(i  ,j,k2),0.)
                cff4=MIN(dRdx(i+1,j,k1),0.)
 
                cff=    cff1*(cff1*dTdr(i,j,k2)-dTdx(i  ,j,k1))
     &                 +cff2*(cff2*dTdr(i,j,k2)-dTdx(i+1,j,k2))
 
     &                 +cff3*(cff3*dTdr(i,j,k2)-dTdx(i  ,j,k2))
     &                 +cff4*(cff4*dTdr(i,j,k2)-dTdx(i+1,j,k1))
 
                cff1=MAX(dRde(i,j  ,k1),0.)
                cff2=MAX(dRde(i,j+1,k2),0.)
 
                cff3=MIN(dRde(i,j  ,k2),0.)
                cff4=MIN(dRde(i,j+1,k1),0.)
 
                cff=cff+cff1*(cff1*dTdr(i,j,k2)-dTde(i,j  ,k1))
     &                 +cff2*(cff2*dTdr(i,j,k2)-dTde(i,j+1,k2))
 
     &                 +cff3*(cff3*dTdr(i,j,k2)-dTde(i,j  ,k2))
     &                 +cff4*(cff4*dTdr(i,j,k2)-dTde(i,j+1,k1))
 
                FC(i,j,k2)=0.5*cff*diff2(i,j,itrc)*FC(i,j,k2)
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) + dt*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
        endif  ! <-- k.gt.0
      enddo   ! --> k
#endif /* T_DIF2 */
 
 
#ifdef TS_DIF4
# ifdef EW_PERIODIC
#  define I_EXT_RANGE istr-1,iend+1
# else
#  define I_EXT_RANGE max(istr-1,1),min(iend+1,Lm)
# endif
# ifdef NS_PERIODIC
#  define J_EXT_RANGE jstr-1,jend+1
# else
#  define J_EXT_RANGE max(jstr-1,1),min(jend+1,Mm)
# endif
!
!  Add in horizontal biharmonic diffusion along rotated surfaces.
!  The biharmonic operator is computed by applying the rotated
!  Laplacian operator twice.
!
!                                            !!! WARNING: RECURSIVE
!  The FIRST rotated Laplacian operator:           BLOCKING SEQUENCE
!--------------------------------------------------------------------
!                   !  Vertical    dTdx,dTde(:,:,k2) k+1   rho-points
      k2=1          !  placement     dTdr,FC(:,:,k2) k+1/2   W-points
      do k=0,N,+1   !  of fields:  dTdx,dTde(:,:,k1) k     pho-points
       k1=k2        !                dTdr,FC(:,:,k1) k-1/2   W-points
       k2=3-k1
        if (k.lt.N) then
          do j=J_EXT_RANGE
            do i=I_EXT_RANGE+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dRdx(i,j,k2)=cff*(rho(i,j,k+1)-rho(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-
     &                               t(i-1,j,k+1,nrhs,itrc))
            enddo
          enddo
          do j=J_EXT_RANGE+1
            do i=I_EXT_RANGE
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dRde(i,j,k2)=cff*(rho(i,j,k+1)-rho(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-
     &                               t(i,j-1,k+1,nrhs,itrc))
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=-1+J_EXT_RANGE+1
            do i=-1+I_EXT_RANGE+1
              FC(i,j,k2)=0.0
              dTdr(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=-1+J_EXT_RANGE+1
            do i=-1+I_EXT_RANGE+1
#if defined MAX_SLOPE
              cff=MAX( rho(i,j,k)-rho(i,j,k+1),
     &                 1.E-14,
     &                 0.25*MAX_SLOPE*(z_r(i,j,k+1)-z_r(i,j,k))
     &                     *sqrt( dRdx(i,j,k2)**2+dRdx(i+1,j,k2)**2
     &                           +dRdx(i,j,k1)**2+dRdx(i+1,j,k1)**2
     &                           +dRde(i,j,k2)**2+dRde(i,j+1,k2)**2
     &                           +dRde(i,j,k1)**2+dRde(i,j+1,k1)**2)
     &                                                             )
              cff1=-1.0/cff
#elif defined MIN_STRAT
              cff1=-1.0/MAX( rho(i,j,k)-rho(i,j,k+1),
     &                       MIN_STRAT*(z_r(i,j,k+1)-z_r(i,j,k)))
#elif defined EPS
              cff1=-1.0/MAX( rho(i,j,k)-rho(i,j,k+1), EPS )
#endif
 
              FC(i,j,k2)=cff1*(z_r(i,j,k+1)-z_r(i,j,k))
              dTdr(i,j,k2)=cff1*(t(i,j,k+1,nrhs,itrc)
     &                             -t(i,j,k,nrhs,itrc))
            enddo
          enddo
        endif
 
        if (k.gt.0) then
          do j=J_EXT_RANGE
            do i=I_EXT_RANGE+1
              FX(i,j)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             MAX(dRdx(i,j,k1),0.)*(dTdr(i-1,j,k1)+dTdr(i,j,k2))
     &            +MIN(dRdx(i,j,k1),0.)*(dTdr(i-1,j,k2)+dTdr(i,j,k1))
     &                                                             ))
            enddo
          enddo
 
# ifdef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=J_EXT_RANGE
              .....
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=J_EXT_RANGE
              .....
            enddo
          endif
# endif
 
          do j=J_EXT_RANGE+1
            do i=I_EXT_RANGE
              FE(i,j)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             MAX(dRde(i,j,k1),0.)*(dTdr(i,j-1,k1)+dTdr(i,j,k2))
     &            +MIN(dRde(i,j,k1),0.)*(dTdr(i,j-1,k2)+dTdr(i,j,k1))
     &                                                             ))
            enddo
          enddo
# ifdef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=I_EXT_RANGE
              .....
            enddo
          endif
          if (NORTHERN_EDGE) then
            do i=I_EXT_RANGE
               ....
            enddo
          endif
# endif
 
          if (k.lt.N) then
            do j=J_EXT_RANGE
              do i=I_EXT_RANGE
c                cff1=MAX(dRdx(i,j,k1)+dRdx(i+1,j,k2),0.)
c                cff2=MIN(dRdx(i,j,k2)+dRdx(i+1,j,k1),0.)
c
c                cff3=MAX(dRde(i,j,k1)+dRde(i,j+1,k2),0.)
c                cff4=MIN(dRde(i,j,k2)+dRde(i,j+1,k1),0.)
c
c                FC(i,j,k2)=FC(i,j,k2)*0.25*(
c     &                  dTdr(i,j,k2)*((cff1+cff2)**2+(cff3+cff4)**2)
c
c     &                           -cff1*(dTdx(i,j,k1)+dTdx(i+1,j,k2))
c     &                           -cff2*(dTdx(i,j,k2)+dTdx(i+1,j,k1))
c
c     &                           -cff3*(dTde(i,j,k1)+dTde(i,j+1,k2))
c     &                           -cff4*(dTde(i,j,k2)+dTde(i,j+1,k1))
c     &                                                             )
 
 
                cff1=MAX(dRdx(i  ,j,k1),0.)
                cff2=MAX(dRdx(i+1,j,k2),0.)
 
                cff3=MIN(dRdx(i  ,j,k2),0.)
                cff4=MIN(dRdx(i+1,j,k1),0.)
 
                cff=    cff1*(cff1*dTdr(i,j,k2)-dTdx(i  ,j,k1))
     &                 +cff2*(cff2*dTdr(i,j,k2)-dTdx(i+1,j,k2))
 
     &                 +cff3*(cff3*dTdr(i,j,k2)-dTdx(i  ,j,k2))
     &                 +cff4*(cff4*dTdr(i,j,k2)-dTdx(i+1,j,k1))
 
                cff1=MAX(dRde(i,j  ,k1),0.)
                cff2=MAX(dRde(i,j+1,k2),0.)
 
                cff3=MIN(dRde(i,j  ,k2),0.)
                cff4=MIN(dRde(i,j+1,k1),0.)
 
                cff=cff+cff1*(cff1*dTdr(i,j,k2)-dTde(i,j  ,k1))
     &                 +cff2*(cff2*dTdr(i,j,k2)-dTde(i,j+1,k2))
 
     &                 +cff3*(cff3*dTdr(i,j,k2)-dTde(i,j  ,k2))
     &                 +cff4*(cff4*dTdr(i,j,k2)-dTde(i,j+1,k1))
 
                FC(i,j,k2)=0.5*cff*FC(i,j,k2)
              enddo
            enddo
          endif
          do j=J_EXT_RANGE
            do i=I_EXT_RANGE
              LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                       +FE(i,j+1)-FE(i,j))
     &                     +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
        endif   ! <-- k.gt.0
      enddo   ! --> k
!
!  Apply boundary conditions to the Laplacian
!  (for the cases other than periodic: closed or gradient)
!
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do k=1,N
            do j=J_EXT_RANGE
#  ifdef WESTERN_WALL
              LapT(istr-1,j,k)=0.
#  else
              LapT(istr-1,j,k)=LapT(istr,j,k)
#  endif
            enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do k=1,N
            do j=J_EXT_RANGE
#  ifdef EASTERN_WALL
              LapT(iend+1,j,k)=0.
#  else
              LapT(iend+1,j,k)=LapT(iend,j,k)
#  endif
            enddo
          enddo
        endif
# endif /* !EW_PERIODIC */
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do k=1,N
            do i=I_EXT_RANGE
#  ifdef SOUTHERN_WALL
              LapT(i,jstr-1,k)=0.
#  else
              LapT(i,jstr-1,k)=LapT(i,jstr,k)
#  endif
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do k=1,N
            do i=I_EXT_RANGE
#  ifdef NORTHERN_WALL
              LapT(i,jend+1,k)=0.
#  else
              LapT(i,jend+1,k)=LapT(i,jend,k)
#  endif
            enddo
          enddo
        endif
# endif /* !NS_PERIODIC */
# undef I_EXT_RANGE
# undef J_EXT_RANGE
!                                          !!! WARNING: RECURSIVE
!  The SECOND rotated Laplacian operator:        BLOCKING SEQUENCE
!--------------------------------------------------------------------
!                   !  Vertical    dTdx,dTde(:,:,k2) k+1   rho-points
      k2=1          !  placement     dTdr,FC(:,:,k2) k+1/2   W-points
      do k=0,N,+1   !  of fields:  dTdx,dTde(:,:,k1) k     pho-points
       k1=k2        !                dTdr,FC(:,:,k1) k-1/2   W-points
       k2=3-k1
 
        if (k.lt.N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dRdx(i,j,k2)=cff*(rho(i,j,k+1) -rho(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i-1,j,k+1))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dRde(i,j,k2)=cff*(rho(i,j,k+1) -rho(i,j-1,k+1))
              dTde(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i,j-1,k+1))
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              FC(i,j,k2)=0.0
              dTdr(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
#if defined MAX_SLOPE
              cff=MAX( rho(i,j,k)-rho(i,j,k+1),
     &                 1.E-14,
     &                 0.25*MAX_SLOPE*(z_r(i,j,k+1)-z_r(i,j,k))
     &                     *sqrt( dRdx(i,j,k2)**2+dRdx(i+1,j,k2)**2
     &                           +dRdx(i,j,k1)**2+dRdx(i+1,j,k1)**2
     &                           +dRde(i,j,k2)**2+dRde(i,j+1,k2)**2
     &                           +dRde(i,j,k1)**2+dRde(i,j+1,k1)**2)
     &                                                             )
              cff1=-1.0/cff
#elif defined MIN_STRAT
              cff1=-1.0/MAX( rho(i,j,k)-rho(i,j,k+1),
     &                       MIN_STRAT*(z_r(i,j,k+1)-z_r(i,j,k)))
#elif defined EPS
              cff1=-1.0/MAX( rho(i,j,k)-rho(i,j,k+1), EPS )
#endif
 
              FC(i,j,k2)=cff1*(z_r(i,j,k+1)-z_r(i,j,k))
              dTdr(i,j,k2)=cff1*(LapT(i,j,k+1)-LapT(i,j,k))
            enddo
          enddo
        endif
 
        if (k.gt.0) then
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=-0.25*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             MAX(dRdx(i,j,k1),0.)*(dTdr(i-1,j,k1)+dTdr(i,j,k2))
     &            +MIN(dRdx(i,j,k1),0.)*(dTdr(i-1,j,k2)+dTdr(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=-0.25*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                     *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             MAX(dRde(i,j,k1),0.)*(dTdr(i,j-1,k1)+dTdr(i,j,k2))
     &            +MIN(dRde(i,j,k1),0.)*(dTdr(i,j-1,k2)+dTdr(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k.lt.N) then
            do j=jstr,jend
              do i=istr,iend
c                cff1=MAX(dRdx(i,j,k1)+dRdx(i+1,j,k2),0.)
c                cff2=MIN(dRdx(i,j,k2)+dRdx(i+1,j,k1),0.)
c
c                cff3=MAX(dRde(i,j,k1)+dRde(i,j+1,k2),0.)
c                cff4=MIN(dRde(i,j,k2)+dRde(i,j+1,k1),0.)
c
c                FC(i,j,k2)=-0.25*FC(i,j,k2)*diff4(i,j,itrc)*(
c     &                             dTdr(i,j,k2)*( (cff1+cff2)**2
c     &                                           +(cff3+cff4)**2)
c
c     &                           -cff1*(dTdx(i,j,k1)+dTdx(i+1,j,k2))
c     &                           -cff2*(dTdx(i,j,k2)+dTdx(i+1,j,k1))
c
c     &                           -cff3*(dTde(i,j,k1)+dTde(i,j+1,k2))
c     &                           -cff4*(dTde(i,j,k2)+dTde(i,j+1,k1))
c     &                                                             )
 
                cff1=MAX(dRdx(i  ,j,k1),0.)
                cff2=MAX(dRdx(i+1,j,k2),0.)
 
                cff3=MIN(dRdx(i  ,j,k2),0.)
                cff4=MIN(dRdx(i+1,j,k1),0.)
 
                cff=    cff1*(cff1*dTdr(i,j,k2)-dTdx(i  ,j,k1))
     &                 +cff2*(cff2*dTdr(i,j,k2)-dTdx(i+1,j,k2))
 
     &                 +cff3*(cff3*dTdr(i,j,k2)-dTdx(i  ,j,k2))
     &                 +cff4*(cff4*dTdr(i,j,k2)-dTdx(i+1,j,k1))
 
                cff1=MAX(dRde(i,j  ,k1),0.)
                cff2=MAX(dRde(i,j+1,k2),0.)
 
                cff3=MIN(dRde(i,j  ,k2),0.)
                cff4=MIN(dRde(i,j+1,k1),0.)
 
                cff=cff+cff1*(cff1*dTdr(i,j,k2)-dTde(i,j  ,k1))
     &                 +cff2*(cff2*dTdr(i,j,k2)-dTde(i,j+1,k2))
 
     &                 +cff3*(cff3*dTdr(i,j,k2)-dTde(i,j  ,k2))
     &                 +cff4*(cff4*dTdr(i,j,k2)-dTde(i,j+1,k1))
 
                FC(i,j,k2)=-0.5*cff*diff4(i,j,itrc)*FC(i,j,k2)
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) + dt*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
        endif  ! <-- k.gt.0
      enddo   ! --> k
#  endif /* T_DIF4 */

# ifdef EXCHANGE
        call exchange_tile (istr,iend,jstr,jend,
     &                      t(START_2D_ARRAY,1,nnew,itrc),N)
# endif
      return
      end
 
