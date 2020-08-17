#include "cppdefs.h"
#ifdef AVERAGES
 
      subroutine set_avg (tile)
      implicit none
      integer tile
# include "param.h"
# ifdef SOLVE3D
#  include "ncvars.h"
#  include "work.h"
#  include "private_scratch.h"
# endif
# include "compute_tile_bounds.h"
# ifdef SOLVE3D
      if (wrtavg(indxW)) then
        call wvlcty_tile (istr,iend,jstr,jend, work, A2d(1,1),
     &                                   A2d(1,2), A2d(1,3))
      endif
# endif
      call set_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_avg_tile (istr,iend,jstr,jend)

! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of the syncronization issues, a delayed-mode procedure is
! used for averaging.  This implies that all fields to be averaged are
! sampled during the next time step, rather than at the end of the time
! step when they are computed.
!
! Although this algorithm results in somewhat awkward controlling logic
! it has the advantage that that all fields to be sampled correspond to
! exactly the same time, which is time step "n".  Particularly, this is
! done this way because vertical velocity corresponding to the newly
! computed horizontal velocities becomes available only during the
! following time step. The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the average arrays: when  mod(iic-1,navg).eq.1  the target arrays
! are set to the first contribution; (2) accumulation of averaged data,
! when mod(iic-1,navg).gt.1; and (3) adding the last contribution and
! scaling.

      implicit none
      integer istr,iend,jstr,jend, i,j
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "mixing.h"
# include "averages.h"
# ifdef SOLVE3D
#  if defined SEDIMENT_BIOLOGY && !defined OUTPUT_AVG_PHYS_BIO
#   include "bgcflux.h"
#  endif
#  if defined BIOLOGY_BEC && !defined OUTPUT_AVG_PHYS_BIO
#   include "ecosys_bec.h"
#   include "bec_avg.h"
#  endif
# ifdef BULK_DIAGS
#  include "forces.h"
# endif
# include "work.h"
# endif

# include "compute_auxiliary_bounds.h"

      if (iic.gt.ntsavg) then
        cff=1./dble(navg) !<-- scaling
        if (navg .eq. 1 .or. mod(iic-ntsavg,navg).eq.1) then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &               'set_avg :: started averaging at iic=', iic,
     &                          'ntsavg=', ntsavg, 'navg=', navg
          endif
          cff1=0.   !--> initialize
        else
          cff1=1.   !--> keep adding
        endif

        if (ZEROTH_TILE) time_avg=cff1*time_avg + cff*time

        if (wrtavg(indxZ)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              zeta_avg(i,j)=cff1*zeta_avg(i,j) + cff*zeta(i,j,knew)
            enddo
          enddo
        endif
        if (wrtavg(indxUb)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              ubar_avg(i,j)=cff1*ubar_avg(i,j) + cff*ubar(i,j,knew)
            enddo
          enddo
        endif
        if (wrtavg(indxVb)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              vbar_avg(i,j)=cff1*vbar_avg(i,j) +cff*vbar(i,j,knew)
            enddo
          enddo
        endif
# ifdef SOLVE3D
        if (wrtavg(indxU)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                u_avg(i,j,k)=cff1*u_avg(i,j,k) +cff*u(i,j,k,nstp)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxV)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                v_avg(i,j,k)=cff1*v_avg(i,j,k) + cff*v(i,j,k,nstp)
              enddo
            enddo
          enddo
        endif
# if defined OUTPUT_AVG_PHYS_BIO
          do itrc=1,1+ntrc_salt
# else
        do itrc=1,NT
# endif
          if (wrtavg(indxT+itrc-1)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  t_avg(i,j,k,itrc)=cff1*t_avg(i,j,k,itrc)
     &                              +cff*t(i,j,k,nstp,itrc)
                enddo
              enddo
            enddo
          endif
        enddo
#  if defined SEDIMENT_BIOLOGY && !defined OUTPUT_AVG_PHYS_BIO
          do itrc=1,NT_sed
             if (wrtavg(indxSedFirst+itrc-1)) then
                do j=jstrR,jendR
                   do i=istrR,iendR
                      t_sed_avg(i,j,itrc)=cff1*t_sed_avg(i,j,itrc)
     &                              +cff*t_sed(i,j,itrc)
                   enddo
                enddo
             endif
          enddo
#  endif /* SEDIMENT_BIOLOGY && !OUTPUT_AVG_PHYS_BIO */
#  if ((defined BIOLOGY_NPZDOC && defined CARBON) || defined BIOLOGY_BEC) && ! defined OUTPUT_AVG_PHYS_BIO
          do j=jstrR,jendR
             do i=istrR,iendR
                pH_avg(i,j) = cff1*pH_avg(i,j) + cff*pH_hist(i,j)
                pCO2_avg(i,j) = cff1*pCO2_avg(i,j) +cff*pCO2sw(i,j)
                pCO2air_avg(i,j) =cff1*pCO2air_avg(i,j)+cff*pCO2air(i,j)
             end do
          end do
#  endif /* BIOLOGY_NPZDOC ... */
#  if (defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC) && ! defined OUTPUT_AVG_PHYS_BIO
          do j=jstrR,jendR
             do i=istrR,iendR
                PARinc_avg(i,j) = cff1*PARinc_avg(i,j) +cff*PARinc(i,j)
             end do
          end do
          do k=1,N
            do j=jstrR,jendR
               do i=istrR,iendR
                  PAR_avg(i,j,k) = cff1*PAR_avg(i,j,k) +cff*PAR(i,j,k)
               end do
            end do
          end do
#   ifdef OUTPUT_SPECDOM
          do j=jstrR,jendR
             do i=istrR,iendR
                dom_sp_sfc_avg(i,j) = cff1*dom_sp_sfc_avg(i,j) +
     &               cff*dom_sp_sfc(i,j)
             end do
          end do
          do j=jstrR,jendR
             do i=istrR,iendR
                dom_diat_sfc_avg(i,j) = cff1*dom_diat_sfc_avg(i,j) +
     &               cff*dom_diat_sfc(i,j)
             end do
          end do
          do j=jstrR,jendR
             do i=istrR,iendR
                dom_diaz_sfc_avg(i,j) = cff1*dom_diaz_sfc_avg(i,j) +
     &               cff*dom_diaz_sfc(i,j)
             end do
          end do
          do j=jstrR,jendR
             do i=istrR,iendR
                dom_sp_int_avg(i,j) = cff1*dom_sp_int_avg(i,j) +
     &               cff*dom_sp_int(i,j)
             end do
          end do
          do j=jstrR,jendR
             do i=istrR,iendR
                dom_diat_int_avg(i,j) = cff1*dom_diat_int_avg(i,j) +
     &               cff*dom_diat_int(i,j)
             end do
          end do
          do j=jstrR,jendR
             do i=istrR,iendR
                dom_diaz_int_avg(i,j) = cff1*dom_diaz_int_avg(i,j) +
     &               cff*dom_diaz_int(i,j)
             end do
          end do
#   endif /* OUTPUT_SPECDOM */
#  endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC */
        if (wrtavg(indxR)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
#  ifdef SPLIT_EOS
                rho_avg(i,j,k)=cff1*rho_avg(i,j,k) +cff*( rho1(i,j,k)
     &                                        -qp1(i,j,k)*z_r(i,j,k) )
#  else
                rho_avg(i,j,k)=cff1*rho_avg(i,j,k) +cff*rho(i,j,k)
#  endif
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxO)) then
          if (mod(iic-ntsavg,navg).eq.1) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=We(i,j,k)+Wi(i,j,k)
                enddo
              enddo
            enddo
          elseif (mod(iic-ntsavg,navg).gt.1) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=w_avg(i,j,k) +We(i,j,k)+Wi(i,j,k)
                enddo
              enddo
            enddo
          elseif (mod(iic-ntsavg,navg).eq.0) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=cff*( w_avg(i,j,k)+We(i,j,k)+Wi(i,j,k)
     &                                              )*pm(i,j)*pn(i,j)
                enddo
              enddo
            enddo
          endif
        endif
        if (wrtavg(indxW)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                wvl_avg(i,j,k)=cff1*wvl_avg(i,j,k) +cff*work(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkv)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akv_avg(i,j,k)=cff1*akv_avg(i,j,k) +cff*Akv(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkt)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akt_avg(i,j,k)=cff1*akt_avg(i,j,k)+cff*Akt(i,j,k,itemp)
              enddo
            enddo
          enddo
        endif
#  ifdef SALINITY
        if (wrtavg(indxAks)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                aks_avg(i,j,k)=cff1*aks_avg(i,j,k)+cff*Akt(i,j,k,isalt)
              enddo
            enddo
          enddo
        endif
#  endif
#  ifdef LMD_KPP
        if (wrtavg(indxHbls)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              hbl_avg(i,j)=cff1*hbl_avg(i,j) +cff*hbls(i,j,nstp)
            enddo
          enddo
        endif
#  endif
#  ifdef LMD_BKPP
        if (wrtavg(indxHbbl)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              hbbl_avg(i,j)=cff1*hbbl_avg(i,j) +cff*hbbls(i,j,nstp)
            enddo
          enddo
        endif
#  endif
!#  ifdef BULK_DIAGS
!          do j=jstrR,jendR
!            do i=istrR,iendR
!              sustr_blk_avg(i,j)=cff1*sustr_blk_avg(i,j) +cff*sustr_blk(i,j)
!              svstr_blk_avg(i,j)=cff1*svstr_blk_avg(i,j) +cff*svstr_blk(i,j)
!              shflx_net_avg(i,j)=cff1*shflx_net_avg(i,j) +cff*shflx_net(i,j)
!              shflx_lat_avg(i,j)=cff1*shflx_lat_avg(i,j) +cff*shflx_lat(i,j)
!              shflx_sen_avg(i,j)=cff1*shflx_sen_avg(i,j) +cff*shflx_sen(i,j)
!              shflx_rad_avg(i,j)=cff1*shflx_rad_avg(i,j) +cff*shflx_rad(i,j)
!              swflx_emp_avg(i,j)=cff1*swflx_emp_avg(i,j) +cff*swflx_emp(i,j)
!              shflx_wwk_avg(i,j)=cff1*shflx_wwk_avg(i,j) +cff*shflx_wwk(i,j)
!            enddo
!          enddo
!#  endif

# endif /* SOLVE3D */

        if (mod(iic-ntsavg,navg).eq.0) then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &              'set_avg :: finished averaging at iic=', iic,
     &                         'ntsavg=', ntsavg,  'navg=', navg
          endif
        endif
      endif       !<-- iic.gt.ntsavg
      return
      end
#else
      subroutine set_avg_empty
      end
#endif /* AVERAGES */
