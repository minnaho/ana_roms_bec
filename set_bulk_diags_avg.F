#include "cppdefs.h"
#if defined BULK_DIAGS && defined AVERAGES

      subroutine set_bulk_diags_avg(tile)

      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_bulk_diags_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_bulk_diags_avg_tile(Istr,Iend,Jstr,Jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Thought this algorithm results in somewhat ackwad controlling
! logic it has the advantage that that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(iic-nts_bulk_diags_avg,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(iic-nts_bulk_diags_avg,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
!      USE scalars
!      USE grid
!      USE averages
!      USE ocean2d
!      USE physflux
!      USE ncscrum

      implicit none
      integer Istr,Iend,Jstr,Jend, i, j, k, iflux, iflux_out
      real cff, cff1
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "averages.h"
# include "forces.h"
!
#include "compute_auxiliary_bounds.h"
!
! calculation of averaged fluxes will only be performed if n_bulk_diags_avg
! is a positive number
      if (iic.gt.nts_bulk_diags_avg) then
        cff=1./dble(n_bulk_diags_avg) !<-- scaling

         if (n_bulk_diags_avg .eq. 1 .or. 
     &        mod(iic-nts_bulk_diags_avg,n_bulk_diags_avg).eq.1) then
           if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &               'set_bdiags_avg :: started averaging at iic=', iic,
     &                          'ntsavg=', ntsavg, 'navg=', navg
           endif
           cff1=0.   !--> initialize
        else
          cff1=1.   !--> keep adding
        endif

        if (ZEROTH_TILE) time_bulk_diags_avg=cff1*time_bulk_diags_avg + cff*time


            do j = jstrR, jendR
               do i = istrR, iendR
                  zeta_bulk_diags_avg(i,j) = 
     &                         cff1*zeta_bulk_diags_avg(i,j) + 
     &                         cff*zeta(i,j,knew)
                  sustr_blk_avg(i,j)=cff1*sustr_blk_avg(i,j) + 
     &                         cff*sustr_blk(i,j)
                  svstr_blk_avg(i,j)=cff1*svstr_blk_avg(i,j) +
     &                         cff*svstr_blk(i,j)
                  shflx_net_avg(i,j)=cff1*shflx_net_avg(i,j) + 
     &                         cff*shflx_net(i,j)
                  shflx_lat_avg(i,j)=cff1*shflx_lat_avg(i,j) + 
     &                         cff*shflx_lat(i,j)
                  shflx_sen_avg(i,j)=cff1*shflx_sen_avg(i,j) + 
     &                         cff*shflx_sen(i,j)
                  shflx_rad_avg(i,j)=cff1*shflx_rad_avg(i,j) + 
     &                         cff*shflx_rad(i,j)
                  swflx_emp_avg(i,j)=cff1*swflx_emp_avg(i,j) + 
     &                         cff*swflx_emp(i,j)
                  shflx_wwk_avg(i,j)=cff1*shflx_wwk_avg(i,j) +
     &                         cff*shflx_wwk(i,j)
                  surf_u_avg(i,j)=cff1*surf_u_avg(i,j) +
     &                         cff*surf_u(i,j)
                  surf_v_avg(i,j)=cff1*surf_v_avg(i,j) +
     &                         cff*surf_v(i,j)
            enddo
          enddo

      endif ! iic.gt.nts_bulk_diags_avg
      return
      end
#else /* BULK_DIAGS*/
      subroutine set_bulk_diags_avg_empty
      end
#endif /* BULK_DIAGS */

