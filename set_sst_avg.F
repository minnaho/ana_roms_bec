!
!
#include "cppdefs.h"
#if defined SOLVE3D && defined OUTPUT_SST_AVG
      subroutine set_sst_avg (tile)
      implicit none
      integer tile
# include "param.h"
#  include "work.h"
#  include "ncvars.h"
# include "compute_tile_bounds.h"

      call set_sst_avg_tile (Istr,Iend,Jstr,Jend)

      return
      end

      subroutine set_sst_avg_tile (Istr,Iend,Jstr,Jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Though this algorithm results in somewhat ackward controlling
! logic it has the advantage that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(ilc-1,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(ilc-1,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
# include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j, ilc
      real cff
      integer itrc,k
#  include "work.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "averages.h"
# include "mixing.h"
# include "forces.h"
!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
      if (ilc.gt.nts_sst_avg) then 
        if (mod(ilc-1,n_sst_avg).eq.1) then
           do j=JstrR,JendR
              do i=IstrR,IendR
                 sst_avg(i,j)=t(i,j,N,nstp,1)
              enddo
           enddo
        elseif (mod(ilc-1,n_sst_avg).gt.1) then
           do j=JstrR,JendR
              do i=IstrR,IendR
                 sst_avg(i,j)=sst_avg(i,j)+t(i,j,N,nstp,1)
              enddo
           enddo
        elseif (mod(ilc-1,n_sst_avg).eq.0) then
          cff=1./float(n_sst_avg)
          if (ZEROTH_TILE) time_sst_avg=time_sst_avg+float(n_sst_avg)*dt
          do j=JstrR,JendR
             do i=IstrR,IendR
                sst_avg(i,j)=cff*( sst_avg(i,j)
     &                           + t(i,j,N,nstp,1))
             enddo
          enddo
        endif
      endif
      return
      end
#else
      subroutine set_sst_avg_empty
      end
#endif /* OUTPUT_SST_AVG */

