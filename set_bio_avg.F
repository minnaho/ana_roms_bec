#include "cppdefs.h"
#if defined AVERAGES && defined OUTPUT_AVG_PHYS_BIO
      subroutine set_bio_avg (tile)
      implicit none
      integer tile
# include "param.h"
# ifdef SOLVE3D
#  include "ncvars.h"
#  include "work.h"
#  include "private_scratch.h"
# endif
# include "compute_tile_bounds.h"

      call set_bio_avg_tile (Istr,Iend,Jstr,Jend)

      return
      end
 
      subroutine set_bio_avg_tile (istr,iend,jstr,jend)

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
! of the average arrays: when  mod(iic-1,n_bio_avg).eq.1  the target arrays
! are set to the first contribution; (2) accumulation of averaged data,
! when mod(iic-1,n_bio_avg).gt.1; and (3) adding the last contribution and
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
# include "mixing.h"
# include "averages.h"
# ifdef SOLVE3D
#  if defined SEDIMENT_BIOLOGY
#   include "bgcflux.h"
#  endif
#  if defined BIOLOGY_BEC
#   include "ecosys_bec.h"
#   include "bec_avg.h"
#  endif
!
# include "compute_auxiliary_bounds.h"

      if (iic.gt.nts_bio_avg) then
        cff=1./dble(n_bio_avg) !<-- scaling
        if (n_bio_avg .eq. 1 .or. mod(iic-nts_bio_avg,n_bio_avg).eq.1) 
     &       then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &               'set_avg :: started averaging at iic=', iic,
     &                          'nts_bio_avg=', nts_bio_avg, 
     &            'n_bio_avg=', n_bio_avg
          endif
          cff1=0.   !--> initialize
        else
          cff1=1.   !--> keep adding
        endif

        if (ZEROTH_TILE) time_bio_avg=cff1*time_bio_avg + cff*time

        if (wrtavg(indxZ)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              zeta_bio_avg(i,j)=cff1*zeta_bio_avg(i,j) + 
     &              cff*zeta(i,j,knew)
            enddo
          enddo
        endif
        do itrc=ntrc_salt+1,NT
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
#  if defined SEDIMENT_BIOLOGY
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
#  endif /* SEDIMENT_BIOLOGY */
#  if (defined BIOLOGY_NPZDOC && defined CARBON) || defined BIOLOGY_BEC
          do j=jstrR,jendR
             do i=istrR,iendR
                pH_avg(i,j) = cff1*pH_avg(i,j) + cff*pH_hist(i,j)
                pCO2_avg(i,j) = cff1*pCO2_avg(i,j) +cff*pCO2sw(i,j)
                pCO2air_avg(i,j) =cff1*pCO2air_avg(i,j)+cff*pCO2air(i,j)
             end do
          end do
#  endif /* BIOLOGY_NPZDOC ... */
#  if defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
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
# endif /* SOLVE3D */

        if (mod(iic-nts_bio_avg,n_bio_avg).eq.0) then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &              'set_bio_avg :: finished averaging at iic=', iic,
     &                         'nts_bio_avg=', nts_bio_avg,  
     &            'n_bio_avg=', n_bio_avg
          endif
        endif
      endif       !<-- iic.gt.nts_bio_avg
      return
      end
#else
      subroutine set_avg_empty
      end
#endif /* AVERAGES */
