!   $Log: set_phys_flux_avg.F,v $
!   Revision 1.3  2002/04/19 00:52:55  frenzel
!   Added code for the analysis of mixing and nudging fluxes.
!
!   Revision 1.2  2002/02/06 19:00:04  frenzel
!   Correction: Previously, averaging was not correct if the averaging was
!   performed every time step.
!
!   Revision 1.1  2002/01/16 19:58:53  frenzel
!   New routines related to flux analysis.
!
!
#include "cppdefs.h"
#if defined PHYS_FLUX_ANALYSIS && defined AVERAGES

      subroutine set_phys_flux_avg(tile)

      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_phys_flux_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_phys_flux_avg_tile(Istr,Iend,Jstr,Jend)
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
! of the averages arrays: when mod(iic-nts_phys_flux_avg,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(iic-nts_phys_flux_avg,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
!      USE scalars
!      USE grid
!      USE averages
!      USE ocean2d
!      USE physflux
!      USE ncscrum

#ifdef PFA_BY_TRACER
      use phys_flux
#endif
      implicit none
      integer Istr,Iend,Jstr,Jend, i, j, k, iflux, iflux_out
      real cff, cff1
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
!# include "averages.h"
# include "physflux.h"
# include "ncroms_physflux.h"
!
#include "compute_auxiliary_bounds.h"
!
! calculation of averaged fluxes will only be performed if n_phys_flux_avg
! is a positive number
      if (n_phys_flux_avg < 1) return

      if (iic.gt.nts_phys_flux_avg) then
         if (n_phys_flux_avg .eq. 1 .or. 
     &        mod(iic-nts_phys_flux_avg,n_phys_flux_avg).eq.1) then
            cff =1.0
            cff1=0.0
            if (ZEROTH_TILE) then
               time_phys_flux_avg=time
                  mpi_master_only write(*,'(6x,A,I11,2(X,A,I5))')
     &                 'SET_PHYS_FLUX_AVG -- Started averaging at iic=',
     &                 iic,' : nts_phys_flux_avg=',nts_phys_flux_avg,
     &                 'n_phys_flux_avg=',n_phys_flux_avg
               endif
            elseif (mod(iic-nts_phys_flux_avg,n_phys_flux_avg).gt.1)then
               cff =1.0
               cff1=1.0
               if (ZEROTH_TILE) time_phys_flux_avg=
     &              time_phys_flux_avg+time
            elseif (mod(iic-nts_phys_flux_avg,n_phys_flux_avg).eq.0)then
               cff=1./float(n_phys_flux_avg)
               cff1=1.0
               if (ZEROTH_TILE) then
                  time_phys_flux_avg=cff*(time_phys_flux_avg+time)
                  mpi_master_only write(*,'(6x,A,I11,2(X,A,I5))')
     &                'SET_PHYS_FLUX_AVG -- Finished averaging at iic=',
     &                 iic,' : nts_phys_flux_avg=',nts_phys_flux_avg,
     &                 'n_phys_flux_avg=',n_phys_flux_avg
               endif
            endif


            do j = jstrR, jendR
               do i = istrR, iendR
                  zeta_phys_flux_avg(i,j) =  
     &                 cff * ( cff1*zeta_phys_flux_avg(i,j) +
     &                 zeta(i,j,knew) )
               end do
            end do
            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 1, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        HorXAdvFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*
     &                       HorXAdvFlux_avg(i,j,k,iflux_out) +
     &                       HorXAdvFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do ! iflux
            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 1, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                         HorYAdvFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*
     &                       HorYAdvFlux_avg(i,j,k,iflux_out) +
     &                       HorYAdvFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do ! iflux
            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 0, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        VertAdvFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*
     &                       VertAdvFlux_avg(i,j,k,iflux_out) +
     &                       VertAdvFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do ! iflux

#ifdef HOR_DIFF_ANALYSIS
            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 1, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        HorXMixFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*HorXMixFlux_avg(i,j,k,iflux_out) +
     &                       HorXMixFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do

            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 1, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        HorYMixFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*HorYMixFlux_avg(i,j,k,iflux_out) +
     &                       HorYMixFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do
            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 1, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        VertMixFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*VertMixFlux_avg(i,j,k,iflux_out) +
     &                       VertMixFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do
#endif
#ifdef VERT_DIFF_ANALYSIS
            do iflux = 1, NT_PFA
# ifdef PFA_BY_TRACER
             if (wrt_pfa(iflux)) then
               iflux_out = pfa_out(iflux)
# else 
               iflux_out = iflux   
# endif /* PFA_BY_TRACER */ 
               do k = 0, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        VertDiffFlux_avg(i,j,k,iflux_out) = 
     &                       cff * ( cff1*
     &                       VertDiffFlux_avg(i,j,k,iflux_out)+
     &                       VertDiffFlux(i,j,k,iflux_out) )
                     end do
                  end do
               end do
# ifdef PFA_BY_TRACER
             endif
# endif /* PFA_BY_TRACER */ 
            end do ! iflux
#endif
# ifdef FULL_PHYS_FLUX_ANALYSIS
            do iflux = 1, NT_PFA
               do j = jstrR, jendR
                  do i = istrR, iendR 
                     TopFlux_avg(i,j,iflux) = 
     &                    cff * ( cff1*TopFlux_avg(i,j,iflux) + 
     &                    TopFlux(i,j,iflux) )
                  end do
               end do
            end do
#  ifdef BOTTOM_FLUX
            do iflux = 1, NT_PFA
               do j = jstrR, jendR
                  do i = istrR, iendR 
                     BottomFlux_avg(i,j,iflux) = 
     &                    cff * ( cff1*BottomFlux_avg(i,j,iflux) + 
     &                    BottomFlux(i,j,iflux) )
                  end do
               end do
            end do
#  endif /* BOTTOM_FLUX */
            do iflux = 1, NT_PFA
               do k = 1, N
                  do j = jstrR, jendR
                     do i = istrR, iendR 
                        NudgingFlux_avg(i,j,k,iflux) = 
     &                       cff * ( cff1*NudgingFlux_avg(i,j,k,iflux) +
     &                       NudgingFlux(i,j,k,iflux) )
                     end do
                  end do
               end do
            end do
# endif /* FULL_PHYS_FLUX_ANALYSIS */
      endif ! iic.gt.nts_phys_flux_avg
      return
      end
#else /* PHYS_FLUX_ANALYSIS && AVERAGES */
      subroutine set_phys_flux_avg_empty
      end
#endif /* PHYS_FLUX_ANALYSIS && AVERAGES */

