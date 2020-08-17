#include "cppdefs.h"

#ifdef PHYS_FLUX_ANALYSIS 
! this routine MUST be called after allocate_physflux_arrays if
! PFA_BY_TRACER is defined!
      subroutine init_arrays_physflux (tile)
      implicit none
      integer tile, i,j
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"

#ifdef NOTUSED_SOFAR
#ifdef SOLVE3D
      do j=1,6           ! Initialize (first touch) private
        do i=1,N3d       ! scratch arrays in parallel by each
          A3d(i,j)=0.    ! thread.
        enddo
      enddo
      do i=1,N2d
        iA2d(i)=0
      enddo
#endif
      do j=1,32
        do i=1,N2d
          A2d(i,j)=0.
        enddo
      enddo
#endif
      call init_arrays_physflux_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine init_arrays_physflux_tile (istr,iend,jstr,jend)

! Initialize (first touch) globally accessable arrays to zeros. 
! Because of the "first touch"
! default data distribution policy, this operation actually performs
! distribution of the shared arrays accross the cluster, unless
! another distribution policy is specified to override the default.

#ifdef PFA_BY_TRACER
      use phys_flux
#endif
      implicit none
      integer istr,iend,jstr,jend, i,j,k,itrc
      real, parameter :: init=0.    !!!!  0xFFFA5A5A ==> NaN
#define ALL_DATA
#include "param.h"
#include "scalars.h"
#ifndef PFA_BY_TRACER
# include "physflux.h"
#endif
#ifdef PRINT_TILE_RANGES
# ifdef MPI
#  include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank, ierr
# endif
#endif

#include "compute_extended_bounds.h"

#ifdef PRINT_TILE_RANGES 
# ifdef MPI
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, ocean_grid_comm, status, ierr)
      endif
      i=mynode
# else
      i=proc(2)
# endif
!      write(*,'(I4/2(6x,A6,I3,3x,A6,I3))') i, 'istr =',istr,
!     &        'iend =',iend,   'jstr =',jstr, 'jend =',jend
!      write(*,'(4x,2(6x,A6,I3,3x,A6,I3)/)')   'istrR=',istrR,
!     &        'iendR=',iendR, 'jstrR=',jstrR, 'jendR=',jendR
# ifdef MPI
      if (mynode .lt. NNODES) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, ocean_grid_comm,  ierr)
      endif
# endif
#endif

! set all fluxes to 0
      HorXAdvFlux = 0.0
      HorYAdvFlux = 0.0
      VertAdvFlux = 0.0
#  ifdef HOR_DIFF_ANALYSIS
      HorXMixFlux = 0.0
      HorYMixFlux = 0.0
      VertMixFlux = 0.0
#  endif
#  ifdef VERT_DIFF_ANALYSIS
      VertDiffFlux = 0.0
#  endif
# ifdef FULL_PHYS_FLUX_ANALYSIS
error not yet coded
# endif

# ifdef AVERAGES
      HorXAdvFlux_avg = 0.0
      HorYAdvFlux_avg = 0.0
      VertAdvFlux_avg = 0.0
#  ifdef HOR_DIFF_ANALYSIS
      HorXMixFlux_avg = 0.0
      HorYMixFlux_avg = 0.0
      VertMixFlux_avg = 0.0
#  endif
#  ifdef VERT_DIFF_ANALYSIS
      VertDiffFlux_avg = 0.0
#  endif
# endif /* AVERAGES */

      return
      end
#else /* PHYS_FLUX_ANALYSIS */
      subroutine init_arrays_physflux_empty
      end
#endif /* PHYS_FLUX_ANALYSIS */ 
