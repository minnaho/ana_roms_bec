#include "cppdefs.h"
#if defined SOLVE3D && defined PHYS_FLUX_ANALYSIS && defined PFA_BY_TRACER
      module phys_flux
! allocatable arrays can not be in common blocks
! use this module so that the physical fluxes can be used 
! by different subroutines
      implicit none
      save

      real, allocatable, dimension(:,:,:,:) :: HorXAdvFlux
      real, allocatable, dimension(:,:,:,:) :: HorYAdvFlux
      real, allocatable, dimension(:,:,:,:) :: VertAdvFlux

#  ifdef HOR_DIFF_ANALYSIS
      real, allocatable, dimension(:,:,:,:) :: HorXMixFlux
      real, allocatable, dimension(:,:,:,:) :: HorYMixFlux
      real, allocatable, dimension(:,:,:,:) :: VertMixFlux
#  endif

#  ifdef VERT_DIFF_ANALYSIS
      real, allocatable, dimension(:,:,:,:) :: VertDiffFlux
#  endif
  
# ifdef AVERAGES
      real, allocatable, dimension(:,:,:,:) :: HorXAdvFlux_avg
      real, allocatable, dimension(:,:,:,:) :: HorYAdvFlux_avg
      real, allocatable, dimension(:,:,:,:) :: VertAdvFlux_avg

#   ifdef HOR_DIFF_ANALYSIS
      real, allocatable, dimension(:,:,:,:) :: HorXMixFlux_avg
      real, allocatable, dimension(:,:,:,:) :: HorYMixFlux_avg
      real, allocatable, dimension(:,:,:,:) :: VertMixFlux_avg
#   endif
#   ifdef VERT_DIFF_ANALYSIS
      real, allocatable, dimension(:,:,:,:) :: VertDiffFlux_avg
#   endif
# endif /* AVERAGES */
      contains
      subroutine allocate_physflux_arrays()
      implicit none
#include "param.h"
#include "ncvars.h"
#include "physflux.h"
#include "ncroms_physflux.h"

      integer itrc, status, count
      integer NT_PFA_out

      count = 0
      do itrc = 1, NT
         if (wrt_pfa(itrc)) then
            count = count + 1
            pfa_out(itrc) = count
         end if
      end do
      NT_PFA_out = count
      mpi_master_only print *, 'ALLOCATING phys flux arrays with NT=',
     &     NT_PFA_out

      allocate (HorXAdvFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,NT_PFA_out))
      allocate (HorYAdvFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,NT_PFA_out))
      allocate (VertAdvFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,0:N,
     &     NT_PFA_out))
#  ifdef HOR_DIFF_ANALYSIS
      allocate (HorXMixFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
      allocate (HorYMixFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
      allocate (VertMixFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
#  endif
#  ifdef VERT_DIFF_ANALYSIS
      allocate (VertDiffFlux(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,0:N,
     &     NT_PFA_out))
#  endif

#  ifdef AVERAGES
      allocate (HorXAdvFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
      allocate (HorYAdvFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
      allocate (VertAdvFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,0:N,
     &     NT_PFA_out))
#   ifdef HOR_DIFF_ANALYSIS
      allocate (HorXMixFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
      allocate (HorYMixFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
      allocate (VertMixFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,N,
     &     NT_PFA_out))
#   endif
#   ifdef VERT_DIFF_ANALYSIS
      allocate (VertDiffFlux_avg(GLOBAL_I_ARRAY,GLOBAL_J_ARRAY,0:N,
     &     NT_PFA_out))
#   endif
#  endif /* AVERAGES */

      end subroutine allocate_physflux_arrays
      end module phys_flux
#else
      subroutine phys_flux_mod_empty
      end subroutine phys_flux_mod_empty
#endif /* SOLVE3D && PHYS_FLUX_ANALYSIS && PFA_BY_TRACER */
