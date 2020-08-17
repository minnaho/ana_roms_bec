/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ====
*/

# define FLAT_BOTTOM

/*
   Main switch starts here: model configuration choice.
*/

#if defined FLAT_BOTTOM           /* Flat bottom with idealized stratification*/
                                  /* and open boundaries*/

# define SOLVE3D
# define UV_COR
# define UV_ADV
# define CURVGRID
c# define SPHERICAL
# define MASKING
# define MASK_LAND_DATA
# define VERT_COORD_TYPE_SM09

c# define EXACT_RESTART
# define AVERAGES
                        ! Equation of State
# define SALINITY
# define NONLIN_EOS
# define SPLIT_EOS


# undef WIND_STRESS_AMP

                       ! Lateral Mixing
# define UV_VIS2
# undef VIS_GRID

# undef ADV_ISONEUTRAL
# define TS_DIF2
# undef DIF_GRID
                       ! Vertical Mixing
# define LMD_MIXING
c# define LMD_KPP
c# define LMD_BKPP
c# define LMD_RIMIX
c# define LMD_CONVEC
c# define LMD_NONLOCAL
c# undef LMD_DDMIX

c# define BULK_FLUX   % FOR FBTEST

# ifdef BULK_FLUX
#  define DAILYPAR_BEC
c#  define SIMPLE_BULK
#  define COUPLED_SURF_CURR
c#  define COUPLED_SURF_CURR_ONLY_DRAG
c#  define COUPLED_SURF_CURR_ONLY_REL
#  define WND_AT_RHO_POINTS
#  define BULK_DIAGS
#  undef DIURNAL_SRFLUX
#  define SALINITY_RESTORING
# else
#  undef QCORRECTION
#  undef DIURNAL_SRFLUX
#  define STR_AT_RHO_POINTS
# endif

# define ANA_SMFLUX    
# define ANA_STFLUX    
# define ANA_SSFLUX    
# define ANA_SRFLUX    
# define ANA_INITIAL   
# define ANA_GRID

                       ! Open Boundary Conditions
# define OBC_WEST    % FOR FBTEST
# define OBC_SOUTH  % FOR FBTEST
# define OBC_NORTH  % FOR FBTEST
# define OBC_EAST
# define OBC_M2ORLANSKI
c# define NO_NUDGING_EAST
c# define OBC_M2FLATHER   % FOR FBTEST
# define OBC_M3ORLANSKI   % FOR FBTEST
# define OBC_TORLANSKI   % FOR FBTEST

# undef PASSIVE_TRACER

             /* point source */
# define PSOURCE
!# define Q_PSOURCE 
# define ANA_PSOURCE
# define VFLX_CORR

# define LARGE_HIS_FILES

# define BIOLOGY_BEC

# ifdef BIOLOGY_BEC
# define BEC_LIANG
# define NCYCLE_ANOXIC
# endif

 
#include "set_global_definitions.h"
 
#endif
