/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ====
*/

#undef BALTIC
#undef BASIN           /* Big Bad Basin */
#undef CANBAS2
#undef CANYON_A        /* Canyon_A */
#undef CANYON_B        /* Canyon_B */
#undef DAMEE_B         /* North Atlantic DAMEE, Big Domain */
#undef EKMAN_SPIRAL
#undef GRAV_ADJ        /* Graviational Adjustment */
#undef ISWAKE         /* Island Wake Problem */
#undef NJ_BIGHT        /* New Jersey Bight Application */
#undef PACIFIC         /* North-Equatorial Pacific Application */
#undef PACIFIC_2D     /* Pacific Tsunami model */
#undef OVERFLOW        /* Graviational/Overflow */
#undef SEAMOUNT        /* Seamount */
#undef SOLITON         /* Equatorial Rossby Wave */
#undef RIVER          /* River runoff test problem */
#undef UPWELLING       /* Upwelling */
#undef USWEST          /* US West Coast Application */
#undef WAVE_RAD        /* A test for wave radiation boundaries */

!#define USWC10KM_2D
!#define PACIFIC_25KM
!# define USWC42
# define USSWC1
/*
    Embedded (nested) grid configuration segment
*/

c--#ifndef MAX_GRID_LEVEL
c--# define MAX_GRID_LEVEL 2
c--# include "redefs.X"
c--#endif


/*
   Main switch starts here: model configuration choice.
*/

#if defined BASIN    /* Big Bad Basin Configuration */
# define SOLVE3D
 
# define UV_ADV
# define UV_COR
 
# undef  SALINITY
# undef  NONLIN_EOS
 
# undef  TS_DIF2
# undef  TS_DIF4
 
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# undef  ANA_VMIX
 
#elif defined CANYON_A      /* Canyon A Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define TS_DIF2
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
 
#elif defined CANYON_B      /* Canyon B Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define TS_DIF2
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_VMIX
 
#elif defined DAMEE_B               /* North Atlantic DAMEE configuration */
# define SOLVE3D
c--# define AVERAGES
# define UV_COR
# define UV_ADV
!                       Equation of State
# define SALINITY
# define NONLIN_EOS
# define SPLIT_EOS
!                       Forcing and Climatology
# define TCLIMATOLOGY
# define TNUDGING
# define QCORRECTION
# define SFLX_CORR
!                       Lateral Mixing
# define VIS_GRID
# define DIF_GRID
!                       Vertical Mixing
# define LMD_MIXING
#  undef LMD_RIMIX
#  define LMD_CONVEC
c--#  define LMD_DDMIX
 
c--#  define LMD_KPP
c--#  define LMD_NONLOCAL

!       Grid Configuration and Boundaries

# define CURVGRID
# define SPHERICAL
# define MASKING
# define EASTERN_WALL
# define WESTERN_WALL
# define SOUTHERN_WALL
# define NORTHERN_WALL
 
 
c--# define REST_STATE_TEST    /* Rest-state unforced problem */
# ifdef REST_STATE_TEST     /* (pressure gradient error test) */
#  define ANA_INITIAL
#  define NONLIN_EOS
#  undef SPLIT_EOS
#  define SALINITY
#  define ANA_SMFLUX
#  define ANA_SSFLUX
#  define ANA_STFLUX
#  define ANA_SST
#  define ANA_SRFLUX
#  undef TCLIMATOLOGY
#  undef TNUDGING
#  undef QCORRECTION
#  undef LMD_MIXING
#  define UV_VIS2
# endif


#elif defined EKMAN_SPIRAL     /* Ekman Spiral Test Problem */
# define SOLVE3D
# define UV_ADV
# define UV_COR

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX

# define EW_PERIODIC
# define NS_PERIODIC

#define STATIONS

#elif defined GRAV_ADJ     /* Gravitational Adjustment */
# define SOLVE3D
 
# define UV_ADV
# undef UV_COR
# define UV_VIS2
 
# define TS_DIF2
 
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX

c--# define OBC_WEST
c--# define OBC_EAST 
c--# define OBC_M2ORLANSKI
c--# define OBC_M3ORLANSKI
c--# define OBC_TORLANSKI
!
#elif defined NJ_BIGHT        /* New Jersey Bight Configuration */
# define SOLVE3D
 
# define UV_ADV
# define UV_COR
# define UV_VIS2
 
# define SALINITY
 
# define TS_DIF2
 
# define CURVGRID
# define SPHERICAL
 
# define OBC_EAST
# define OBC_NORTH
# define OBC_SOUTH
# define OBC_FSORLANSKI
# define OBC_M2ORLANSKI
# define OBC_M3ORLANSKI
# define LMD_MIXING
# undef LMD_RIMIX
# define LMD_KPP
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_STFLUX

#elif defined CHABU
# define SOLVE3D
# define UV_COR
c---# define NON_TRADITIONAL
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING
# define MASK_LAND_DATA
# define VERT_COORD_TYPE_SM09

# define EXACT_RESTART
# define AVERAGES
                        ! Equation of State
# define SALINITY
# define NONLIN_EOS
# define SPLIT_EOS

# define DIURNAL_SRFLUX
# define QCORRECTION
!# define SFLX_CORR
!# define SSS_dQdSST    ! dQdSSSt must be read together with
                       ! SSS or not (in this case with SST)

# undef WIND_STRESS_AMP

                       ! Lateral Mixing
# define UV_VIS2
# undef VIS_GRID

# undef ADV_ISONEUTRAL
# define TS_DIF2
# undef DIF_GRID
                       ! Vertical Mixing
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP

# undef LMD_RIMIX
# define LMD_CONVEC
# undef LMD_NONLOCAL
# undef LMD_DDMIX
                       ! Open Boundary Conditions
# define OBC_EAST
# define OBC_WEST
# define OBC_SOUTH
# define OBC_NORTH

# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY 
# define T_FRC_BRY

# define SPONGE


!!hf 
!#elif defined PACIFIC_25KM   /* North-Equatorial Pacific Configuration */
!#elif defined USW_2D || defined USWC10KM_2D
#elif defined USWC42
c# define DIAGNOSTICS_EDDY
c#define BIO_1ST_USTREAM_TEST
c#define MYMIXING

# define SOLVE3D
# define UV_COR
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING
# define MASK_LAND_DATA
# define VERT_COORD_TYPE_SM09

# define EXACT_RESTART
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
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
# undef LMD_DDMIX

#define BULK_FLUX

# ifdef BULK_FLUX
# define DAILYPAR_BEC
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

                       ! Open Boundary Conditions
# define OBC_WEST
# define OBC_SOUTH
# define OBC_NORTH
# define OBC_EAST
# define NO_NUDGING_EAST
# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY 
# define T_FRC_BRY

# define SPONGE

# undef PASSIVE_TRACER
# define BIOLOGY_BEC

#ifdef BIOLOGY_BEC
! this affects parameters for both diaz and phaeo
#define BEC_LIANG
#ifdef BEC_LIANG
!# define PHAEOCYSTIS
!# define UMAX35
!# define UMAX51
!# define PCHL25
#endif

#define SEDIMENT_BIOLOGY

# define BENTHIC_FE_FLUX_SEVERMANN
!# define REDUCED_DIAT_MORT2
!# define QUARTER_DIAT_MORT2
!#define NP16_SP_DIAT
!#define NO_SP_GROWTH
!#define NO_DIAZ_GROWTH

# define OUTPUT_SPECDOM

# ifdef EXPLICIT_SINKING
#  define SEDIMENT_BIOLOGY
# endif
/* this is for the Columbia River, not the outside boundaries */
!#define TCLIM_LOCAL
!#define ANA_TCLIMA
!#define FIXCONC_COLUMBIA

!# define DMS_CYCLE
# define NCYCLE_ANOXIC
!# define NO_PELAGIC_DENIT
!# define ALT_SED_DENIT
!# define N15_CYCLE
# undef BIOLOGY_NPZDOC
# ifdef TCLIMATOLOGY
#  define MULT_CLIM_FILES
# endif
# ifdef DIURNAL_SRFLUX
#  define DAILYPAR_BEC
# endif
# define READ_VARIABLE_ATM_PCO2
!# define VARIABLE_ATM_PCO2
# define VFLX_CORR
# define BGC_FLUX_ANALYSIS
# define LARGE_BGC_FILES
# define PHYS_FLUX_ANALYSIS
# define LARGE_PHYS_FILES
c#define VERT_DIFF_ANALYSIS
c#define HOR_DIFF_ANALYSIS
c#define PFA_BY_TRACER

c#define OUTPUT_AVG_PHYS_BIO
#endif /* BIOLOGY_BEC */


#define LARGE_AVG_FILES
#define LARGE_HIS_FILES
#define LARGE_RST_FILES

c====================================================c
#elif defined USSWC1               /* US SOUTHWEST COAST DX=1KM Configuration */
c# define DIAGNOSTICS_EDDY
#define BIO_1ST_USTREAM_TEST
c#define MYMIXING

# define SOLVE3D
# define UV_COR
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING
# define MASK_LAND_DATA
# define VERT_COORD_TYPE_SM09

# define EXACT_RESTART
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
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
# undef LMD_DDMIX

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

                       ! Open Boundary Conditions
c# define OBC_WEST    % FOR FBTEST
c# define OBC_SOUTH  % FOR FBTEST
c# define OBC_NORTH  % FOR FBTEST
c# undef  OBC_EAST
c# define NO_NUDGING_EAST
c# define OBC_M2FLATHER   % FOR FBTEST
c# define OBC_M3ORLANSKI   % FOR FBTEST
c# define OBC_TORLANSKI   % FOR FBTEST

c# define Z_FRC_BRY   % FOR FBTEST
c# define M2_FRC_BRY    % FOR FBTEST
c# define M3_FRC_BRY     % FOR FBTEST
c# define T_FRC_BRY    % FOR FBTEST

c# define SPONGE    % FOR FBTEST

# undef PASSIVE_TRACER
# define BIOLOGY_BEC   

# ifdef BIOLOGY_BEC
! this affects parameters for both diaz and phaeo
#  define BEC_LIANG
#  ifdef BEC_LIANG
!# define PHAEOCYSTIS
!# define UMAX35
!# define UMAX51
!# define PCHL25
#  endif

c# define SEDIMENT_BIOLOGY    % FOR FBTEST

c# define BENTHIC_FE_FLUX_SEVERMANN    % FOR FBTEST
!# define REDUCED_DIAT_MORT2
!# define QUARTER_DIAT_MORT2
!#define NP16_SP_DIAT
!#define NO_SP_GROWTH
!#define NO_DIAZ_GROWTH

# define OUTPUT_SPECDOM

# ifdef EXPLICIT_SINKING
#  define SEDIMENT_BIOLOGY
# endif
/* this is for the Columbia River, not the outside boundaries */
!#define TCLIM_LOCAL
!#define ANA_TCLIMA
!#define FIXCONC_COLUMBIA

!# define DMS_CYCLE
# define NCYCLE_ANOXIC
!# define NO_PELAGIC_DENIT
!# define ALT_SED_DENIT
!# define N15_CYCLE
# undef BIOLOGY_NPZDOC
# ifdef TCLIMATOLOGY
#  define MULT_CLIM_FILES
# endif
# ifdef DIURNAL_SRFLUX
#  define DAILYPAR_BEC
# endif
!# define READ_VARIABLE_ATM_PCO2     % FOR FBTEST
!# define VARIABLE_ATM_PCO2

             /* point source */
# define PSOURCE
# define Q_PSOURCE

# define VFLX_CORR
# define BGC_FLUX_ANALYSIS
# define LARGE_BGC_FILES
# define PHYS_FLUX_ANALYSIS
# define LARGE_PHYS_FILES
# define VERT_DIFF_ANALYSIS
c#define HOR_DIFF_ANALYSIS
# define PFA_BY_TRACER

c#define OUTPUT_AVG_PHYS_BIO
# endif /* BIOLOGY_BEC */

c# define TIDES      % FOR FBTEST
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
#  define TIDERAMP
# endif



# define LARGE_AVG_FILES
# define LARGE_HIS_FILES
# define LARGE_RST_FILES
c=======================================================c

#elif defined PACIFIC_2D   /* Pacific Tsynami model */
# define UV_COR
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING
# undef VIS_GRID
# define UV_VIS2
# define ANA_SMFLUX
# define ANA_INITIAL

# define OBC_WEST
# define OBC_SOUTH
# define OBC_M2FLATHER
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
c--# define OBC_M2ORLANSKI
c--# define OBC_VOLCONS
# define SPONGE

#elif defined OVERFLOW      /* Gravitational Overflow */
# define SOLVE3D
 
# define UV_ADV
# define UV_VIS2
 
# define TS_DIF2
 
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
 
#elif defined SEAMOUNT     /* Seamount Configuration */
# define SOLVE3D
# define EXACT_RESTART

# define UV_ADV
# define UV_COR
c--# define UV_VIS2
 
# undef TS_DIF2
# undef  TS_DIF4
 
c--# define  NONLIN_EOS
# undef SALINITY
 
# define EW_PERIODIC
# define NS_PERIODIC
 
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX


#elif defined ISWAKE     /* Island Wake Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR

#define SPLIT_EOS

# define UV_VIS2
c---# define SPONGE
# define LINEAR_DRAG_ONLY

# define OBC_WEST
# define OBC_EAST
c--# define OBC_SOUTH
c--# define OBC_NORTH

c--# define OBC_M2ORLANSKI
c--# define OBC_M3ORLANSKI
c--# define OBC_TORLANSKI
c--# define OBC_M2SPECIFIED
# define OBC_M3SPECIFIED
# define OBC_TSPECIFIED

#define OBC_M2FLATHER
 
#define EXACT_RESTART

# define ANA_BRY
# define T_FRC_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY

c--# define AVERAGES

# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

# define MASK_LAND_DATA

! Vertical mixing: nothing is defined here => use
! externally supplied constant background value.

#elif defined SOLITON    /* Equatorial Rossby Soliton */
# undef  SOLVE3D
 
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
 
# define UV_COR
# define UV_ADV
# undef UV_VIS2
# define EW_PERIODIC
c--# define NS_PERIODIC

c--# define OBC_WEST
c--# define OBC_EAST
c--# define OBC_NORTH
c--# define OBC_SOUTH
c--# define OBC_M2ORLANSKI
 
#elif defined RIVER     /* River run-off test problem */
# define SOLVE3D
# define UV_ADV
# define UV_COR

# define NONLIN_EOS
# define SPLIT_EOS
# define SALINITY

# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

# define PSOURCE
# define ANA_PSOURCE

# define EASTERN_WALL
# define WESTERN_WALL
# define NORTHERN_WALL
# define OBC_SOUTH
# define OBC_TORLANSKI
# define OBC_M2ORLANSKI
# define OBC_M3ORLANSKI
 
#elif defined UPWELLING     /* Upwelling Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR

c--# define ADV_ISONEUTRAL

# undef NONLIN_EOS
# define SALINITY

# define EW_PERIODIC
 
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# undef ANA_VMIX

# ifndef ANA_VMIX
#  define LMD_MIXING
#  define LMD_KPP
#  undef LMD_BKPP
c--# define LMD_RIMIX
c--# define LMD_CONVEC
c--# define LMD_NONLOCAL
#endif

#elif defined CANBAS2   /* Canary Basin model */

# define SOLVE3D
# define UV_ADV
# define UV_COR

# define NONLIN_EOS
# define SPLIT_EOS
# define SALINITY

# define EXACT_RESTART
# define AVERAGES

# define UV_VIS2
# define TS_DIF2

# define QCORRECTION

# define LMD_MIXING
# define LMD_KPP
# undef LMD_RIMIX
# define LMD_CONVEC
# undef LMD_NONLOCAL

# define CURVGRID
# define SPHERICAL
# define MASKING

# define OBC_WEST
# define OBC_EAST
# define OBC_NORTH
# define OBC_SOUTH

# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY

#elif defined BALTIC  /* Baltic Sea model */

# define SOLVE3D
# define UV_ADV
# define UV_COR

# define NONLIN_EOS
# define SPLIT_EOS
# define SALINITY

# define EXACT_RESTART
# undef AVERAGES

# define UV_VIS2
# define TS_DIF2


# define LMD_MIXING
# define LMD_KPP
# undef LMD_RIMIX
# define LMD_CONVEC
# undef LMD_NONLOCAL

# define CURVGRID
# define SPHERICAL
# define MASKING

# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

 
#elif defined USWEST   /* US West Coast Configuration */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                       /* Equation of State */
# define NONLIN_EOS
# define SPLIT_EOS
# define SALINITY

                       /* Forcing and Climatology */
# define QCORRECTION
# define SPONGE
# define SFLX_CORR

                      /* Lateral viscosity/mixing  */
# define UV_VIS2
# undef VIS_GRID
# define TS_DIF2
# undef DIF_GRID
                      /* Vertical Mixing */
# define LMD_MIXING
# undef LMD_RIMIX
# define LMD_CONVEC
# undef LMD_DDMIX
# define LMD_KPP
# undef LMD_NONLOCAL
                      /* Grid Configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
                      /* Open Boundary Conditions */
# undef OBC_EAST
c# define OBC_WEST
c# define OBC_NORTH
c# define OBC_SOUTH
c--> # define OBC_VOLCONS
c--> # define OBC_FLUX_CORR

# undef OBC_M2ORLANSKI
# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

#define BRY
#ifdef BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# undef TNUDGING
# undef M3NUDGING
# undef M2NUDGING
# undef CLIMAT_TS_MIXH
#else
# define TNUDGING
# define M3NUDGING
# define M2NUDGING
#endif

# undef OBC_TSPECIFIED
# undef OBC_M2SPECIFIED
# undef OBC_M3SPECIFIED

# define EXACT_RESTART
# define AVERAGES



#elif defined WAVE_RAD
# undef  SOLVE3D
# define UV_COR
# define UV_ADV
# undef UV_VIS2

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX

c# define OBC_WEST
c# define OBC_EAST
c# define OBC_NORTH
c# define OBC_SOUTH
c--# define OBC_M2ORLANSKI
# define OBC_M2FLATHER
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY

#endif
 
#include "set_global_definitions.h"
 
! a few checks
#ifdef FULL_PHYS_FLUX_ANALYSIS
# ifdef PFA_BY_TRACER
error this is not yet implemented
# endif
#endif
