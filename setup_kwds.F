#include "cppdefs.h"

      subroutine setup_kwds (ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new keywords were !!!!!!
!!!!!! introduced into "read_inp.F". To create or refresh this !!!!!!
!!!!!! file use compile and execute "checkkwds.F" as an        !!!!!!
!!!!!! independent program, or use commands "make checkkwds"   !!!!!!
!!!!!! or "make depend".                                       !!!!!!

      implicit none
      integer ierr, is,ie
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      do is=1,max_opt_size
        kwds(is:is)=' '
      enddo
      is=1
#if defined BULK_DIAGS || defined READ_ATM_PCO2 || defined READ_VARIABLE_ATM_PCO
#endif
#ifndef NONLIN_EOS
#endif
#ifdef FLOATS
#endif
#ifdef PSOURCE
#endif
#ifdef BGC_FLUX_ANALYSIS
#endif
#ifdef PHYS_FLUX_ANALYSIS
#endif
#ifdef WKB_WWAVE
#endif
#ifdef MPI
# ifdef HOFFMAN2
# else
# endif
#endif
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
# ifdef LMD_KPP
# endif
#endif
#if defined SOLITON
#elif defined SEAMOUNT
#elif defined WAVE_RAD
#elif defined ISWAKE
#elif defined UPWELLING
#elif defined GRAV_ADJ
#elif defined IBERIA
#elif defined RIVER
#elif defined PACIFIC
#elif defined PACIFIC_2D
#elif defined CANBAS2
#elif defined SPIRAL
#elif defined PACBIG
#elif defined NWPAC
#elif defined SMSC
#elif defined SMSC_SUB
#elif defined USWEST
#elif defined SCB_TEST
#elif defined DAMEE_B
#else
#endif
#ifdef MPI
#else
#endif
#ifdef AVERAGES
#endif
#ifdef STATIONS
#endif
#ifdef DIAGNOSTICS_EDDY
# ifdef AVERAGES
# endif
#endif
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='title'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='time_stepping'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='S-coord'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='rho0'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifndef NONLIN_EOS
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='lin_rho_eos'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef SALINITY
#  endif
#  ifdef SALINITY
#  endif
# endif
#endif /* SOLVE3D */
#ifdef UV_VIS2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='lateral_visc'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if !defined EW_PERIODIC || !defined NS_PERIODIC
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='gamma2'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='tracer_diff2'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bottom_drag'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='vertical_mixing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
# ifdef MY25_MIXING
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='MY_bak_mixing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef Q_DIF2
#  endif
#  ifdef Q_DIF4
#  endif
# endif
#endif /* SOLVE3D */
#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING \
  || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M2NUDGING \
                       || defined M3NUDGING  || defined WKB_FRC_BRY
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='nudg_cof'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
#endif
#ifdef SPONGE
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='v_sponge'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_FLUX_CORR
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='tau_obc'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifndef ANA_GRID
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='grid'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='initial'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef ANA_INITIAL
#endif
#if defined MPI && defined PARALLEL_FILES
#endif
#ifdef ANA_INITIAL
#endif
#if !defined ANA_SMFLUX  ||  defined SOLVE3D  &&\
   ( !defined ANA_STFLUX  ||\
     (defined SG_BBL96    && !defined ANA_BSEDIM) ||\
     (defined SG_BBL96    && !defined ANA_WWAVE)  ||\
     (defined QCORRECTION && !defined ANA_SST)    ||\
     (defined SALINITY    && !defined ANA_SSFLUX) ||\
     (defined LMD_KPP     && !defined ANA_SRFLUX))
# if !defined BULK_FLUX
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='forcing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
# endif /* !BULK_FLUX */
# if defined BULK_FLUX
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bulk_forcing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
# endif /* BULK_FLUX */
#endif
#if defined PSOURCE && !defined ANA_PSOURCE
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='point_source'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
#if defined SOLVE3D &&( defined TCLIMATOLOGY ||( defined TNUDGING &&\
                        defined T_FRC_BRY )) && !defined ANA_TCLIMA
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='climatology'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
#if defined T_FRC_BRY || defined M2_FRC_BRY || defined M3_FRC_BRY \
                      || defined Z_FRC_BRY  || defined WKB_FRC_BRY
# ifndef ANA_BRY
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='boundary'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
# endif
#endif
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='restart'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='history'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='primary_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
# ifdef SEDIMENT_BIOLOGY
# endif /* SEDIMENT_BIOLOGY */
#endif
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
# ifdef SEDIMENT_BIOLOGY
# endif /* SEDIMENT_BIOLOGY */
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='auxiliary_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SALINITY
# endif
# ifdef LMD_KPP
# endif
# ifdef LMD_BKPP
# endif
# ifdef SALINITY
# endif
# ifdef LMD_KPP
# endif
# ifdef LMD_BKPP
# endif
# ifdef SALINITY
# endif
# ifdef LMD_KPP
# endif
# ifdef LMD_BKPP
# endif
#endif /* SOLVE3D */
#ifdef AVERAGES
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='primary_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
#  ifdef SEDIMENT_BIOLOGY
#  endif /* SEDIMENT_BIOLOGY */
# endif
# ifdef SOLVE3D
# endif
#ifdef SOLVE3D
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='auxiliary_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
# ifdef OUTPUT_AVG_PHYS_BIO
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='averages_bio'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif /* OUTPUT_AVG_PHYS_BIO */
#ifdef OUTPUT_SST_AVG
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='SST_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif /* OUTPUT_SST_AVG */
# endif /* SOLVE3D */
#endif /* AVERAGES */
#ifdef BGC_FLUX_ANALYSIS
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bgc_flux_histories'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bgc_flux_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
# endif /* AVERAGES */
#endif /* BGC_FLUX_ANALYSIS */
#ifdef PHYS_FLUX_ANALYSIS
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='phys_flux_histories'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='phys_flux_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
#  ifdef PFA_BY_TRACER
#  endif
#  ifdef PFA_BY_TRACER
#  endif
# endif /* AVERAGES */
#endif /* PHYS_FLUX_ANALYSIS */
#ifdef BULK_DIAGS
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bulk_diags_histories'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bulk_diags_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
# endif /* AVERAGES */
#endif /* BULK_DIAGS */
#ifdef STATIONS
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='stations'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='primary_stations'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
#ifdef SOLVE3D
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='auxiliary_stations'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
# endif /* SOLVE3D */
#endif /* STATIONS */
#ifdef FLOATS
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='floats'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='float_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef AGRIF
# endif
# ifdef AGRIF
# endif
#endif /* FLOATS */
#if defined BIOLOGY_BEC || defined BIOLOGY_NPZDOC
#ifdef READ_ATM_PCO2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='pCO2_atm'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#endif /* READ_ATM_PCO2 */
# ifdef VARIABLE_ANN_ATM_PCO2
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='start_year'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='futr_scen'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif /* VARIABLE_ANN_ATM_PCO2 */
# ifdef READ_VARIABLE_ATM_PCO2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='pCO2_atm_file'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif /* READ_VARIABLE_ATM_PCO2 */
#endif /* BIOLOGY_BEC || BIOLOGY_NPZDOC */
#if defined DIAGNOSTICS_EDDY
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='diags_eddy'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='diags_eddy_avg'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
      ie=is +25
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='diags_eddy_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +25
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='diags_eddy_average_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_EDDY */
#ifdef MPI
#endif
      return
  99  mpi_master_only write(*,'(/1x,2A/25x,A/)')
     &  '### ERROR: setup_kwds :: Insufficient size of string "kwds" ',
     &  'in file "strings.h".',  'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
