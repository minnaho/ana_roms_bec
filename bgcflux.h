! This header file contains all variables for the 
! computation of biogeochemical fluxes.
! It must always be included in the biology and some other routines.
! The variables for the netcdf output are in ncroms_bgcflux.h, which
! needs to be included only if BGC_FLUX_ANALYSIS is defined.

#if defined SOLVE3D && defined BIOLOGY_NPZDOC

! Fluxes: all bgc fluxes are defined in a way that they will be positive under
! normal circumstances
      ! fluxes of Nitrogen [mmol N m^-2]
      integer NFlux_NewProd, NFlux_RegProd, NFlux_Grazing
      integer NFlux_SlopFeed, NFlux_Zfecp, NFlux_Pmort, NFlux_Zmetab
      integer NFlux_Zexcr, NFlux_ZmortS, NFlux_ZmortL, NFlux_DetCoagP
      integer NFlux_DetCoagD, NFlux_ReminS, NFlux_ReminL, NFlux_Nitrif
      integer NumFluxTermsN, NumFluxTerms

      parameter(NFlux_NewProd = 1)   ! new production (NO3 -> Phyto.)
      parameter(NFlux_RegProd = 2)   ! regenerated prod. (NH4 -> Phyto.)
      parameter(NFlux_Grazing = 3)   ! grazing (Phytoplankton -> Zoo.)
      parameter(NFlux_SlopFeed = 4)  ! sloppy feeding (Phyto. -> SDetN)
      parameter(NFlux_Zfecp = 5)     ! fecal pellets (Phyto. -> LDetN)
      parameter(NFlux_Pmort = 6)     ! Phyto. mortality (Phy. -> SDetN)
      parameter(NFlux_Zmetab = 7)    ! Zoopl. metabolism (Zoo. -> NH4)
      parameter(NFlux_Zexcr = 8)     ! Zoopl. excretion (Zoo. -> NH4)
      parameter(NFlux_ZmortS = 9)     ! Zoopl. mortality (Zoo. -> SDetN)
      parameter(NFlux_ZmortL = 10)    ! Zoopl. mortality (Zoo. -> LDetN)
      parameter(NFlux_DetCoagP = 11)  ! Coagul. of detr. (Phy. -> LDetN) 
      parameter(NFlux_DetCoagD = 12) ! Coagul. of detr. (SDetN -> LDetN)
      parameter(NFlux_ReminS = 13)   ! Remineralization (SDetN -> NH4)
      parameter(NFlux_ReminL = 14)   ! Remineralization (LDetN -> NH4)
      parameter(NFlux_Nitrif = 15)   ! Nitrification (NH4 -> NO3)
      parameter(NumFluxTermsN = NFlux_Nitrif)
#   ifdef CARBON
      integer CFlux_Zresp
      integer CFlux_DetCoagD, CFlux_ReminS, CFlux_ReminL, CFlux_Dissolv
      integer NumFluxTermsC
      ! fluxes of and changes in carbon [mmol C m^-2]
      parameter(CFlux_Zresp = NFlux_Nitrif + 1)      ! Zoopl. respiration (Zoo. -> DIC)
      parameter(CFlux_DetCoagD = CFlux_Zresp + 1)! Coagul. of det. (SDetC -> LDetC)
      parameter(CFlux_ReminS = CFlux_DetCoagD + 1)  ! Remineralization (SDetC -> NH4)
      parameter(CFlux_ReminL = CFlux_ReminS + 1)    ! Remineralization (LDetC -> NH4)
      parameter(CFlux_Dissolv = CFlux_ReminL + 1)   ! Dissolution of CaCO3 (CaCO3 -> DIC)
      parameter(NumFluxTermsC = CFlux_Dissolv - NumFluxTermsN)

      parameter(NumFluxTerms = NumFluxTermsN + NumFluxTermsC)
#   else /* CARBON */
      parameter(NumFluxTerms = NumFluxTermsN)
#   endif /* CARBON */

      integer NumGasExcTerms 
#   ifdef OXYGEN
      ! gas exchange fluxes
      integer OFlux_GasExc
      parameter(OFlux_GasExc = 1)
#    ifdef CARBON
      integer CFlux_GasExc
      parameter(CFlux_GasExc = 2)
      parameter(NumGasExcTerms = 2)
#    else /* CARBON */
      parameter(NumGasExcTerms = 1)
#    endif /* CARBON */
#   else /* OXYGEN */
      parameter(NumGasExcTerms = 0)
#   endif /* OXYGEN */  

      ! vertical sinking fluxes
      integer PFlux_VSink, SDNFlux_VSink, LDNFlux_VSink, NumVSinkTerms
      parameter(PFlux_VSink = 1)    ! phytoplankton
      parameter(SDNFlux_VSink = 2)  ! small detritus N
      parameter(LDNFlux_VSink = 3)  ! large detritus N
#   ifdef CARBON
      integer SDCFlux_VSink, LDCFlux_VSink, CaCO3Flux_VSink
      parameter(SDCFlux_VSink = 4)  ! small detritus C
      parameter(LDCFlux_VSink = 5)  ! large detritus C
      parameter(CaCO3Flux_VSink = 6)! CaCO3
      parameter(NumVSinkTerms = 6)
#   else /* CARBON */
      parameter(NumVSinkTerms = 3)
#   endif /* CARBON */
#   ifdef SEDIMENT_BIOLOGY
      integer NFlux_ReminSed, NumSedFluxTerms
      parameter(NFlux_ReminSed = 1) ! Remineralization (Sed. OrgN -> NH4)
#     ifdef CARBON
      integer CFlux_ReminSed, CFlux_DissolvSed
      parameter(CFlux_ReminSed = NFlux_ReminSed + 1) ! Remineralization (Sed. OrgC -> DIC)
      parameter(CFlux_DissolvSed = CFlux_ReminSed + 1) ! Dissolution of SedCaCO3
      parameter(NumSedFluxTerms = CFlux_DissolvSed)
#     else /* CARBON */
      parameter(NumSedFluxTerms = NFlux_ReminSed)
#    endif /* CARBON */
#   endif /* SEDIMENT_BIOLOGY */

! The PARinc value represents the incoming radiation. PAR values are
! at rho-points, which is the value used for determining new and
! regenerated production. All values are averaged over the number of
! biological time steps.
      real PAR(GLOBAL_2D_ARRAY,N)
      real PARinc(GLOBAL_2D_ARRAY)
      common /ocean_par/PAR,PARinc

#ifdef OXYGEN
      real u10(GLOBAL_2D_ARRAY)
      real Kv_O2(GLOBAL_2D_ARRAY)
      real O2satu(GLOBAL_2D_ARRAY)
      common /gasexc_o2/ u10, Kv_O2, O2satu
#  ifdef CARBON
      real Kv_CO2(GLOBAL_2D_ARRAY)
      real CO2sol(GLOBAL_2D_ARRAY)
      !DL: use pH_hist instead of pHsw:
      real pH_hist(GLOBAL_2D_ARRAY)
      real pCO2sw(GLOBAL_2D_ARRAY)
      real pCO2air(GLOBAL_2D_ARRAY)
      common /gasexc_co2/Kv_CO2, CO2sol, pH_hist, pCO2sw, pCO2air
#  endif /* CARBON */
#endif /* OXYGEN */

# ifdef SEDIMENT_BIOLOGY
! These components represent tracer concentrations
! in the sediment. The unit is mmol m-2.
      real t_sed(GLOBAL_2D_ARRAY,NT_sed)
      common /ocean_t_sed/t_sed
# endif /* SEDIMENT_BIOLOGY */

      integer bio_count
      common /communicators_npzdoc/ bio_count

      real*QUAD global_sum(0:2*NT+1)
#  if defined SALINITY && defined VFLX_CORR
     &     , global_srf_sum(0:NT)
#  endif /* SALINITY && VFLX_CORR */
      common /communicators_rq_npzdoc/ global_sum
# if defined SALINITY && defined VFLX_CORR
     &     , global_srf_sum
# endif /* SALINITY && VFLX_CORR */

# ifdef BGC_FLUX_ANALYSIS
#  ifdef SEDIMENT_BIOLOGY
      real SedFlux(GLOBAL_2D_ARRAY,NumSedFluxTerms)
      common /ocean_bgc_sedflux/SedFlux
#  endif /* SEDIMENT_BIOLOGY */

      real Flux(GLOBAL_2D_ARRAY,N,NumFluxTerms)
      common /ocean_bgc_flux/Flux

#  ifdef OXYGEN
      real GasExcFlux(GLOBAL_2D_ARRAY,NumGasExcTerms)
      common /ocean_bgc_gasexcflux/GasExcFlux
#  endif /* OXYGEN */

      real VSinkFlux(GLOBAL_2D_ARRAY,0:N,NumVSinkTerms)
      common /ocean_bgc_vsinkflux/VSinkFlux

# endif /* BGC_FLUX_ANALYSIS */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ifdef AVERAGES

      real PAR_avg(GLOBAL_2D_ARRAY,N)
      real PARinc_avg(GLOBAL_2D_ARRAY)
      common /avg_par/ PAR_avg,PARinc_avg

#  ifdef BGC_FLUX_ANALYSIS
      real PAR_flux_avg(GLOBAL_2D_ARRAY,N)
      real PARinc_flux_avg(GLOBAL_2D_ARRAY)
      common /avg_par_flux/ PAR_flux_avg,PARinc_flux_avg

      real zeta_bgc_flux_avg(GLOBAL_2D_ARRAY)
      common /ocean_zeta_bgc_flux_avg/zeta_bgc_flux_avg
#  endif /* BGC_FLUX_ANALYSIS */

#  ifdef OXYGEN
      real u10_avg(GLOBAL_2D_ARRAY)
      real Kv_O2_avg(GLOBAL_2D_ARRAY)
      real O2satu_avg(GLOBAL_2D_ARRAY)
      common /gasexc_o2_avg/ u10_avg, Kv_O2_avg, O2satu_avg
#   ifdef CARBON
      real Kv_CO2_avg(GLOBAL_2D_ARRAY)
      real CO2sol_avg(GLOBAL_2D_ARRAY)
      real pco2_avg(GLOBAL_2D_ARRAY)
      real pCO2air_avg(GLOBAL_2D_ARRAY)
      real pH_avg(GLOBAL_2D_ARRAY)
      common /gasexc_co2_avg/ Kv_CO2_avg, CO2sol_avg, pco2_avg, 
     &     pCO2air_avg, pH_avg
#   endif /* CARBON */
#  endif /* OXYGEN */

# ifdef SEDIMENT_BIOLOGY
      real t_sed_avg(GLOBAL_2D_ARRAY,NT_sed)
      common /ocean_t_sed_avg/t_sed_avg
# endif /* SEDIMENT_BIOLOGY */

#  ifdef BGC_FLUX_ANALYSIS
#   ifdef SEDIMENT_BIOLOGY
      real SedFlux_avg(GLOBAL_2D_ARRAY,NumSedFluxTerms)
      common /ocean_bgc_sedflux_avg/SedFlux_avg
#   endif /* SEDIMENT_BIOLOGY */

      real Flux_avg(GLOBAL_2D_ARRAY,N,NumFluxTerms)
      common /bgc_flux_avg/Flux_avg

#   ifdef OXYGEN
      real GasExcFlux_avg(GLOBAL_2D_ARRAY,NumGasExcTerms)
      common /bgc_gasexcflux_avg/GasExcFlux_avg
#   endif /* OXYGEN */

      real VSinkFlux_avg(GLOBAL_2D_ARRAY,0:N,NumVSinkTerms)
      common /bgc_vsinkflux_avg/VSinkFlux_avg
#  endif /* BGC_FLUX_ANALYSIS */
# endif /* AVERAGES */
#endif /* SOLVE3D && BIOLOGY_NPZDOC */
