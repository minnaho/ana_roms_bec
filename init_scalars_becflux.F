#include "cppdefs.h"
#if defined BIOLOGY_BEC && defined BGC_FLUX_ANALYSIS
       subroutine init_scalars_becflux()
       implicit none

#include "param.h"
#include "ecosys_bec.h"
#include "ncroms_bgcflux.h"

       integer::ind

!
! Initialize netCDF files IDs to closed status.
!
       ncid_bgc_flux_his=-1
       nrec_bgc_flux_his=0
#ifdef AVERAGES
       ncid_bgc_flux_avg=-1
       nrec_bgc_flux_avg=0
#endif



#ifdef REDUCED_BGC_OUTPUT
                write_bgcflux = .true.
#endif

       ind=1
       vname_bgcflux(1,ind)='WS'
       vname_bgcflux(2,ind)='Wind speed'
       vname_bgcflux(3,ind)='m/s '

       ind=ind+1
       vname_bgcflux(1,ind)='XKW'
       vname_bgcflux(2,ind)='XKW_AVG'
       vname_bgcflux(3,ind)='m/s '


       ind=ind+1
       vname_bgcflux(1,ind)='ATM_PRESS'
       vname_bgcflux(2,ind)='Atmospheric pressure'
       vname_bgcflux(3,ind)='atm '


       ind=ind+1
       vname_bgcflux(1,ind)='SCHMIDT_O2'
       vname_bgcflux(2,ind)='Schmidt number for O2'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif


       ind=ind+1
       vname_bgcflux(1,ind)='O2SAT'
       vname_bgcflux(2,ind)='O2 saturation'
       vname_bgcflux(3,ind)='mmol/m3 '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif



       ind=ind+1
       vname_bgcflux(1,ind)='FG_O2'
       vname_bgcflux(2,ind)='Air-sea flux of O2'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif


# ifdef DMS_CYCLE
       ind=ind+1
       vname_bgcflux(1,ind)='FG_DMS'
       vname_bgcflux(2,ind)='Air-sea flux of DMS'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
       ind=ind+1
       vname_bgcflux(1,ind)='FG_N2O'
       vname_bgcflux(2,ind)='Air-sea flux of excess N2O'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='FG_N2'
       vname_bgcflux(2,ind)='Air-sea flux of excess N2'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Sed_denitr'
       vname_bgcflux(2,ind)='Sediment denitrification'
       vname_bgcflux(3,ind)='mmol/m3/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#  ifdef N15_CYCLE
       ind=ind+1
       vname_bgcflux(1,ind)='FG_N2O_15N'
       vname_bgcflux(2,ind)='Air-sea flux of excess N2O_15N'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='FG_N2_15N'
       vname_bgcflux(2,ind)='Air-sea flux of excess N2_15N'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Sed_denitr_15N'
       vname_bgcflux(2,ind)='Sediment denitrification of 15N'
       vname_bgcflux(3,ind)='mmol/m3/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#  endif
# endif /* NCYCLE_ANOXIC */

# ifdef NOT_USED___NEEDS_RECODING_IF_USED
       ind=ind+1
       vname_bgcflux(1,ind)='pH'
       vname_bgcflux(2,ind)='pH value'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='pCO2'
       vname_bgcflux(2,ind)='Surface water pCO2'
       vname_bgcflux(3,ind)='ppmv '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif
 
        ind=ind+1
       vname_bgcflux(1,ind)='DpCO2'
       vname_bgcflux(2,ind)='DPCO2'
       vname_bgcflux(3,ind)='ppmv '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# endif /* NOT_USED___NEEDS_RECODING_IF_USED */

       ind=ind+1
       vname_bgcflux(1,ind)='SCHMIDT_CO2'
       vname_bgcflux(2,ind)='Schmidt number for CO2'
       vname_bgcflux(3,ind)=' '

#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='CO2STAR'
       vname_bgcflux(2,ind)='CO2STAR_AVG'
       vname_bgcflux(3,ind)='mmol/m3 '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif


       ind=ind+1
       vname_bgcflux(1,ind)='DCO2STAR'
       vname_bgcflux(2,ind)='DCO2STAR_AVG'
       vname_bgcflux(3,ind)='mmol/m3'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='FG_CO2'
       vname_bgcflux(2,ind)='Air-sea flux of CO2'
       vname_bgcflux(3,ind)='mmol/m2/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='IRON_FLUX'
       vname_bgcflux(2,ind)='Iron surface flux'
       vname_bgcflux(3,ind)='mmol/m2/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PARinc'
       vname_bgcflux(2,ind)=
     &      'Inc. Photosynthetically available radiation'
       vname_bgcflux(3,ind)='W/m2'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#ifdef SEDIMENT_BIOLOGY
       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_POC'
       vname_bgcflux(2,ind)=
     &      'Flux of POC into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_CaCO3'
       vname_bgcflux(2,ind)=
     &      'Flux of CaCO3 into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_Si'
       vname_bgcflux(2,ind)=
     &      'Flux of silicate into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# ifdef N15_CYCLE
       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_PON15'
       vname_bgcflux(2,ind)=
     &      'Flux of PON15 into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# endif
#endif /* SEDIMENT_BIOLOGY */

!!!!!!!!!!!!!!!!!!!!!!
! 3D fluxes begin here

       ind=ind+1
       vname_bgcflux(1,ind)='PO4_RESTORE'
       vname_bgcflux(2,ind)='PO4 restoring flux'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='NO3_RESTORE'
       vname_bgcflux(2,ind)='NO3 restoring flux'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SIO3_RESTORE'
       vname_bgcflux(2,ind)='SiO3 restoring flux'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PAR'
       vname_bgcflux(2,ind)='Photosynthetically available radiation'
       vname_bgcflux(3,ind)='W/m2'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PO4STAR'
       vname_bgcflux(2,ind)='PO4STAR_AVG'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='POC_FLUX_IN'
       vname_bgcflux(2,ind)='POC_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='POC_PROD'
       vname_bgcflux(2,ind)='POC_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='POC_REMIN'
       vname_bgcflux(2,ind)='POC remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_FLUX_IN'
       vname_bgcflux(2,ind)='CACO3_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_PROD'
       vname_bgcflux(2,ind)='CaCO3 production'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_REMIN'
       vname_bgcflux(2,ind)='CaCO3 remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SIO2_FLUX_IN'
       vname_bgcflux(2,ind)='SIO2_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SIO2_PROD'
       vname_bgcflux(2,ind)='SiO2 production'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SIO2_REMIN'
       vname_bgcflux(2,ind)='SiO2 remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DUST_FLUX_IN'
       vname_bgcflux(2,ind)='DUST_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DUST_REMIN'
       vname_bgcflux(2,ind)='Dust remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='P_IRON_FLUX_IN'
       vname_bgcflux(2,ind)='P_IRON_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='P_IRON_PROD'
       vname_bgcflux(2,ind)='P_IRON_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='P_IRON_REMIN'
       vname_bgcflux(2,ind)='P_IRON remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_SP'
       vname_bgcflux(2,ind)='Grazing of small phytoplankton'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_DIAT'
       vname_bgcflux(2,ind)='Grazing of diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_TOT'
       vname_bgcflux(2,ind)='Total grazing flux'
       vname_bgcflux(3,ind)='mmol/m3/s'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SP_LOSS'
       vname_bgcflux(2,ind)='SP_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_LOSS'
       vname_bgcflux(2,ind)='DIAT_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='ZOO_LOSS'
       vname_bgcflux(2,ind)='ZOO_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SP_AGG'
       vname_bgcflux(2,ind)='Aggregation of small phytoplankton'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_AGG'
       vname_bgcflux(2,ind)='Aggregation of diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_SP'
       vname_bgcflux(2,ind)='PHOTOC_SP_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='f_ratio_SP'
       vname_bgcflux(2,ind)='f-ratio for small phytoplankton'
       vname_bgcflux(3,ind)='-'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_DIAT'
       vname_bgcflux(2,ind)='PHOTOC_DIAT_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='f_ratio_DIAT'
       vname_bgcflux(2,ind)='f-ratio for diatoms'
       vname_bgcflux(3,ind)='-'
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='TOT_PROD'
       vname_bgcflux(2,ind)='Net primary production'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='no3_v_sp'
       vname_bgcflux(2,ind)='NO3 uptake by SP'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='nh4_v_sp'
       vname_bgcflux(2,ind)='NH4 uptake by SP'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='no3_v_diat'
       vname_bgcflux(2,ind)='NO3 uptake by diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='nh4_v_diat'
       vname_bgcflux(2,ind)='NH4 uptake by diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# ifdef PHAEOCYSTIS
       ind=ind+1
       vname_bgcflux(1,ind)='no3_v_phaeo'
       vname_bgcflux(2,ind)='NO3 uptake by phaeocystis'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='nh4_v_phaeo'
       vname_bgcflux(2,ind)='NH4 uptake by phaeocystis'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# endif /* PHAEOCYSTIS */

       ind=ind+1
       vname_bgcflux(1,ind)='DOC_PROD'
       vname_bgcflux(2,ind)='Production of DOC'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DOC_REMIN'
       vname_bgcflux(2,ind)='Remineralization of DOD'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='FE_SCAVENGE'
       vname_bgcflux(2,ind)='Fe scavenging'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SP_N_LIM'
       vname_bgcflux(2,ind)='N limitation of small phytoplankton'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SP_FE_LIM'
       vname_bgcflux(2,ind)='Fe limitation of small phytoplankton'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SP_PO4_LIM'
       vname_bgcflux(2,ind)='PO4 limitation of small phytoplankton'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='SP_LIGHT_LIM'
       vname_bgcflux(2,ind)='Light limitation of small phytoplankton'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_N_LIM'
       vname_bgcflux(2,ind)='N limitation of diatoms'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_FE_LIM'
       vname_bgcflux(2,ind)='Fe limitation of diatoms'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_PO4_LIM'
       vname_bgcflux(2,ind)='PO4 limitation of diatoms'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_SIO3_LIM'
       vname_bgcflux(2,ind)='SiO3 limitation of diatoms'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_LIGHT_LIM'
       vname_bgcflux(2,ind)='Light limitation of diatoms'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_FORM'
       vname_bgcflux(2,ind)='CACO3_FORM_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
#ifdef PHAEOCYSTIS
       vname_bgcflux(1,ind)='PHAE_N_LIM'
       vname_bgcflux(2,ind)='PHAE_N_LIM_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_PHAE'
       vname_bgcflux(2,ind)='GRAZE_PHAE_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHAE_LOSS'
       vname_bgcflux(2,ind)='PHAE_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_PHAE'
       vname_bgcflux(2,ind)='PHOTOC_PHAE_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s   '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHAE_P_LIM'
       vname_bgcflux(2,ind)='PO4 limitation of Phaeocystis'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHAE_FE_LIM'
       vname_bgcflux(2,ind)='Fe limitation of Phaeocystis'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHAE_LIGHT_LIM'
       vname_bgcflux(2,ind)='Light limitation of Phaeocystis'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#else /* PHAEOCYSTIS */
       vname_bgcflux(1,ind)='DIAZ_NFIX'
       vname_bgcflux(2,ind)='DIAZ_NFIX_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s   '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_DIAZ'
       vname_bgcflux(2,ind)='Grazing of diazotrophs'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_LOSS_AGG'
       vname_bgcflux(2,ind)='DIAZ_LOSS_AGG_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_DIAZ'
       vname_bgcflux(2,ind)='PHOTOC_DIAZ_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_P_LIM'
       vname_bgcflux(2,ind)='PO4 limitation of diazotrophs'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_FE_LIM'
       vname_bgcflux(2,ind)='Fe limitation of diazotrophs'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_LIGHT_LIM'
       vname_bgcflux(2,ind)='Light limitation of diazotrophs'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#endif /* PHAEOCYSTIS */

       ind=ind+1
       vname_bgcflux(1,ind)='FE_SCAVENGE_RATE'
       vname_bgcflux(2,ind)='Iron scavenging rate'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DON_PROD'
       vname_bgcflux(2,ind)='Production of DON'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DON_REMIN'
       vname_bgcflux(2,ind)='Remineralization of DON'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DOFE_PROD'
       vname_bgcflux(2,ind)='Production of DOFe'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DOFE_REMIN'
       vname_bgcflux(2,ind)='Remineralization of DOFe'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DOP_PROD'
       vname_bgcflux(2,ind)='Production of DOP'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DOP_REMIN'
       vname_bgcflux(2,ind)='Remineralization of DOP'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='BSI_FORM'
       vname_bgcflux(2,ind)='BSI_FORM_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
# ifdef PHAEOCYSTIS
       vname_bgcflux(1,ind)='PHOTOFE_PHAE'
       vname_bgcflux(2,ind)='PHOTOFE_PHAE_AVG'
# else
       vname_bgcflux(1,ind)='PHOTOFE_DIAZ'
       vname_bgcflux(2,ind)='PHOTOFE_DIAZ_AVG'
# endif
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOFE_DIAT'
       vname_bgcflux(2,ind)='PHOTOFE_DIAT_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '

#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOFE_SP'
       vname_bgcflux(2,ind)='PHOTOFE_SP_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif
  
       ind=ind+1
       vname_bgcflux(1,ind)='NITRIF'
       vname_bgcflux(2,ind)='Nitrification'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='J_O2'
       vname_bgcflux(2,ind)='Oxygen sources - sinks'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# ifdef DMS_CYCLE
       ind=ind+1
       vname_bgcflux(1,ind)='DMSP_prod'
       vname_bgcflux(2,ind)='DMSP production'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DMSP_uptake'
       vname_bgcflux(2,ind)='DMSP uptake by phytoplankton'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DMSP_lysis'
       vname_bgcflux(2,ind)='DMSP lysis'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DMS_prod'
       vname_bgcflux(2,ind)='DMS production'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif


       ind=ind+1
       vname_bgcflux(1,ind)='DMS_photlys'
       vname_bgcflux(2,ind)='DMS photolysis'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif


# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
       ind=ind+1
       vname_bgcflux(1,ind)='NH4_oxid'
       vname_bgcflux(2,ind)='NH4 oxidation'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Denitrif'
       vname_bgcflux(2,ind)='Denitrification'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='N2O_prod'
       vname_bgcflux(2,ind)='N2O production'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='N2_prod'
       vname_bgcflux(2,ind)='N2 production'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

# endif /* NCYCLE_ANOXIC */

# ifdef N15_CYCLE
       ind=ind+1
       vname_bgcflux(1,ind)='NO3_v_sp_15N'
       vname_bgcflux(2,ind)='NO3_15N uptake by SP'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='NH4_v_sp_15N'
       vname_bgcflux(2,ind)='NH4_15N uptake by SP'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='NO3_v_diat_15N'
       vname_bgcflux(2,ind)='NO3_15N uptake by diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='NH4_v_diat_15N'
       vname_bgcflux(2,ind)='NH4_15N uptake by diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Graze_sp_15N'
       vname_bgcflux(2,ind)='Grazing of SP_15N'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Graze_diat_15N'
       vname_bgcflux(2,ind)='Grazing of DIAT_15N'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
# ifdef PHAEOCYSTIS
       vname_bgcflux(1,ind)='Graze_phae_15N'
       vname_bgcflux(2,ind)='Grazing of PHAE_15N'
# else
       vname_bgcflux(1,ind)='Graze_diaz_15N'
       vname_bgcflux(2,ind)='Grazing of DIAZ_15N'
# endif
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Sp_loss_15N'
       vname_bgcflux(2,ind)='SP_15N loss'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Diat_loss_15N'
       vname_bgcflux(2,ind)='DIAT_15N loss'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
# ifdef PHAEOCYSTIS
       vname_bgcflux(1,ind)='Phae_loss_15N'
       vname_bgcflux(2,ind)='PHAE_15N loss'
# else
       vname_bgcflux(1,ind)='Diaz_loss_15N'
       vname_bgcflux(2,ind)='DIAZ_15N loss'
# endif
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Zoo_loss_15N'
       vname_bgcflux(2,ind)='ZOO_15N loss'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Sp_agg_15N'
       vname_bgcflux(2,ind)='SP_15N aggregation'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Diat_agg_15N'
       vname_bgcflux(2,ind)='DIAT_15N aggregation'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
# ifdef PHAEOCYSTIS
       vname_bgcflux(1,ind)='Phaeo_15N_uptake'
       vname_bgcflux(2,ind)='Phae_15N N uptake'
# else
       vname_bgcflux(1,ind)='Diaz_Nfix_15N'
       vname_bgcflux(2,ind)='DIAZ_15N N fixation'
# endif
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='PON_remin_15N'
       vname_bgcflux(2,ind)='PON_15N N remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DON_prod_15N'
       vname_bgcflux(2,ind)='DON_15N production'
       vname_bgcflux(3,ind)=' '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='DON_remin_15N'
       vname_bgcflux(2,ind)='DON_15N remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Nitrif_15N'
       vname_bgcflux(2,ind)='Nitrification of 15N'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#  ifdef NCYCLE_ANOXIC
       ind=ind+1
       vname_bgcflux(1,ind)='NH4_oxid_15N'
       vname_bgcflux(2,ind)='NH4 oxidation of 15N'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='Denitrif_15N'
       vname_bgcflux(2,ind)='Denitrification of 15N'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='N2O_prod_15N'
       vname_bgcflux(2,ind)='N2O production of 15N'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

       ind=ind+1
       vname_bgcflux(1,ind)='N2_prod_15N'
       vname_bgcflux(2,ind)='N2 production of 15N'
       vname_bgcflux(3,ind)='mmol/m3/s  '
#ifdef REDUCED_BGC_OUTPUT
               write_bgcflux(ind) = .false.
#endif

#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */

       end subroutine init_scalars_becflux

#else /* BIOLOGY_BEC && BGC_FLUX_ANALYSIS */
      subroutine init_scalars_becflux_empty()
      return
      end
#endif  /* BIOLOGY_BEC && BGC_FLUX_ANALYSIS */
 
