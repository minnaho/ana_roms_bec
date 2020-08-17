#ifdef BIOLOGY_BEC
! these arrays are needed even without flux analysis
       real,dimension(GLOBAL_2D_ARRAY):: 
     &     PH_AVG, pCO2_AVG, pCO2air_AVG, PARinc_avg
        real,dimension(GLOBAL_2D_ARRAY,N)::  
     &    PAR_avg
        common /time_avg/
     &    PH_AVG, pCO2_AVG, pCO2air_AVG, PARinc_avg,
     &    PAR_avg
#ifdef OUTPUT_SPECDOM
        real,dimension(GLOBAL_2D_ARRAY) :: 
     &       dom_sp_sfc_avg, dom_diat_sfc_avg, dom_diaz_sfc_avg, 
     &       dom_sp_int_avg, dom_diat_int_avg, dom_diaz_int_avg
        common /specdom_avg/ 
     &       dom_sp_sfc_avg, dom_diat_sfc_avg, dom_diaz_sfc_avg, 
     &       dom_sp_int_avg, dom_diat_int_avg, dom_diaz_int_avg
#endif
       
#if defined EXPLICIT_SINKING || defined SEDIMENT_BIOLOGY
      real t_sed_avg(GLOBAL_2D_ARRAY,NT_sed)
      common /ocean_t_sed_avg/t_sed_avg
#endif

#ifdef BGC_FLUX_ANALYSIS
       real,dimension(GLOBAL_2D_ARRAY):: WS_AVG, XKW_AVG, 
     &     AP_AVG, SCHMIDT_O2_AVG, O2SAT_AVG, FG_O2_AVG, 
     &     SCHMIDT_CO2_AVG, CO2STAR_AVG, DCO2STAR_AVG, 
     &     FG_CO2_AVG, IRON_FLUX_AVG, 
     &     nitrate_FLUX_AVG,ammonium_FLUX_AVG,alkatm_FLUX_AVG,
     &     anthpco2_FLUX_AVG,
     &     PARinc_flux_avg, zeta_bgc_flux_avg
# ifdef DMS_CYCLE
     &     , fg_dms_avg
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
     &     , fg_n2o_avg, fg_n2_avg
#  ifdef N15_CYCLE
     &     , fg_n2o_n15_avg, fg_n2_n15_avg
#  endif
# endif /* NCYCLE_ANOXIC */
        real,dimension(GLOBAL_2D_ARRAY,N)::  
     &    PO4_RESTORE_AVG, NO3_RESTORE_AVG, 
     &    SiO3_RESTORE_AVG, PAR_flux_avg, PO4STAR_AVG, 
     &    POC_FLUX_IN_AVG, POC_PROD_AVG, POC_REMIN_AVG, 
     &    CaCO3_FLUX_IN_AVG, CaCO3_PROD_AVG, 
     &    CaCO3_REMIN_AVG,  SiO2_FLUX_IN_AVG, 
     &    SiO2_PROD_AVG, SiO2_REMIN_AVG, dust_FLUX_IN_AVG, 
     &    dust_REMIN_AVG, P_iron_FLUX_IN_AVG, 
     &    P_iron_PROD_AVG, P_iron_REMIN_AVG, 
     &    graze_sp_AVG, graze_diat_AVG, graze_tot_AVG, 
     &    sp_loss_AVG, diat_loss_AVG, zoo_loss_AVG, 
     &    sp_agg_AVG, diat_agg_AVG, 
     &    photoC_sp_AVG, f_ratio_sp_avg,
     &    photoC_diat_AVG, f_ratio_diat_avg, tot_prod_AVG, 
     &    no3_v_sp_avg, nh4_v_sp_avg, 
     &    no3_v_diat_avg, nh4_v_diat_avg,
#ifdef PHAEOCYSTIS
     &    no3_v_phaeo_avg, nh4_v_phaeo_avg,
#endif
     &    DOC_prod_AVG, DOC_remin_AVG, Fe_scavenge_AVG, 
     &    sp_N_lim_AVG, sp_Fe_lim_AVG, sp_PO4_lim_AVG, 
     &    sp_light_lim_AVG, diat_N_lim_AVG, diat_Fe_lim_AVG, 
     &    diat_PO4_lim_AVG, diat_SiO3_lim_AVG
        real,dimension(GLOBAL_2D_ARRAY,N)::
     &    diat_light_lim_AVG, CaCO3_form_AVG, 
     &    diaz_Nfix_AVG, graze_diaz_AVG, diaz_loss_AVG,
     &     photoC_diaz_AVG, diaz_P_lim_AVG, 
     &    diaz_Fe_lim_AVG, diaz_light_lim_AVG, 
     &     Fe_scavenge_rate_AVG, DON_prod_AVG, 
     &    DON_remin_AVG, DOFe_prod_AVG, 
     &    DOFe_remin_AVG, DOP_prod_AVG, 
     &    DOP_remin_AVG, bSI_form_AVG, 
     &    photoFe_diaz_AVG, photoFe_diat_AVG, 
     &    photoFe_sp_AVG,nitrif_AVG,j_o2_avg
# ifdef DMS_CYCLE
         real,dimension(GLOBAL_2D_ARRAY,N) :: dmsp_prod_avg,
     &       dmsp_p_uptake_avg, dmsp_lysis_avg,
     &       dms_prod_avg, dms_photlys_avg
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
         real,dimension(GLOBAL_2D_ARRAY,N) ::  ammox_avg,
     &        denitr_avg, n2o_prod_avg, n2_prod_avg
         real,dimension(GLOBAL_2D_ARRAY) :: denitr_sed_avg
# endif /* NCYCLE_ANOXIC */
# ifdef N15_CYCLE
         real,dimension(GLOBAL_2D_ARRAY,N) ::  no3_v_sp_n15_avg,
     &        nh4_v_sp_n15_avg,
     &        no3_v_diat_n15_avg, nh4_v_diat_n15_avg,
     &        graze_sp_n15_avg, graze_diat_n15_avg, graze_diaz_n15_avg,
     &        sp_loss_n15_avg, diat_loss_n15_avg, 
     &        diaz_loss_n15_avg, zoo_loss_n15_avg,
     &        sp_agg_n15_avg, diat_agg_n15_avg,
     &        diaz_Nfix_n15_avg, pon_remin_n15_avg,    
     &        don_prod_n15_avg, don_remin_n15_avg,
     &        nitrif_n15_avg
#  ifdef NCYCLE_ANOXIC
     &      , ammox_n15_avg, denitr_n15_avg
     &        , n2o_prod_n15_avg, n2_prod_n15_avg
         real,dimension(GLOBAL_2D_ARRAY) :: denitr_sed_n15_avg
#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */ 

       common /time_avg1/WS_AVG, XKW_AVG, 
     &    AP_AVG, SCHMIDT_O2_AVG, O2SAT_AVG, FG_O2_AVG, 
     &    SCHMIDT_CO2_AVG, CO2STAR_AVG, DCO2STAR_AVG, 
     &    FG_CO2_AVG, IRON_FLUX_AVG, 
     &    nitrate_FLUX_AVG,ammonium_FLUX_AVG,alkatm_FLUX_AVG,
     &    anthpco2_FLUX_AVG,
     &    PARinc_flux_avg, zeta_bgc_flux_avg,
     &    PO4_RESTORE_AVG, NO3_RESTORE_AVG, 
     &    SiO3_RESTORE_AVG, PAR_flux_avg, PO4STAR_AVG, 
     &    POC_FLUX_IN_AVG, POC_PROD_AVG, POC_REMIN_AVG, 
     &    CaCO3_FLUX_IN_AVG, CaCO3_PROD_AVG, 
     &    CaCO3_REMIN_AVG,  SiO2_FLUX_IN_AVG, 
     &    SiO2_PROD_AVG, SiO2_REMIN_AVG, dust_FLUX_IN_AVG, 
     &    dust_REMIN_AVG, P_iron_FLUX_IN_AVG, 
     &    P_iron_PROD_AVG, P_iron_REMIN_AVG, 
     &    graze_sp_AVG, graze_diat_AVG, graze_tot_AVG, 
     &    sp_loss_AVG, diat_loss_AVG, zoo_loss_AVG 
# ifdef DMS_CYCLE
     &     , fg_dms_avg
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
     &     , fg_n2o_avg, fg_n2_avg
#  ifdef N15_CYCLE
     &     , fg_n2o_n15_avg, fg_n2_n15_avg
#  endif
# endif /* NCYCLE_ANOXIC */

       common /time_avg2/
     &    sp_agg_AVG, diat_agg_AVG, 
     &    photoC_sp_AVG, f_ratio_sp_avg,
     &    photoC_diat_AVG, f_ratio_diat_avg, tot_prod_AVG, 
     &    no3_v_sp_avg, nh4_v_sp_avg, 
     &    no3_v_diat_avg, nh4_v_diat_avg,
#ifdef PHAEOCYSTIS
     &    no3_v_phaeo_avg, nh4_v_phaeo_avg,
#endif
     &    DOC_prod_AVG, DOC_remin_AVG, Fe_scavenge_AVG, 
     &    sp_N_lim_AVG, sp_Fe_lim_AVG, sp_PO4_lim_AVG, 
     &    sp_light_lim_AVG, diat_N_lim_AVG, diat_Fe_lim_AVG, 
     &    diat_PO4_lim_AVG, diat_SiO3_lim_AVG, 
     &    diat_light_lim_AVG, CaCO3_form_AVG, 
     &    diaz_Nfix_AVG, graze_diaz_AVG, diaz_loss_AVG,
     &     photoC_diaz_AVG, diaz_P_lim_AVG, 
     &    diaz_Fe_lim_AVG, diaz_light_lim_AVG, 
     &     Fe_scavenge_rate_AVG, DON_prod_AVG, 
     &    DON_remin_AVG, DOFe_prod_AVG, 
     &    DOFe_remin_AVG, DOP_prod_AVG, 
     &    DOP_remin_AVG, bSI_form_AVG, 
     &    photoFe_diaz_AVG, photoFe_diat_AVG, 
     &    photoFe_sp_AVG, nitrif_AVG, j_o2_avg
# ifdef DMS_CYCLE
     &      , dmsp_prod_avg, dmsp_p_uptake_avg, dmsp_lysis_avg, 
     &      dms_prod_avg, dms_photlys_avg
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
     &      , ammox_avg, denitr_avg, n2o_prod_avg, n2_prod_avg, 
     &      denitr_sed_avg
# endif /* NCYCLE_ANOXIC */
# ifdef N15_CYCLE
     &      , no3_v_sp_n15_avg, nh4_v_sp_n15_avg, 
     &      no3_v_diat_n15_avg, nh4_v_diat_n15_avg, 
     &      graze_sp_n15_avg, graze_diat_n15_avg, graze_diaz_n15_avg,
     &      sp_loss_n15_avg, diat_loss_n15_avg, 
     &      diaz_loss_n15_avg, zoo_loss_n15_avg,
     &      sp_agg_n15_avg, diat_agg_n15_avg,
     &      diaz_Nfix_n15_avg, pon_remin_n15_avg,    
     &      don_prod_n15_avg, don_remin_n15_avg, 
     &      nitrif_n15_avg
#  ifdef NCYCLE_ANOXIC
     &      , ammox_n15_avg, denitr_n15_avg, 
     &      n2o_prod_n15_avg, n2_prod_n15_avg, denitr_sed_n15_avg
#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */ 
# ifdef SEDIMENT_BIOLOGY
       real bot_flux_poc_avg(GLOBAL_2D_ARRAY), 
     &      bot_flux_caco3_avg(GLOBAL_2D_ARRAY),
     &      bot_flux_si_avg(GLOBAL_2D_ARRAY)
#  ifdef N15_CYCLE
     &      , bot_flux_pon15_avg(GLOBAL_2D_ARRAY)
#  endif
       common /sed_flux/ bot_flux_poc_avg, bot_flux_caco3_avg,
     &      bot_flux_si_avg
#  ifdef N15_CYCLE
     &      , bot_flux_pon15_avg
#  endif
# endif /* SEDIMENT_BIOLOGY */
#endif /* BGC_FLUX_ANALYSIS */
#endif /* BIOLOGY_BEC */
