#ifdef BIOLOGY_BEC
!zhc: This file includes parameters which were originally defined in ecosys_parms.F90 from WHOI

        real c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1., c0=0.0,c2=2., 
     &  c1000=1000.,p5=0.5,
     &  spd = 86400.0,  dps = c1 / spd ,   
     &  t0_kelvin= 273.16)
       
        real  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C 
        parameter ( parm_Red_D_C_P  = 117.0,
     &  parm_Red_P_C_P  = 117.0,     
     &  parm_Red_D_C_N  = 117.0 / 16.0,  
     &  parm_Red_P_C_N  = 117.0 / 16.0,   
     &  parm_Red_D_C_O2 = 117.0 / 150.0, 
     &  parm_Red_P_C_O2 = 117.0 / 150.0,
     &  parm_Red_Fe_C   = 3.0e-6)

        real, parameter :: Q         = 0.137 ! fixed N/C ratio for all plankton
        real, parameter :: Qinv      = 1.0/0.137 ! fixed C/N ratio for all plankton
#ifdef NP16_SP_DIAT
        real, parameter :: Qp        = 0.0625 * Q ! fixed P/C ratio for all plankton except diazotrophs
        real, parameter :: Qp_diat   = Qp ! fixed P/C ratio for diatoms
#else
        real, parameter :: Qp        = 0.00855 ! fixed P/C ratio for all plankton except diazotrophs
!mod Qp        real, parameter :: Qp        = 0.05 * Q ! fixed P/C ratio for small plankton
        real, parameter :: Qp_diat   = 0.1 * Q ! fixed P/C ratio for diatoms
#endif
        real, parameter :: parm_n2o_nev_a1 = 0.26 ! parameter for nevisson parameterization
        real, parameter :: parm_n2o_nev_a2 =  -0.0004 ! parameter for nevisson parameterization

#ifdef NCYCLE_ANOXIC
        real, parameter :: Qdenit    = 0.9811320755 ! 104:106 N/C ratio for denitrification
#endif
!       parm_Fe_bioavail,        fraction of Fe flux that is bioavailable
!       parm_prod_dissolve,      frac. of prod -> DOC
!       parm_o2_min,             lower limit of O2 for prod  consumption (mmol/m^3)
!       parm_Rain_CaCO3,         Rain ratio for CaCO3
!       parm_Rain_SiO2,          Rain ratio for SiO2
!       parm_kappa_nitrif,       nitrification inverse time constant (1/sec)
!       parm_nitrif_par_lim,     PAR limit for nitrif. (W/m^2)
!       parm_POC_flux_ref,       reference POC flux (mmol C/m^2/sec)
!       parm_rest_prod_tau,      time-scale for restoring prod (sec)
!       parm_rest_prod_z_c,      depth-limit for restoring (m)
!       parm_z_umax_0,           max. zoopl growth rate on sphyto at tref (1/sec)
!       parm_diat_umax_0,        max. zoopl growth rate on diatoms at tref (1/sec)
!       parm_z_mort_0,           zoopl linear mort rate (1/sec)
!       parm_z_mort2_0,          zoopl quad mort rate, higher trophic level grazing (1/sec/((mmol C/m3))
!       parm_sd_remin_0,         small detrital remineralization rate (1/sec)
!       parm_sp_kNO3,            small phyto nitrate uptake half saturation coeff. (mmol N/m3)
!       parm_diat_kNO3,          diatom nitrate uptake half saturation coeff. (mmol N/m3)
!       parm_sp_kNH4,            small phyto ammonium uptake half saturation coeff. (mmol N/m3)
!       parm_diat_kNH4,          diatom ammonium uptake half saturation coeff. (mmol N/m3)
#ifdef PHAEOCYSTIS
!       parm_phaeo_kNO3,         phaeocystis nitrate uptake half saturation coeff. (mmol N/m3)
!       parm_phaeo_kNH4,         phaeocystis ammonium uptake half saturation coeff. (mmol N/m3)
#endif
!       parm_sp_kFe,             small phyto iron uptake half saturation coefficient (mmol Fe/m3)
!       parm_diat_kFe,           diatom iron uptake half saturation coefficient (mmol Fe/m3)
!       parm_diat_kSiO3,         diatom si uptake half saturation coefficient (mmol SiO3/m3)
!       parm_sp_kPO4,            small phyto PO4 uptake (mmol P/m^3)
!       parm_diat_kPO4,          diatom PO4 uptate (mmol P/m^3)
!       parm_z_grz,              grazing coefficient for small phyto (mmol C/m^3)
!       parm_alphaChl,           Chl. specific initial slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
!       parm_labile_ratio,       portion of loss to DOC that routed directly to DIC (non-dimensional)
!       parm_alphaDiaz,          chl. spec. init. slope of P_I curve for diazotrophs
!       parm_diaz_umax_0         max. zoopl growth rate on diazotrophs at tre        

        real  parm_Fe_bioavail,   parm_prod_dissolve,   
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,  
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,    
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,       
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,         
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,         
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,            
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,      
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,       
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
#ifdef PHAEOCYSTIS
     &  ,parm_phaeo_kNO3, parm_phaeo_kNH4  
#endif
# ifdef SEDIMENT_BIOLOGY
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
#  ifdef N15_CYCLE
     &  , t_remin_sed_pon15
#  endif
# endif /* SEDIMENT_BIOLOGY */

       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve, 
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2, 
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref, 
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0, 
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0, 
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe, 
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4, 
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl, 
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
#ifdef PHAEOCYSTIS
     &  ,parm_phaeo_kNO3, parm_phaeo_kNH4
#endif
# ifdef SEDIMENT_BIOLOGY
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
#  ifdef N15_CYCLE
     &  , t_remin_sed_pon15
#  endif
# endif /* SEDIMENT_BIOLOGY */

#ifdef DMS_CYCLE
       real :: parm_I_min_DMS  ! W m-2
       real :: parm_I_max_DMS  ! W m-2
       real :: parm_I_mean_DMS ! W m-2
       real :: parm_inv_I_del_DMS ! 1 / (W m-2)

       real :: parm_r_min_DMS  ! mmol S / mmol N
       real :: parm_r_max_DMS  ! mmol S / mmol N
       real :: parm_r_del_DMS  ! mmol S / mmol N

       real :: parm_gamma1     ! [adim]

       real :: parm_f_lysis_DMSP ! s-1

       real :: parm_gamma_s_max_DMS ! s-1

       real :: parm_k_phot_DMS ! s-1

       real :: parm_alpha1     ! [adim]

       real :: parm_alpha_p    ! [adim]
       common/eco_param_dms/ parm_I_min_DMS, parm_I_max_DMS, 
     &      parm_I_mean_DMS, parm_inv_I_del_DMS,
     &      parm_r_min_DMS, parm_r_max_DMS, parm_r_del_DMS,
     &      parm_gamma1, parm_f_lysis_DMSP,
     &      parm_gamma_s_max_DMS, parm_k_phot_DMS, parm_alpha1,
     &      parm_alpha_p
#endif /* DMS_CYCLE */

#ifdef NCYCLE_ANOXIC
        real parm_kappa_nitrif_no2
        real parm_n2o_prod
        real parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
#endif /* NCYCLE_ANOXIC */

#ifdef N15_CYCLE
        real parm_n15_std_fraction
        real parm_alpha_n15_no3_v_sp
        real parm_alpha_n15_nh4_v_sp
        real parm_alpha_n15_no3_v_diat
        real parm_alpha_n15_nh4_v_diat
# ifdef PHAEOCYSTIS
        real parm_alpha_n15_no3_v_phaeo
        real parm_alpha_n15_nh4_v_phaeo
# endif
        real parm_alpha_n15_nitrif
        real parm_alpha_n15_don_remin
        real parm_alpha_n15_don_prod
        real parm_alpha_n15_graze_sp_zoo
        real parm_alpha_n15_graze_sp_pon
        real parm_alpha_n15_graze_sp_don
        real parm_alpha_n15_graze_sp_nh4
        real parm_alpha_n15_sp_loss_pon
        real parm_alpha_n15_sp_loss_don
        real parm_alpha_n15_sp_loss_nh4
        real parm_alpha_n15_sp_agg
        real parm_alpha_n15_graze_diat_zoo
        real parm_alpha_n15_graze_diat_pon
        real parm_alpha_n15_graze_diat_don
        real parm_alpha_n15_graze_diat_nh4
        real parm_alpha_n15_diat_loss_pon
        real parm_alpha_n15_diat_loss_don
        real parm_alpha_n15_diat_loss_nh4
        real parm_alpha_n15_diat_agg
        real parm_alpha_n15_graze_diaz_zoo
        real parm_alpha_n15_graze_diaz_pon
        real parm_alpha_n15_graze_diaz_don
        real parm_alpha_n15_graze_diaz_nh4
        real parm_alpha_n15_diaz_agg
        real parm_alpha_n15_diaz_loss_don
        real parm_alpha_n15_diaz_loss_nh4
        real parm_alpha_n15_diaz_Nfix
        real parm_alpha_n15_diaz_Nexcrete
        real parm_alpha_n15_zoo_loss_don
        real parm_alpha_n15_zoo_loss_nh4
        real parm_alpha_n15_pon_remin
        real parm_alpha_n15_n2o_gasex
        real parm_alpha_n15_n2_gasex
# ifdef NCYCLE_ANOXIC
        real parm_alpha_n15_ammox
        real parm_alpha_n15_denitr
        real parm_alpha_n15_denitr_sed
        real parm_alpha_n15_n2o_prod
        real parm_alpha_n15_n2_prod
# endif /* NCYCLE_ANOXIC */
        common /eco_param_n15_cycle/ parm_n15_std_fraction,
     &       parm_alpha_n15_no3_v_sp, parm_alpha_n15_nh4_v_sp,
     &       parm_alpha_n15_no3_v_diat, parm_alpha_n15_nh4_v_diat,
# ifdef PHAEOCYSTIS
     &       parm_alpha_n15_no3_v_phaeo, parm_alpha_n15_nh4_v_phaeo,
# endif
     &       parm_alpha_n15_nitrif, parm_alpha_n15_don_remin,
     &       parm_alpha_n15_don_prod, 
     &       parm_alpha_n15_graze_sp_zoo, parm_alpha_n15_graze_sp_pon,
     &       parm_alpha_n15_graze_sp_don, parm_alpha_n15_graze_sp_nh4,
     &       parm_alpha_n15_sp_loss_pon,
     &       parm_alpha_n15_sp_loss_don, parm_alpha_n15_sp_loss_nh4,
     &       parm_alpha_n15_sp_agg,
     &       parm_alpha_n15_graze_diat_zoo, 
     &       parm_alpha_n15_graze_diat_pon,
     &       parm_alpha_n15_graze_diat_don, 
     &       parm_alpha_n15_graze_diat_nh4,
     &       parm_alpha_n15_diat_loss_pon,
     &       parm_alpha_n15_diat_loss_don, 
     &       parm_alpha_n15_diat_loss_nh4,
     &       parm_alpha_n15_diat_agg,
     &       parm_alpha_n15_graze_diaz_zoo, 
     &       parm_alpha_n15_graze_diaz_pon,
     &       parm_alpha_n15_graze_diaz_don, 
     &       parm_alpha_n15_graze_diaz_nh4,
     &       parm_alpha_n15_diaz_agg, 
     &       parm_alpha_n15_diaz_loss_don, 
     &       parm_alpha_n15_diaz_loss_nh4,
     &       parm_alpha_n15_diaz_Nfix,
     &       parm_alpha_n15_diaz_Nexcrete,
     &       parm_alpha_n15_zoo_loss_don, 
     &       parm_alpha_n15_zoo_loss_nh4,
     &       parm_alpha_n15_pon_remin,
     &       parm_alpha_n15_n2o_gasex,
     &       parm_alpha_n15_n2_gasex
# ifdef NCYCLE_ANOXIC
     &       , parm_alpha_n15_ammox, parm_alpha_n15_denitr,
     &       parm_alpha_n15_denitr_sed, 
     &       parm_alpha_n15_n2o_prod, parm_alpha_n15_n2_prod
# endif /* NCYCLE_ANOXIC */
#endif /* N15_CYCLE */
#endif /* BIOLOGY_BEC */
