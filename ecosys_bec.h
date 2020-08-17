#ifdef BIOLOGY_BEC
!zhc:  This file includes variables originally defined in ecosys_mod.F90 from WHOI
     
       real tracer(GLOBAL_I_ARRAY,N,ntrc_bio,2)

        common /tracers/ tracer

# ifdef SEDIMENT_BIOLOGY
        real tracer_sed(GLOBAL_I_ARRAY,NT_sed,2)
        common /tracer_sed/ tracer_sed
# endif


        real ifrac(GLOBAL_2D_ARRAY),
     &    press(GLOBAL_2D_ARRAY)
!    IFRAC  sea ice fraction (non-dimensional)
!    PRESS  sea level atmospheric pressure (Pascals)
        common /fic_ap/ifrac,press

!--------------------------------------------------------------------------
!   variables used for time-averaging
!--------------------------------------------------------------------------

! these arrays are needed even without flux analysis
        real PH_HIST(GLOBAL_2D_ARRAY),
     &   pCO2sw(GLOBAL_2D_ARRAY),
     &   pCO2air(GLOBAL_2D_ARRAY) , 
     &   PAR(GLOBAL_2D_ARRAY,N),
     &   PARinc(GLOBAL_2D_ARRAY)
#ifdef BULK_FLUX
     &   ,PARinc_rst(GLOBAL_2D_ARRAY)
#endif 
       common /time_averaging/
     &       ph_hist,pCO2sw,pCO2air,
     &       PAR,PARinc
#ifdef BULK_FLUX
     &  , PARinc_rst
#endif 
#ifdef OUTPUT_SPECDOM
        real,dimension(GLOBAL_2D_ARRAY) :: dom_sp_sfc, dom_diat_sfc,
     &       dom_diaz_sfc, dom_sp_int, dom_diat_int, dom_diaz_int,
     &       spchl_int, diatchl_int, diazchl_int
        common /specdom/ dom_sp_sfc, dom_diat_sfc,
     &       dom_diaz_sfc, dom_sp_int, dom_diat_int, dom_diaz_int,
     &       spchl_int, diatchl_int, diazchl_int
#endif

#ifdef BGC_FLUX_ANALYSIS
       real WS_HIST(GLOBAL_2D_ARRAY), 
     &   XKW_HIST(GLOBAL_2D_ARRAY), 
     &   AP_HIST(GLOBAL_2D_ARRAY), 
     &   SCHMIDT_O2_HIST(GLOBAL_2D_ARRAY), 
     &   O2SAT_HIST(GLOBAL_2D_ARRAY), 
     &   FG_O2_HIST(GLOBAL_2D_ARRAY), 
     &    SCHMIDT_CO2_HIST(GLOBAL_2D_ARRAY), 
     &   CO2STAR_HIST(GLOBAL_2D_ARRAY),  
     &   DCO2STAR_HIST(GLOBAL_2D_ARRAY), 
     &    FG_CO2_HIST(GLOBAL_2D_ARRAY), 
     &   IRON_FLUX_HIST(GLOBAL_2D_ARRAY),
     &   nitrate_FLUX_HIST(GLOBAL_2D_ARRAY),
     &   ammonium_FLUX_HIST(GLOBAL_2D_ARRAY),
     &   alkatm_FLUX_HIST(GLOBAL_2D_ARRAY),
     &   anthpco2_FLUX_HIST(GLOBAL_2D_ARRAY)
# ifdef DMS_CYCLE
       real fg_dms_hist(GLOBAL_2D_ARRAY)
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
       real fg_n2o_hist(GLOBAL_2D_ARRAY)
       real fg_n2_hist(GLOBAL_2D_ARRAY)
#  ifdef N15_CYCLE
       real fg_n2o_n15_hist(GLOBAL_2D_ARRAY)
       real fg_n2_n15_hist(GLOBAL_2D_ARRAY)
#  endif /* N15_CYCLE */
# endif /* NCYCLE_ANOXIC */
        real 
     &    PO4_RESTORE_HIST(GLOBAL_2D_ARRAY,N), 
     &    NO3_RESTORE_HIST(GLOBAL_2D_ARRAY,N), 
     &    SiO3_RESTORE_HIST(GLOBAL_2D_ARRAY,N), 
     &   PO4STAR_HIST(GLOBAL_2D_ARRAY,N), 
     &    POC_FLUX_IN_HIST(GLOBAL_2D_ARRAY,N), 
     &    POC_PROD_HIST(GLOBAL_2D_ARRAY,N), 
     &   POC_REMIN_HIST(GLOBAL_2D_ARRAY,N), 
     &    CaCO3_FLUX_IN_HIST(GLOBAL_2D_ARRAY,N), 
     &    CaCO3_PROD_HIST(GLOBAL_2D_ARRAY,N), 
     &    CaCO3_REMIN_HIST(GLOBAL_2D_ARRAY,N),  
     &    SiO2_FLUX_IN_HIST(GLOBAL_2D_ARRAY,N), 
     &    SiO2_PROD_HIST(GLOBAL_2D_ARRAY,N), 
     &    SiO2_REMIN_HIST(GLOBAL_2D_ARRAY,N), 
     &    dust_FLUX_IN_HIST(GLOBAL_2D_ARRAY,N)
# ifdef N15_CYCLE
     &   , PON15_PROD_HIST(GLOBAL_2D_ARRAY,N), 
     &   PON15_REMIN_HIST(GLOBAL_2D_ARRAY,N)
# endif        
          real  dust_REMIN_HIST(GLOBAL_2D_ARRAY,N), 
     &    P_iron_FLUX_IN_HIST(GLOBAL_2D_ARRAY,N), 
     &    P_iron_PROD_HIST(GLOBAL_2D_ARRAY,N), 
     &    P_iron_REMIN_HIST(GLOBAL_2D_ARRAY,N), 
     &    graze_sp_HIST(GLOBAL_2D_ARRAY,N), 
     &    graze_diat_HIST(GLOBAL_2D_ARRAY,N), 
     &    graze_tot_HIST(GLOBAL_2D_ARRAY,N), 
     &    sp_loss_HIST(GLOBAL_2D_ARRAY,N), 
     &    diat_loss_HIST(GLOBAL_2D_ARRAY,N), 
     &    zoo_loss_HIST(GLOBAL_2D_ARRAY,N), 
     &    sp_agg_HIST(GLOBAL_2D_ARRAY,N), 
     &    diat_agg_HIST(GLOBAL_2D_ARRAY,N) 

          real  photoC_sp_HIST(GLOBAL_2D_ARRAY,N), 
     &    f_ratio_sp_hist(GLOBAL_2D_ARRAY,N),
     &    photoC_diat_HIST(GLOBAL_2D_ARRAY,N), 
     &    f_ratio_diat_hist(GLOBAL_2D_ARRAY,N),
     &    tot_prod_HIST(GLOBAL_2D_ARRAY,N), 
     &    no3_v_sp_hist(GLOBAL_2D_ARRAY,N), 
     &    nh4_v_sp_hist(GLOBAL_2D_ARRAY,N), 
     &    no3_v_diat_hist(GLOBAL_2D_ARRAY,N), 
     &    nh4_v_diat_hist(GLOBAL_2D_ARRAY,N), 
#ifdef PHAEOCYSTIS
     &    no3_v_phaeo_hist(GLOBAL_2D_ARRAY,N), 
     &    nh4_v_phaeo_hist(GLOBAL_2D_ARRAY,N), 
#endif
     &    DOC_prod_HIST(GLOBAL_2D_ARRAY,N), 
     &    DOC_remin_HIST(GLOBAL_2D_ARRAY,N), 
     &    Fe_scavenge_HIST(GLOBAL_2D_ARRAY,N), 
     &    sp_N_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &    sp_Fe_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &     sp_PO4_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &    sp_light_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &    diat_N_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &     diat_Fe_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &    diat_PO4_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &    diat_SiO3_lim_HIST(GLOBAL_2D_ARRAY,N)
       
         real diat_light_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &   CaCO3_form_HIST(GLOBAL_2D_ARRAY,N), 
     &    diaz_Nfix_HIST(GLOBAL_2D_ARRAY,N), 
     &    graze_diaz_HIST(GLOBAL_2D_ARRAY,N), 
     &    diaz_loss_HIST(GLOBAL_2D_ARRAY,N)
        real photoC_diaz_HIST(GLOBAL_2D_ARRAY,N), 
     &    diaz_P_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &    diaz_Fe_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &     diaz_light_lim_HIST(GLOBAL_2D_ARRAY,N), 
     &     Fe_scavenge_rate_HIST(GLOBAL_2D_ARRAY,N), 
     &    DON_prod_HIST(GLOBAL_2D_ARRAY,N), 
     &    DON_remin_HIST(GLOBAL_2D_ARRAY,N), 
     &    DOFe_prod_HIST(GLOBAL_2D_ARRAY,N)
   
         real  DOFe_remin_HIST(GLOBAL_2D_ARRAY,N), 
     &   DOP_prod_HIST(GLOBAL_2D_ARRAY,N), 
     &    DOP_remin_HIST(GLOBAL_2D_ARRAY,N), 
     &    bSI_form_HIST(GLOBAL_2D_ARRAY,N), 
     &    photoFe_diaz_HIST(GLOBAL_2D_ARRAY,N), 
     &    photoFe_diat_HIST(GLOBAL_2D_ARRAY,N), 
     &    photoFe_sp_HIST(GLOBAL_2D_ARRAY,N),
     &    nitrif_HIST(GLOBAL_2D_ARRAY,N),
     &    j_o2_hist(GLOBAL_2D_ARRAY,N)  

# ifdef DMS_CYCLE
         real dmsp_prod_hist(GLOBAL_2D_ARRAY,N),
     &        dmsp_p_uptake_hist(GLOBAL_2D_ARRAY,N),
     &        dmsp_lysis_hist(GLOBAL_2D_ARRAY,N),
     &        dms_prod_hist(GLOBAL_2D_ARRAY,N),
     &        dms_photlys_hist(GLOBAL_2D_ARRAY,N)
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
         real ammox_hist(GLOBAL_2D_ARRAY,N),
     &        denitr_hist(GLOBAL_2D_ARRAY,N),
     &        n2o_prod_hist(GLOBAL_2D_ARRAY,N),
     &        n2_prod_hist(GLOBAL_2D_ARRAY,N),
     &        denitr_sed_hist(GLOBAL_2D_ARRAY)
# endif /* NCYCLE_ANOXIC */
# ifdef N15_CYCLE
         real no3_v_sp_n15_hist(GLOBAL_2D_ARRAY,N)
         real nh4_v_sp_n15_hist(GLOBAL_2D_ARRAY,N)
         real no3_v_diat_n15_hist(GLOBAL_2D_ARRAY,N)
         real nh4_v_diat_n15_hist(GLOBAL_2D_ARRAY,N)
         real graze_sp_n15_hist(GLOBAL_2D_ARRAY,N)
         real graze_diat_n15_hist(GLOBAL_2D_ARRAY,N)
         real graze_diaz_n15_hist(GLOBAL_2D_ARRAY,N)
         real sp_loss_n15_hist(GLOBAL_2D_ARRAY,N)
         real diat_loss_n15_hist(GLOBAL_2D_ARRAY,N)
         real diaz_loss_n15_hist(GLOBAL_2D_ARRAY,N)
         real zoo_loss_n15_hist(GLOBAL_2D_ARRAY,N)
         real nitrif_n15_hist(GLOBAL_2D_ARRAY,N)
         real sp_agg_n15_hist(GLOBAL_2D_ARRAY,N)
         real diat_agg_n15_hist(GLOBAL_2D_ARRAY,N)
         real diaz_Nfix_n15_hist(GLOBAL_2D_ARRAY,N)
         real pon_remin_n15_hist(GLOBAL_2D_ARRAY,N)
         real don_remin_n15_hist(GLOBAL_2D_ARRAY,N)
         real don_prod_n15_hist(GLOBAL_2D_ARRAY,N)
#  ifdef NCYCLE_ANOXIC
         real ammox_n15_hist(GLOBAL_2D_ARRAY,N),
     &        denitr_n15_hist(GLOBAL_2D_ARRAY,N),
     &        n2o_prod_n15_hist(GLOBAL_2D_ARRAY,N),
     &        n2_prod_n15_hist(GLOBAL_2D_ARRAY,N),
     &        denitr_sed_n15_hist(GLOBAL_2D_ARRAY)
#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */
       common /time_averaging1/WS_HIST, XKW_HIST, 
     &   AP_HIST, SCHMIDT_O2_HIST, O2SAT_HIST, FG_O2_HIST, 
     &    SCHMIDT_CO2_HIST, CO2STAR_HIST, 
     &    DCO2STAR_HIST,
     &    FG_CO2_HIST, IRON_FLUX_HIST,
     &    PO4_RESTORE_HIST, NO3_RESTORE_HIST, 
     &    SiO3_RESTORE_HIST, PO4STAR_HIST, 
     &    POC_FLUX_IN_HIST, POC_PROD_HIST, POC_REMIN_HIST, 
     &    CaCO3_FLUX_IN_HIST, CaCO3_PROD_HIST, 
     &    CaCO3_REMIN_HIST,  SiO2_FLUX_IN_HIST, 
     &    SiO2_PROD_HIST, SiO2_REMIN_HIST, dust_FLUX_IN_HIST, 
# ifdef N15_CYCLE
     &    PON15_PROD_HIST, PON15_REMIN_HIST,
# endif        
     &    dust_REMIN_HIST, P_iron_FLUX_IN_HIST, 
     &    P_iron_PROD_HIST, P_iron_REMIN_HIST, 
     &    graze_sp_HIST, graze_diat_HIST, graze_tot_HIST, 
     &    sp_loss_HIST, diat_loss_HIST, zoo_loss_HIST, 
     &    sp_agg_HIST, diat_agg_HIST, 
     &    photoC_sp_HIST, f_ratio_sp_hist,
     &    photoC_diat_HIST, f_ratio_diat_hist, tot_prod_HIST, 
     &    no3_v_sp_hist, nh4_v_sp_hist, 
     &    no3_v_diat_hist, nh4_v_diat_hist,
#ifdef PHAEOCYSTIS
     &    no3_v_phaeo_hist, nh4_v_phaeo_hist,
#endif
     &    DOC_prod_HIST, DOC_remin_HIST, Fe_scavenge_HIST 
# ifdef DMS_CYCLE
     &        , fg_dms_hist
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
     &        , fg_n2o_hist
     &        , fg_n2_hist
     &        , denitr_sed_hist
#  ifdef N15_CYCLE
     &        , fg_n2o_n15_hist
     &        , fg_n2_n15_hist
     &        , denitr_sed_n15_hist
#  endif /* N15_CYCLE */
# endif /* NCYCLE_ANOXIC */

       common /time_averaging2/
     &    sp_N_lim_HIST, sp_Fe_lim_HIST, sp_PO4_lim_HIST, 
     &    sp_light_lim_HIST, diat_N_lim_HIST, diat_Fe_lim_HIST, 
     &    diat_PO4_lim_HIST, diat_SiO3_lim_HIST,
     &    diat_light_lim_HIST, CaCO3_form_HIST, 
     &    diaz_Nfix_HIST, graze_diaz_HIST, diaz_loss_HIST,
     &     photoC_diaz_HIST, diaz_P_lim_HIST, 
     &    diaz_Fe_lim_HIST, diaz_light_lim_HIST, 
     &     Fe_scavenge_rate_HIST, DON_prod_HIST, 
     &    DON_remin_HIST, DOFe_prod_HIST, 
     &    DOFe_remin_HIST, DOP_prod_HIST, 
     &    DOP_remin_HIST, bSI_form_HIST, 
     &    photoFe_diaz_HIST, photoFe_diat_HIST, photoFe_sp_HIST,
     &    nitrif_HIST, j_o2_hist
# ifdef DMS_CYCLE
     &      , dmsp_prod_hist, dmsp_p_uptake_hist, dmsp_lysis_hist,
     &      dms_prod_hist, dms_photlys_hist
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
     &    , ammox_hist, denitr_hist, n2o_prod_hist, n2_prod_hist
# endif /* NCYCLE_ANOXIC */
# ifdef N15_CYCLE
     &      , no3_v_sp_n15_hist, 
     &      nh4_v_sp_n15_hist,
     &      no3_v_diat_n15_hist,
     &      nh4_v_diat_n15_hist,
     &      graze_sp_n15_hist, graze_diat_n15_hist, graze_diaz_n15_hist,
     &      sp_loss_n15_hist, diat_loss_n15_hist, 
     &      diaz_loss_n15_hist, zoo_loss_n15_hist,
     &      sp_agg_n15_hist, diat_agg_n15_hist,
     &      diaz_Nfix_n15_hist, pon_remin_n15_hist,
     &      nitrif_n15_hist, don_remin_n15_hist, don_prod_n15_hist
#  ifdef NCYCLE_ANOXIC
     &    , ammox_n15_hist, denitr_n15_hist, 
     &      n2o_prod_n15_hist, n2_prod_n15_hist
#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */
# ifdef SEDIMENT_BIOLOGY
       real bot_flux_poc_hist(GLOBAL_2D_ARRAY), 
     &      bot_flux_caco3_hist(GLOBAL_2D_ARRAY),
     &      bot_flux_si_hist(GLOBAL_2D_ARRAY)
#  ifdef N15_CYCLE
     &      , bot_flux_pon15_hist(GLOBAL_2D_ARRAY)
#  endif
       common /sed_flux/ bot_flux_poc_hist, bot_flux_caco3_hist,
     &      bot_flux_si_hist
#  ifdef N15_CYCLE
     &      , bot_flux_pon15_hist
#  endif
# endif /* SEDIMENT_BIOLOGY */
#endif /* BGC_FLUX_ANALYSIS */
!--------------------------------------------------------------------------
!   relative tracer indices
!--------------------------------------------------------------------------

!This index is not same as the index in param.h(roms).
       integer  po4_ind , no3_ind,sio3_ind, nh4_ind,fe_ind,
     & o2_ind, dic_ind,alk_ind,doc_ind,spC_ind,spChl_ind,
     & spCaCO3_ind,diatC_ind,diatChl_ind,zooC_ind,spFe_ind,
     &  diatSi_ind,diatFe_ind,diazC_ind,diazChl_ind, diazFe_ind,
     &  don_ind,dofe_ind,dop_ind

       parameter (po4_ind=1 , no3_ind=2,sio3_ind=3, nh4_ind=4,
     &  fe_ind=5,o2_ind=6, dic_ind=7,alk_ind=8,doc_ind=9,
     &  spC_ind=10,spChl_ind=11, spCaCO3_ind=12,diatC_ind=13,
     &  diatChl_ind=14,zooC_ind=15,spFe_ind=16,
     &  diatSi_ind=17,diatFe_ind=18,diazC_ind=19,
     &  diazChl_ind=20, diazFe_ind=21,
     &  don_ind=22,dofe_ind=23,dop_ind=24)

!  po4_ind          =  1,    dissolved inorganic phosphate
!  no3_ind          =  2,    dissolved inorganic nitrate
!  sio3_ind         =  3,    dissolved inorganic silicate
!  nh4_ind          =  4,    dissolved ammonia
!  fe_ind           =  5,    dissolved inorganic iron
!  o2_ind           =  6,    dissolved oxygen
!  dic_ind          =  7,    dissolved inorganic carbon
!  alk_ind          =  8,    alkalinity
!  doc_ind          =  9,    dissolved organic carbon
!  spC_ind          = 10,    small phytoplankton carbon
!  spChl_ind        = 11,    small phytoplankton chlorophyll
!  spCaCO3_ind      = 12,    small phytoplankton caco3
!  diatC_ind        = 13,    diatom carbon
!  diatChl_ind      = 14,    diatom chlorophyll
!  zooC_ind         = 15,    zooplankton carbon
!  spFe_ind         = 16     small phytoplankton iron
!  diatSi_ind       = 17,    diatom silicon
!  diatFe_ind       = 18,    diatom iron
!  diazC_ind        = 19,    diazotroph carbon
!  diazChl_ind      = 20,    diazotroph Chlorophyll
!  diazFe_ind       = 21,    diazotroph iron
!  don_ind          = 22,    dissolved organic nitrogen
!  dofe_ind         = 23,    dissolved organic iron
!  dop_ind          = 24     dissolved organic phosphorus

# ifdef DMS_CYCLE
       integer, parameter :: dmsp_ind = 25
       integer, parameter :: dms_ind = 26
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
       integer, parameter :: no2_ind = 25
#  ifdef DMS_CYCLE
     &      + 2
#  endif /* DMS_CYCLE */
       integer, parameter :: n2o_ind = 26
#  ifdef DMS_CYCLE
     &      + 2
#  endif /* DMS_CYCLE */
       integer, parameter :: n2_ind = 27
#  ifdef DMS_CYCLE
     &      + 2
#  endif /* DMS_CYCLE */
# endif /* NCYCLE_ANOXIC */
# ifdef N15_CYCLE
       integer, parameter :: no3_15n_ind = 25
#  ifdef DMS_CYCLE
     &      + 2
#  endif /* DMS_CYCLE */
#  ifdef NCYCLE_ANOXIC
     &      + 3
#  endif /* NCYCLE_ANOXIC */
       integer, parameter :: nh4_15n_ind = no3_15n_ind + 1
       integer, parameter :: sp_15n_ind = nh4_15n_ind + 1
       integer, parameter :: diat_15n_ind = sp_15n_ind + 1
       integer, parameter :: zoo_15n_ind = diat_15n_ind + 1
       integer, parameter :: diaz_15n_ind = zoo_15n_ind + 1
       integer, parameter :: don_15n_ind = diaz_15n_ind + 1
       integer, parameter :: no2_15n_ind = don_15n_ind + 1
       integer, parameter :: n2o_15n_ind = no2_15n_ind + 1
       integer, parameter :: n2_15n_ind = n2o_15n_ind + 1
# endif /* N15_CYCLE */

# ifdef SEDIMENT_BIOLOGY
       integer, parameter :: sed_poc_ind = 1
       integer, parameter :: sed_caco3_ind = 2
       integer, parameter :: sed_si_ind = 3
#  ifdef N15_CYCLE
       integer, parameter :: sed_pon15_ind = 4
#  endif
# endif /* SEDIMENT_BIOLOGY */

       logical lsource_sink,lflux_gas_o2, lflux_gas_co2,
     &  liron_flux,ldust_flux,lnitrate_flux,lammonium_flux,lalkatm_flux,lanthpco2_flux
        common /ecoflag/lsource_sink,lflux_gas_o2,lflux_gas_co2,
     &   liron_flux,ldust_flux,lnitrate_flux,lammonium_flux,lalkatm_flux,lanthpco2_flux

!------------------------------------------------------------------------
!   restoring climatologies for nutrients
!------------------------------------------------------------------------

       logical lrest_po4,lrest_no3,lrest_sio3 
!     lrest_po4,    restoring on po4 
!     lrest_no3,    restoring on no3 
!     lrest_sio3     restoring on sio3 

       real po4_clim(GLOBAL_2D_ARRAY,N),
     &   no3_clim(GLOBAL_2D_ARRAY,N),
     &   sio3_clim(GLOBAL_2D_ARRAY,N)
           real nutr_rest_time_inv(N)
! inverse restoring time scale for nutrients (1/secs)
        common /restore_flag/lrest_po4,lrest_no3,lrest_sio3
        common /restore_clim/po4_clim,
     &      no3_clim,sio3_clim,nutr_rest_time_inv


       
#if defined EXPLICIT_SINKING
      real t_sed(GLOBAL_2D_ARRAY,NT_sed)
      common /ocean_t_sed/t_sed        
#else
# ifdef SEDIMENT_BIOLOGY
      real t_sed(GLOBAL_2D_ARRAY,NT_sed)
      common /ocean_t_sed/t_sed
# endif
        real sinking_particle_POC(6,GLOBAL_2D_ARRAY),
     & sinking_particle_P_CaCO3(6,GLOBAL_2D_ARRAY),
     & sinking_particle_P_sio2(6,GLOBAL_2D_ARRAY),
     & sinking_particle_dust(6,GLOBAL_2D_ARRAY),
     & sinking_particle_P_iron(6,GLOBAL_2D_ARRAY)
# ifdef N15_CYCLE
     &       , sinking_particle_PON15(6,GLOBAL_2D_ARRAY)
# endif /* N15_CYCLE */

!        sinking_particle_POC,        base units = mmol C
!        sinking_particle_P_CaCO3,    base units = mmol CaCO3
!        sinking_particle_P_SiO2,     base units = mmol SiO2
!        sinking_particle_dust,       base units = kg
!        sinking_particle_P_iron      base units = mmol Fe

!1 for sflux_in,    incoming flux of soft subclass (base units/m^2/sec)
!2 for hflux_in,    incoming flux of hard subclass (base units/m^2/sec)
!3 for prod,        production term (base units/m^3/sec)
!4 for sflux_out,   outgoing flux of soft subclass (base units/m^2/sec)
!5 for hflux_out,   outgoing flux of hard subclass (base units/m^2/sec)
!6 for remin        remineralization term (base units/m^3/sec)

!           diss,        dissolution length for soft subclass
!           gamma,       fraction of production -> hard subclass
!           mass,         mass of 1e6 base units in kg
!           rhoo         QA mass ratio of POC to this particle class
! 1 for poc, 2 for p_caco3, 3 for p_sio2, 4 for dust, 5 for p_iron 

!WARN        real bot_flux_poc(GLOBAL_I_ARRAY)
!WARN        real bot_flux_caco3(GLOBAL_I_ARRAY)
!WARN        real bot_flux_si(GLOBAL_I_ARRAY)
       real diss(5),gamma(5),mass(5),rhoo(5) 
      common /sinking_part/sinking_particle_POC,
     &  sinking_particle_P_CaCO3,sinking_particle_P_SiO2, 
     &  sinking_particle_dust,sinking_particle_P_iron,
     &  diss,gamma,mass,rhoo
# ifdef N15_CYCLE
     &       , sinking_particle_PON15
# endif /* N15_CYCLE */
!WARN      &      , bot_flux_poc, bot_flux_caco3, bot_flux_si
#endif /* EXPLICIT_SINKING */

        logical landmask(GLOBAL_2D_ARRAY)
        common /calcation/landmask
#endif /* BIOLOGY_BEC */
