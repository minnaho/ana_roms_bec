!zhc: initialize ecosystem tracer , grid information, parameters ,variables

#include "cppdefs.h"
#ifdef BIOLOGY_BEC
      subroutine ecosys_init()
      implicit none
 
#include "param.h" 
 !need for im_eco,jm_eco,tracer_cnt_eco
#include "ocean3d.h"  
!need for tracer (t in roms),zt,dz
#include "grid.h" 
 !need for landmask form rmask(roms)
#include "param_bec.h"
#include "ecosys_bec.h"
# ifdef AVERAGES
#include "bec_avg.h"
# endif

!since tracer structure is different to t(roms),first, assign to tracer
      integer::i,j,k,m,nn
      
      sinking_particle_POC=0.0
      sinking_particle_P_CaCO3=0.0
      sinking_particle_P_SiO2=0.0
      sinking_particle_dust=0.0
      sinking_particle_P_iron =0.0
         
       parm_Fe_bioavail    = 0.02
       parm_prod_dissolve  = 0.67
       parm_o2_min         = 1.0
       parm_Rain_CaCO3     = 0.07
       parm_Rain_SiO2      = 0.03
# ifdef NCYCLE_ANOXIC
       parm_kappa_nitrif_no2 = 0.33333 * dps
       parm_n2o_prod       = 0.95
       parm_n2_prod        = 0.33333 * dps
!       parm_n2_prod        = 0.033333 * dps
# endif
       parm_kappa_nitrif   = 0.06 * dps       ! (= 1/( days))
       parm_nitrif_par_lim = 4.0
! HF: changed from 6.6e-5 on Oct 3, 2005, following Keith Moore's changes
       parm_POC_flux_ref   = 2.e-3           ! (mmol C/m^2/sec) - WAS 6.6e-3 nmol C/cm^2/sec
       parm_rest_prod_tau  = 30.0 * dps       ! (= 30 days)
       parm_rest_prod_z_c  = 75.              ! (m)
       parm_z_umax_0       = 2.5 * dps
#ifdef MOORE_GBC_2004
       parm_diat_umax_0    = 2.0 * dps     
#else
       parm_diat_umax_0    = 1.95 * dps ! standard value    
!       parm_diat_umax_0    = 1.035 * dps     
!       parm_diat_umax_0    = 4.14 * dps     
#endif
       parm_z_mort_0       = 0.08 * dps
! HF: changed from 0.46 on Oct 3, 2005, following Keith Moore's changes
#ifdef MOORE_GBC_2004
       parm_z_mort2_0      = 0.46 * dps
#else
       parm_z_mort2_0      = 0.42 * dps
#endif
       parm_sd_remin_0     = 0.01 * dps       ! (= 1/(100 days))
       parm_sp_kNO3        = 0.5
       parm_diat_kNO3      = 2.5
#ifdef MOORE_GBC_2004
       parm_sp_kNH4        = 0.005
#else
       parm_sp_kNH4        = 0.01
#endif
! HF: changed from 0.06 on Oct 3, 2005, following Keith Moore's changes
       parm_diat_kNH4      = 0.1
#ifdef PHAEOCYSTIS
! HF: parameters from Wang & Moore, 2011 (Table 1)
       parm_phaeo_kNO3     = 1.5
       parm_phaeo_kNH4     = 0.06
!12km       parm_phaeo_kNO3     = 1.16667 
!12km      parm_phaeo_kNH4     = 0.04
#endif
       parm_sp_kFe         = 0.035e-3          ! mmol Fe/m3
! HF: changed from 0.16e-3 on Oct 3, 2005, following Keith Moore's changes
#ifdef MOORE_GBC_2004
       parm_diat_kFe       = 0.16e-3          ! mmol Fe/m3
#else
       parm_diat_kFe       = 0.08e-3          ! mmol Fe/m3
#endif
       parm_diat_kSiO3     = 1.0
#ifdef MOORE_GBC_2004
       parm_sp_kPO4        = 0.0003125
#else
       parm_sp_kPO4        = 0.01
#endif
       parm_diat_kPO4      = 0.1
       parm_z_grz          = 1.05
! HF: changed from "0.25 * dps" on Oct 3, 2005, following Keith Moore's changes
#ifdef MOORE_GBC_2004
       parm_alphaChl       = 0.25 * dps
#else
       parm_alphaChl       = 0.3 * dps
#endif
       parm_labile_ratio   = 0.70
#ifdef MOORE_GBC_2004
       parm_alphaDiaz      = 0.028 * dps
       parm_diaz_umax_0    = 1.2 * dps
#else
#ifdef PHAEOCYSTIS
       parm_alphaDiaz      = 0.63 * dps ! Wang & Moore, 2011 (Table 1)
!12km       parm_alphaDiaz      = 0.3 * dps ! same as diat. 
# ifdef UMAX35
       parm_diaz_umax_0    = 3.5 * dps
# elif defined UMAX51
       parm_diaz_umax_0    = 5.1 * dps  ! Wang & Moore, 2011 (Table 1)
# else
       parm_diaz_umax_0    = 1.95 * dps
# endif
#else /* PHAEOCYSTIS */
! HF: changed from "0.025*dps" on Oct 3, 2005, following Keith Moore's changes
       parm_alphaDiaz      = 0.036 * dps
       parm_diaz_umax_0    = 1.2 * dps
#endif /* PHAEOCYSTIS */
#endif

# ifdef SEDIMENT_BIOLOGY
       t_remin_sed_poc   = 0.003 * dps  ! /d
       t_remin_sed_caco3 = 0.003 * dps  ! /d
       t_remin_sed_si    = 0.003 * dps  ! /d
#  ifdef N15_CYCLE
       t_remin_sed_pon15 = t_remin_sed_poc ! PON is coupled to POC
#  endif
# endif /* SEDIMENT_BIOLOGY */


        gQsi_0        = 0.137
        gQsi_coef     = 2.5
        gQsi_max      = 0.685
         
        lrest_po4      = .FALSE.
        lrest_no3      = .FALSE.
        lrest_sio3     = .FALSE.

        lsource_sink    = .TRUE.
        lflux_gas_o2    = .TRUE.  !read gas flux
        lflux_gas_co2   = .TRUE.
        liron_flux      =.TRUE.
        ldust_flux      =.TRUE.

        lnitrate_flux      =.TRUE.
        lammonium_flux      =.TRUE.
        lalkatm_flux      =.TRUE.
        lanthpco2_flux      =.TRUE.

         where(rmask ==1)
              landmask=.TRUE.
            elsewhere 
               landmask=.FALSE.
            endwhere

#ifdef DMS_CYCLE
       parm_I_min_DMS = 25.0  ! W m-2
       parm_I_max_DMS = 150.0 ! W m-2
       parm_I_mean_DMS = 0.5 * (parm_I_min_DMS + parm_I_max_DMS)
       parm_inv_I_del_DMS = 1.0 / (parm_I_max_DMS - parm_I_min_DMS)

       parm_r_min_DMS = 0.04  ! mmol S / mmol N
       parm_r_max_DMS = 0.22  ! mmol S / mmol N
       parm_r_del_DMS = parm_r_max_DMS - parm_r_min_DMS

       parm_gamma1 = 0.05     ! [adim]

       parm_f_lysis_DMSP = 0.01 * dps    ! s-1

       parm_gamma_s_max_DMS = 0.25 * dps ! s-1

       parm_k_phot_DMS = 0.15 * dps ! s-1

       parm_alpha1 = 0.7      ! [adim]

       parm_alpha_p = 0.1     ! [adim]
#endif /* DMS_CYCLE */

#ifdef N15_CYCLE
!       parm_n15_std_fraction = 0.00366
       parm_n15_std_fraction = 1.0
!
       parm_alpha_n15_no3_v_sp = 0.995
       parm_alpha_n15_nh4_v_sp = 0.995
       parm_alpha_n15_no3_v_diat = 0.995
       parm_alpha_n15_nh4_v_diat = 0.995
# ifdef PHAEOCYSTIS
       parm_alpha_n15_no3_v_phaeo = 0.995
       parm_alpha_n15_nh4_v_phaeo = 0.995
# endif
!       parm_alpha_n15_no3_v_sp = 1.0
!       parm_alpha_n15_nh4_v_sp = 1.0
!       parm_alpha_n15_no3_v_diat = 1.0
!       parm_alpha_n15_nh4_v_diat = 1.0
       parm_alpha_n15_nitrif = 1.0
       parm_alpha_n15_don_remin = 1.0
       parm_alpha_n15_don_prod = 1.0
       parm_alpha_n15_graze_sp_zoo = 1.0
       parm_alpha_n15_graze_sp_pon = 1.0
       parm_alpha_n15_graze_sp_don = 1.0
       parm_alpha_n15_graze_sp_nh4 = 1.0
       parm_alpha_n15_sp_loss_pon = 1.0
       parm_alpha_n15_sp_loss_don = 1.0
       parm_alpha_n15_sp_loss_nh4 = 1.0
       parm_alpha_n15_sp_agg = 1.0
       parm_alpha_n15_graze_diat_zoo = 1.0
       parm_alpha_n15_graze_diat_pon = 1.0
       parm_alpha_n15_graze_diat_don = 1.0
       parm_alpha_n15_graze_diat_nh4 = 1.0
       parm_alpha_n15_diat_loss_pon = 1.0
       parm_alpha_n15_diat_loss_don = 1.0
       parm_alpha_n15_diat_loss_nh4 = 1.0
       parm_alpha_n15_diat_agg = 1.0
       parm_alpha_n15_graze_diaz_zoo = 1.0
       parm_alpha_n15_graze_diaz_pon = 1.0
       parm_alpha_n15_graze_diaz_don = 1.0
       parm_alpha_n15_graze_diaz_nh4 = 1.0
       parm_alpha_n15_diaz_agg       = 1.0
       parm_alpha_n15_diaz_loss_don = 1.0
       parm_alpha_n15_diaz_loss_nh4 = 1.0
       parm_alpha_n15_diaz_Nfix = 1.0
       parm_alpha_n15_diaz_Nexcrete = parm_alpha_n15_diaz_Nfix
       parm_alpha_n15_zoo_loss_don = 1.0
       parm_alpha_n15_zoo_loss_nh4 = 0.995
!       parm_alpha_n15_zoo_loss_nh4 = 1.0
       parm_alpha_n15_pon_remin = 1.0
# ifdef NCYCLE_ANOXIC
       parm_alpha_n15_ammox = 1.0
       parm_alpha_n15_denitr = 0.975
!       parm_alpha_n15_denitr = 1.0
       parm_alpha_n15_denitr_sed = 1.0
       parm_alpha_n15_n2o_prod = 1.0
       parm_alpha_n15_n2_prod = 1.0
# endif /* NCYCLE_ANOXIC */
#endif /* N15_CYCLE */

    !------------------------------------------------------------------------
    ! initialize tracers and apply land mask to tracers (suppose done by roms
    !------------------------------------------------------------------------
        

    !------------------------------------------------------------------------
    !   read gas flux forcing (if required)(fice_file,xkw_file and ap_file
    !   dont need right now 
    !------------------------------------------------------------------------
         ifrac=0.0
         press=1.


    !-------------------------------------------------------------------------
    !   initialize restoring timescale (if required)
    !------------------------------------------------------------------------
      IF (lrest_po4 .OR. lrest_no3 .OR. lrest_sio3) THEN
!zhc: all flags are false, no action right now
      END IF
    !-----------------------------------------------------------------------
    !   load restoring fields (if required)
    !-----------------------------------------------------------------------

    !------------------------------------------------------------------------
    !   load dust flux fields (if required)dust_flux,iron_flux
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !  If co2 gas flux is to be computed, initialize PH_HIST. The way in which
    !   this is done is based upon aux_restart_option.
    !------------------------------------------------------------------------
!!HF      IF (lflux_gas_co2) THEN
!!HF         ph_hist=0.0            !zhc:initial ph with zero
!!HF      endif

       pH_hist              = c0    
       pCO2sw               = c0    
       pCO2air              = c0
       PAR                  = c0    
       PARinc               = c0    
# ifdef BGC_FLUX_ANALYSIS
       WS_HIST              = c0   
       XKW_HIST             = c0    
       AP_HIST              = c0    
       SCHMIDT_O2_HIST      = c0    
       O2SAT_HIST           = c0    
       FG_O2_HIST           = c0    
       SCHMIDT_CO2_HIST     = c0    
# ifdef DMS_CYCLE
       fg_dms_hist          = c0
# endif /* DMS_CYCLE */ 
# ifdef NCYCLE_ANOXIC
       fg_n2o_hist          = c0
       fg_n2_hist           = c0
#  ifdef N15_CYCLE
       fg_n2o_n15_hist      = c0
       fg_n2_n15_hist       = c0
#  endif
# endif /* NCYCLE_ANOXIC */ 
       CO2STAR_HIST         = c0    
       DCO2STAR_HIST        = c0    
       FG_CO2_HIST          = c0    
       IRON_FLUX_HIST       = c0    
       nitrate_FLUX_HIST    = c0
       ammonium_FLUX_HIST   = c0
       alkatm_FLUX_HIST     = c0
       anthpco2_FLUX_HIST   = c0
       PO4_RESTORE_HIST     = c0    
       NO3_RESTORE_HIST     = c0    
       SiO3_RESTORE_HIST    = c0    
       PO4STAR_HIST         = c0    
       POC_FLUX_IN_HIST     = c0    
       POC_PROD_HIST        = c0    
       POC_REMIN_HIST       = c0    
       CaCO3_FLUX_IN_HIST   = c0    
       CaCO3_PROD_HIST      = c0    
       CaCO3_REMIN_HIST     = c0    
       SiO2_FLUX_IN_HIST    = c0    
       SiO2_PROD_HIST       = c0    
       SiO2_REMIN_HIST      = c0    
       dust_FLUX_IN_HIST    = c0    
       dust_REMIN_HIST      = c0    
       P_iron_FLUX_IN_HIST  = c0   
       P_iron_PROD_HIST     = c0    
       P_iron_REMIN_HIST    = c0    
       graze_sp_HIST        = c0    
       graze_diat_HIST      = c0    
       graze_tot_HIST       = c0    
       sp_loss_HIST         = c0    
       diat_loss_HIST       = c0    
       zoo_loss_HIST        = c0    
       sp_agg_HIST          = c0    
       diat_agg_HIST        = c0    
       photoC_sp_HIST       = c0    
       photoC_diat_HIST     = c0    
       tot_prod_HIST        = c0    
       no3_v_sp_HIST        = c0    
       nh4_v_sp_HIST        = c0    
       no3_v_diat_HIST      = c0    
       nh4_v_diat_HIST      = c0
# ifdef PHAEOCYSTIS
       no3_v_phaeo_HIST     = c0    
       nh4_v_phaeo_HIST     = c0    
# endif
       DOC_prod_HIST        = c0    
       DOC_remin_HIST       = c0    
       Fe_scavenge_HIST     = c0    
       sp_N_lim_HIST        = c0    
       sp_Fe_lim_HIST       = c0    
       sp_PO4_lim_HIST      = c0    
       sp_light_lim_HIST    = c0    
       diat_N_lim_HIST      = c0    
       diat_Fe_lim_HIST     = c0    
       diat_PO4_lim_HIST    = c0    
       diat_SiO3_lim_HIST   = c0    
       diat_light_lim_HIST  = c0    
       CaCO3_form_HIST      = c0    
       diaz_Nfix_HIST       = c0    
       graze_diaz_HIST      = c0    
       diaz_loss_HIST       = c0    
       photoC_diaz_HIST     = c0    
       diaz_P_lim_HIST      = c0    
       diaz_Fe_lim_HIST     = c0    
       diaz_light_lim_HIST  = c0    
       Fe_scavenge_rate_HIST = c0   
       DON_prod_HIST     = c0    
       DON_remin_HIST    = c0    
       DOFe_prod_HIST    = c0    
       DOFe_remin_HIST   = c0    
       DOP_prod_HIST     = c0    
       DOP_remin_HIST    = c0    
       bSi_form_HIST     = c0    
       photoFe_diat_HIST    = c0    
       photoFe_diaz_HIST    = c0    
       photoFe_sp_HIST      = c0
       nitrif_hist       = c0
# ifdef DMS_CYCLE
       dmsp_prod_hist    = c0
       dmsp_lysis_hist   = c0
       dms_prod_hist     = c0
       dms_photlys_hist  = c0
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
       ammox_hist = c0
       denitr_hist = c0
       n2o_prod_hist = c0
       n2_prod_hist = c0
# endif /* NCYCLE_ANOXIC */
# endif /* BGC_FLUX_ANALYSIS */

# ifdef AVERAGES
       pH_avg              = c0    
       pCO2_avg          = c0    
       pCO2air_avg         = c0
       PAR_avg             = c0    
       PARinc_avg          = c0    
# ifdef BGC_FLUX_ANALYSIS
       WS_AVG              = c0   
       XKW_AVG             = c0    
       AP_AVG              = c0    
       SCHMIDT_O2_AVG      = c0    
       O2SAT_AVG           = c0    
       FG_O2_AVG           = c0    
       SCHMIDT_CO2_AVG     = c0    
# ifdef DMS_CYCLE
       fg_dms_avg          = c0
# endif /* DMS_CYCLE */ 
# ifdef NCYCLE_ANOXIC
       fg_n2o_avg          = c0
       fg_n2_avg           = c0
#  ifdef N15_CYCLE
       fg_n2o_n15_avg      = c0
       fg_n2_n15_avg       = c0
#  endif
# endif /* NCYCLE_ANOXIC */ 
       CO2STAR_AVG         = c0    
       DCO2STAR_AVG        = c0    
       FG_CO2_AVG          = c0    
       IRON_FLUX_AVG       = c0    
       nitrate_FLUX_AVG    = c0
       ammonium_FLUX_AVG   = c0
       alkatm_FLUX_AVG     = c0
       anthpco2_FLUX_AVG   = c0
       PO4_RESTORE_AVG     = c0    
       NO3_RESTORE_AVG     = c0    
       SiO3_RESTORE_AVG    = c0    
       PO4STAR_AVG         = c0    
       POC_FLUX_IN_AVG     = c0    
       POC_PROD_AVG        = c0    
       POC_REMIN_AVG       = c0    
       CaCO3_FLUX_IN_AVG   = c0    
       CaCO3_PROD_AVG      = c0    
       CaCO3_REMIN_AVG     = c0    
       SiO2_FLUX_IN_AVG    = c0    
       SiO2_PROD_AVG       = c0    
       SiO2_REMIN_AVG      = c0    
       dust_FLUX_IN_AVG    = c0    
       dust_REMIN_AVG      = c0    
       P_iron_FLUX_IN_AVG  = c0    
       P_iron_PROD_AVG     = c0    
       P_iron_REMIN_AVG    = c0    
       graze_sp_AVG        = c0    
       graze_diat_AVG      = c0    
       graze_tot_AVG       = c0    
       sp_loss_AVG         = c0    
       diat_loss_AVG       = c0    
       zoo_loss_AVG        = c0    
       sp_agg_AVG          = c0    
       diat_agg_AVG        = c0    
       photoC_sp_AVG       = c0    
       photoC_diat_AVG     = c0    
       tot_prod_AVG        = c0    
       no3_v_sp_AVG        = c0    
       nh4_v_sp_AVG        = c0    
       no3_v_diat_AVG      = c0    
       nh4_v_diat_AVG      = c0
# ifdef PHAEOCYSTIS
       no3_v_phaeo_AVG     = c0    
       nh4_v_phaeo_AVG     = c0    
# endif
       DOC_prod_AVG        = c0    
       DOC_remin_AVG       = c0    
       Fe_scavenge_AVG     = c0    
       sp_N_lim_AVG        = c0    
       sp_Fe_lim_AVG       = c0    
       sp_PO4_lim_AVG      = c0    
       sp_light_lim_AVG    = c0    
       diat_N_lim_AVG      = c0    
       diat_Fe_lim_AVG     = c0    
       diat_PO4_lim_AVG    = c0    
       diat_SiO3_lim_AVG   = c0    
       diat_light_lim_AVG  = c0    
       CaCO3_form_AVG      = c0    
       diaz_Nfix_AVG       = c0    
       graze_diaz_AVG      = c0    
       diaz_loss_AVG       = c0    
       photoC_diaz_AVG     = c0    
       diaz_P_lim_AVG      = c0    
       diaz_Fe_lim_AVG     = c0    
       diaz_light_lim_AVG  = c0    
       Fe_scavenge_rate_AVG = c0   
       DON_prod_AVG      = c0    
       DON_remin_AVG     = c0    
       DOFe_prod_AVG     = c0    
       DOFe_remin_AVG    = c0    
       DOP_prod_AVG      = c0    
       DOP_remin_AVG     = c0    
       bSi_form_AVG      = c0    
       photoFe_diat_AVG    = c0    
       photoFe_diaz_AVG    = c0    
       photoFe_sp_AVG      = c0
       nitrif_avg       = c0
# ifdef DMS_CYCLE
       dmsp_prod_avg    = c0
       dmsp_lysis_avg   = c0
       dms_prod_avg     = c0
       dms_photlys_avg  = c0
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
       ammox_avg = c0
       denitr_avg = c0
       n2o_prod_avg = c0
       n2_prod_avg = c0
# endif /* NCYCLE_ANOXIC */
# endif /* BGC_FLUX_ANALYSIS */
# endif /* AVERAGES */

# ifdef BGC_FLUX_ANALYSIS
       call init_scalars_becflux()
# endif /* BGC_FLUX_ANALYSIS */
#else
       subroutine ecosys_init_empty()
#endif  /*BIOLOGY_BEC*/

       return 
       end
