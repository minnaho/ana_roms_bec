#include "cppdefs.h"
#if defined BIOLOGY_BEC && defined BGC_FLUX_ANALYSIS && defined AVERAGES
      subroutine set_bec_flux_avg(tile)
      implicit none
      integer tile
# include "param.h"

# include "compute_tile_bounds.h"
      call set_bec_flux_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

         subroutine set_bec_flux_avg_tile (Istr,Iend,Jstr,Jend)
         implicit none
# include "param.h"
# include "ecosys_bec.h"
# include "bec_avg.h"
# include "ncroms_bgcflux.h"
# include "scalars.h"
# include "ocean2d.h"
         integer istr,jstr,iend,jend
         real cff, cff1

# include "compute_auxiliary_bounds.h"
! calculation of averaged fluxes will only be performed if n_bgc_flux_avg
! is a positive number
         if (n_bgc_flux_avg < 1) return

         if (iic.gt.nts_bgc_flux_avg) then
            if (n_bgc_flux_avg .eq. 1 .or.
     &           mod(iic-nts_bgc_flux_avg,n_bgc_flux_avg).eq.1) then
               cff =1.0
               cff1=0.0
               if (ZEROTH_TILE) then
                  time_bgc_flux_avg=time
                  mpi_master_only write(*,'(6x,A,I11,2(X,A,I7))')
     &                 'SET_BGC_FLUX_AVG -- Started averaging at iic=',
     &                 iic,' : nts_bgc_flux_avg=',nts_bgc_flux_avg,
     &                 'n_bgc_flux_avg=',n_bgc_flux_avg
               endif
            elseif (mod(iic-nts_bgc_flux_avg,n_bgc_flux_avg).gt.1) then
               cff =1.0
               cff1=1.0
               if (ZEROTH_TILE) time_bgc_flux_avg=time_bgc_flux_avg+time
            elseif (mod(iic-nts_bgc_flux_avg,n_bgc_flux_avg).eq.0) then
               cff=1./float(n_bgc_flux_avg)
               cff1=1.0
               if (ZEROTH_TILE) then
                  time_bgc_flux_avg=cff*(time_bgc_flux_avg+time)
                  mpi_master_only write(*,'(6x,A,I11,2(X,A,I7))')
     &                 'SET_BGC_FLUX_AVG -- Finished averaging at iic=',
     &                 iic,' : nts_bgc_flux_avg=',nts_bgc_flux_avg,
     &                 'n_bgc_flux_avg=',n_bgc_flux_avg
               endif
            endif

! 2D variables
            zeta_bgc_flux_avg(istrR:iendR,jstrR:jendR) = 
     &        cff * ( cff1*zeta_bgc_flux_avg(istrR:iendR,jstrR:jendR) +
     &        zeta(istrR:iendR,jstrR:jendR,knew) )
            WS_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*WS_AVG(istrR:iendR,jstrR:jendR) 
     &           + WS_hist(istrR:iendR,jstrR:jendR) )
            XKW_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*XKW_AVG(istrR:iendR,jstrR:jendR) 
     &           + XKW_hist(istrR:iendR,jstrR:jendR) )
            ap_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*ap_AVG(istrR:iendR,jstrR:jendR) 
     &           + ap_hist(istrR:iendR,jstrR:jendR) )
            SCHMIDT_O2_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*SCHMIDT_O2_AVG(istrR:iendR,jstrR:jendR) 
     &           + SCHMIDT_O2_hist(istrR:iendR,jstrR:jendR) )
            O2SAT_AVG(istrR:iendR,jstrR:jendR)= 
     &           cff * ( cff1*O2SAT_AVG(istrR:iendR,jstrR:jendR) 
     &           +O2SAT_hist(istrR:iendR,jstrR:jendR) )
            FG_O2_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_O2_AVG(istrR:iendR,jstrR:jendR) 
     &           +FG_O2_hist(istrR:iendR,jstrR:jendR) )
# ifdef DMS_CYCLE
            FG_DMS_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_DMS_avg(istrR:iendR,jstrR:jendR)
     &           + FG_DMS_hist(istrR:iendR,jstrR:jendR) )
# endif /* DMS_CYCLE */
# ifdef NCYCLE_ANOXIC
            FG_N2O_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_N2O_avg(istrR:iendR,jstrR:jendR)
     &           + FG_N2O_hist(istrR:iendR,jstrR:jendR) )
            FG_N2_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_N2_avg(istrr:iendR,jstrR:jendR)
     &           + FG_N2_hist(istrR:iendR,jstrR:jendR) )
#  ifdef N15_CYCLE
            FG_N2O_N15_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_N2O_N15_avg(istrR:iendR,jstrR:jendR)
     &           + FG_N2O_N15_hist(istrR:iendR,jstrR:jendR) )
            FG_N2_N15_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_N2_N15_avg(istrR:iendR,jstrR:jendR)
     &           + FG_N2_N15_hist(istrR:iendR,jstrR:jendR) )
#  endif /* N15_CYCLE */
# endif /* NCYCLE_ANOXIC */
            SCHMIDT_CO2_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*SCHMIDT_CO2_AVG(istrR:iendR,jstrR:jendR) 
     &           + SCHMIDT_CO2_hist(istrR:iendR,jstrR:jendR) )
            CO2STAR_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*CO2STAR_AVG(istrR:iendR,jstrR:jendR) 
     &           + CO2STAR_hist(istrR:iendR,jstrR:jendR) )
            DCO2STAR_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*DCO2STAR_AVG(istrR:iendR,jstrR:jendR) 
     &           + DCO2STAR_hist(istrR:iendR,jstrR:jendR) )
            FG_CO2_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*FG_CO2_AVG(istrR:iendR,jstrR:jendR) 
     &           + FG_CO2_hist(istrR:iendR,jstrR:jendR) )
            IRON_FLUX_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*IRON_FLUX_AVG(istrR:iendR,jstrR:jendR) 
     &           + IRON_FLUX_hist(istrR:iendR,jstrR:jendR) )

            nitrate_FLUX_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*nitrate_FLUX_AVG(istrR:iendR,jstrR:jendR)
     &           + nitrate_FLUX_hist(istrR:iendR,jstrR:jendR) )

            ammonium_FLUX_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*ammonium_FLUX_AVG(istrR:iendR,jstrR:jendR)
     &           + ammonium_FLUX_hist(istrR:iendR,jstrR:jendR) )

            alkatm_FLUX_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*alkatm_FLUX_AVG(istrR:iendR,jstrR:jendR)
     &           + alkatm_FLUX_hist(istrR:iendR,jstrR:jendR) )

            anthpco2_FLUX_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*anthpco2_FLUX_AVG(istrR:iendR,jstrR:jendR)
     &           + anthpco2_FLUX_hist(istrR:iendR,jstrR:jendR) )

            PARinc_flux_AVG(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*PARinc_flux_AVG(istrR:iendR,jstrR:jendR) 
     &           + PARinc(istrR:iendR,jstrR:jendR) )
# ifdef SEDIMENT_BIOLOGY
            bot_flux_poc_avg(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*bot_flux_poc_avg(istrR:iendR,jstrR:jendR) 
     &           + bot_flux_poc_hist(istrR:iendR,jstrR:jendR) )

            bot_flux_caco3_avg(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*bot_flux_caco3_avg(istrR:iendR,jstrR:jendR) 
     &           + bot_flux_caco3_hist(istrR:iendR,jstrR:jendR) )

            bot_flux_si_avg(istrR:iendR,jstrR:jendR) = 
     &           cff * ( cff1*bot_flux_si_avg(istrR:iendR,jstrR:jendR) 
     &           + bot_flux_si_hist(istrR:iendR,jstrR:jendR) )
#  ifdef N15_CYCLE
            bot_flux_pon15_avg(istrR:iendR,jstrR:jendR) = 
     &           cff * (cff1*bot_flux_pon15_avg(istrR:iendR,jstrR:jendR) 
     &           + bot_flux_pon15_hist(istrR:iendR,jstrR:jendR) )
#  endif
# endif /* SEDIMENT_BIOLOGY */
! 3D variables
            PO4_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * ( cff1*PO4_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + PO4_RESTORE_hist(istrR:iendR,jstrR:jendR,:) )
            NO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * ( cff1*NO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:)  
     &           + NO3_RESTORE_hist(istrR:iendR,jstrR:jendR,:) )
            SiO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * (cff1*SiO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + SiO3_RESTORE_hist(istrR:iendR,jstrR:jendR,:) )
            PAR_flux_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * (cff1*PAR_flux_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + PAR(istrR:iendR,jstrR:jendR,:) )
            PO4STAR_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * (cff1*PO4STAR_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + PO4STAR_hist(istrR:iendR,jstrR:jendR,:) )
            POC_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * (cff1*POC_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)  
     &           + POC_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            POC_PROD_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * (cff1*POC_PROD_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + POC_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            POC_REMIN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff * (cff1*POC_REMIN_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + POC_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*CaCO3_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + CaCO3_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_PROD_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*CaCO3_PROD_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + CaCO3_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_REMIN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*CaCO3_REMIN_AVG(istrR:iendR,jstrR:jendR,:) +
     &           CaCO3_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            SiO2_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*SiO2_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)+
     &           SiO2_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            SiO2_PROD_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*SiO2_PROD_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + SiO2_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            SiO2_REMIN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*SiO2_REMIN_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + SiO2_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            dust_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*dust_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + dust_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            dust_REMIN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff *(cff1*dust_REMIN_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +dust_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            P_iron_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*P_iron_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) 
     &          +P_iron_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            P_iron_PROD_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*P_iron_PROD_AVG(istrR:iendR,jstrR:jendR,:) 
     &             + P_iron_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            P_iron_REMIN_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*P_iron_REMIN_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + P_iron_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            graze_sp_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*graze_sp_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + graze_sp_hist(istrR:iendR,jstrR:jendR,:) )
            graze_diat_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*graze_diat_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + graze_diat_hist(istrR:iendR,jstrR:jendR,:) )
            graze_tot_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*graze_tot_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + graze_tot_hist(istrR:iendR,jstrR:jendR,:) )
            sp_loss_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*sp_loss_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +sp_loss_hist(istrR:iendR,jstrR:jendR,:) )
            diat_loss_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diat_loss_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diat_loss_hist(istrR:iendR,jstrR:jendR,:) )
            zoo_loss_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*zoo_loss_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + zoo_loss_hist(istrR:iendR,jstrR:jendR,:) )
            sp_agg_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*sp_agg_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +sp_agg_hist(istrR:iendR,jstrR:jendR,:) )
            diat_agg_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diat_agg_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diat_agg_hist(istrR:iendR,jstrR:jendR,:) )
            photoC_sp_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*photoC_sp_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + photoC_sp_hist(istrR:iendR,jstrR:jendR,:) )
            f_ratio_sp_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*f_ratio_sp_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + f_ratio_sp_hist(istrR:iendR,jstrR:jendR,:) )
            photoC_diat_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*photoC_diat_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + photoC_diat_hist(istrR:iendR,jstrR:jendR,:) )
            f_ratio_diat_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*f_ratio_diat_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + f_ratio_diat_hist(istrR:iendR,jstrR:jendR,:) )
            tot_prod_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*tot_prod_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +tot_prod_hist(istrR:iendR,jstrR:jendR,:) )
            no3_v_sp_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*no3_v_sp_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +no3_v_sp_hist(istrR:iendR,jstrR:jendR,:) )
            nh4_v_sp_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*nh4_v_sp_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +nh4_v_sp_hist(istrR:iendR,jstrR:jendR,:) )
            no3_v_diat_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*no3_v_diat_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +no3_v_diat_hist(istrR:iendR,jstrR:jendR,:) )
            nh4_v_diat_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*nh4_v_diat_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +nh4_v_diat_hist(istrR:iendR,jstrR:jendR,:) )
# ifdef PHAEOCYSTIS
            no3_v_phaeo_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*no3_v_phaeo_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +no3_v_phaeo_hist(istrR:iendR,jstrR:jendR,:) )
            nh4_v_phaeo_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*nh4_v_phaeo_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +nh4_v_phaeo_hist(istrR:iendR,jstrR:jendR,:) )
# endif /* PHAEOCYSTIS */
            DOC_prod_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*DOC_prod_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + DOC_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DOC_remin_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*DOC_remin_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +DOC_remin_hist(istrR:iendR,jstrR:jendR,:) )
            Fe_scavenge_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*Fe_scavenge_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + Fe_scavenge_hist(istrR:iendR,jstrR:jendR,:) )
            sp_N_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*sp_N_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + sp_N_lim_hist(istrR:iendR,jstrR:jendR,:) )
            sp_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*sp_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + sp_Fe_lim_hist(istrR:iendR,jstrR:jendR,:) )
            sp_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*sp_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + sp_PO4_lim_hist(istrR:iendR,jstrR:jendR,:) )
            sp_light_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_light_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + sp_light_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_N_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diat_N_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diat_N_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diat_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diat_Fe_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diat_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diat_PO4_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_SiO3_lim_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*diat_SiO3_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + diat_SiO3_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_light_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diat_light_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + diat_light_lim_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_form_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*CaCO3_form_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +CaCO3_form_hist(istrR:iendR,jstrR:jendR,:) )
#ifdef PHAEOCYSTIS
! the "Nfix" variable is used for N limitation instead
#endif
            diaz_Nfix_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diaz_Nfix_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diaz_Nfix_hist(istrR:iendR,jstrR:jendR,:) )
            graze_diaz_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*graze_diaz_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + graze_diaz_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_loss_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diaz_loss_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + diaz_loss_hist(istrR:iendR,jstrR:jendR,:) )
            photoC_diaz_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoC_diaz_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + photoC_diaz_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_P_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diaz_P_lim _AVG(istrR:iendR,jstrR:jendR,:) 
     &           + diaz_P_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diaz_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           +diaz_Fe_lim _hist(istrR:iendR,jstrR:jendR,:) )
            diaz_light_lim_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*diaz_light_lim_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + diaz_light_lim_hist(istrR:iendR,jstrR:jendR,:) )
            Fe_scavenge_rate_AVG(istrR:iendR,jstrR:jendR,:) = 
     &         cff*(cff1*Fe_scavenge_rate_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + Fe_scavenge_rate_hist(istrR:iendR,jstrR:jendR,:) )
            DON_prod_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*DON_prod_AVG(istrR:iendR,jstrR:jendR,:) + 
     &           DON_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DON_remin_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*DON_remin_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + DON_remin_hist(istrR:iendR,jstrR:jendR,:) )
            DOFe_prod_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*DOFe_prod_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DOFe_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DOFe_remin_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*DOFe_remin_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + DOFe_remin_hist(istrR:iendR,jstrR:jendR,:) )
            DOP_prod_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOP_prod_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + DOP_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DOP_remin_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*DOP_remin_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + DOP_remin_hist(istrR:iendR,jstrR:jendR,:) )
            bSI_form_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*bSI_form_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + bSI_form_hist(istrR:iendR,jstrR:jendR,:) )
            photoFe_diaz_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*photoFe_diaz_AVG(istrR:iendR,jstrR:jendR,:)
     &           + photoFe_diaz_hist(istrR:iendR,jstrR:jendR,:) )
            photoFe_diat_AVG(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*photoFe_diat_AVG(istrR:iendR,jstrR:jendR,:) +  
     &           photoFe_diat_hist(istrR:iendR,jstrR:jendR,:) )
            photoFe_sp_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*photoFe_sp_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + photoFe_sp_hist(istrR:iendR,jstrR:jendR,:) )
            nitrif_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*nitrif_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + nitrif_hist(istrR:iendR,jstrR:jendR,:) )
            j_o2_AVG(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*j_o2_AVG(istrR:iendR,jstrR:jendR,:) 
     &           + j_o2_hist(istrR:iendR,jstrR:jendR,:) )
# ifdef DMS_CYCLE
            dmsp_prod_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*dmsp_prod_avg(istrR:iendR,jstrR:jendR,:)
     &           + dmsp_prod_hist(istrR:iendR,jstrR:jendR,:) )
            dmsp_p_uptake_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*dmsp_p_uptake_avg(istrR:iendR,jstrR:jendR,:)
     &           + dmsp_p_uptake_hist(istrR:iendR,jstrR:jendR,:) )
            dmsp_lysis_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*dmsp_lysis_avg(istrR:iendR,jstrR:jendR,:)
     &           + dmsp_lysis_hist(istrR:iendR,jstrR:jendR,:) )
            dms_prod_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*dms_prod_avg(istrR:iendR,jstrR:jendR,:)
     &           + dms_prod_hist(istrR:iendR,jstrR:jendR,:) )
            dms_photlys_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*dms_photlys_avg(istrR:iendR,jstrR:jendR,:)
     &           + dms_photlys_hist(istrR:iendR,jstrR:jendR,:) )
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
            ammox_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*ammox_avg(istrR:iendR,jstrR:jendR,:)
     &           + ammox_hist(istrR:iendR,jstrR:jendR,:) )

            denitr_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*denitr_avg(istrR:iendR,jstrR:jendR,:)
     &           + denitr_hist(istrR:iendR,jstrR:jendR,:) )

            n2o_prod_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*n2o_prod_avg(istrR:iendR,jstrR:jendR,:)
     &           + n2o_prod_hist(istrR:iendR,jstrR:jendR,:) )
            n2_prod_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*n2_prod_avg(istrR:iendR,jstrR:jendR,:)
     &           + n2_prod_hist(istrR:iendR,jstrR:jendR,:) )

            denitr_sed_avg(istrR:iendR,jstrR:jendR) = 
     &           cff*(cff1*denitr_sed_avg(istrR:iendR,jstrR:jendR)
     &           + denitr_sed_hist(istrR:iendR,jstrR:jendR) )
# endif /* NCYCLE_ANOXIC */

# ifdef N15_CYCLE
            no3_v_sp_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*no3_v_sp_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + no3_v_sp_n15_hist(istrR:iendR,jstrR:jendR,:) )

            nh4_v_sp_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*nh4_v_sp_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + nh4_v_sp_n15_hist(istrR:iendR,jstrR:jendR,:) )

            no3_v_diat_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*no3_v_diat_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + no3_v_diat_n15_hist(istrR:iendR,jstrR:jendR,:) )

            nh4_v_diat_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*nh4_v_diat_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + nh4_v_diat_n15_hist(istrR:iendR,jstrR:jendR,:) )

            graze_sp_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*graze_sp_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + graze_sp_n15_hist(istrR:iendR,jstrR:jendR,:) )

            graze_diat_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*graze_diat_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + graze_diat_n15_hist(istrR:iendR,jstrR:jendR,:) )

            graze_diaz_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*graze_diaz_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + graze_diaz_n15_hist(istrR:iendR,jstrR:jendR,:) )

            sp_loss_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*sp_loss_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + sp_loss_n15_hist(istrR:iendR,jstrR:jendR,:) )

            diat_loss_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*diat_loss_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + diat_loss_n15_hist(istrR:iendR,jstrR:jendR,:) )

            diaz_loss_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*diaz_loss_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + diaz_loss_n15_hist(istrR:iendR,jstrR:jendR,:) )

            zoo_loss_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*zoo_loss_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + zoo_loss_n15_hist(istrR:iendR,jstrR:jendR,:) )

            sp_agg_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*sp_agg_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + sp_agg_n15_hist(istrR:iendR,jstrR:jendR,:) )

            diat_agg_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*diat_agg_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + diat_agg_n15_hist(istrR:iendR,jstrR:jendR,:) )

            diaz_Nfix_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*diaz_Nfix_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + diaz_Nfix_n15_hist(istrR:iendR,jstrR:jendR,:) )

            pon_remin_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*pon_remin_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + pon_remin_n15_hist(istrR:iendR,jstrR:jendR,:) )

            don_prod_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*don_prod_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + don_prod_n15_hist(istrR:iendR,jstrR:jendR,:) )

            don_remin_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*don_remin_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + don_remin_n15_hist(istrR:iendR,jstrR:jendR,:) )

            nitrif_n15_avg(istrR:iendR,jstrR:jendR,:) =  
     &           cff*(cff1*nitrif_n15_avg(istrR:iendR,jstrR:jendR,:) 
     &           + nitrif_n15_hist(istrR:iendR,jstrR:jendR,:) )
#  ifdef NCYCLE_ANOXIC
            ammox_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*ammox_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + ammox_n15_hist(istrR:iendR,jstrR:jendR,:) )

            denitr_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*denitr_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + denitr_n15_hist(istrR:iendR,jstrR:jendR,:) )

            n2_prod_n15_avg(istrR:iendR,jstrR:jendR,:) = 
     &           cff*(cff1*n2_prod_n15_avg(istrR:iendR,jstrR:jendR,:)
     &           + n2_prod_n15_hist(istrR:iendR,jstrR:jendR,:) )

            denitr_sed_n15_avg(istrR:iendR,jstrR:jendR) = 
     &           cff*(cff1*denitr_sed_n15_avg(istrR:iendR,jstrR:jendR)
     &           + denitr_sed_n15_hist(istrR:iendR,jstrR:jendR) )
#  endif /* NCYCLE_ANOXIC */
# endif /* N15_CYCLE */
         end if ! iic.gt.nts_bgc_flux_avg

         return
         end

#else
      subroutine set_bec_flux_avg_empty
      end
#endif /* BIOLOGY_BEC && BGC_FLUX_ANALYSIS */
