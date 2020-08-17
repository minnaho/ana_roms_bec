/* This is include file "averages.h": time-averaged fields
 for output purposes:
*/
#ifdef AVERAGES
      real zeta_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE zeta_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real ubar_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  ubar_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real vbar_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  vbar_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      common /avg_zeta/zeta_avg /avg_ubar/ubar_avg
     &                          /avg_vbar/vbar_avg
# ifdef OUTPUT_AVG_PHYS_BIO
      real zeta_bio_avg(GLOBAL_2D_ARRAY)
      common /avg_zeta_bio/zeta_bio_avg
# endif /* OUTPUT_AVG_PHYS_BIO */
# ifdef SOLVE3D
      real u_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE u_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real v_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE v_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real t_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE t_avg(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real rho_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE rho_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real w_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE w_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_u/u_avg /avg_v/v_avg /avg_t/t_avg
     &                /avg_rho/rho_avg /avg_w/w_avg
#  ifdef OUTPUT_SST_AVG
      real sst_avg(GLOBAL_2D_ARRAY)
      common /avg_sst/ sst_avg
#  endif
      real wvl_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE wvl_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_wvl/wvl_avg
      real akv_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE akv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real akt_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE akt_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_akv/akv_avg /avg_akt/akt_avg
#  ifdef SALINITY
      real aks_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE aks_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_aks/aks_avg
#  endif
#  ifdef LMD_KPP
      real hbl_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE hbl_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      common /avg_hbl/hbl_avg
#  endif
#  ifdef LMD_BKPP
      real hbbl_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE hbbl_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_hbbl/hbbl_avg
#  endif
#  ifdef BULK_DIAGS
      real sustr_blk_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE sustr_blk_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real svstr_blk_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE svstr_blk_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real shflx_net_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_net_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real shflx_lat_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_lat_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real shflx_sen_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_sen_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real shflx_rad_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_rad_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real swflx_emp_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE swflx_emp_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real shflx_wwk_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_wwk_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real surf_u_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE surf_u_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real surf_v_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE surf_v_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      common /avg_sustr_blk/sustr_blk_avg
      common /avg_svstr_blk/svstr_blk_avg
      common /avg_shflx_net/shflx_net_avg
      common /avg_shflx_lat/shflx_lat_avg
      common /avg_shflx_sen/shflx_sen_avg
      common /avg_shflx_rad/shflx_rad_avg
      common /avg_swflx_emp/swflx_emp_avg
      common /avg_shflx_wwk/shflx_wwk_avg
      common /avg_surf_u/surf_u_avg
      common /avg_surf_v/surf_v_avg
#  endif
              
# endif
#endif /* AVERAGES */
 
