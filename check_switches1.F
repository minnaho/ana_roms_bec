#include "cppdefs.h"


      subroutine check_switches1 (ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new CPP-switches  !!!!!!
!!!!!! were introduced into "cppdefs.h".  NO ACTION IS NEEDED  !!!!!!
!!!!!! if changes in "cppdefs.h" are limited to activation or  !!!!!!
!!!!!! deactivation of previously known switches.              !!!!!!
!!!!!! To refresh this file compile and execute "cppcheck.F"   !!!!!!
!!!!!! as an independent program, or use commands              !!!!!!
!!!!!! "make checkdefs" or "make depend".                      !!!!!!
!!!!!! Number of Configuration Choices:                      1 !!!!!!
!!!!!! Total number of CPP-switches:                        49 !!!!!!

      implicit none
      integer ierr, is,ie, iexample
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      do is=1,max_opt_size
        cpps(is:is)=' '
      enddo
      iexample=0
!
! cppdefs.h                       
!
      is=1
      ie=is+10
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='<cppdefs.h>'
      cpps(ie+1:ie+1)=' '

#ifdef FLAT_BOTTOM
      iexample=iexample+1
      is=ie+2
      ie=is+10
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='FLAT_BOTTOM'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SOLVE3D
      is=ie+2
      ie=is+6
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='SOLVE3D'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_COR
      is=ie+2
      ie=is+5
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='UV_COR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_ADV
      is=ie+2
      ie=is+5
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='UV_ADV'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CURVGRID
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='CURVGRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef MASKING
      is=ie+2
      ie=is+6
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='MASKING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef MASK_LAND_DATA
      is=ie+2
      ie=is+13
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='MASK_LAND_DATA'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VERT_COORD_TYPE_SM09
      is=ie+2
      ie=is+19
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='VERT_COORD_TYPE_SM09'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef AVERAGES
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='AVERAGES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SALINITY
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='SALINITY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NONLIN_EOS
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='NONLIN_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPLIT_EOS
      is=ie+2
      ie=is+8
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='SPLIT_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef WIND_STRESS_AMP
      is=ie+2
      ie=is+14
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='WIND_STRESS_AMP'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_VIS2
      is=ie+2
      ie=is+6
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='UV_VIS2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VIS_GRID
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='VIS_GRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ADV_ISONEUTRAL
      is=ie+2
      ie=is+13
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ADV_ISONEUTRAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TS_DIF2
      is=ie+2
      ie=is+6
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='TS_DIF2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DIF_GRID
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='DIF_GRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_MIXING
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='LMD_MIXING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BULK_FLUX
      is=ie+2
      ie=is+8
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='BULK_FLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DAILYPAR_BEC
      is=ie+2
      ie=is+11
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='DAILYPAR_BEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef COUPLED_SURF_CURR
      is=ie+2
      ie=is+16
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='COUPLED_SURF_CURR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef WND_AT_RHO_POINTS
      is=ie+2
      ie=is+16
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='WND_AT_RHO_POINTS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BULK_DIAGS
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='BULK_DIAGS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DIURNAL_SRFLUX
      is=ie+2
      ie=is+13
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='DIURNAL_SRFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SALINITY_RESTORING
      is=ie+2
      ie=is+17
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='SALINITY_RESTORING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef QCORRECTION
      is=ie+2
      ie=is+10
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='QCORRECTION'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef STR_AT_RHO_POINTS
      is=ie+2
      ie=is+16
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='STR_AT_RHO_POINTS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SMFLUX
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_SMFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_STFLUX
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_STFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SSFLUX
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_SSFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SRFLUX
      is=ie+2
      ie=is+9
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_SRFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_INITIAL
      is=ie+2
      ie=is+10
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_INITIAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_GRID
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_GRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_WEST
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_WEST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_SOUTH
      is=ie+2
      ie=is+8
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_SOUTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_NORTH
      is=ie+2
      ie=is+8
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_NORTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_EAST
      is=ie+2
      ie=is+7
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_EAST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M2ORLANSKI
      is=ie+2
      ie=is+13
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_M2ORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M3ORLANSKI
      is=ie+2
      ie=is+13
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_M3ORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_TORLANSKI
      is=ie+2
      ie=is+12
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='OBC_TORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef PASSIVE_TRACER
      is=ie+2
      ie=is+13
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='PASSIVE_TRACER'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef PSOURCE
      is=ie+2
      ie=is+6
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='PSOURCE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_PSOURCE
      is=ie+2
      ie=is+10
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='ANA_PSOURCE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VFLX_CORR
      is=ie+2
      ie=is+8
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='VFLX_CORR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LARGE_HIS_FILES
      is=ie+2
      ie=is+14
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='LARGE_HIS_FILES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BIOLOGY_BEC
      is=ie+2
      ie=is+10
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='BIOLOGY_BEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BEC_LIANG
      is=ie+2
      ie=is+8
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='BEC_LIANG'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NCYCLE_ANOXIC
      is=ie+2
      ie=is+12
      if (ie.ge.max_opt_size) goto 99
      cpps(is:ie)='NCYCLE_ANOXIC'
      cpps(ie+1:ie+1)=' '
#endif
      if (iexample.eq.0) then
        mpi_master_only write(*,'(/1x,A)')
     & '### ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
        mpi_master_only write(*,'(/1x,2A/)')
     &  '### ERROR: More than one main configuration ',
     &             'switch is defined in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  mpi_master_only write(*,'(/1x,2A/14x,A)')
     &  '### ERROR: checkdefs :: Insufficient size of string "cpps" ',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
