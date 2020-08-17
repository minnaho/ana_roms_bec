#include "cppdefs.h"

      subroutine closecdf       ! Close output netCDF files.
      implicit none
      integer ierr
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"
#ifdef FLOATS
# include "ncvars_floats.h"
#endif
 
      if (ncidrst.ne.-1) ierr=nf_close(ncidrst)
      if (ncidhis.ne.-1) ierr=nf_close(ncidhis)
#ifdef AVERAGES
      if (ncidavg.ne.-1) ierr=nf_close(ncidavg)
#endif
#ifdef FLOATS
      if (ncidflt.ne.-1) ierr=nf_close(ncidflt)
#endif
#if defined DIAGNOSTICS_EDDY
      if (nciddiags_eddy.ne.-1) ierr=nf_close(nciddiags_eddy)
# ifdef AVERAGES
      if (nciddiags_eddy_avg.ne.-1) ierr=nf_close(nciddiags_eddy_avg)
# endif
#endif

# ifdef MPI_SILENT_MODE
      if (mynode.eq.0) then
# endif
        write(*,'(/1x,A,6(1x,A,I5)/)') 'MAIN: DONE - Records written:',
     &                     'history:', nrechis ,  'restart:', nrecrst
#ifdef AVERAGES
     &                                         ,  'averages:', nrecavg
#endif
#ifdef FLOATS
     &                                         ,   'floats:', nrecflt
#endif
#ifdef STATIONS
     &                                         ,  'stations:', nrecstn
#endif
     &                                             MYID

        if (may_day_flag.eq.1) then
          write(*,'(/1x,A/)') 'MAIN: Abnormal termination: BLOWUP'
        elseif (may_day_flag.eq.2) then
          write(*,'(/1x,2A/)')  '### ERROR: Cause of termination: ',
     &                                            'netCDF input.'
        elseif (may_day_flag.eq.3) then
          write(*,'(/1x,2A/)')  '### ERROR: Cause of termination: ',
     &                                           'netCDF output.'
        elseif (may_day_flag.eq.4) then
          write(*,'(/1x,A/)') '### ERROR: Cannot open netCDF file.'
        elseif (may_day_flag.eq.5) then
          write(*,'(/1x,A/)') '### ERROR: Error in input script.'
        elseif (may_day_flag.eq.6) then
          write(*,'(/1x,2A/)') '### ERROR: read_inp :: An input ',
     &                                   'file is not available.'
        elseif (may_day_flag.eq.7) then
          write(*,'(/1x,A/)') '### ERROR: KPP algoritm failure.'
        elseif (may_day_flag.eq.8) then
          write(*,'(/1x,A/)') '### ERROR: Wrong choice of tiling.'
        elseif (may_day_flag.eq.9) then
          write(*,'(/1x,2A/)') '### ERROR: Wrong choice of ',
     &                               'S-coordinate parameters.'
        endif
# ifdef MPI_SILENT_MODE
      endif
# endif
      return
      end
