#include "cppdefs.h"

#if defined BIOLOGY_BEC && defined BULK_FLUX
      subroutine det_srflx_dailyavg_tile(Istr,Iend,Jstr,Jend)
      implicit none
#include "param.h"
#include "scalars.h"
#include "forces.h"
#include "param_bec.h"
#include "ecosys_bec.h"

      integer Istr,Iend,Jstr,Jend
      integer i,j,k,l
      real frac_iptr

! first determine if the array has been completely filled with values
      do j = Jstr, Jend
         do i = Istr, Iend
            if (iptr_srflx_day_set(i,j) .eq. 0) then
               if (PARinc_rst_read .eq. 1) then
! array has not yet been fully set up
! if PAR was read in from a restart file, 
! use those values after rescaling
                  srflx_dailyavg(i,j) = PARinc_rst(i,j) / (0.45*rho0*Cp)
               else
! average over the values that have been filled in so far
                  srflx_dailyavg(i,j) = c0
! it will be 0 during the initial call to ecosys_bec; avoid /0
                  if (iptr_srflx_day(i,j) .gt. 0) then
                     frac_iptr = 1.0 / real(iptr_srflx_day(i,j))
                  else
                     frac_iptr = 1.0
                  end if
                  do l = 1, iptr_srflx_day(i,j)
                     srflx_dailyavg(i,j) = srflx_dailyavg(i,j) + 
     &                    srflx_day(i,j,l)
                  end do ! l
                  srflx_dailyavg(i,j) = srflx_dailyavg(i,j) * frac_iptr
               end if
            else 
! array has been fully set up, compute average
               srflx_dailyavg(i,j) = c0
               do l = 1, num_srflx_day
                  srflx_dailyavg(i,j) = srflx_dailyavg(i,j) + 
     &                 srflx_day(i,j,l)
               end do ! l
               srflx_dailyavg(i,j) = srflx_dailyavg(i,j) * srflx_frac
            end if
         end do ! i
      end do ! j
      return
      end subroutine det_srflx_dailyavg_tile
#else  /* BIOLOGY_BEC && BULK_FLUX */
      subroutine det_srflx_dailyavg_empty
      end subroutine det_srflx_dailyavg_empty
#endif /* BIOLOGY_BEC && BULK_FLUX */
