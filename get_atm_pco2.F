#include "cppdefs.h"

#ifdef READ_VARIABLE_ATM_PCO2
      function get_atm_pco2(tdays)
      implicit none
# include "param.h"
# include "forces.h"
# include "netcdf.inc"

      real get_atm_pco2, tdays, cff
      integer, parameter :: MAX_NTIME=840 ! 70 years
      real(kind=4) buffer(MAX_NTIME)  ! always float!
      integer ncid, ierr, lstr, lvar, lenstr, pco2_id, pco2_tid
      integer pco2_tdim, ntime, i
      character(len=9), parameter :: pco2_tvar='pco2_time'
      character(len=4), parameter :: pco2_var='pco2'
      integer start(1), count(1)

      lstr=lenstr(pco2_atm_file)
      ierr=nf_open(pco2_atm_file(1:lstr), nf_nowrite, ncid)
      if (ierr .ne. nf_noerr) then
         stop 'could not open pco2 file' ! WARN not the best error handling
      end if

      ierr=nf_inq_dimid(ncid, pco2_tvar,pco2_tdim)
      if (ierr .ne. nf_noerr) then
         stop 'could not get tdim in pco2 file' ! WARN not the best error handling
      end if

      ! inquire length of time dimension
      ierr=nf_inq_dimlen(ncid, pco2_tdim, ntime)
      if (ierr .ne. nf_noerr) then
         stop 'could not read tdim in pco2 file' ! WARN not the best error handling
      end if
      lvar=lenstr(pco2_var)
      ierr=nf_inq_varid (ncid, pco2_var(1:lvar), pco2_id)
      if (ierr .ne. nf_noerr) then
         stop 'could not read pco2 var from file' ! WARN not the best error handling
      end if
      ierr=nf_inq_varid (ncid, 'pco2_time', pco2_tid)
      if (ierr .ne. nf_noerr) then
         stop 'could not read pco2 time from file' ! WARN not the best error handling
      end if

      ! arrays are so small, it's easiest to read them all at 
      ! always using float for both variables
      start(1) = 1
      count(1) = ntime
      ierr=nf_get_vara_real(ncid,pco2_tid,start,count,buffer)
      if (ierr .ne. nf_noerr) then
         stop 'could not read pco2 time values from file' ! WARN not the best error handling
      end if

      i = 2
      do while (buffer(i) .lt. tdays)
!DEBUG    print *,'buffer,tdays:',buffer(i),tdays
         i = i+1
      end do
      cff = ( (tdays-buffer(i-1)) / (buffer(i) - buffer(i-1)) )
!DEBUG print *,'tdays,time,cff:',tdays,buffer(i-1:i),cff

      start(1) = i-1
      count(1) = 2
      ierr=nf_get_vara_real(ncid,pco2_id,start,count,buffer)
      if (ierr .ne. nf_noerr) then
         stop 'could not read pco2 values from file' ! WARN not the best error handling
      end if

      get_atm_pco2 = cff * (buffer(2) - buffer(1)) + buffer(1)

!DEBUG      print *,'pco2:',buffer(1:2),get_atm_pco2

      ierr = nf_close(ncid)

      end function get_atm_pco2
#else /* READ_VARIABLE_ATM_PCO2 */
      subroutine get_atm_pco2_empty
      end subroutine get_atm_pco2_empty
#endif /* READ_VARIABLE_ATM_PCO2 */
