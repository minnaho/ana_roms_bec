      program nsub

! Takes one argument which is expected to be in MxN format where M,N
! are two integer numbers separated by lowercase "x" (e.g., 4x32) and
! adjusts parameter settings in "param.h" as
!
!              NSUB_X=M,  NSUB_E=M,
!
! affecting only lines not commented out by FORTRAN comment rules
! (i.e, with "C", "c" in first position, or "!" in front), however,
! if there are different settings controlled by CPP-branches, ALL
! non-commented lines with NSUB_X=..., NSUB_E=... will be replaced.
!
! Basically this is automatic search-and-replace string editing.
! This code can to work within an automatic script to facilitate
! compiling multiple executables, when it is not practical to do it
! manually, e.g.,
!
!        #/bin/tcsh
!        foreach CASE (2x2 2x4 2x8 2x16 3x4 3x8 3x16 3x32 3x38 3x46)
!          nsub $CASE
!          make clean ; make ; mv roms roms.$CASE
!        end
!
! Useful for cache performance testing or whatever other purposes.
! Activating TEST_ONLY CPP-switch prevents overwriting "param.h" --
! a safety feature -- in this case the output remains in tmpfile
! with name specified as parameter below.

#define TEST_ONLY
      implicit none
      character(len=7), parameter :: paramfile='param.h',
     &                                 tmpfile='nsb.tmp'
      integer, parameter :: max_str=80
      character(len=max_str) str
      character(len=4) nsub_x,nsub_e
      integer ierr, i,j,ixe, is,ie, ipass, lstr, iargc,lenstr

      lstr=iargc()
      if (lstr.gt.0) then
        do i=1,max_str
          str(i:i)=' '
        enddo
        call getarg (1,str)
        lstr=lenstr(str)
        j=0
        ierr=0
        do i=1,lstr
          if (j.eq.0 .and. str(i:i).eq.'x') then
            j=i
          elseif (str(i:i).lt.'0' .or. str(i:i).gt.'9') then
            ierr=ierr+1
          endif
        enddo
        if (j.gt.1 .and. ierr.eq.0) then
          nsub_x=str(1:j-1)
          nsub_e=str(j+1:lstr)

          open(unit=11, file=paramfile, status='old')
          open(unit=12, file=tmpfile)
          do while (ierr.eq.0)
            do i=1,max_str
              str(i:i)=' '    !<-- empty-out string before reading
            enddo                  
            read(11,'(A)',iostat=ierr) str
            lstr=max_str
            do while(str(lstr:lstr).eq.' ' .and. lstr.gt.1)
              lstr=lstr-1
            enddo
            if (str(lstr:lstr).eq.' ') lstr=0  !<-- blank string 
            if (lstr.gt.0) then
              i=1
              do while(str(i:i).eq.' ' .and. i.lt.lstr)
                i=i+1
              enddo
              if (str(i:i).ne.'!' .and. str(1:1).ne.'c' .and.
     &                                  str(1:1).ne.'C') then
                do ipass=1,2
                  j=index(str(i:lstr),'NSUB_')
                  if (j.gt.0) then
                    i=i+j+4
                    if (str(i:i).eq.'X') then
                      ixe=1 
                    elseif (str(i:i).eq.'E') then
                      ixe=2 
                    else
                      ixe=0
                    endif
                    if (ixe.gt.0) then

! Search for "=" sign followed by an integer number...

                      i=i+1
                      j=index(str(i:lstr),'=')
                      if (j.gt.0) then
                        i=i+j

                        do while(str(i:i).eq.' ' .and. i.lt.lstr)
                          i=i+1
                        enddo
                        is=0
                        ie=0
                        if (str(i:i).ge.'0' .and. str(i:i).le.'9') then
                          is=i
                          i=i+1
                          do while(i.lt.lstr .and. str(i:i).ge.'0'
     &                                       .and. str(i:i).le.'9')
                            i=i+1
                          enddo
                          if (str(i:i).ge.'0'.and.str(i:i).le.'9') then
                            ie=i
                          else
                            ie=i-1
                          endif

! At this moment str(is:ie) contains the old parameter setting (an
! integer number) which needs to be substituted with the new one: 

                          if (ixe.eq.1) then
                            j=lenstr(nsub_x)
                            str=str(1:is-1)/ /nsub_x(1:j)/
     &                                       /str(ie+1:lstr)
                          elseif (ixe.eq.2) then 
                            j=lenstr(nsub_e)
                            str=str(1:is-1)/ /nsub_e(1:j)/ 
     &                                       /str(ie+1:lstr)
                          endif
                          lstr=lstr +j-(ie-is+1) !<-- adjust length
                        endif
                      endif
                    endif
                  endif
                enddo  !<-- ipass, keep "i" to roll-over for next ipass
              endif
              write(12,'(A)') str(1:lstr)
            else                              ! write an empty line if
              if (ierr.eq.0) write(12,'(A)')  ! lstr.eq.0 except at the
            endif           !<-- lstr.gt.0    ! very end of the file 

          enddo
          close(unit=11)
          close(unit=12)

#ifndef TEST_ONLY
          str='mv '/ /tmpfile/ /' '/ /paramfile
          lstr=lenstr(str)
c**       write(*,*) '''', str(1:lstr), ''''
          call system(str(1:lstr))
#endif
        elseif (j.eq.0) then
          write(*,'(/4A/)') '### ERROR: There is no ''x'' within the ',
     &                        'given argument ''',  str(1:lstr), '''.'
        elseif (ierr.ne.0) then
          write(*,'(/4A/)')  '### ERROR: An illegal character within ',
     &                         'the argument ''',  str(1:lstr),  '''.'
        endif
      else
        write(*,'(/1x,2A,2/16x,A,2/1x,2A/1x,2A/)')    'This operator ',
     &          'should be invoked as',  'nsub MxN',  'where M and N ',
     &          'are two integer numbers separated by ''x'';', 'nsub ',
     &          'edits "param.h" to set NSUB_X=M, NSUB_E=N there.'
      endif
      stop
      end

#include "lenstr.F"
