#ifdef REPLACE 
      program mreplace
#else
      program mgrep
#endif

! Either a basic tool to rename objects (usually variable names and CPP
! macros) consistently throughout a set of multiple files, or to detect
! the presence on named object within a file or files.  The usage is as
! follows:
!              mgrep name file1 [file2 ...]
! and
!              mreplace old_name new_name file1 [file2 ...]
!
!
! where "name" (also "old_name") is the name of the presumably existing
! object to be detected (replaced by "new_name"); all the remaining
! arguments "file1", "file2" ... are the file names (it is acceptable
! to use wildcards to form the list).

! Operation "mgrep" is very similar to standard UNIX grep command,
! except that instead matching character string patterns, it responds
! to the names following syntax conventions of  Fortran and C, i.e.,
! words containing only upper- and lower- case letters A-Z, digits 0-9,
! and underscore. [Very close to the behavior of grep with -w flag.]
 
! This program works as follows: it sequentially reads the files
! line-by-line and searches for word matching "old_name".  Matching
! means that the entire word is to be matched, not just a sub-sequence
! of characters within it (e.g., "var_1" and "var_1x" are considered
! different names). The search algorithm is the same for "mgrep" and
! "mreplace", hence most of the code is the same; "mreplace" writes
! a temporal file made of the lines with "old_name" replaced by
! "new_name", if matches found, or the original lines. Before start
! processing an individual file it sets a flag "change_file" to .false.
! and then sets it to .true. when match is found. If at the end of
! file processing the flag is .true., the original file is replaced
! by temporal file using system call; Otherwise the temporal file is
! simply discarded, and the original file remains untouched (i.e.,
! retaining its date and attributes).


      implicit none
      integer, parameter :: max_fname_size=80, max_string_size=512,
     &                            max_wsize=32, maxwords=32, iin=11
      logical change_file, switch
      integer nargs, num_changed_files, is(maxwords), ie(maxwords),
     &        nwrds, n,i,m, lold, lfnm, lstr, ierr, iargc, lenstr
      character(len=max_wsize) old_name
      character(len=max_fname_size) fname
      character(len=2*max_string_size) string
      
#ifdef REPLACE
      integer, parameter :: iout=12
      character(len=17), parameter :: tmpfile='/tmp/mreplace.tmp'
      
      character(len=max_wsize) new_name
      integer lnew
#else
      logical change_string
#endif
      num_changed_files=0

      nargs=iargc()
#ifdef REPLACE
      if (nargs.gt.2) then
#else
      if (nargs.gt.1) then
#endif
        call getarg(1,old_name)
        lold=lenstr(old_name)
#ifdef REPLACE
        call getarg(2,new_name)
        lnew=lenstr(new_name)
        do n=3,nargs
#else
        do n=2,nargs
#endif
          do i=1,max_fname_size
            fname(i:i)=' '
          enddo
          call getarg(n,fname)
          lfnm=lenstr(fname)
          open(unit=iin,  file=fname(1:lfnm),   form='formatted',
     &                                status='old', iostat=ierr)
          if (ierr.ne.0) then
            write(*,'(/1x,3A/)') '### ERROR: cannot open file ''',
     &                                       fname(1:lfnm), '''.'
            stop
          endif

#ifdef REPLACE
          open(unit=iout, file=tmpfile, form='formatted')
#endif
          change_file=.false.
          do while (ierr.eq.0)
            do i=1,max_string_size
              string(i:i)=' '
            enddo
            read(iin,'(A)',iostat=ierr) string
            lstr=max_string_size
            do while(string(lstr:lstr).eq.' ' .and. lstr.gt.1)
              lstr=lstr-1
            enddo
            if (string(lstr:lstr).eq.' ') lstr=lstr-1
            if (lstr.gt.0) then
              nwrds=0
              switch=.true.
#ifndef REPLACE
              change_string=.false.
#endif
              do i=1,lstr
                if ((string(i:i).ge.'A' .and. string(i:i).le.'Z') .or.
     &              (string(i:i).ge.'a' .and. string(i:i).le.'z')  .or.
     &              (string(i:i).ge.'0' .and. string(i:i).le.'9')  .or.
     &                                        string(i:i).eq.'_') then
                  if (switch) then
                    switch=.false.        ! Logical variable "switch"
                    nwrds=nwrds+1         ! goes T -> F when entering
                    is(nwrds)=i           ! a legal Fortran name, i.e.
                  endif                   ! a word made of permitted
                  ie(nwrds)=i             ! characters (upper- and
                else                      ! lowercase letters, digits,
                  switch=.true.           ! and underscores).
                endif                     ! It changes back F -> T
              enddo                       ! upon exiting such word.

              do m=nwrds,1,-1
                if (ie(m)-is(m)+1.eq.lold) then
                  if (string(is(m):ie(m)) .eq. old_name(1:lold)) then
#ifdef REPLACE
                    if (lnew.ne.lold) then
                      string(is(m)+lnew:lstr+lnew-lold)
     &                            =string(ie(m)+1:lstr)
                      lstr=lstr+lnew-lold
                    endif
                    string(is(m):is(m)+lnew-1)=new_name(1:lnew)
#else
                    change_string=.true.
#endif
                    change_file=.true.
                  endif
                endif
              enddo
#ifdef REPLACE
              if (string(lstr:lstr).eq.' ') lstr=lstr-1
              write(iout,'(A)',iostat=ierr) string(1:lstr)
            else
              if (ierr.eq.0) write(iout,*)  !<-- to keep empty line as
#else
              if (change_string) write(*,'(3A)') fname(1:lfnm), ':',
     &                                           string(1:lstr)
#endif
            endif    !<-- lstr.gt.0          !   in the original file
          enddo    !<-- while(ierr.eq.0)
          close(iin)
          if (change_file) then
            num_changed_files=num_changed_files+1
#ifdef REPLACE
            close(iout)
            string='/bin/mv '/ /tmpfile/ /' '/ /fname(1:lfnm)
            lstr=lenstr(string)
c**         write(*,*) '''', string(1:lstr), ''''
            call system(string(1:lstr))
            write(*,*) fname(1:lfnm), ': renamed ''', old_name(1:lold),
     &                            ''' into ''', new_name(1:lnew), ''''
#else
c**         write(*,'(3(1x,A))') old_name(1:lold), 'found in',
c**  &                                          fname(1:lfnm)
#endif
#ifdef REPLACE
          else
            close(iout, status='delete')
# ifdef VERBOSE
            write(*,*) fname(1:lfnm), ': no changes'
# endif
#endif
          endif
        enddo   !<-- n=3,nargs, loop over arguments
#ifndef VERBOSE
        if (num_changed_files.eq.0) then
          write(*,'(1x,A,I4,1x,2A)') 'Processed', nargs-2, 'files. ',
     &                          'No matches found. No changes made.'
        endif
#endif
      else
#ifdef REPLACE
        write(*,'(/1x,3A/1x,2A,2/12x,A/2(/1x,2A)/1x,A/)')  'This is ',
     &   '"mreplace": a tool to rename objects (variable names, ',
     &   'macros, etc)', 'consistently throughout a set of multiple ',
     &   'files. The usage should be:',
     &              'mreplace old_name new_name file1 [file2 ...]',
     &   'where "old_name" is the name of an object to be renamed ',
     &   'into "new_name";',  '"file1",  "file2" ... are the names ',
     &   'of the files to be operated on', '(wildcards are accepted).'
#else
        write(*,'(/1x,3A/1x,2A,2/12x,A,2/8(1x,2A/))')  'This is ',
     &   '"mgrep": a tool to identify the presence of named objects ',
     &  '(variable', 'names, macros, etc) within a file or multiple ',
     &   'files. The usage should be:',
     &              'mgrep name file1 [file2 ...]',
     &   'where "name" is the name of an object to be detected ',
     &   '"file1", "file2"... are',  'the name[s] of the file[s] ',
     &   'to be checked (wildcards are accepted within the',
     &   'list). This operation does not make any changes within ',
     &   'the files and produces', 'output similar to that of the ',
     &   'standard UNIX grep command.  However, unlike the',
     &   'standard "grep" which matches character string patterns, ',
     &   '"mgrep" responds to',  '"legal" variable names, as they ',
     &   'are defined under Fortran of C syntax rules,',
     &   'e.g., if a file contains objects named ABC, BC, BCD, ',
     &   'and ABCDE, then command', '"grep BC filename" will ',
     &   'highlight all of them while mgrep highlights only BC.'
#endif
      endif
      stop
      end
 
#include "lenstr.F"
