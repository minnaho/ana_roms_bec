      program ncextract

! Purpose: to exstract a subset of declarations from "netcdf.inc",
! which are relevant to ROMS code. To do so, first search through the
! whole ROMS code (that is all *.F and *.h files) and build catalog
! of all "nf_" words which occur in the code. Once the catalog is
! built, it read the standard "netcdf.inc" file line-by-line and save
! only lines, which contains nf_words which are present in the
! catalog.

! NOTE: This measure is provided purely for convenience just to speed
! up things a little bit by shortening the "netcdf.inc" include file
! (typically from 1300+ lines of the standard "netcdf.inc" file down
! to 140 lines, which leads to more compact and easier to search .f
! files).  It does not affect the result of compilation of the entire
! model. Some manual editing of the resultant "netcdf.inc.ext" file
! may ne needed in order to eliminate conflicts during subsequent
! compilation --- typically eliminate usage of few obsolecent ncxxxxx
! parameters, which are provided by UNIDATA for compatibility with
! older versions of netCDF library.

      implicit none
      integer max_string, max_size, max_found, input, iout
      parameter (max_string=72,  max_size=32,  max_found=512,
     &                                     input=11, iout=12)
      character string*(max_string), buff*(max_string),
     &  quote*1, double_quote*1, nfword(max_found)*(max_size)
      integer nlen(max_found), nfound, iargc, last_arg, n,
     &    case_fold, i,j, ie,is, istrt,length, line,iocheck
      logical not_end_of_file, lswtch, lprint
 
      nfound=0
      quote=char(39)
      double_quote=char(34)
      case_fold=ichar('a')-ichar('A')
 
      last_arg=iargc()
      if (last_arg.gt.0) then
        do n=1,last_arg+1
          if (n.ne.last_arg+1) then
            call getarg(n, string)
          else
            string='/usr/local/include/netcdf.inc'
            open (unit=iout, file='netcdf.inc.ext',form='formatted')
          endif
          open(unit=input, file=string, form='formatted')
          line=0
          not_end_of_file=.true.
          do while (not_end_of_file)
            line=line+1
            string(1:1)=' '
            read(input,'(A)',iostat=iocheck,end=1) string
            goto 2
   1        not_end_of_file=.false.
   2        length=0
            lswtch=.false.
            do i=1,max_string,+1
              if (string(i:i).eq.quote) then
                lswtch=.not.lswtch
              elseif (string(i:i).eq.'!' .and. .not.lswtch) then
                goto 3
              endif
              if (string(i:i).ne.' ') length=i
            enddo
   3        continue
            do while (length.gt.0 .and. string(length:length).eq.' ')
              length=length-1
            enddo
            istrt=1
            do while (istrt.lt.length.and.string(istrt:istrt).eq.' ')
              istrt=istrt+1
            enddo
!
! Save the non-modified copu of the string into buffer (only for
! "netcdf.inc" file, then to simplify comparison and extraction fold
! cases to lower case and suppress all symbols except letters, digits
! and underscore.
!
            if (n.eq.last_arg+1) buff=string(1:length)
            do i=istrt,length
              if (string(i:i).ge.'A' .and. string(i:i).le.'Z') then
                string(i:i)=char(ichar(string(i:i))+case_fold)
              elseif (string(i:i).ge.'a'.and.string(i:i).le.'z') then
 
              elseif (string(i:i).ge.'0'.and.string(i:i).le.'9') then
 
              elseif (string(i:i).eq.'_') then
 
              else
                string(i:i)=' '
              endif
            enddo
!
! Search for nf_words inside the string; once found, check whether it
! is already in the catalog, and if it is not, put it in (all passes,
! except the last, where input file is "netcdf.inc") during the last
! pass, if, as opposite, nf_word is in catalog, set switch to print
! the whole non-modified version of the line into the new file.
!
            lprint=.false.
            i=istrt
            do while (i.lt.length)
              if (string(i:i+2).eq.'nf_') then
                is=i
              else
                is=0
              endif
              i=i+3
              do while (string(i:i).ne.' ')
                i=i+1
              enddo
              if (is.gt.0) then
                ie=i-1
                lswtch=.true.
                do j=1,nfound
                  if (ie-is+1.eq.nlen(j)) then
                    if (string(is:ie).eq.nfword(j)(1:nlen(j))) then
                      lswtch=.false.
                    endif
                  endif
                enddo
                if (lswtch .and. n.ne.last_arg+1) then
                  nfound=nfound+1
                  nlen(nfound)=ie-is+1
                  nfword(nfound)=string(is:ie)
                  write(*,'(3A)') '''',
     &              nfword(nfound)(1:nlen(nfound)), ''''
                elseif (.not.lswtch .and. n.eq.last_arg+1) then
                  lprint=.true.
                endif
              endif
              do while (string(i:i).eq.' ' .and. i.lt.length)
                i=i+1
              enddo
            enddo
            if (lprint) write(iout,'(A)') buff(1:length)
          enddo
          close (unit=input)
        enddo
        close (unit=iout)
      else
        write(*,'(/1x,2A/, /8x,A/)') '#### ERROR: No filenames ',
     &          'are specified. Usage of ncextract should be',
     &          'ncextract *.F *.h'
      endif
      stop
      end
 
 
 
