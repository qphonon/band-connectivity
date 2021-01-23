! 
!     Copyright (C) 2021  Chee Kwan Gan (ihpcganck@gmail.com)
!     Copyright (C) 2020  Chee Kwan Gan (ihpcganck@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
! 

module commod
  use sgsym
  use derived_constants
  implicit none

  interface fargv
    module procedure getarg_str,getarg_int,getarg_double
  end interface

contains

  subroutine getarg_str(n,str)
    integer :: n
    character(len=*) :: str
    call get_command_argument(n,str)
  end subroutine getarg_str

  subroutine getarg_int(n,k)
    integer :: n, k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(I10)') k
  end subroutine getarg_int

  subroutine getarg_double(n,k)
    integer :: n
    real(double) :: k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(D22.15)') k
  end subroutine getarg_double


  subroutine inversemapping(n,fdindex,bdindex)
    integer :: n,i,ind
    integer :: fdindex(n)
    integer :: bdindex(n)
    do i = 1, n
      ind = fdindex(i)
      bdindex(ind) = i
    enddo
  end subroutine inversemapping

  function num2str(v,nspaces) result (str)
    integer :: tmpv,i,v,nspaces,digit,ind
    character(len=nspaces) :: str
    tmpv = v
    if(v < 0) then
      write(*,'(A)') 'err: v is negative.'
      stop 1
    endif
    do i = 1, nspaces
      digit = modulo(tmpv,10)
      tmpv = tmpv/10
      ind = nspaces-i+1
      str(ind:ind) = DIGITARR(digit)
    enddo
  end function num2str

  function str2N(str) result(n)
    character(len=*) :: str
    integer :: n
    read(str,'(I10)') n
  end function str2N

  function N2str(i) result (str)
    integer :: i,j,tmpv,intarr(20),nd,digit
    character(len=20) :: str

    if(i < 0) then
      tmpv = -i
    else
      tmpv = I
    endif

    nd = 0
    do
      digit = modulo(tmpv,10)
      nd = nd + 1
      intarr(nd) = digit
      tmpv = tmpv/10
      if(tmpv == 0) exit
    enddo
    str = ' '
    do j = 1, nd
      str(j:j) = digitarr(intarr(nd-j+1))
    enddo
    if(i < 0) then
      nd = nd + 1
      do j = nd,2,-1
        str(j:j) = str(j-1:j-1)
      enddo
      str(1:1) = '-'
    endif
  end function N2str

  subroutine fargn(n)
    integer :: n
    integer :: m

    m = command_argument_count()
    if(n /= m) then
      write(*,'(A,I4,A)') 'we need ',n,' arguments'
      write(*,'(A)') 'err: check number of arguments'
      stop 1
    endif
  end subroutine fargn

  subroutine lowercase(string)
    integer :: i,j
    character(len=*)  :: string
    do i = 1,len(string)
      j = index(upper,string(i:i))
      if(j /= 0) string(i:i) = lower(j:j)
    enddo
  end subroutine lowercase

  function num_of_strings(tmp)
    character(len=*) :: tmp
    integer :: num_of_strings
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100

    character(len=100) :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)

      if(io /= 0) then
        exit
      endif
    enddo
    num_of_strings = i-1
    if(num_of_strings == max_num_cols) then

    endif
  end function num_of_strings

  function num_of_integers(tmp)
    character(len=*) :: tmp
    integer :: num_of_integers
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100
    integer :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)
      if(io /= 0) then
        exit
      endif
    enddo
    num_of_integers = i-1
    if(num_of_integers == max_num_cols) then
      write(*,*) 'count_int_number= ',num_of_integers
      write(*,*) 'err: we may need to increase max_num_cols.'
      stop 1
    endif
  end function num_of_integers

  subroutine read_gendata(dstype,fileu,filename,nametag,v,dopt,vdefault,vec)
    integer,optional :: vec
    integer :: dopt
    integer :: dstype,fileu
    type(dtype) :: v,vdefault
    character(len=*) :: nametag,filename
    logical :: done
    integer :: ios,s
    character(len=DSL) :: str1,str2
    logical :: found
    integer :: i,j
    character(len=DSL) :: tmpstr

    if(present(vec)) then
      if(vec > DTYPEMAXLEN) then
        write(*,*) 'vec= ',vec
        write(*,*) 'DTYPEMAXLEN= ',DTYPEMAXLEN
        write(*,'(A)') 'err: vec > DTYPEMAXLEN.'
        stop 1
      endif
    endif

    found = .false.
    v%li = .false.
    v%lr = .false.
    v%ls = .false.
    if(dstype == itype) then
      v%li = .true.
    else if(dstype == rtype) then
      v%lr = .true.
    else if(dstype == stype) then
      v%ls = .true.
    endif

    write(*,'(A)') 'Search keyword=|'//trim(nametag)//'| in '//trim(filename)
    open(fileu,file=trim(filename),status='old',action='read')

    done = .false.
    do
      if(done) exit
      read(fileu,DSF,iostat=ios) str1
      if(ios /= 0) then
        done = .true.
      else
        call onespace(DSL,str1,str2)

        s = scan(str2,' ')
        if(s > 0) then
          if(str2(1:s-1) == '#keyword') then
            str2 = str2(s+1:DSL)
            s = scan(str2,' ')
            if(str2(1:s-1) == trim(nametag)) then
              done = .true.
              found = .true.
              str2 = str2(s+1:DSL)
              if(v%li) then

                write(*,'(A)') 'str2 is |'//trim(str2)//'|'
                if(present(vec)) then
                  read(str2,*) v%ivec(1:vec)
                  tmpstr = N2str(v%ivec(1))
                  do j = 2, vec
                    tmpstr = trim(tmpstr)//','//trim(N2str(v%ivec(j)))
                  enddo
                  write(*,'(A)') 'Integer variable '//trim(nametag)//' is |'//trim(tmpstr)//'|'
                else
                  read(str2,*) v%i
                  write(*,'(A)') 'Integer variable '//trim(nametag)//' is '//trim(N2str(v%i))
                endif
              else if(v%lr) then

                write(*,'(A)') 'str2 is |'//trim(str2)//'|'
                if(present(vec)) then
                  read(str2,*) v%rvec(1:vec)
                  write(*,*) 'Real variable '//trim(nametag)//' is ',v%rvec(1:vec)
                else
                  read(str2,*) v%r
                  write(*,*) 'Real variable '//trim(nametag)//' is ',v%r
                endif
              else if(v%ls) then

                write(*,'(A)') 'str2 is |'//trim(str2)//'|'
                if(present(vec)) then
                  read(str2,*) v%svec(1:vec)
                  write(*,*) 'String variable '//trim(nametag)//' is:'
                  do i = 1, vec
                    write(*,'(A)') trim(v%svec(i))
                  enddo
                else
                  read(str2,*) v%s
                  write(*,'(A)') 'String variable '//trim(nametag)//' is '//trim(v%s)
                endif
              endif
            endif
          else

          endif
        else
          done = .true.
        endif
      endif
    enddo
    close(fileu)

    if(.not. found .and. dopt == NONDEFAULTOPT) then
      write(*,*) 'cannot find the keyword='//trim(nametag)
      write(*,'(A)') 'err: cannot find the the keyword.'
      stop 1
    else if(.not. found .and. dopt == DEFAULTOPT) then
      write(*,'(/A)') 'Could not find '//trim(nametag)
      if(v%li) then
        if(present(vec)) then
          v%ivec(1:vec) = vdefault%ivec(1:vec)
          write(*,*) 'Warning: assume default value of ',v%ivec(1:vec)
        else
          v%i = vdefault%i
          write(*,*) 'Warning: assume default value of ',v%i
        endif
      else if(v%lr) then
        if(present(vec)) then
          v%rvec(1:vec) = vdefault%rvec(1:vec)
          write(*,*) 'Warning: assume default value of ',v%rvec(1:vec)
        else
          v%r = vdefault%r
          write(*,*) 'Warning: assume default value of ', v%r
        endif
      else if(v%ls) then
        if(present(vec)) then
          v%svec(1:vec) = vdefault%svec(1:vec)
          write(*,*) 'Warning: assume default value of ',v%svec(1:vec)
        else
          v%s = vdefault%s
          write(*,'(A)') 'Warning: assume default value of '//trim(v%s)
        endif
      endif
    endif
  end subroutine read_gendata

  subroutine getlenang(a,L)
    real(double) :: a(3,3),L(6)
    real(double) :: sintheta,costheta

    L(1) = vecmag3(a(1:3,1))
    L(2) = vecmag3(a(1:3,2))
    L(3) = vecmag3(a(1:3,3))

    costheta = dotprod3(a(1:3,2),a(1:3,3))/(L(2)*L(3))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(4) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,3),a(1:3,1))/(L(3)*L(1))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(5) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,1),a(1:3,2))/(L(1)*L(2))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(6) = atan2(sintheta,costheta)*rad2deg
  end subroutine getlenang

  subroutine real2recip(reallatt,reciplatt)
    real(double) :: reallatt(3,3),reciplatt(3,3)
    real(double) :: vol,tmp2(3,3)

    vol = det3(reallatt(1,1))
    if(vol < 0) then
      write(*,*) 'vol = ',vol
      write(*,'(A)') 'err: vol is negative. Check reallatt.'
      stop 1
    endif

    tmp2 = inv3x3(reallatt(1,1))

    reciplatt(1:3,1) = tmp2(1,1:3)
    reciplatt(1:3,2) = tmp2(2,1:3)
    reciplatt(1:3,3) = tmp2(3,1:3)
  end subroutine real2recip

  function crossprod(A,B)  result(C)
    real(double) :: A(1:3),B(1:3),C(1:3)
    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)
  end function crossprod

  function dotprod3(A,B)
    real(double) :: A(1:3),B(1:3),dotprod3
    dotprod3 = dotprodn(3,a(1),b(1))
  end function dotprod3

  function dotprodn(n,A,B)
    integer :: n,i
    real(double) :: A(n),B(n),dotprodn
    dotprodn = zero
    do i = 1, n
      dotprodn = dotprodn + A(i)*B(i)
    enddo
  end function dotprodn

  function tripprod(A)
    real(double) :: A(3,3),tripprod
    tripprod = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end function tripprod

  function det3(a)
    real(double) :: det3,a(3,3)
    det3 = tripprod(a)
  end function det3

  subroutine frac2abs(reallatt,FracC,AbsC)
    real(double) :: FracC(1:3),AbsC(1:3),reallatt(3,3)
    integer :: i,j
    do i = 1, 3
      AbsC(i) = 0
      do j = 1, 3
        AbsC(i) = AbsC(i) + reallatt(i,j)*FracC(j)
      enddo
    enddo
  end subroutine frac2abs

  subroutine abs2frac(reciplatt,AbsC,FracC)
    real(double) :: AbsC(3),FracC(3),reciplatt(3,3)
    integer :: i
    do i = 1, 3
      FracC(i) = dotprod3(AbsC(1:3),reciplatt(1:3,i))
    enddo
  end subroutine abs2frac

  function vecmag3(x)
    real(double) :: vecmag3,x(3)
    vecmag3 = sqrt(dotprod3(x,x))
  end function vecmag3

  function inv3x3(a) result(inv)
    real(double) :: a(3,3),inv(3,3)
    real(double) :: det,absdet

    inv(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
    inv(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
    inv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
    inv(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
    inv(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
    inv(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))
    inv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    inv(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
    inv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    det = a(1,1)*inv(1,1)+a(1,2)*inv(2,1)+a(1,3)*inv(3,1)
    absdet = abs(det)

    if(absdet < 1.0d-13) then
      write(*,*) 'absdet = ',absdet
      write(*,'(A)') 'err: det is too small in inv3x3.'
      stop 1
    endif
    inv = inv/det
  end function inv3x3

  subroutine transposen(n,A,B)
    integer :: i,j,n
    real(double) :: A(n,n),B(n,n)
    do i = 1, n
      do j = 1, n
        B(i,j) = A(j,i)
      enddo
    enddo
  end subroutine transposen

  subroutine matvecn(n,a,b,c)
    integer :: j,k,n
    real(double) :: a(n,n),b(n),c(n)
    do j = 1, n
      c(j) = 0D0
      do k = 1,n
        c(j) = c(j) + a(j,k)*b(k)
      enddo
    enddo
  end subroutine matvecn

  subroutine matvec3(a,b,c)
    real(double) :: a(3,3),b(3),c(3)
    call matvecn(3,a(1,1),b(1),c(1))
  end subroutine matvec3

  subroutine modulo1(x1)
    real(double) :: x1
    integer :: xint
    real(double) :: pluseps=1.0d-9

    if(x1 < 0.0d0) then

      xint = int(x1)
      x1 = x1 - xint + 1d0
      if(x1 == 1.0d0) x1 = 0d0
    else if(x1 >= 1.0d0) then

      xint = int(x1)
      x1 = x1 - xint
    endif

    if(x1 >= one .or. x1 < zero) then
      write(*,*) 'x1 is ',x1
      write(*,'(A)') 'err: must stop, x1 must be in [0,1) '
      stop 1
    endif

    if(abs(x1 - one) < pluseps) then

      x1 = zero
    endif
  end subroutine modulo1

  subroutine supercell_2_vasp(unitn,filename,atmtype,zord,s,writeforce,ATMSswitch)
    integer,optional :: writeforce
    integer,optional :: ATMSswitch
    integer :: unitn
    character(len=*) :: filename
    type(supercell) :: s
    integer :: atmtype,zord(atmtype)
    integer :: totat,i,ind,sizei,j,n
    integer :: groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer,allocatable :: fmap(:)
    character(len=DSL) :: commentline
    real(double),parameter :: maxlen=1000d0
    real(double) :: len6(6),frac(3),cartesian(3),num3(3)
    integer,allocatable :: localbmap(:)
    integer :: k

    n = s%n
    allocate(fmap(n))
    allocate(localbmap(n))

    ind = 0
    totat = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
          outcell(ind)%force(1:3) = s%at(j)%force(1:3)
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_vasp'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_vasp: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo
    open(unit=unitn,file=trim(filename),status='replace')

    if(present(ATMSswitch)) then
      if(ATMSswitch == 1) then
        commentline = 'ATMS: '//trim(N2str(n))//' '//trim(N2str(atmtype))
        do i = 1, atmtype

          commentline = trim(commentline)//' '//trim(N2str(groupsize(i)))//' '//trim(Cats(zord(i)))
        enddo
        commentline = trim(commentline)//' SG '//trim(s%sg_label)
      endif
    else
      commentline = ''

      if(trim(s%sg_label) /= '' .and. trim(s%sg_label) /= 'undefined') then
        commentline = 'SG '//trim(s%sg_label)
      endif
    endif

    write(unitn,'(A)') trim(commentline)

    write(unitn,'(F6.3)') 1.0d0
    write(unitn,'(3F18.12)') s%a(1:3,1)
    write(unitn,'(3F18.12)') s%a(1:3,2)
    write(unitn,'(3F18.12)') s%a(1:3,3)

    commentline=''
    do i = 1, atmtype
      commentline = trim(commentline)//' '//trim(CATS(zord(i)))
    enddo
    write(unitn,'(A)') trim(commentline)
    write(unitn,'(10I4)') (groupsize(i),i=1,atmtype)
    if(trim(s%fracabs) == 'frac') then
      write(unitn,'(A)') 'Direct'
    else if(trim(s%fracabs) == 'abs') then
      write(unitn,'(A)') 'Cartesian'
    else
      write(*,'(A)') 'err: Must be either frac or abs.'
      stop 1
    endif

    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
        frac(1:3) = outcell(ind)%p(1:3)

        call frac2abs(s%a,frac(1),cartesian(1))

        if(trim(s%fracabs) == 'frac') then
          num3(1:3) = frac(1:3)

          do k = 1, 3
            call modulo1(num3(k))
          enddo
        else if(trim(s%fracabs) == 'abs') then
          num3(1:3) = cartesian(1:3)
        endif
        if(.not. present(writeforce)) then
          write(unitn,'(3F18.12)') num3(1:3)
        else
          write(unitn,'(3F18.12,A,3F18.12)') num3(1:3),' | ',outcell(ind)%force(1:3)
        endif
      enddo
    enddo
    close(unitn)

    call getlenang(s%a,len6(1))
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_vasp

  subroutine onespace(ns,comment1,comment)
    integer :: ns
    character(len=ns) :: comment1,comment
    integer :: i,j

    comment = 's'

    i = 1
    j = 1

    loop2: do
      loop1: do
        if(comment1(i:i) /= ' ') exit loop1
        i = i + 1
        if(i == ns) exit loop2
      enddo loop1

      loop3: do
        if(comment1(i:i) == ' ') exit loop3
        comment(j:j) = comment1(i:i)
        j = j + 1
        i = i + 1
      enddo loop3
      comment(j:j) = ' '
      j = j + 1
    enddo loop2
  end subroutine onespace

  subroutine assign_str_fr_struct_out(inu,s)
    integer :: inu
    type(supercell) :: s
    integer :: i

    do i = 1, 3
      read(inu,*) s%a(1:3,i)
    enddo
    read(inu,*) s%n
    do i = 1, s%n
      read(inu,*) s%at(i)%gr, s%at(i)%z,s%at(i)%f(1:3)
    enddo
    call real2recip(s%a,s%b)
    close(inu)
  end subroutine assign_str_fr_struct_out

  subroutine assign_str_fr_xsf(controlu,sc)
    integer :: i,controlu
    type(supercell) :: sc
    character(LEN=DSL) :: tmp

    read(controlu,*) tmp
    if(trim(tmp) /= 'CRYSTAL') then
      write(*,'(A)') 'err: Keyword CRYSTAL not found.'
      stop 1
    endif
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMVEC') then
      write(*,'(A)') 'err: Keyword PRIMVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    read(controlu,*) tmp
    if(trim(tmp) /= 'CONVVEC') then
      write(*,'(A)') 'err: Keyword CONVVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    call real2recip(sc%a,sc%b)
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMCOORD') then
      write(*,'(A)') 'err: Keyword PRIMCOORD not found.'
      stop 1
    endif
    read(controlu,*) sc%n
    do i = 1, sc%n
      read(controlu,*) sc%at(i)%z, sc%at(i)%ac
      call abs2frac(sc%b,sc%at(i)%ac,sc%at(i)%f)
    enddo
    sc%sg_label = '1-P1'
  end subroutine assign_str_fr_xsf

  subroutine assign_str_fr_poscar(controlu,sc)
    integer :: controlu
    type(supercell) :: sc
    real(double) :: totalscaling,vol
    integer :: i,j,s,nsp,v,z,ind
    character(LEN=DSL) :: comment1,comment2,tmp1,tmp2,valuestr,targetlabel
    real(double) :: pos(3)
    integer,allocatable :: numarr(:)
    integer :: totn,k,labellen,SGnumber,loc
    logical :: done,found,puredigits
    integer :: nvalid
    real(double) :: density

    read(controlu,DSF) comment1

    call onespace(DSL,comment1,comment2)

    sc%commentline = trim(comment2)
    write(*,'(A,A,A)') 'POSCAR comment is [',trim(comment2),'].'
    s = scan(comment2,' ')

    if(comment2(1:s-1) == 'ATMS:') then
      tmp1= comment2(s+1:DSL)

      read(tmp1,*) v
      sc%n = v
      if(sc%n < 1) then
        write(*,*) 'sc%n = ',sc%n
        write(*,'(A)') 'err: check sc%n in assign_str_fr_poscar'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)

      read(tmp1,*) v
      sc%nsp = v
      if(sc%nsp < 1) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,'(A)') 'err: check sc%nsp in assign_str_fr_poscar.'
        stop 1
      elseif(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp,zordmax=',sc%nsp,zordmax
        write(*,'(A)') 'err: check sc%nsp, i.e., the number of species/elements'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)
      nsp = v
      if(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,*) 'zordmax = ',zordmax
        write(*,'(A)') 'err: out of bound.'
        stop 1
      endif

      do i = 1, nsp

        read(tmp1,*) v

        sc%n_per_species(i) = v

        s = scan(tmp1,' ')
        tmp1 = tmp1(s+1:DSL)

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        found = .false.
        do j = 1, zordmax
          if(found) exit
          if(trim(valuestr) == trim(Ats(j)) .or. trim(valuestr) == trim(CAts(j)) ) then
            z = j
            found = .true.
          endif
        enddo
        if(.not. found) then
          write(*,*) 'valuestr = ',trim(valuestr)
          write(*,'(A)') 'err: this valuestr is not found in the periodic table.'
          stop 1
        endif
        sc%zarr(i) = z
        tmp1 = tmp1(s+1:DSL)
      enddo

      ind = 0
      do j = 1, sc%nsp
        ind = ind + sc%n_per_species(j)
      enddo

      if(ind /= sc%n) then
        write(*,*) 'ind,sc%n=',ind,sc%n
        write(*,'(A)') 'err: missing atoms in the assign_str_fr_poscar.'
        stop 1
      endif

      sc%sg_label = ''
      if(trim(tmp1) /= '') then

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        if(trim(valuestr) /= 'SG') then
          write(*,*) 'valuestr is '//trim(valuestr)
          write(*,*) 'The first optional parameter must be SG.'
          write(*,'(A)') 'err: Check poscar.'
          stop 1
        endif
        tmp1 = tmp1(s+1:DSL)
        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        sc%sg_label = trim(valuestr)
      else
        write(*,*) 'No SG label'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      allocate(numarr(sc%nsp))

      read(controlu,DSF) tmp1
      write(*,*) 'tmp1 is .'//trim(tmp1)//'.'
      call onespace(DSL,tmp1,tmp2)
      s = scan(tmp2,'0123456789')

      if(s == 0) then

        write(*,'(A)') 'We are going to process the species line |'//trim(tmp2)//'|'

        nsp = num_of_strings(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        do i = 1, sc%nsp
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)
          write(*,*) 'valuestr is '//trim(valuestr)
          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              z = j
              if(z /= sc%zarr(i)) then
                write(*,*) 'z, sc%zarr(i) = ',z,sc%zarr(i)
                write(*,'(A)') 'err: Inconsistency.'
                stop 1
              endif
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
        enddo

        read(controlu,DSF) tmp2
        nsp = num_of_integers(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif

        read(tmp2,*) numarr(1:sc%nsp)
      else

        nsp = num_of_integers(tmp1)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        read(tmp1,*) numarr(1:sc%nsp)
      endif

      do i = 1, sc%nsp
        if(numarr(i) /= sc%n_per_species(i)) then
          write(*,*) 'species: i=',i
          write(*,*) 'numarr(i) = '//trim(N2str(numarr(i)))//', sc%n_per_species(i) = '//trim(N2str(sc%n_per_species(i)))
          write(*,'(A)') 'err: numarr(i) and sc%n_per_species(i) are not consistent'
          stop 1
        endif
      enddo
      deallocate(numarr)

    else
      if(comment2(1:s-1) == 'SG') then

        tmp2 = comment2(s+1:DSL)
        targetlabel = trim(tmp2)

        write(*,*) 'targetlabel = '//trim(targetlabel)

        puredigits = .true.
        labellen = len(trim(targetlabel))
        do i = 1, labellen
          if(.not. puredigits) then
            exit
          endif
         loc = scan(targetlabel(i:i),'0123456789')
         if(loc /= 1) then
           puredigits = .false.
         endif
        enddo

        if(puredigits)  then
          SGnumber = str2N(trim(targetlabel))
          write(*,*) 'SGnumber is ',SGnumber
          if(SGnumber >= 1 .and. SGnumber <= 230) then
            targetlabel = trim(SGbase(1,SGnumber)%sglabel)
            write(*,*) 'Fully expand the default SG label to the full SG label:'
            write(*,*) 'SGnumber= ',SGnumber
            write(*,*) 'Full targetlabel= ',trim(targetlabel)
          else
            write(*,*) 'err: must be between 1 and 230 inclusive. impossible. '
            stop 1
          endif
        endif
        sc%sg_label = trim(targetlabel)
        write(*,'(A)') 'Possibly adjusted sg_label is '//trim(sc%sg_label)
      else
        write(*,'(A)') 'WARNING WARNING: No ATMS:, no SG, hence we assume SG 1'
        sc%sg_label = '1-P1'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      read(controlu,DSF) tmp1
      call onespace(DSL,tmp1,tmp2)

      s = scan(tmp2,'0123456789')

      if(s == 0) then

        nsp = 0
        done = .false.
        do
          if(done) exit
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)

          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              nsp = nsp + 1
              sc%zarr(nsp) = j
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
          if(len(trim(tmp2)) == 0) then
            write(*,'(A)') 'nsp = '//trim(N2str(nsp))
            sc%nsp = nsp
            done = .true.
          endif
        enddo

        read(controlu,DSF) tmp1

        nvalid = num_of_integers(tmp1)
        if(nsp /= nvalid) then
          write(*,*) 'nsp,nvalid=',nsp,nvalid
          write(*,*) 'err: inconsistency in species line and species line.'
          stop 1
        endif
        read(tmp1,*) sc%n_per_species(1:nsp)

      else
        write(*,'(A)') 'err: Since this is without ATMs, we must have species line.'
        stop 1
      endif
    endif

    totn = 0
    do j = 1, sc%nsp
      v = sc%n_per_species(j)
      do k = 1, v
        totn = totn + 1
        sc%at(totn)%z = sc%zarr(j)
        sc%at(totn)%mass = massofa(sc%zarr(j))
      enddo
    enddo
    write(*,'(A)') 'totn = '//trim(N2str(totn))
    sc%n = totn

    sc%a = sc%a*totalscaling

    call real2recip(sc%a,sc%b)

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Real latt:', sc%a(1:3,1),'    | ',vecmag3(sc%a(1:3,1)),sc%a(1:3,2),'    | ',vecmag3(sc%a(1:3,2)),sc%a(1:3,3),'    | ',vecmag3(sc%a(1:3,3))
    vol = det3(sc%a(1,1))
    if(vol .le. 0) then
      write(*,'(A)') 'err: vol is not positive.'
      stop 1
    endif

    sc%vol = vol

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Reciprocal latt:', sc%b(1:3,1),'    | ',vecmag3(sc%b(1:3,1)),sc%b(1:3,2),'    | ',vecmag3(sc%b(1:3,2)),sc%b(1:3,3),'    | ',vecmag3(sc%b(1:3,3))

    read(controlu,*) tmp1
    if(trim(tmp1) == 'Selective') then
      read(controlu,*) tmp1
    endif
    if(trim(tmp1) /= 'Direct' .and. trim(tmp1) /= 'direct' .and. trim(tmp1) /= 'Cartesian' ) then
      write(*,'(A)') 'err: dummy should be Direct, direct, or Cartesian'
      stop 1
    endif

    do i = 1, sc%n
      sc%at(i)%force(1:3) = (/zero,zero,zero/)
    enddo
    do i = 1, sc%n
      if(trim(tmp1) == 'Direct' .or. trim(tmp1) == 'direct') then
        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) sc%at(i)%f(1:3)

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      else if(trim(tmp1) == 'Cartesian') then

        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) pos(1:3)

        call abs2frac(sc%b(1,1),pos(1),sc%at(i)%f(1))

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      endif
    enddo
    call getlenang(sc%a,sc%la(1))
    write(*,'(A)') 'la(1:6)='
    do i = 1, 6
      write(*,'(3F15.10)') sc%la(i)
    enddo
    write(*,*) 'Volume of the supercell= ',vol
    density = crys_density(sc)
    sc%density = density
    write(*,*) 'Density is ',density,' g/cm^3'

  end subroutine assign_str_fr_poscar

  subroutine assign_str_fr_fdf(inu,s)
    integer :: inu,natom,nspecies,sp,i,j
    type(supercell) :: s
    character(len=DSL) :: str1,str2,str3
    real(double) :: x(3)
    character(len=DSL) :: icf

    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfAtoms') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofAtoms'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') natom
      write(*,*) 'natom = ',natom
      s%n = natom
    endif
    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfSpecies') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') nspecies
      write(*,*) 'nspecies = ',nspecies
      s%nsp = nspecies
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, nspecies
      read(inu,*) j,sp
      if(j/=i) then
        write(*,*) 'j,i=',j,i
        write(*,*) 'sp = ',sp
        write(*,'(A)') 'err: something wrong with the species numbering ?'
        stop 1
      else
        s%zarr(i) = sp
        write(*,*) 's%zarr(',i,')=',sp
      endif
    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2,str3
    if(trim(str1) /= 'LatticeConstant') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be LatticeConstant'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str2) /= '1.0') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be 1.0'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str3) /= 'Ang') then
      write(*,*) 'str3 = ',trim(str3)
      write(*,*) 'but str3 must be Ang'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be  LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) s%a(1:3,1)
    read(inu,*) s%a(1:3,2)
    read(inu,*) s%a(1:3,3)
    write(*,*) 's%a is '
    write(*,*) s%a(1:3,1)
    write(*,*) s%a(1:3,2)
    write(*,*) s%a(1:3,3)
    call real2recip(s%a,s%b)
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) == 'Fractional') then
      icf = 'frac'
    else if(trim(str2) == 'NotScaledCartesianAng') then
      icf = 'abs'
    else
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be either Fractional or NotScaledCartesianAng'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, natom
      read(inu,*) x(1:3),sp
      if(trim(icf) == 'frac') then
        s%at(i)%f(1:3) = x(1:3)
      else if(trim(icf) == 'abs') then
        call abs2frac(s%b,x,s%at(i)%f)
      else
        write(*,'(A)') 'err: frac and abs problem.'
        stop 1
      endif
      s%at(i)%z = s%zarr(sp)

    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
  end subroutine assign_str_fr_fdf

  function massofa(z)
    integer :: z
    real(double) :: massofa
    if(z < 0 .or. z > zordmax) then
      write(*,*) 'z = ',z
      write(*,'(A)') 'massofa not defined yet.'
      write(*,*) 'err: massofa problem'
      stop 1
    endif
    massofa = MASSOFATOM(z)
  end function massofa

  function crys_density(s)
    type(supercell) :: s
    real(double) :: mass,crys_density
    integer :: i

    mass = 0.0d0
    do i = 1, s%n
      mass = mass + massofatom(s%at(i)%z)
    enddo
    crys_density = (mass*AMU/(det3(s%a(1:3,1:3))*1.d-30))*(1.0d3/1.0d6)
  end function crys_density

  subroutine dbl_sort(n,v,ind,kflag)
    integer,parameter :: maxrecur=50
    integer :: i,n,segment,tind,ind(n),eindex(maxrecur),bindex(maxrecur),bi,ei,kflag,rightindex,leftindex
    real(double) :: v(n),vref,t
    logical :: foundright,foundleft,cross,inseq

    if(kflag /= 1 .and. kflag /= -1) then
      write(*,*) 'kflag in dbl_sort = ',kflag
      write(*,'(A)') 'err: kflag should be 1 (ascending) or -1 (descending).'
      stop 1
    endif

    if(n < 1) then
      write(*,*) 'dbl_sort, n = ',n
      write(*,'(A)') 'err: wrong array size for sorting.'
      stop 1
    endif
    if(n == 1) then
      return
    endif

    inseq = .true.
    i = 2
    do
      if(i > n) then
        exit
      endif
      if(.not. inseq) then
        exit
      endif
      if(kflag == 1) then
        inseq = v(i) >= v(i-1)
      else if(kflag == -1) then
        inseq = v(i) <= v(i-1)
      endif
      i = i + 1
    enddo
    if(inseq) then
     return
    endif

    if(kflag == -1) then
      do i = 1, n
        v(i) = -v(i)
      enddo
    endif

    segment = 1
    bindex(1) = 1
    eindex(1) = n
    do

      if(segment == 0) exit

      bi = bindex(segment)
      ei = eindex(segment)
      vref = v(bi)
      rightindex = ei
      leftindex = bi+1
      cross = .false.
      do
        if(cross) then
          exit
        endif
        foundright = .false.
        do
          if(foundright) exit
          if(v(rightindex) >= vref) then
            rightindex = rightindex - 1
            if(rightindex < bi) then
              foundright = .true.
            endif
          else
            foundright = .true.
          endif
        enddo
        foundleft = .false.
        do
          if(foundleft) exit
          if(v(leftindex) <= vref) then
            leftindex = leftindex + 1
            if(leftindex > ei) then
              foundleft = .true.
            endif
          else
            foundleft = .true.
          endif
        enddo
        if(leftindex > rightindex) then
          cross = .true.

          if(rightindex < bi) then

            segment = segment - 1

            if(ei-bi > 1) then
              segment = segment + 1
              bindex(segment) = bi+1
              eindex(segment) = ei
            endif
          else

            if(v(rightindex) < vref) then
              v(bi) = v(rightindex)
              v(rightindex) = vref
              tind = ind(bi)
              ind(bi) = ind(rightindex)
              ind(rightindex) = tind
            endif

            segment = segment - 1

            if(rightindex - bi > 1) then
              segment = segment + 1
              bindex(segment) = bi
              eindex(segment) = rightindex - 1
            endif
            if(ei-rightindex > 1) then
              segment = segment + 1
              if(segment > maxrecur) then
                write(*,*) 'maxrecur = ',maxrecur
                write(*,'(A)') 'maxrecur is not large enough'
              endif
              bindex(segment) = rightindex + 1
              eindex(segment) = ei
            endif
          endif
        else
          t = v(leftindex)
          v(leftindex) = v(rightindex)
          v(rightindex) = t
          tind = ind(leftindex)
          ind(leftindex) = ind(rightindex)
          ind(rightindex) = tind
          leftindex = leftindex+1
          rightindex = rightindex-1
        endif
      enddo
    enddo
    do i = 2, n
      if(v(i) < v(i-1)) then
        write(*,*) 'v(i),v(i-1)=',v(i),v(i-1)
        write(*,'(A)') 'err: v(i) must be greater or equal to v(i-1).'
        stop 1
      endif
    enddo
    if(kflag == -1) then
      do i = 1, n
        v(i) = -v(i)
      enddo
    endif
  end subroutine dbl_sort

  subroutine read_struc(inputu,inputfile,filestem,inputformat,s)
    integer :: inputu,length,ind
    character(len=*) :: inputfile,filestem,inputformat
    type(supercell) :: s

    length = len_trim(inputfile)
    ind = index(trim(inputfile),'.',BACK=.TRUE.)
    filestem=inputfile(1:ind-1)

    inputformat=inputfile(ind+1:length)

    open(inputu,file=trim(inputfile),status='old',action='read')

    if(trim(inputformat)=='vasp' .or. trim(inputformat) == 'VASP') then
      call assign_str_fr_poscar(inputu,s)

    else if(trim(inputformat)=='STRUCT_OUT') then
      call assign_str_fr_struct_out(inputu,s)
    else if(trim(inputformat)=='fdf') then
      call assign_str_fr_fdf(inputu,s)
    else if(trim(inputformat)=='xsf') then
      call assign_str_fr_xsf(inputu,s)
    else
      write(*,*)
      write(*,*) 'WARNING: accepted format are vasp,arc/car,STRUCT_OUT,gjf,fdf'
      write(*,*) 'but inputformat is ',trim(inputformat)
      write(*,*) 'unrecognized input file format.'
      write(*,'(A)') 'err: check input file format.'
      stop 1
    endif
    close(inputu)

  end subroutine read_struc

  subroutine assign_m12_tvec(tvec,m12,ns1,s1,s2,distance_thres)
    real(double) :: tvec(3)
    integer :: ns1,m12(ns1)
    type(supercell) :: s1,s2
    real(double) ::distance_thres
    integer :: i,j,ii,jj,kk
    real(double) :: v1(3),v2(3),v3(3),newf(3)
    real(double) :: mindis,dis
    m12(1:ns1) = -1

    do i = 1, ns1
      call frac2abs(s1%a(1,1),s1%at(i)%f(1),v1(1))

      v1(1:3) = v1(1:3) + tvec(1:3)

      mindis = 1d100
      do j = 1, s2%n

        do kk = -1, 1
          do jj = -1, 1
            do ii = -1, 1
              newf = s2%at(j)%f + (/ii,jj,kk/)
              call frac2abs(s2%a(1,1),newf(1),v2(1))
              v3 = v2-v1
              dis = vecmag3(v3(1))
              if(dis < mindis) then
                mindis = dis
              endif
              if(dis < distance_thres) then
                m12(i) = j
                exit
              endif
            enddo
          enddo
        enddo
      enddo

      if(m12(i) == -1) then
        write(*,*) 'mindis = ',mindis
        write(*,*) 'atm index=', i, ', cannot be found in supercell.poscar.'
        write(*,'(A)') 'err: serious mapping error in assign_m12'
        stop 1
      endif
    enddo

  end subroutine assign_m12_tvec

  subroutine assign_m21(s2,s1,m21,ns2)
    type(supercell) :: s1,s2
    integer :: ns2,ns1
    integer :: i,j,m21(ns2)
    real(double) :: v2(3),v1(3),diffv(3),a1(3,3),inva1(3,3)
    real(double) :: sol(3),devia(3)

    m21(1:ns2) = -1

    a1(1:3,1:3) = s1%a(1:3,1:3)

    inva1 = inv3x3(a1)

    ns1 = s1%n

    do i = 1, ns2
      call frac2abs(s2%a,s2%at(i)%f,v2)
      do j = 1, ns1
        call frac2abs(s1%a,s1%at(j)%f,v1)

        diffv = v2-v1
        call matvecn(3,inva1,diffv,sol)

        devia(1:3) = abs(sol(1:3)-nint(sol(1:3)))

        if(vecmag3(devia(1)) < 1.d-10) then
          m21(i) = j

        endif
      enddo
      if(m21(i) == -1) then
        write(*,*) 'i = ',i
        write(*,*) 'cannot find an equivalent atom in the pPOSCAR'
        write(*,'(A)') 'err: serious translational error.'
        stop 1
      endif
    enddo
  end subroutine assign_m21

  subroutine least_sq_fit(ndeg,ndata,offset,xval,yval,s,M,MT,MTM,MTy,sol,ipiv)
    integer :: ndeg,ndata,ndeg1
    real(double) :: offset
    real(double) :: xval(ndata),yval(ndata)
    real(double) :: M(ndata,0:ndeg),MT(0:ndeg,ndata),MTM(0:ndeg,0:ndeg),MTy(0:ndeg)
    real(double) :: sol(0:ndeg,1)
    real(double) :: s(0:ndeg)
    integer :: ipiv(0:ndeg)
    integer :: i,j,info
    real(double) :: alpha,beta

    do i = 1, ndata
      do j = 0, ndeg
        if(j == 0) then
          M(i,j) = 1.0D0
        else
          M(i,j) = (xval(i)-offset)**(j*1.0D0)
        endif
      enddo
    enddo

    do i = 0, ndeg
      do j = 1, ndata
        MT(i,j) = M(j,i)
      enddo
    enddo
    ndeg1=ndeg+1
    alpha = One
    beta = Zero
    call dgemm('N','N',ndeg1,ndeg1,ndata,1.0D0,MT(0,1),ndeg1,M,ndata,0.0D0,MTM(0,0),ndeg1)
    alpha = One
    beta = Zero
    call dgemv('N',ndeg1,ndata,alpha,MT(0,1),ndeg1,yval(1),1,beta,MTy(0),1)
    sol(:,1) = MTy(:)

    CALL dgesv(ndeg1,1,MTM(0,0),ndeg1,ipiv(0),sol(0,1),ndeg1,info)
    if(info /= 0) then
      write(*,*) 'In least_sq_fit, info = ',info
      write(*,'(A)') 'err: info not zero.'
      stop 1
    endif
    s(:) = sol(:,1)
    do i = 0, ndeg

    enddo
  end subroutine least_sq_fit

  subroutine fx_fpx_least_sq_fit(x,a,ndeg,offset,fx,fpx)
    real(double) :: fx,fpx,x,a(0:ndeg),offset
    integer :: i,ndeg
    fx = zero
    fpx = zero
    do i = 0, ndeg
      fx = fx + a(i)*(x-offset)**(i*one)
      if(i > 0) then
        fpx = fpx + i*a(i)*(x-offset)**(i-one)
      endif
    enddo
  end subroutine fx_fpx_least_sq_fit

  subroutine find_ellipse_area2(zv,area,r1,r2)
    complex(double) :: zv(3)
    real(double) :: area,r1,r2
    real(double) :: area2,Rz(3),Iz(3),crossp(3)

    real(double) :: a1,b1,a2,b2,a3,b3
    real(double) :: A,B,C,Ap,Cp
    real(double) :: R,S
    real(double) :: konst

    a1 = real(zv(1))
    b1 = aimag(zv(1))

    a2 = real(zv(2))
    b2 = aimag(zv(2))
    a3 = real(zv(3))
    b3 = aimag(zv(3))

    A = a1*a1 + a2*a2 + a3*a3
    B = b1*b1 + b2*b2 + b3*b3
    C = -(a1*b1 + a2*b2 + a3*b3)

    Ap = (A-B)/2d0
    Cp = C

    S = (A+B)/2d0

    R = sqrt(Ap*Ap + Cp*Cp)

    r1 = sqrt(S+R)

    if(R > S) then

      if(abs(R-S) < 1d-15) then
        r2 = zero
      else
        write(*,*) 'too large a difference.'
        write(*,'(A)') 'err: Give up.'
        stop 1
      endif
    else
      r2 = sqrt(S-R)
    endif

    konst = pi

    area = konst*r1*r2

    Rz(1:3) = real(zv(1:3))
    Iz(1:3) = aimag(zv(1:3))

    crossp = crossprod(rz(1),iz(1))
    area2 = konst*vecmag3(crossp(1))

    if( abs(area2 - area) > 1d-7) then

      write(*,*) 'area2,area=',area2,area

      write(*,'(A)') 'err: Too big difference from two appraoches of obtaining areas.'
      stop 1
    else

    endif
  end subroutine find_ellipse_area2

  subroutine form_groups(nk,rsq,ig,ng,degentol)
    integer :: nk
    real(double) :: rsq(nk)
    integer :: ig(nk+1),ng,i,subdim
    real(double) :: degentol
    real(double) :: rg
    integer :: tsum

    rg = rsq(1)
    ng = 1
    ig(1) = 1
    do i = 2, nk
      if(abs(rsq(i)-rg) > degentol) then
        rg = rsq(i)
        ng = ng + 1
        ig(ng) = i
      endif
    enddo
    ig(ng+1) = nk+1

    tsum = 0

    do i = 1, ng
      subdim = ig(i+1)-ig(i)
      tsum = tsum + subdim

    enddo
    if(tsum /= nk) then
      write(*,*) 'tsum,nk=',tsum,nk
    endif
  end subroutine form_groups

  subroutine pt_shiftedfreq(aev,pev,zfreq,pfreq,nstages,minstage,degentol,nb,D1,D2,deltaD,tmpv1,tmpv2,tmpv3,T1,T2,w1,w2,indarr,CWORK,LCWORK,RWORK,LRWORK,Di1,Di2)

    integer :: nb,nstages,stage,LCWORK,LRWORK,minstage
    integer :: ng,m,k,dimen,mm,info,indarr(nb)
    integer :: dp,ind,ind2,kk,indw2,recordkk,tind

    integer,allocatable :: localmap1(:),localmap2(:),ig(:),gmap(:,:)
    integer,allocatable :: ipiv(:),store(:)

    integer,parameter :: ndata = 3
    integer,parameter :: ndeg=1

    complex(double) :: CWORK(LCWORK),D1(nb,nb),D2(nb,nb),deltaD(nb,nb)
    complex(double) :: tmpv1(nb,nb),tmpv2(nb,nb)
    complex(double) :: tmpv3(nb,nb)
    complex(double) :: T1(nb,nb),T2(nb,nb)
    complex(double) :: Di1(nb,nb),Di2(nb,nb)
    complex(double) :: Calpha,Cbeta

    real(double) :: pev(nb,0:nstages),aev(nb,0:nstages)
    real(double) :: offset,newv,mindis,dis,newvp,freqsq
    real(double) :: RWORK(LRWORK),zfreq(nb),pfreq(nb),w1(nb),w2(nb),freq,degentol

    real(double),allocatable :: tmpnu(:),deltawsqr(:),predictE(:),lambda(:),xval(:),yval(:)
    real(double),allocatable :: Mv(:,:),MT(:,:),MTM(:,:),MTy(:),sol(:,:),s(:)
    real(double),allocatable :: freqsqarr(:,:),sortarr(:),sev(:,:)
    integer :: debug

    allocate(xval(ndata))
    allocate(yval(ndata))
    allocate(MV(ndata,0:ndeg))
    allocate(MT(0:ndeg,ndata))
    allocate(MTM(0:ndeg,0:ndeg))
    allocate(MTy(0:ndeg))
    allocate(sol(0:ndeg,1))
    allocate(s(0:ndeg))
    allocate(ipiv(0:ndeg))

    allocate(store(nb))

    allocate(lambda(0:nstages))
    allocate(freqsqarr(nb,0:nstages))

    allocate(gmap(nb,0:nstages))
    allocate(ig(nb+1))
    allocate(deltawsqr(nb))
    allocate(predictE(nb))
    allocate(localmap1(nb))
    allocate(localmap2(nb))

    allocate(sortarr(nb))
    allocate(sev(nb,0:nstages))
    allocate(tmpnu(nb))

    if(nstages < 1) then
      write(*,*) 'nstages =',nstages
      stop 1
    endif

    do stage = 0, nstages
      lambda(stage) = stage*1.0D0/(nstages*1.0D0)
    enddo

    tmpv1(1:nb,1:nb) = D1(1:nb,1:nb)
    call zheev('V','U',nb,tmpv1(1,1),nb,w1(1),CWORK(1),LCWORK,RWORK(1),info)
    if(info /= 0) then
      write(*,*) 'info = ',info
      write(*,'(A)') 'err: zheev error.'
      stop 1
    endif
    freqsqarr(1:nb,0) = w1(1:nb)

    do k = 1, nb
      gmap(k,0) = k
    enddo

    aev(1:nb,0) = w1(1:nb)
    pev(1:nb,0) = w1(1:nb)

    do k = 1, nb
      freqsq = w1(k)
      freq = sqrt(abs(freqsq))
      if(freqsq < zero) then
        tmpnu(k) = -freq*Rydoverh2invcm
      else
        tmpnu(k) = freq*Rydoverh2invcm
      endif
    enddo

    debug = 0
    if(debug == 1) then
      write(*,'(A)') 'eigenvalues (in 1/cm) of zeroth structure, w1, are:'
      call print_arr_in_n_columns(tmpnu(1),nb,6)
    endif

    w2(1:nb) = w1(1:nb)
    tmpv2(1:nb,1:nb) = tmpv1(1:nb,1:nb)
    Di2(1:nb,1:nb) = D1(1:nb,1:nb)

    do k = 1, nb
      do m = 1, nb
        deltaD(m,k) = (D2(m,k) - D1(m,k))/(nstages*one)
      enddo
    enddo

    do stage = 1, nstages

      w1(1:nb) = w2(1:nb)
      tmpv1(1:nb,1:nb) = tmpv2(1:nb,1:nb)
      Di1(1:nb,1:nb) = Di2(1:nb,1:nb)

      Di2(1:nb,1:nb) = D1(1:nb,1:nb) + dcmplx(lambda(stage),0)*(D2(1:nb,1:nb) - D1(1:nb,1:nb))
      tmpv2(1:nb,1:nb) = Di2(1:nb,1:nb)

      call zheev('V','U',nb,tmpv2(1,1),nb,w2(1),CWORK(1),LCWORK,RWORK(1),info)
      if(info /= 0) then
        write(*,*) 'info = ',info
        write(*,'(A)') 'err: zheev error.'
        stop 1
      endif
      freqsqarr(1:nb,stage) = w2(1:nb)

      predictE(1:nb) = zero

      if(stage < minstage) then

        call form_groups(nb,w1(1),ig(1),ng,degentol)

        do k = 1, ng
          dimen = ig(k+1) - ig(k)

          do m = 1, dimen
            do mm = 1, nb
              tmpv3(mm,m) = tmpv1(mm,ig(k)+(m-1))
            enddo
          enddo

          T1(1:nb,1:dimen) = CZero
          Calpha = COne
          Cbeta = CZero
          call zgemm('N','N',nb,dimen,nb,Calpha,deltaD(1,1),nb,tmpv3(1,1),nb,Cbeta,T1(1,1),nb)
          T2(1:dimen,1:dimen) = CZero
          Calpha = COne
          Cbeta = CZero
          call zgemm('C','N',dimen,dimen,nb,Calpha,tmpv3(1,1),nb,T1(1,1),nb,Cbeta,T2(1,1),nb)

          call zheev('N','U',dimen,T2(1,1),nb,deltawsqr(1),CWORK(1),LCWORK,RWORK(1),info)
          if(info /= 0) then
            write(*,*) 'info = ',info
            write(*,'(A)') 'err: zheev error.'
            stop 1
          endif

          do m = 1, dimen
            predictE(ig(k)+(m-1)) = w1(ig(k)+(m-1)) + deltawsqr(m)
          enddo
        enddo

        do k = 1, nb
          indarr(k) = k
        enddo
        call dbl_sort(nb,predictE(1),indarr(1),1)

        call inversemapping(nb,indarr(1),localmap1(1))

      else

        do k = 1, nb
          store(k) = k
        enddo

        do k = 1, nb

          do dp = 1, ndata
            xval(dp) = lambda(stage-ndata+dp-1)
            yval(dp) = aev(k,stage-ndata+dp-1)
          enddo
          offset = lambda(stage)
          call least_sq_fit(ndeg,ndata,offset,xval(1),yval(1),s(0),Mv(1,0),MT(0,1),MTM(0,0),MTy(0),sol(0,1),ipiv(0))

          call fx_fpx_least_sq_fit(lambda(stage),s(0),ndeg,offset,newv,newvp)

          mindis = 1.0d100
          recordkk = -1
          do kk = nb, k, -1
            indw2 = store(kk)
            dis = abs(w2(indw2)-newv)
            if(dis < mindis) then

              mindis = dis
              recordkk = kk
            endif
          enddo

          if(recordkk == -1) then
            write(*,*) 'Bad recordkk= ',recordkk
            stop 1
          endif
          localmap2(k) = store(recordkk)

          tind = store(k)
          store(k) = store(recordkk)
          store(recordkk) = tind

        enddo

      endif

      do k = 1, nb
        ind = gmap(k,stage-1)
        if(stage < minstage) then
          ind2 = localmap1(ind)
          gmap(k,stage) = ind2
        elseif(stage >= minstage) then
          gmap(k,stage) = localmap2(k)
        endif
      enddo

      do k = 1, nb
        pev(k,stage) = predictE(gmap(k,stage))
        aev(k,stage) = w2(gmap(k,stage))
      enddo
    enddo

    do stage = nstages-ndata, 0, -1

      do k = 1, nb
        store(k) = k
      enddo

      do k = 1, nb
        do dp = 1, ndata
          xval(dp) = lambda(stage+dp)
          yval(dp) = aev(k,stage+dp)
        enddo
        offset = lambda(stage)
        call least_sq_fit(ndeg,ndata,offset,xval(1),yval(1),s(0),Mv(1,0),MT(0,1),MTM(0,0),MTy(0),sol(0,1),ipiv(0))
        call fx_fpx_least_sq_fit(lambda(stage),s(0),ndeg,offset,newv,newvp)
        mindis = 1.0d100

        do kk = nb, k, -1
          indw2 = store(kk)
          dis = abs(freqsqarr(indw2,stage)- newv)
          if(dis < mindis) then
            mindis = dis
            recordkk = kk
          endif
        enddo

        aev(k,stage) = freqsqarr(  store(recordkk), stage )

        tind = store(k)
        store(k) = store(recordkk)
        store(recordkk) = tind

      enddo
    enddo

    do k = 1, nb
      indarr(k) = k
      sortarr(k) = aev(k,0)
    enddo
    call dbl_sort(nb,sortarr(1),indarr(1),1)
    do k = 1, nb
      do stage = 0, nstages
        sev(k,stage) = aev(indarr(k),stage)
      enddo
    enddo

    do k = 1, nb
      do stage = 0, nstages
        aev(k,stage) = sev(k,stage)
      enddo
    enddo

    do k = 1, nb
      freqsq = aev(k,0)
      freq = sqrt(abs(freqsq))
      if(freqsq < zero) then
        zfreq(k) = -freq*Rydoverh2invcm
      else
        zfreq(k) = freq*Rydoverh2invcm
      endif
    enddo

    do k = 1, nb
      freqsq = aev(k,nstages)
      freq = sqrt(abs(freqsq))
      if(freqsq < zero) then
        pfreq(k) = -freq*Rydoverh2invcm
      else
        pfreq(k) = freq*Rydoverh2invcm
      endif
    enddo

    deallocate(xval)
    deallocate(yval)
    deallocate(Mv)
    deallocate(MT)
    deallocate(MTM)
    deallocate(MTy)
    deallocate(s)
    deallocate(sol)
    deallocate(ipiv)

    deallocate(store)

    deallocate(lambda)
    deallocate(freqsqarr)

    deallocate(gmap)
    deallocate(ig)
    deallocate(deltawsqr)
    deallocate(predictE)
    deallocate(localmap1)
    deallocate(localmap2)

    deallocate(sortarr)
    deallocate(sev)
    deallocate(tmpnu)

  end subroutine pt_shiftedfreq

  subroutine print_arr_in_n_columns(r,nb,n)
    integer :: nb,n
    real(double) :: r(nb)
    integer :: k,rem

    if(n < 1) then
      write(*,*) 'In print_arr_in_n_columns: n=',n
      write(*,*) 'n must be greater than 1'
      stop 1
    endif
    do k = 1, nb/n
      write(*,'(50E18.8)') r( (k-1)*n+1: n*k)
    enddo

    rem = modulo(nb,n)
    if(rem > 0) then
      write(*,'(50E18.8)') r( (nb/n)*n+1: (nb/n)*n+rem)
    endif
  end subroutine print_arr_in_n_columns
end module commod
