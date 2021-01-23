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

module qmod
  use commod
  implicit none

  character(len=*),parameter :: parfile                     ='qphonon.par'
  character(len=*),parameter :: qpathfile                   ='qpath.dat'
  character(len=*),parameter :: pposcar                     ='p.vasp'
  character(len=*),parameter :: sposcar                     ='supercell.vasp'
  character(len=*),parameter :: visposcar                   ='vis.vasp'
  character(len=*),parameter :: forcesdat                   ='forces.dat'
  character(len=*),parameter :: vasp_qpoints                ='vasp_line_mode_kpoints.dat'
  character(len=*),parameter :: freqgpf                     ='unsorted-freq-gives-GP.dat'

  character(len=*),parameter :: poscardir                   ='dir-poscars'

  character(len=*),parameter :: evecdir                     ='dir-eigenvectors'
  character(len=*),parameter :: ellipse_areaf               ='all-bs-ellipsearea.dat'
  character(len=*),parameter :: allqf                       ='allq.dyn'

  character(len=*),parameter :: qevalf                      ='qevals.dat'
  character(len=*),parameter :: qevecf                      ='qevec.dat'

  character(len=*),parameter :: fc_cutoff_file              ='fconst-cutoff-radius.dat'
  character(len=*),parameter :: calfconstsfile              ='cal-fconsts.dat'
  character(len=*),parameter :: indexf                      ='list-of-ids'
  character(len=*),parameter :: rawdf                       ='raw.dat'
  character(len=*),parameter :: gvbsf                       ='gv-bs.dat'
  character(len=*),parameter :: bcbsf                       ='bc-bs.dat'
  character(len=*),parameter :: pbcbsf                      ='pbc-bs.dat'

  integer,parameter :: inu           =11
  integer,parameter :: forceu        =12
  integer,parameter :: bsu           =13
  integer,parameter :: bsklistu      =14
  integer,parameter :: visu          =15
  integer,parameter :: qtravu        =16
  integer,parameter :: kinu          =17
  integer,parameter :: calfcu        =18
  integer,parameter :: zstaru        =19
  integer,parameter :: ofconstu      =20
  integer,parameter :: nfconstu      =21
  integer,parameter :: radu          =22
  integer,parameter :: deltaUu       =23
  integer,parameter :: evecu         =24
  integer,parameter :: evalqu        =25
  integer,parameter :: evecqu        =26
  integer,parameter :: vaspu         =27
  integer,parameter :: ellipse_areau =28
  integer,parameter :: allqu         =29
  integer,parameter :: posevecu      =30
  integer,parameter :: xcrysu        =31
  integer,parameter :: indexu        =32
  integer,parameter :: rawdu         =33
  integer,parameter :: freqgpu       =34
  integer,parameter :: gvbsu         =35
  integer,parameter :: bcbsu         =36
  integer,parameter :: pbcbsu        =37

  integer,parameter :: verbose = 1
  integer,parameter :: MODE_INDEX_LEN=4
  integer,parameter :: nrange=5
  integer,parameter :: ngmax = 200

  real(double),parameter :: distance_thres=1.0d-8

  type controlt
    character(len=DSL) :: title
    character(len=DSL) :: dir
    real(double) :: delta
    integer :: diffscheme
    integer :: ASR
    integer :: iradcut
    integer :: dynmat
    real(double) :: special_kpt(3)
    integer :: BS
    integer :: evec
    integer :: IndexFC(2)
    integer :: fcmulti
    integer :: storedyn
    integer :: frequnit

    integer :: gv
    real(double) :: gvdelta
    integer :: gvscheme

    integer :: GP
    integer :: GPscheme
    character(len=DSL) :: positivedir
    character(len=DSL) :: negativedir
    real(double) :: GPeps
    real(double) :: GPscalefactor

    integer :: dynnstages
    real(double) :: dyndegentol

    integer :: bc
    real(double) :: bcdelta

    integer :: bcsubn
    integer :: bcminstage
    integer :: bcndata
    integer :: bcndeg
    real(double) :: bcdegentol
  end type controlt

  type eigencheckert
    real(double) :: amp_evec
    real(double) :: amp_deltaU
    real(double) :: amp_vis
    integer :: n_deltaU
    integer :: n_vis
  end type eigencheckert

  type bst
    integer :: nsegments
    integer :: npoints
    real(double),allocatable :: bvec(:,:),evec(:,:)
    real(double),allocatable :: q(:,:,:)
    complex(double),allocatable :: freq(:,:,:)
    real(double),allocatable :: dis(:,:)
    real(double),allocatable :: r(:)
  end type bst

  type bstype
    integer :: ns
    integer :: np
    integer :: nb
    real(double),allocatable :: omega(:,:,:)
  end type bstype

  type gvtype
    integer :: ns
    integer :: np
    integer :: nb
    real(double),allocatable :: velocity(:,:,:,:)
  end type gvtype

  type lotot
    integer :: loto
    integer :: ns1
    integer :: ns2
    real(double) :: epsil(3,3)
    real(double),allocatable :: zeu(:,:,:)
    real(double),allocatable :: corfc(:,:,:,:)
  end type lotot

  type unitt
    real(double) :: length
    real(double) :: energy
    real(double) :: mass
    real(double) :: freqfac
    real(double) :: omega2nuRydinvh
  end type unitt

contains

  subroutine FiniteCrystalStat(s1,ns1,Infinite_NumGr_At,Infinite_MemNum_GrAt,Infinite_Dis_GrAt,s2,ns2,Finite_NumGr_At,Finite_MemNum_GrAt,Finite_Dis_GrAt,CompleteShellInd)
    type(supercell) :: s1,s2
    integer :: ns1,ns2,i,j
    integer :: Infinite_NumGr_At(ns1),Infinite_MemNum_GrAt(ngmax,ns1)
    real(double) :: dis,Infinite_Dis_GrAt(ngmax,ns1)
    real(double) :: Finite_Dis_GrAt(ns2,ns1)
    real(double) :: frac1(3),frac2(3),pos1(3),pos2(3)
    integer :: n1,n2,n3
    real(double),allocatable :: disarr(:)
    real(double) :: mindis,RG
    integer,allocatable :: indarr(:)
    integer :: Finite_NumGr_At(ns1)
    integer :: CompleteShellInd(ns1)
    integer :: Finite_MemNum_GrAt(ns2,ns1)
    integer :: NG,shellind
    integer,allocatable :: bindG(:)

    allocate(disarr(ns2))
    allocate(indarr(ns2))
    allocate(bindG(ns2+1))
    write(*,*) 'ns2 = ',ns2

    do i = 1, ns1
      frac1 = s1%at(i)%f
      call frac2abs(s1%a,frac1,pos1)
      do j = 1, ns2
        mindis = 1d100
        do n1 = -1, 1
          do n2 = -1, 1
            do n3 = -1, 1
              frac2 = (/n1*1d0,n2*1d0,n3*1d0/) + s2%at(j)%f
              call frac2abs(s2%a,frac2,pos2)
              dis = vecmag3(pos2-pos1)
              if(dis < mindis) then
                mindis = dis
              endif
            enddo
          enddo
        enddo
        disarr(j) = mindis
      enddo

      do j = 1, ns2
        indarr(j) = j
      enddo
      call dbl_sort(ns2,disarr(1),indarr(1),1)
      RG = disarr(1)
      NG = 1
      Finite_Dis_GrAt(NG,i) = RG
      bindG(1) = 1
      do j = 2, ns2
        if(abs(disarr(j) - RG) > 1d-8) then
          RG = disarr(j)
          NG = NG + 1
          Finite_Dis_GrAt(NG,i) = RG
          bindG(NG) = j
        endif
      enddo
      Finite_NumGr_At(i) = NG
      NG = NG + 1
      bindG(NG) = ns2+1
      do j = 1, Finite_NumGr_At(i)
        Finite_MemNum_GrAt(j,i) = bindG(j+1)-bindG(j)
      enddo
      write(*,'(A,I5,I5,I5,I5)') 'Dealing with statistic of atom,Finite_NumGr_At(i),Infinite_NumGr_At(i),ngmax=',i,Finite_NumGr_At(i),Infinite_NumGr_At(i),ngmax
      do j = 1, min(Finite_NumGr_At(i),Infinite_NumGr_At(i))

      enddo

      shellind = 0
      do
        shellind = shellind + 1
        if(shellind > Finite_NumGr_At(i)) then
          write(*,*) 'shellind,Finite_NumGr_At(i)= ',shellind,Finite_NumGr_At(i)
          exit
        endif
        if(Finite_MemNum_GrAt(shellind,i) == Infinite_MemNum_GrAt(shellind,i) .and. abs(Finite_Dis_GrAt(shellind,i)-Infinite_Dis_GrAt(shellind,i)) < 1d-8) then
          CompleteShellInd(i) = shellind
        else
          exit
        endif
      enddo
      write(*,*) 'atomi,Complete shell index is ',i,CompleteShellInd(i)
      write(*,*) ' '
    enddo
    deallocate(disarr)
    deallocate(indarr)
    deallocate(bindG)
  end subroutine FiniteCrystalStat

  subroutine InfiniteCrystalStat(s1,ns1,ng_ns1,numG,disG)
    type(supercell) :: s1
    real(double) :: frac1(3),frac2(3),pos1(3),pos2(3),distance
    integer :: ns1,i,j,n1,n2,n3,ind
    real(double),allocatable :: dis(:)
    integer,allocatable :: indarr(:),bindG(:)
    integer :: ng_ns1(ns1),numG(ngmax,ns1)
    integer :: newNG,NG,onedim,threedim,arrsize
    real(double) :: RG
    real(double) :: disG(ngmax,ns1)

    onedim = 2*nrange + 1
    threedim = onedim * onedim * onedim
    arrsize = ns1*threedim
    write(*,*) 'onedim,threedim,arrsize=',onedim,threedim,arrsize
    allocate(dis(arrsize))
    allocate(indarr(arrsize))
    allocate(bindG(ngmax))
    do i = 1, ns1
      frac1 = s1%at(i)%f
      call frac2abs(s1%a,frac1,pos1)
      ind = 0
      do n1 = -nrange,nrange
        do n2 = -nrange,nrange
          do n3 = -nrange,nrange
            do j = 1, ns1
              frac2 = s1%at(j)%f + (/n1*1d0,n2*1d0,n3*1d0/)
              call frac2abs(s1%a,frac2,pos2)
              distance = vecmag3(pos2-pos1)
              ind = ind+1
              dis(ind) = distance
            enddo
          enddo
        enddo
      enddo
      if(ind /= arrsize) then
        write(*,'(A)') 'err: array size problem'
        stop 1
      endif

      do  j = 1, arrsize
        indarr(j) = j
      enddo
      call dbl_sort(arrsize,dis(1),indarr(1),1)

      RG = dis(1)
      NG = 1
      disG(NG,i) = RG

      numG(NG,i) = 1
      do j = 2, arrsize
        if(abs(dis(j)-RG) > 1d-8) then
          RG = dis(j)
          newNG = NG + 1
          if(newNG > ngmax) then
            exit
          endif
          if(NG > ngmax) then

          endif
          NG = NG + 1
          disG(NG,i) = RG

          numG(NG,i) = 1
        else

          numG(NG,i) = numG(NG,i) + 1
        endif
      enddo
      ng_ns1(i) = NG
    enddo
    do i = 1, ns1
      do j = 1, ng_ns1(i)

      enddo
    enddo
    deallocate(dis)
    deallocate(indarr)
    deallocate(bindG)
  end subroutine InfiniteCrystalStat

  subroutine get_phys_unit(lengthstr,energystr,massstr,phyunit)
    type(unitt) :: phyunit
    character(len=DSL) :: lengthstr,energystr,massstr
    character(len=DSL) :: dimen,str

    write(*,*) 'lengthstr,energystr,massstr='//trim(lengthstr)//' '//trim(energystr)//' '//trim(massstr)
    call lowercase(lengthstr)

    select case(trim(lengthstr))
    case("ang")
      phyunit%length = Angstrom
    case("bohr")
      phyunit%length = Bohr
    case default
      dimen=lengthstr(1:2)
      write(*,*) 'dimen is '//trim(dimen)
      if(trim(dimen) == 'l=') then
        str = lengthstr(3:DSL)
        read(str,*) phyunit%length
      else
        write(*,*) 'supported unit of length is either Angstrom or Bohr or L=xxx'
        write(*,*) 'but the unit of length is ',trim(lengthstr)
        write(*,'(A)') 'err: unit of length '
        stop 1
      endif
    end select
    write(*,*) 'phyunit%length= ',phyunit%length

    call lowercase(energystr)
    select case(trim(energystr))
    case("ev")
      phyunit%energy = eV
    case("ryd")
      phyunit%energy = Rydberg
    case default
      dimen=energystr(1:2)
      write(*,*) 'dimen is '//trim(dimen)
      if(trim(dimen) == 'e=') then
        str = energystr(3:DSL)
        read(str,*) phyunit%energy
        write(*,*) phyunit%energy
      else
        write(*,*) 'supported unit of energy is either eV or Ryd or E=xxx'
        write(*,*) 'but the unit of energy is ',trim(energystr)
        write(*,'(A)') 'err: unit of energy'
        stop 1
      endif
    end select
    write(*,*) phyunit%energy

    call lowercase(massstr)
    select case(trim(massstr))
    case("amu")
      phyunit%mass = AMU
    case default
      dimen = massstr(1:2)
      if(trim(dimen) == 'm=') then
        str = massstr(3:DSL)
        read(str,*) phyunit%mass
      else
        write(*,*) 'supported unit of mass is AMU or M=xxx'
        write(*,*) 'but the unit of mass is ',trim(massstr)
        write(*,'(A)') 'err: unit of mass'
        stop 1
      endif
    end select
    write(*,*) phyunit%mass
  end subroutine get_phys_unit

  subroutine output_bs(bs,nb,dir)
    type(bst) :: bs
    integer :: i,j,k,kptind,nb,npoints,nsegments
    character(len=*) :: dir
    character(len=DSL) :: tbEIGENVALfile
    character(len=DSL) :: tbEIGENVALfile_klist
    integer,allocatable :: indarr(:)

    allocate(indarr(nb))

    npoints = bs%npoints
    nsegments = bs%nsegments

    write(*,*) 'output_bs: nsegments = ',nsegments, ', npoints =',npoints
    tbEIGENVALfile=trim(dir)//'/tbEIGENVAL'
    tbEIGENVALfile_klist=trim(dir)//'/tbEIGENVAL-klist-only'
    open(unit=bsu,file=trim(tbEIGENVALfile),status='replace')
    open(unit=bsklistu,file=trim(tbEIGENVALfile_klist),status='replace')
    write(bsu,'(A)') ' dummy1'
    write(bsu,'(A)') ' dummy2'
    write(bsu,'(A)') ' dummy3'
    write(bsu,'(A)') ' dummy4'
    write(bsu,'(A)') ' dummy5'
    write(bsu,'(1I4,2I10)') 0, npoints*nsegments, nb
    write(bsu,*) ''
    kptind = 0
    write(bsklistu,'(I10)') nsegments*npoints
    do i = 1, nsegments
      do j = 1, npoints
        kptind = kptind + 1
        write(bsu,'(3F20.10)') bs%q(1:3,j,i)
        write(bsklistu,'(3F20.10)') bs%q(1:3,j,i)

        do k = 1, nb
          bs%r(k) = real(bs%freq(k,j,i))
          indarr(k) = k
        enddo
        call dbl_sort(nb,bs%r(1),indarr(1),1)
        do k = 1, nb
          write(bsu,'(I4,2F20.10)') k, bs%freq(indarr(k),j,i)
        enddo
        write(bsu,*) ''
      enddo
    enddo
    close(bsu)
    close(bsklistu)
    write(*,'(A)') 'file = '//trim(tbEIGENVALfile)//' is closed'
    write(*,'(A)') 'file = '//trim(tbEIGENVALfile_klist)//' is closed'
    deallocate(indarr)
  end subroutine output_bs

  subroutine minv3(v1,v2,v3,s2)
    real(double) :: v(3),mindis,vdiff(3),v1(3),v2(3),v3(3),supercellv(3)
    type(supercell) :: s2
    integer :: i,j,k
    real(double) :: dis

    vdiff = v2 - v1
    mindis = 1d100
    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          supercellv = i*s2%a(:,1) + j*s2%a(:,2) + k*s2%a(:,3)
          v = vdiff + supercellv
          dis = vecmag3(v)
          if(dis < mindis) then
            mindis = dis
            v3 = v
          endif
        enddo
      enddo
    enddo
  end subroutine minv3

  subroutine minv3_multiplicity(v1,v2,v3,s2,storedv3,multiplicity)
    real(double) :: v(3),mindis,vdiff(3),v1(3),v2(3),v3(3),supercellv(3)
    type(supercell) :: s2
    integer :: i,j,k
    real(double) :: dis
    real(double) :: storedv3(3,50)
    integer :: multiplicity

    vdiff = v2 - v1
    mindis = 1.0d100
    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          supercellv = i*s2%a(:,1) + j*s2%a(:,2) + k*s2%a(:,3)
          v = vdiff + supercellv
          dis = vecmag3(v)
          if(abs(dis-mindis) < 1.0d-6) then
            multiplicity = multiplicity + 1
            v3 = v
            storedv3(1:3,multiplicity) = v3(1:3)
          else if(dis < mindis) then
            mindis = dis
            v3 = v
            multiplicity = 1
            storedv3(1:3,multiplicity) = v3(1:3)
          endif
        enddo
      enddo
    enddo

    do i = 1, multiplicity

    enddo
  end subroutine minv3_multiplicity

  subroutine STV(tn,v3,i,j,s1)
    real(double) :: v3(3),tn(3),d1(3),d2(3),d3(3)
    integer :: i,j
    type(supercell) :: s1
    real(double) :: A1(3,3),invA1(3,3),sol(3)
    integer :: ints1,ints2,ints3
    real(double) :: soltol

    call frac2abs(s1%a,s1%at(i)%f,d1)

    call frac2abs(s1%a,s1%at(j)%f,d2)

    d3 = d2-d1
    tn = v3-d3

    A1(1:3,1:3) = s1%a(1:3,1:3)

    invA1 = inv3x3(A1)
    call matvecn(3,invA1,tn,sol)

    if(verbose == 5) then

      ints1=nint(sol(1))
      ints2=nint(sol(2))
      ints3=nint(sol(3))

      soltol=1d-5
      if(abs(ints1-sol(1))>soltol .or. abs(ints2-sol(2)) > soltol .or. abs(ints3-sol(3))> soltol) then
        write(*,*) 'ints1,int2,int3= ',ints1,ints2,ints3
        write(*,*) 'sol(1:3) = ',sol(1:3)
        write(*,'(A)') 'err: integer ints1 problem.'
        stop 1
      endif
    endif
  end subroutine STV

  subroutine generate_perturb_struc(control,m12,ns1,s1,s2)
    type(controlt) :: control
    type(supercell) :: s1,s2
    integer :: ns1,m12(ns1)
    type(supercell) :: ps
    integer :: i,j,k,m,ns2,pm,indexins2
    real(double) :: delta,ac(3)
    character(len=MODE_INDEX_LEN) :: str
    integer :: run_ind,atmtype,zord(zordmax)
    character(len=DSL) :: absposcarfile

    delta = control%delta

    pm = 2
    ps%n = s2%n
    ps%a = s2%a
    ps%b = s2%b
    ns2 = s2%n
    do i = 1, ns2
      ps%at(i)%f = s2%at(i)%f
      ps%at(i)%z = s2%at(i)%z
    enddo

    run_ind = 0
    str(1:MODE_INDEX_LEN) = num2str(run_ind,MODE_INDEX_LEN)

    open(unit=indexu,file=trim(control%dir)//'/'//trim(poscardir)//'/'//trim(indexf),status='replace')
    write(indexu,'(A)') trim(str)

    atmtype = s1%nsp
    zord(1:atmtype) = s1%zarr(1:atmtype)
    ps%fracabs='frac'
    ps%sg_label=''

    absposcarfile = trim(control%dir)//'/'//trim(poscardir)//'/'//trim(str)//'.vasp'
    call supercell_2_vasp(vaspu,trim(absposcarfile),atmtype,zord(1),ps)

    run_ind = run_ind + 1
    do i = 1, ns1
      indexins2 = m12(i)
      do j = 1, 3
        do k = 1, pm

          do m = 1, ns2
            ps%at(m)%f = s2%at(m)%f
            ps%at(m)%z = s2%at(m)%z
            call frac2abs(ps%a,ps%at(m)%f,ps%at(m)%ac)
          enddo
          ac = ps%at(indexins2)%ac
          ac(j) = ac(j) + (-1)**(k-1)*delta
          ps%at(indexins2)%ac = ac
          call abs2frac(ps%b,ps%at(indexins2)%ac,ps%at(indexins2)%f)
          str(1:MODE_INDEX_LEN) = num2str(run_ind,MODE_INDEX_LEN)
          write(indexu,'(A)') trim(str)

          atmtype = s1%nsp
          zord(1:atmtype) = s1%zarr(1:atmtype)

          ps%fracabs='frac'
          absposcarfile = trim(control%dir)//'/'//trim(poscardir)//'/'//trim(str)//'.vasp'
          call supercell_2_vasp(vaspu,trim(absposcarfile),atmtype,zord(1),ps)
          run_ind = run_ind + 1
        enddo
      enddo
    enddo
    close(indexu)
  end subroutine generate_perturb_struc

  subroutine read_allforces(ns1,ns2,zf,nzf,dir)
    integer :: ns1,ns2
    real(double) :: zf(3,ns2)
    real(double) :: nzf(3,ns2,3,ns1,2)
    integer :: i,j,k,m,pm
    character(len=*) :: dir
    character(len=DSL) :: forcesfile

    real(double),allocatable :: sortf(:)
    integer,allocatable :: intarr(:)
    integer :: minprint

    forcesfile=trim(dir)//'/'//trim(forcesdat)
    write(*,'(A)') 'Try to open file ='//trim(forcesfile)
    open(unit=forceu,file=trim(forcesfile),action='read')

    allocate(sortf(ns2),intarr(ns2))

    do i = 1, ns2
      read(forceu,*) zf(1:3,i)
      intarr(i) = i
      sortf(i) = vecmag3(zf(1:3,i))
    enddo

    call dbl_sort(ns2,sortf(1),intarr(1),-1)
    minprint = 5
    if(ns2 < minprint) then
      minprint = ns2
    endif
    write(*,*) 'Residual force on atoms (sorted in descending order) :'
    do i = 1, minprint
      write(*,*) 'i,forcemag=',intarr(i),sortf(i)
    enddo

    pm = 2

    do i = 1, ns1
      do j = 1, 3
        do m = 1, pm
          do k = 1, ns2

            read(forceu,*) nzf(1:3,k,j,i,m)
          enddo
        enddo
      enddo
    enddo
    close(forceu)
    write(*,'(A)') 'file='//trim(forcesfile)//' is now closed'
    deallocate(sortf,intarr)
  end subroutine read_allforces

  subroutine cal_fconst(s1,s2,control,Finite_NumGr_At,Finite_Dis_GrAt,CompleteShellInd,m12,ns1,ns2,fconst,zf,nzf)
    type(supercell) :: s1,s2
    type(controlt) :: control
    integer :: ns1,ns2
    real(double) :: fconst(3,ns2,3,ns1),zf(3,ns2),nzf(3,ns2,3,ns1,2)
    integer :: Finite_NumGr_At(ns1)
    real(double) :: Finite_Dis_GrAt(ns2,ns1)
    integer :: CompleteShellInd(ns1)
    real(double) :: delta
    integer :: i,k,m,n,j
    integer :: m12(ns1)
    real(double) :: sumv,diagv,offdiagsum
    real(double),allocatable :: radius(:)
    real(double) :: v1(3),v2(3),v3(3),mindis
    integer :: ind,int1

    delta = control%delta

    Finite_NumGr_At(1) = Finite_NumGr_At(1)

    allocate(radius(ns1))
    write(*,'(A)') 'Print the values of cutoff radii in cal_fconst..'
    if(control%iradcut == 1) then
      radius(1:ns1) = 1d100
      write(*,'(A)') 'All interactions are used. radius(:) = 1d100 '
    else if(control%iradcut == 2) then
      do i = 1, ns1
        ind = CompleteShellInd(i)
        radius(i) = 0.5d0* ( Finite_Dis_GrAt(ind+1,i) + Finite_Dis_GrAt(ind,i) )
        write(*,*) 'atomi,radius=',i,radius(i)
      enddo
    else if(control%iradcut == 3) then
       open(unit=radu,file=trim(control%dir)//'/'//trim(fc_cutoff_file),action='read')
       read(radu,*) int1
       if(int1 /= ns1) then
         write(*,*) 'int1,ns1=',int1,ns1
         write(*,'(A)') 'err: int1 must be equal to ns1'
         stop 1
       else
         do i = 1, ns1
           read(radu,*) int1,radius(i)
           if(int1 /= i) then
             write(*,*) 'int1,i=',int1,i
             write(*,'(A)') 'err: int1 must be equal to i'
             stop 1
           endif
           write(*,*) 'atomi,radius=',i,radius(i)
         enddo
       endif
       close(radu)
    endif

    open(unit=ofconstu,file=trim(control%dir)//'/distance-give-old-fconst.dat',status='replace')
    write(ofconstu,'(A,I8,I3)') '# ', control%IndexFC(1), control%IndexFC(2)
    open(unit=nfconstu,file=trim(control%dir)//'/distance-give-new-fconst.dat',status='replace')
    write(nfconstu,'(A,I8,I3)') '# ', control%IndexFC(1), control%IndexFC(2)

    do i = 1, ns1
      do k = 1, 3
        do m = 1, ns2
          do n = 1, 3
              if(control%diffscheme == 1) then
                fconst(n,m,k,i) = (nzf(n,m,k,i,1)-nzf(n,m,k,i,2))/(2d0*delta)
              else if(control%diffscheme == 2) then
                fconst(n,m,k,i) = (nzf(n,m,k,i,1) - zf(n,m))/delta
              else if(control%diffscheme == 3) then
                fconst(n,m,k,i) = (zf(n,m) - nzf(n,m,k,i,2))/delta
              endif
              fconst(n,m,k,i) = -fconst(n,m,k,i)

            call frac2abs(s1%a,s1%at(i)%f,v1)
            call frac2abs(s2%a,s2%at(m)%f,v2)
            call minv3(v1,v2,v3,s2)
            mindis = vecmag3(v3(1))
            if(i == control%IndexFC(1) .and. k == control%IndexFC(2)) then
              write(ofconstu,*) mindis,fconst(n,m,k,i)
            endif
            if(mindis > radius(i)) then

              fconst(n,m,k,i) = zero
            endif
            if(i == control%IndexFC(1) .and. k == control%IndexFC(2)) then
              write(nfconstu,*) mindis,fconst(n,m,k,i)
            endif
          enddo
        enddo
      enddo
    enddo
    close(unit=ofconstu)
    close(unit=nfconstu)

    if(control%ASR == 0) then
      write(*,'(A)') 'No Newton''s third law: do nothing to the force constants'
    elseif(control%ASR == 1) then
      write(*,'(A)') 'Imposing the diagonal-element-wise Newton''s third law on the force constants'
      write(*,'(A)') 'You may want to set the switch to 2 to double check'
      do i = 1, ns1
        j = m12(i)
        do k = 1, 3

          sumv = 0
          do m = 1, ns2
            do n = 1, 3
              sumv = sumv + fconst(n,m,k,i)
            enddo
          enddo

          diagv = fconst(k,j,k,i)
          offdiagsum = sumv - diagv

          fconst(k,j,k,i) = -offdiagsum
        enddo
      enddo
    elseif(control%ASR == 2) then
      write(*,'(A)') 'Imposing the (more stringent/accurate) xyz-component-wise Newton''s third law on each force constant row'
      do i = 1, ns1
        j = m12(i)
        do k = 1, 3
          do n = 1, 3
            sumv = 0
            do m = 1, ns2
              sumv = sumv + fconst(n,m,k,i)
            enddo

            diagv = fconst(n,j,k,i)
            offdiagsum = sumv - diagv

            fconst(n,j,k,i) = -offdiagsum
          enddo
        enddo
      enddo
    endif
    deallocate(radius)
  end subroutine cal_fconst

  subroutine read_epsil_and_Born_effective_charges(zstarfile,loto)
    integer :: i,j,ns1
    type(lotot) :: loto
    character(len=*) :: zstarfile

    open(unit=zstaru,file=trim(zstarfile),action='read')
    read(zstaru,*)
    do i = 1, 3
      read(zstaru,*) loto%epsil(i,1:3)
    enddo
    read(zstaru,*)
    ns1 = loto%ns1
    do i = 1, ns1
      do j = 1, 3
        read(zstaru,*) loto%zeu(j,1:3,i)
        write(*,*) 'loto%zeu(j,1:3,i) = ',loto%zeu(j,1:3,i)
      enddo
    enddo
    close(unit=zstaru)
  end subroutine read_epsil_and_Born_effective_charges

  subroutine cal_dynmat(control,phyunit,fracqvec,fconst,loto,mass,dym,nb,s1,s2,m21,ns1,ns2)
    type(controlt) :: control
    type(unitt) :: phyunit
    integer :: ii,i,j,k,m,nb,ns1,ns2
    integer :: m21(ns2)
    type(supercell) :: s1,s2
    complex(double) :: dym(nb,nb),cphase
    real(double) :: v1(3),v2(3),v3(3),tn(3)
    real(double) :: mass(ns1)
    real(double) :: fconst(3,ns2,3,ns1)
    type(lotot) :: loto
    real(double) :: fracqvec(3),qvec(3),qdir(3)
    real(double) :: nsc_real
    integer :: nsc
    real(double) :: epsil(3,3),qeq,zag(3),zbg(3)

    integer :: n_a,n_b,a_blk,b_blk
    real(double) :: volume,correction,eprimesq
    real(double) :: storedv3(3,50)
    integer :: multiplicity,mm

    call real2recip(s1%a(1,1), s1%b(1,1))
    qvec(1:3) = two*pi*  ( fracqvec(1)*s1%b(1:3,1) + fracqvec(2)*s1%b(1:3,2) + fracqvec(3)*s1%b(1:3,3))

    do i = 1, nb
      do j = 1, nb
        dym(j,i) = zero
      enddo
    enddo

    do i = 1, ns1
      do j = 1, 3
        do k = 1, ns2
          do m = 1, 3
            loto%corfc(m,k,j,i) = fconst(m,k,j,i)
          enddo
        enddo
      enddo
    enddo

    if(loto%loto == 1) then
      volume  = det3(s1%a)

      eprimesq = (EChg**2/(4d0*pi*eps0))*(1.0d0/phyunit%Energy)*(1.0d0/phyunit%Length)
      nsc_real = ns2/(1d0*ns1)
      nsc = int( nsc_real )
      if(abs(nsc_real - nsc) > 1.0d-8) then
        write(*,*) 'ns1,ns2=',ns2,ns1
        write(*,*) 'nsc_real,nsc=',nsc_real,nsc
        write(*,'(A)') 'err: ns2 must be an integer multiple of ns1'
        stop 1
      endif
      qdir(1:3) = qvec(1:3)
      epsil(1:3,1:3) = loto%epsil(1:3,1:3)
      qeq = (qdir(1)*(epsil(1,1)*qdir(1)+epsil(1,2)*qdir(2)+epsil(1,3)*qdir(3))+    &
            qdir(2)*(epsil(2,1)*qdir(1)+epsil(2,2)*qdir(2)+epsil(2,3)*qdir(3))+    &
            qdir(3)*(epsil(3,1)*qdir(1)+epsil(3,2)*qdir(2)+epsil(3,3)*qdir(3)))
      if(qeq < 1.0d-12) then
        write(*,*) 'qeq has a very small magnitude'
        write(*,*) 'A direction for q was not specified:'
        write(*,*) 'TO-LO splitting will be absent.'
        write(*,*) 'qvec = ',qvec
      else
        do n_a = 1, ns1
          a_blk = n_a
          do n_b = 1, ns2
            b_blk = m21(n_b)
            if(b_blk < 1 .or. b_blk > s1%n) then
              write(*,*) 'b_blk = ',b_blk
              write(*,'(A)') 'err: invalid b_blk in cal_dynmat'
              stop 1
            endif
            do i = 1, 3
              zag(i) = qdir(1)*loto%zeu(1,i,a_blk) +  qdir(2)*loto%zeu(2,i,a_blk) +  qdir(3)*loto%zeu(3,i,a_blk)
              zbg(i) = qdir(1)*loto%zeu(1,i,b_blk) +  qdir(2)*loto%zeu(2,i,b_blk) +  qdir(3)*loto%zeu(3,i,b_blk)
            enddo
            do i = 1, 3
              do j = 1, 3
                correction =  4.0d0*pi*eprimesq*zag(i)*zbg(j)/qeq/volume
                loto%corfc(j,n_b,i,n_a) = loto%corfc(j,n_b,i,n_a) + correction/(nsc*1d0)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    do i = 1, ns1
      do j = 1, 3
        do k = 1, ns2

          call frac2abs(s1%a,s1%at(i)%f,v1)
          call frac2abs(s2%a,s2%at(k)%f,v2)

          call minv3(v1,v2,v3,s2)

          call minv3_multiplicity(v1,v2,v3,s2,storedv3(1,1),multiplicity)

          ii = m21(k)
          if(ii < 1 .or. ii > s1%n) then
            write(*,*) 'ii = ',ii
            write(*,'(A)') 'err: invalid ii in cal_dynmat.'
            stop 1
          endif

          if(control%fcmulti == 0) then
            call STV(tn,v3,i,ii,s1)
            do m = 1, 3
              cphase=exp(-Imag*dotprod3(tn,qvec))

              dym((i-1)*3+j,(ii-1)*3+m)=dym((i-1)*3+j,(ii-1)*3+m) + (1d0/sqrt(mass(i)*mass(ii)))*loto%corfc(m,k,j,i)*cphase
            enddo
          else
            do mm = 1, multiplicity
              v3 = storedv3(1:3,mm)
              call STV(tn,v3,i,ii,s1)
              do m = 1, 3
                cphase=exp(-Imag*dotprod3(tn,qvec))

                dym((i-1)*3+j,(ii-1)*3+m)=dym((i-1)*3+j,(ii-1)*3+m) + (1d0/sqrt(mass(i)*mass(ii)))*loto%corfc(m,k,j,i)*cphase/(one*multiplicity)
              enddo
            enddo
           endif
        enddo
      enddo
    enddo
  end subroutine cal_dynmat

  subroutine readrawd(dunit,absrawdf,control,phyunit,ns1,ns2,m21,mass,fconst,loto,s1,s2,bs)
    integer :: dunit
    character(len=*) :: absrawdf
    type(controlt) :: control
    type(unitt) :: phyunit
    integer :: ns1,ns2
    integer :: m21(ns2)
    real(double) :: mass(ns1)
    real(double) :: fconst(3,ns2,3,ns1)
    type(lotot) :: loto
    type(supercell) :: s1,s2
    type(bst) :: bs
    integer :: i,j,k,m
    character(len=DSL) :: tmpstr
    integer :: ns1_check,ns2_check

    open(unit=dunit,file=trim(absrawdf),action='read')

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'control') then
      write(*,*) 'err: read control'
      stop 1
    endif
    read(dunit,*) control%dynmat

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'phyunit') then
      write(*,*) 'err: read phyunit'
      stop 1
    endif
    read(dunit,*) phyunit%freqfac

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'ns1') then
      write(*,*) 'err: read ns1'
      stop 1
    endif
    read(dunit,*) ns1_check
    if(ns1_check /= ns1) then
      write(*,*) 'err: ns1_check,ns1= ',ns1_check,ns1
      stop 1
    endif

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'ns2') then
      write(*,*) 'err: read ns2'
      stop 1
    endif
    read(dunit,*) ns2_check
    if(ns2_check /= ns2) then
      write(*,*) 'err: ns2_check,ns2= ',ns2_check,ns2
      stop 1
    endif

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'm21') then
      write(*,*) 'err: read m21'
      stop 1
    endif
    do i = 1, ns2
      read(dunit,*) m21(i)
    enddo

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'mass') then
      write(*,*) 'err: read mass'
      stop 1
    endif
    do i = 1, ns1
      read(dunit,*) mass(i)
    enddo

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'fconst') then
      write(*,*) 'err: read fconst'
      stop 1
    endif
    do i = 1, ns1
      do j = 1, 3
        do k = 1, ns2
          do m = 1, 3
            read(dunit,*) fconst(m,k,j,i)
          enddo
        enddo
      enddo
    enddo

    read(dunit,'(A)') tmpstr
    if(trim(tmpstr) /= 'loto') then
      write(*,*) 'err: read loto'
      stop 1
    endif
    read(dunit,*) loto%loto
    if(loto%loto == 1) then
      read(dunit,*) loto%ns1
      read(dunit,*) loto%ns2
      do i = 1, 3
        do j = 1, 3
          read(dunit,*) loto%epsil(j,i)
        enddo
      enddo
      do i = 1, ns1
        do j = 1, 3
          read(dunit,*) loto%zeu(j,1:3,i)
        enddo
      enddo
    endif

    call rawreadsupercell(dunit,s1)
    call rawreadsupercell(dunit,s2)
    call rawreadbs(dunit,bs)
    close(dunit)

  end subroutine readrawd

  subroutine writerawd(dunit,absrawdf,control,phyunit,ns1,ns2,m21,mass,fconst,loto,s1,s2,bs)
    integer :: dunit
    character(len=*) :: absrawdf
    type(controlt) :: control
    type(unitt) :: phyunit
    integer :: ns1,ns2
    integer :: m21(ns2)
    real(double) :: mass(ns1)
    real(double) :: fconst(3,ns2,3,ns1)
    type(lotot) :: loto
    type(supercell) :: s1,s2
    type(bst) :: bs
    integer :: i,j,k,m

    open(unit=dunit,file=trim(absrawdf),status='replace')
    write(dunit,'(A)') 'control'
    write(dunit,*) control%dynmat
    write(dunit,'(A)') 'phyunit'
    write(dunit,*) phyunit%freqfac
    write(dunit,'(A)') 'ns1'
    write(dunit,*) ns1
    write(dunit,'(A)') 'ns2'
    write(dunit,*) ns2
    write(dunit,'(A)') 'm21'
    do i = 1, ns2
      write(dunit,*) m21(i)
    enddo
    write(dunit,'(A)') 'mass'
    do i = 1, ns1
      write(dunit,*) mass(i)
    enddo
    write(dunit,'(A)') 'fconst'
    do i = 1, ns1
      do j = 1, 3
        do k = 1, ns2
          do m = 1, 3
            write(dunit,*) fconst(m,k,j,i)
          enddo
        enddo
      enddo
    enddo
    write(dunit,'(A)') 'loto'
    write(dunit,*) loto%loto
    if(loto%loto == 1) then
      write(dunit,*) loto%ns1
      write(dunit,*) loto%ns2
      do i = 1, 3
        do j = 1, 3
          write(dunit,*) loto%epsil(j,i)
        enddo
      enddo
      do i = 1, ns1
        do j = 1, 3
          write(dunit,*) loto%zeu(j,1:3,i)
        enddo
      enddo
    endif
    call rawwritesupercell(dunit,s1)
    call rawwritesupercell(dunit,s2)
    call rawwritebs(dunit,bs)

    close(dunit)
  end subroutine writerawd

  subroutine rawreadbs(dunit,bs)
    integer :: dunit
    type(bst) :: bs
    integer :: i
    character(len=DSL) :: tmpstr

    read(dunit,*) tmpstr
    if(trim(tmpstr) /= 'bs-bvec-evec') then
      write(*,*) 'err: bs-bvec-evec'
      stop 1
    endif
    read(dunit,*) bs%nsegments
    read(dunit,*) bs%npoints
    do i = 1, bs%nsegments
      read(dunit,*) bs%bvec(1:3,i)
    enddo
    do i = 1, bs%nsegments
      read(dunit,*) bs%evec(1:3,i)
    enddo

  end subroutine rawreadbs

  subroutine rawwritebs(dunit,bs)
    integer :: dunit
    type(bst) :: bs
    integer :: i

    write(dunit,'(A)') 'bs-bvec-evec'
    write(dunit,*) bs%nsegments
    write(dunit,*) bs%npoints
    do i = 1, bs%nsegments
      write(dunit,*) bs%bvec(1:3,i)
    enddo
    do i = 1, bs%nsegments
      write(dunit,*) bs%evec(1:3,i)
    enddo

  end subroutine rawwritebs

  subroutine rawreadsupercell(dunit,s)
    integer :: dunit
    type(supercell) :: s
    integer :: i
    character(len=DSL) :: tmpstr

    read(dunit,*) tmpstr
    if(trim(tmpstr) /= 'supercell') then
      write(*,*) 'err: supercell'
      stop 1
    endif
    read(dunit,*) s%la(1:3)
    read(dunit,*) s%la(4:6)
    read(dunit,*) s%a(1:3,1)
    read(dunit,*) s%a(1:3,2)
    read(dunit,*) s%a(1:3,3)
    read(dunit,*) s%b(1:3,1)
    read(dunit,*) s%b(1:3,2)
    read(dunit,*) s%b(1:3,3)
    read(dunit,*) s%n
    read(dunit,*) s%nsp
    do i = 1, s%nsp
      read(dunit,*) s%n_per_species(i)
    enddo
    do i = 1, s%nsp
      read(dunit,*) s%zarr(i)
    enddo

    do i = 1, s%n
      read(dunit,*) s%at(i)%z
      read(dunit,*) s%at(i)%mass
      read(dunit,*) s%at(i)%f(1:3)
      read(dunit,*) s%at(i)%ac(1:3)
    enddo
  end subroutine rawreadsupercell

  subroutine rawwritesupercell(dunit,s)
    integer :: dunit
    type(supercell) :: s
    integer :: i

    write(dunit,'(A)') 'supercell'
    write(dunit,*) s%la(1:3)
    write(dunit,*) s%la(4:6)
    write(dunit,*) s%a(1:3,1)
    write(dunit,*) s%a(1:3,2)
    write(dunit,*) s%a(1:3,3)
    write(dunit,*) s%b(1:3,1)
    write(dunit,*) s%b(1:3,2)
    write(dunit,*) s%b(1:3,3)
    write(dunit,*) s%n
    write(dunit,*) s%nsp
    do i = 1, s%nsp
      write(dunit,*) s%n_per_species(i)
    enddo
    do i = 1, s%nsp
      write(dunit,*) s%zarr(i)
    enddo

    do i = 1, s%n
      write(dunit,*) s%at(i)%z
      write(dunit,*) s%at(i)%mass
      write(dunit,*) s%at(i)%f(1:3)
      write(dunit,*) s%at(i)%ac(1:3)
    enddo
  end subroutine rawwritesupercell

  subroutine bandconnect(control,phyunit,m21,mass,fconst,loto,ns1,ns2,nb,s1,s2,bs)
    type(controlt) :: control
    type(unitt) :: phyunit
    integer :: ns1,ns2,nb
    type(bst) :: bs
    real(double) :: bq(3),eq(3)
    type(supercell) :: s1,s2
    real(double) :: fconst(3,ns2,3,ns1)
    integer :: m21(ns2)
    real(double) :: mass(ns1)
    type(lotot) :: loto
    character(len=DSL) :: absbcbsf
    character(len=DSL) :: abspbcbsf
    integer :: LCWORK,LRWORK
    integer :: ind,kk,mm,i,k,nsegments
    integer :: intn,minn,maxn,minstage,ind2,tind,ndeg,ndata
    integer :: tmpn,stage,nstages,INFO,debug
    integer :: ng,m,dimen,dp,indw2,recordkk,ii,jj

    real(double) :: distold,dist,disti,freq,degentol
    real(double) :: qvec(3),Cb(3),Ce(3),Deltav(3),distance,realn
    real(double) :: offset,del(3),newv,mindis,dis,newvp,freqsq

    complex(double) :: Calpha,Cbeta,cmplv

    integer,allocatable :: indarr(:),store(:),segn(:),ipiv(:)
    integer,allocatable :: localmap1(:),localmap2(:),ig(:),gmap(:,:)

    real(double),allocatable :: RWORK(:)
    real(double),allocatable :: pev(:,:)
    real(double),allocatable :: aev(:,:)
    real(double),allocatable :: aevinvcm(:,:)
    real(double),allocatable :: sev(:,:)
    real(double),allocatable :: freqarr(:)
    real(double),allocatable :: MV(:,:),MT(:,:),MTM(:,:),MTy(:),sol(:,:),s(:),xval(:)
    real(double),allocatable :: lambda(:)
    real(double),allocatable :: freqsqarr(:,:),sortarr(:)
    real(double),allocatable :: tmpnu(:),deltawsqr(:),predictE(:),yval(:)
    real(double),allocatable :: w1(:),w2(:),w3(:)

    complex(double),allocatable :: CWORK(:)
    complex(double),allocatable :: Dm1(:,:),Dm2(:,:),deltaD(:,:)
    complex(double),allocatable :: tmpv1(:,:),tmpv2(:,:),tmpv3(:,:)
    complex(double),allocatable :: T1(:,:)
    complex(double),allocatable :: T2(:,:)

    integer :: rk,prefer_choice
    integer,allocatable :: prefer(:)

    write(*,*) 'bandconnect: nb= ',nb

    if(control%dynmat == 3) then
      write(*,*) 'dynmat= ',control%dynmat
      write(*,*) 'err: We only handle dynmat=1 and dynmat=2 only.'
      stop 1
    endif

    minstage = control%bcminstage
    if(control%bcminstage < 1) then
      write(*,*) 'err: bcminstage = ',control%bcminstage
      stop 1
    endif
    write(*,*) 'minstage= ',minstage

    if(control%bcsubn < 1) then
      write(*,*) 'err: bcsubn = ',control%bcsubn
      stop 1
    endif
    write(*,*) 'bcsubn = ',control%bcsubn

    LRWORK = nb*3
    LCWORK = nb*2
    allocate(RWORK(LRWORK))
    allocate(CWORK(LCWORK))

    nsegments = bs%nsegments
    allocate(segn(nsegments))

    degentol = control%bcdegentol
    write(*,*) 'degentol = ',degentol

    ndeg = control%bcndeg
    write(*,*) 'ndeg = ',ndeg

    minn = 10**5
    maxn = -10**5
    do i = 1, nsegments
      bq(1:3) = bs%bvec(1:3,i)
      eq(1:3) = bs%evec(1:3,i)

      call matvec3(s1%b(1,1),bq(1),Cb(1))
      call matvec3(s1%b(1,1),eq(1),Ce(1))
      Cb(1:3) = twopi*Cb(1:3)
      Ce(1:3) = twopi*Ce(1:3)
      write(*,'(A,3F10.5)') 'Cb(1:3)= ',Cb(1:3)
      write(*,'(A,3F10.5)') 'Ce(1:3)= ',Ce(1:3)

      Deltav(1:3) = Ce(1:3)-Cb(1:3)
      distance = vecmag3(deltav(1))
      realn = distance/control%bcdelta
      tmpn = int(realn+1)

      intn = tmpn*control%bcsubn
      segn(i) = intn
      write(*,*) 'i,realn,tmpn,intn= ',i,realn,tmpn,intn
      write(*,*) '|deltaq|= ',distance/(intn*one)
      if(intn > maxn) then
        maxn = intn
      endif
      if(intn < minn) then
        minn = intn
      endif
    enddo
    write(*,*) 'minn,maxn= ',minn,maxn

    allocate(aev(nb,0:maxn))
    allocate(pev(nb,0:maxn))
    allocate(aevinvcm(nb,0:maxn))
    allocate(gmap(nb,0:maxn))
    allocate(freqsqarr(nb,0:maxn))
    allocate(lambda(0:maxn))
    allocate(sev(nb,0:maxn))

    ndata = control%bcndata
    write(*,*) 'ndata = ',ndata

    allocate(xval(ndata))
    allocate(yval(ndata))

    allocate(MV(ndata,0:ndeg))
    allocate(MT(0:ndeg,ndata))
    allocate(MTM(0:ndeg,0:ndeg))
    allocate(MTy(0:ndeg))
    allocate(s(0:ndeg))
    allocate(sol(0:ndeg,1))
    allocate(ipiv(0:ndeg))

    allocate(store(nb))
    allocate(freqarr(nb))

    allocate(Dm1(nb,nb))
    allocate(Dm2(nb,nb))
    allocate(DeltaD(nb,nb))
    allocate(tmpv1(nb,nb))
    allocate(tmpv2(nb,nb))
    allocate(tmpv3(nb,nb))
    allocate(T1(nb,nb))
    allocate(T2(nb,nb))
    allocate(w1(nb))
    allocate(w2(nb))
    allocate(w3(nb))
    allocate(indarr(nb))
    allocate(tmpnu(nb))

    allocate(ig(nb+1))
    allocate(deltawsqr(nb))
    allocate(predictE(nb))
    allocate(localmap1(nb))
    allocate(localmap2(nb))

    allocate(sortarr(nb))

    allocate(prefer(nb))

    write(*,*)
    write(*,*) 'Start to do real stuff.'

    absbcbsf = trim(control%dir)//'/'//trim(bcbsf)
    write(*,'(A)') 'absbcbsf is '//trim(absbcbsf)
    abspbcbsf = trim(control%dir)//'/'//trim(pbcbsf)
    write(*,'(A)') 'abspbcbsf is '//trim(abspbcbsf)
    open(bcbsu,file=trim(absbcbsf),status='replace')
    open(pbcbsu,file=trim(abspbcbsf),status='replace')

    ind = 0
    distold = zero

    do i = 1, nsegments

      nstages = segn(i)

      call frac2abs(s1%b(1,1),bs%bvec(1,i),bq(1))
      call frac2abs(s1%b(1,1),bs%evec(1,i),eq(1))
      dist = vecmag3(eq(1:3)-bq(1:3))
      write(*,*) 'i, delta=',i,dist/nstages

      bq(1:3) = bs%bvec(1:3,i)
      eq(1:3) = bs%evec(1:3,i)

      del(1:3) = (eq(1:3)-bq(1:3))/(nstages*one)

      write(*,*) 'i,nstages= ',i,nstages
      do stage = 0, nstages
        lambda(stage) = stage*1.0D0/(nstages*1.0D0)
      enddo

      qvec = bq(1:3)
      call cal_dynmat(control,phyunit,qvec(1),fconst(1,1,1,1),loto,mass(1),Dm1(1,1),nb,s1,s2,m21(1),ns1,ns2)

      if(control%dynmat == 1 .or. control%dynmat == 2) then
        if(control%dynmat == 2) then
          do ii = 1, nb
            do jj = 1, ii
              cmplv = (Dm1(ii,jj) + conjg(Dm1(jj,ii)))/2.0d0
              Dm1(ii,jj) = cmplv
              Dm1(jj,ii) = conjg(cmplv)
            enddo
          enddo
        endif
        do ii = 1, nb
          do jj = 1, nb
            Dm1(jj,ii) = Dm1(jj,ii)*phyunit%omega2nuRydinvh**2
          enddo
        enddo
      endif

      tmpv1(1:nb,1:nb) = Dm1(1:nb,1:nb)
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

      Dm2(1:nb,1:nb) = Dm1(1:nb,1:nb)
      w2(1:nb) = w1(1:nb)
      tmpv2(1:nb,1:nb) = tmpv1(1:nb,1:nb)

      do stage = 1, nstages

        Dm1(1:nb,1:nb) = Dm2(1:nb,1:nb)
        w1(1:nb) = w2(1:nb)
        tmpv1(1:nb,1:nb) = tmpv2(1:nb,1:nb)

        qvec(1:3) = qvec(1:3) + del(1:3)
        call cal_dynmat(control,phyunit,qvec(1),fconst(1,1,1,1),loto,mass(1),Dm2(1,1),nb,s1,s2,m21(1),ns1,ns2)

        if(control%dynmat == 1 .or. control%dynmat == 2) then
          if(control%dynmat == 2) then
            do ii = 1, nb
              do jj = 1, ii
                cmplv = (Dm2(ii,jj) + conjg(Dm2(jj,ii)))/2.0d0
                Dm2(ii,jj) = cmplv
                Dm2(jj,ii) = conjg(cmplv)
              enddo
            enddo
          endif

          do ii = 1, nb
            do jj = 1, nb
              Dm2(jj,ii) = Dm2(jj,ii)*phyunit%omega2nuRydinvh**2
            enddo
          enddo
        endif

        tmpv2(1:nb,1:nb) = Dm2(1:nb,1:nb)
        call zheev('V','U',nb,tmpv2(1,1),nb,w2(1),CWORK(1),LCWORK,RWORK(1),info)
        if(info /= 0) then
          write(*,*) 'info = ',info
          write(*,'(A)') 'err: zheev error.'
          stop 1
        endif
        freqsqarr(1:nb,stage) = w2(1:nb)

        predictE(1:nb) = zero
        sortarr(1:nb) = zero

        if(stage < minstage) then

          write(*,*) 'stage = ',stage

          write(*,*) 'degentol= ',degentol

          if(stage == 4702) then
            write(*,*) 'w1: '
            write(*,'(500F20.10)') w1(1:nb)
          endif

          do k = 1, nb
            freqsq = w1(k)
            freq = sqrt(abs(freqsq))
            if(freqsq < zero) then
              tmpnu(k) = -freq*Rydoverh2invcm
            else
              tmpnu(k) = freq*Rydoverh2invcm
            endif
          enddo

          call form_groups(nb,tmpnu(1),ig(1),ng,degentol)

          deltaD(1:nb,1:nb) = Dm2(1:nb,1:nb) - Dm1(1:nb,1:nb)

          do k = 1, ng
            dimen = ig(k+1) - ig(k)

            if(dimen > 1) then
              write(*,*) 'dimen= ',dimen
            endif

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
            sortarr(k) = predictE(k)
          enddo

          call dbl_sort(nb,sortarr(1),indarr(1),1)
          call inversemapping(nb,indarr(1),localmap1(1))

        else

          do k = 1, nb
            store(k) = k
          enddo

          do k = 1, nb

            do dp = 1, ndata

              ind = stage-ndata+dp-1

              if(ind < 0 .or. ind > nstages) then
                write(*,*) 'err: array out of bound: ind= ',ind
                stop 1
              endif
              xval(dp) = lambda(ind)
              yval(dp) = aev(k,ind)
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
          ind = gmap(k,stage)

          pev(k,stage) = sortarr(ind)
          aev(k,stage) = w2(ind)
        enddo
      enddo

      do stage = nstages-ndata, 0, -1
        write(*,*) 'second pass: stage= ',stage

        do k = 1, nb
          store(k) = k
        enddo

        prefer_choice=3
        if(prefer_choice == 1) then
          do k = 1, nb
            prefer(k) = k
          enddo
        else if(prefer_choice == 2) then
          do k = 1, nb
            prefer(k) = nb-k+1
          enddo
        else if(prefer_choice == 3) then
          do k = 1, nb
            indarr(k) = k
          enddo
          do k = 1, nb
            sortarr(k) = aev(k,nstages)
          enddo

          call dbl_sort(nb,sortarr(1),indarr(1),1)
          do k = 1, nb
            prefer(k) = indarr(k)
          enddo
        endif

        do k = 1, nb
          rk = prefer(k)
          do dp = 1, ndata
            ind = stage+dp
            xval(dp) = lambda(ind)

            yval(dp) = aev(rk,ind)
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

          aev(rk,stage) = freqsqarr(  store(recordkk), stage )

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
        do stage = 0, nstages
          sev(k,stage) = pev(indarr(k),stage)
        enddo
      enddo

      do k = 1, nb
        do stage = 0, nstages
          pev(k,stage) = sev(k,stage)
        enddo
      enddo

      do stage = 0, nstages
        disti = stage*dist/(nstages*one)
        do k = 1, nb
          freqsq = aev(k,stage)
          freq = sqrt(abs(freqsq))
          if(freqsq < zero) then
            tmpnu(k) = -freq*Rydoverh2invcm
          else
            tmpnu(k) = freq*Rydoverh2invcm
          endif
        enddo
        write(bcbsu,'(F12.6,500F20.10)') distold + disti,tmpnu(1:nb)

      enddo

      do stage = 0, nstages
        disti = stage*dist/(nstages*one)
        do k = 1, nb
          freqsq = pev(k,stage)
          freq = sqrt(abs(freqsq))
          if(freqsq < zero) then
            tmpnu(k) = -freq*Rydoverh2invcm
          else
            tmpnu(k) = freq*Rydoverh2invcm
          endif
        enddo
        write(pbcbsu,'(F12.6,500F20.10)') distold + disti,tmpnu(1:nb)
      enddo

      distold = distold + dist

    enddo

    close(bcbsu)
    close(pbcbsu)

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

  end subroutine bandconnect

  subroutine cal_band_struct(control,phyunit,m21,dym,mass,fconst,loto,ns1,ns2,nb,s1,s2,bs)
    type(controlt) :: control
    type(unitt) :: phyunit
    integer :: ns1,ns2,nb
    type(bst) :: bs
    real(double) :: bq(3),eq(3),qvec(3),qvec_2(3),qvec_3(3)
    type(supercell) :: s1,s2
    real(double) :: fconst(3,ns2,3,ns1)
    integer :: m21(ns2)
    real(double) :: mass(ns1)
    type(lotot) :: loto
    complex(double) :: dym(nb,nb)

    complex(double),allocatable :: omega(:)
    integer :: LCWORK,LRWORK,LDVL,LDVR,INFO

    complex(double),allocatable :: CW(:)
    complex(double),allocatable :: VL(:,:),VR(:,:),CWORK(:)
    real(double),allocatable :: RWORK(:)

    integer :: ii,jj
    complex(double) :: cmplv
    integer :: totalp,ind
    integer :: kk,mm
    integer :: i,j,k,nsegments,npoints
    complex(double) :: zv(3)
    real(double) :: area,totalarea
    real(double) :: fracvec(3)
    real(double) :: distold,dist,disti,r1,r2
    real(double) :: freq,romegasq
    real(double),allocatable :: pev(:,:)
    real(double),allocatable :: aev(:,:)
    real(double),allocatable :: aevinvcm(:,:)
    integer :: nstages
    real(double) :: degentol
    complex(double),allocatable :: D1(:,:),D2(:,:),D3(:,:),Di1(:,:),Di2(:,:),deltaD(:,:)
    complex(double),allocatable :: tmpv1(:,:),tmpv2(:,:),tmpv3(:,:)
    complex(double),allocatable :: T1(:,:)
    complex(double),allocatable :: T2(:,:)
    real(double),allocatable :: w1(:),w2(:),w3(:)
    integer,allocatable :: indarr(:)
    real(double),allocatable :: zfreq(:)
    real(double),allocatable :: pfreq(:)
    real(double),allocatable :: mfreq(:)
    character(len=DSL) :: absrawdf

    type(bst) :: bs_2
    integer :: ns1_2,ns2_2
    real(double),allocatable :: fconst_2(:,:,:,:)
    type(controlt) :: control_2
    type(unitt) :: phyunit_2
    integer,allocatable :: m21_2(:)
    real(double),allocatable :: mass_2(:)
    type(lotot) :: loto_2
    type(supercell) :: s1_2
    type(supercell) :: s2_2
    complex(double),allocatable :: dym_2(:,:)

    type(bst) :: bs_3
    integer :: ns1_3,ns2_3
    real(double),allocatable :: fconst_3(:,:,:,:)
    type(controlt) :: control_3
    type(unitt) :: phyunit_3
    integer,allocatable :: m21_3(:)
    real(double),allocatable :: mass_3(:)
    type(lotot) :: loto_3
    type(supercell) :: s1_3
    type(supercell) :: s2_3
    complex(double),allocatable :: dym_3(:,:)
    real(double) :: diffv(3)

    type(bstype) :: gru_bs
    type(gvtype) :: groupv
    real(double) :: GPeps,GPscalefactor
    real(double) :: Binv(3,3),Cartdq(3),dq(3)
    real(double) :: Atrans(3,3)
    real(double) :: gvdelta,gv
    integer :: cartesian_comp
    real(double),allocatable :: gparr(:)
    integer,parameter :: minstage = 4

    allocate(omega(nb))
    LRWORK = nb*3
    LCWORK = nb*2
    LDVL=nb
    LDVR=nb
    allocate(CW(nb))
    allocate(RWORK(LRWORK))
    allocate(CWORK(LCWORK))
    allocate(VL(nb,nb),VR(nb,nb))

    nsegments = bs%nsegments
    npoints = bs%npoints
    totalp = nsegments*npoints

    gru_bs%ns = nsegments
    gru_bs%np = npoints
    gru_bs%nb = nb
    allocate(gru_bs%omega(nb,npoints,nsegments))

    groupv%ns = nsegments
    groupv%np = npoints
    groupv%nb = nb
    allocate(groupv%velocity(3,nb,npoints,nsegments))
    allocate(gparr(nb))

    distold = zero
    do i = 1, nsegments
      call frac2abs(s1%b(1,1),bs%bvec(1,i),bq(1))
      call frac2abs(s1%b(1,1),bs%evec(1,i),eq(1))
      dist = vecmag3(eq(1:3) - bq(1:3))
      do j = 1, npoints
        disti = (j-one)*dist/(npoints-one)
        bs%dis(j,i) = distold + disti
      enddo
      distold = distold + dist
    enddo

    call transposen(3,s1%a(1,1),ATrans(1,1))
    Binv(1:3,1:3) = ATrans(1:3,1:3)/twopi
    write(*,*) 'Binv is ',Binv(1:3,1:3)

    absrawdf=trim(control%dir)//'/'//trim(rawdf)
    write(*,'(A)') 'absrawdf is '//trim(absrawdf)
    call writerawd(rawdu,absrawdf,control,phyunit,ns1,ns2,m21(1),mass(1),fconst(1,1,1,1),loto,s1,s2,bs)

    if(control%gv == 1) then

      call readrawd(rawdu,absrawdf,control,phyunit,ns1,ns2,m21(1),mass(1),fconst(1,1,1,1),loto,s1,s2,bs)
    endif

    if(control%GP == 1 .or. control%gv == 1) then

      write(*,*) 'ns1,ns2,nsegments,npoints,nb=',ns1,ns2,nsegments,npoints,nb
      allocate(fconst_2(3,ns2,3,ns1))
      allocate(bs_2%bvec(3,nsegments))
      allocate(bs_2%evec(3,nsegments))
      allocate(bs_2%q(3,npoints,nsegments))
      allocate(bs_2%freq(nb,npoints,nsegments))
      allocate(bs_2%dis(npoints,nsegments))
      allocate(bs_2%r(nb))
      allocate(m21_2(ns2))
      allocate(mass_2(ns1))
      allocate(loto_2%corfc(3,ns2,3,ns1))
      allocate(loto_2%zeu(3,3,ns1))
      allocate(dym_2(nb,nb))

      allocate(fconst_3(3,ns2,3,ns1))
      allocate(bs_3%bvec(3,nsegments))
      allocate(bs_3%evec(3,nsegments))
      allocate(bs_3%q(3,npoints,nsegments))
      allocate(bs_3%freq(nb,npoints,nsegments))
      allocate(bs_3%dis(npoints,nsegments))
      allocate(bs_3%r(nb))
      allocate(m21_3(ns2))
      allocate(mass_3(ns1))
      allocate(loto_3%corfc(3,ns2,3,ns1))
      allocate(loto_3%zeu(3,3,ns1))
      allocate(dym_3(nb,nb))
    endif
    if(control%GP == 1) then
      write(*,*) 'read positive epsilon'
      absrawdf=trim(control%dir)//'/'//trim(control%positivedir)//'/raw.dat'
      write(*,'(A)') 'absrawdf is '//trim(absrawdf)
      ns1_2 = ns1
      ns2_2 = ns2
      call readrawd(rawdu,absrawdf,control_2,phyunit_2,ns1_2,ns2_2,m21_2(1),mass_2(1),fconst_2(1,1,1,1),loto_2,s1_2,s2_2,bs_2)

      if(ns1_2 /= ns1) then
        write(*,*) 'err: ns1'
        stop 1
      endif
      if(ns2_2 /= ns2) then
        write(*,*) 'err: ns2'
        stop 1
      endif
      if(bs_2%npoints /= bs%npoints) then
        write(*,*) 'err: bs%npoints'
        stop 1
      endif
      if(bs_2%nsegments /= bs%nsegments) then
        write(*,*) 'err: bs%npoints'
        stop 1
      endif
      if(control_2%dynmat /= control%dynmat) then
        write(*,*) 'err: dynmat'
        stop 1
      endif
      do i = 1, bs_2%nsegments
        diffv(1:3) = bs_2%bvec(1:3,i) - bs%bvec(1:3,i)
        if(vecmag3(diffv) > 1.0d-5) then
          write(*,*) 'warning: bs_2%bvec(1:3,i)=',bs_2%bvec(1:3,i)
          write(*,*) '           bs%bvec(1:3,i)=',bs%bvec(1:3,i)
        endif
      enddo
      do i = 1, bs_2%nsegments
        diffv(1:3) = bs_2%evec(1:3,i) - bs%evec(1:3,i)
        if(vecmag3(diffv) > 1.0d-5) then
          write(*,*) 'warning: bs_2%evec(1:3,i)=',bs_2%evec(1:3,i)
          write(*,*) '           bs%evec(1:3,i)=',bs%evec(1:3,i)
        endif
      enddo
      if(loto_2%loto /= loto%loto) then
        write(*,*) 'err: loto'
        stop 1
      endif
      do i = 1, ns1_2
        if(abs(mass_2(i) - mass(i)) > 1.0d-5) then
          write(*,*) 'err: mass'
          stop 1
        endif
      enddo

      write(*,*) 'read negative epsilon'
      absrawdf=trim(control%dir)//'/'//trim(control%negativedir)//'/raw.dat'
      write(*,'(A)') 'absrawdf is '//trim(absrawdf)

      ns1_3 = ns1
      ns2_3 = ns2
      call readrawd(rawdu,absrawdf,control_3,phyunit_3,ns1_3,ns2_3,m21_3(1),mass_3(1),fconst_3(1,1,1,1),loto_3,s1_3,s2_3,bs_3)

      if(ns1_3 /= ns1) then
        write(*,*) 'err: ns1'
        stop 1
      endif
      if(ns2_3 /= ns2) then
        write(*,*) 'err: ns2'
        stop 1
      endif
      if(bs_3%npoints /= bs%npoints) then
        write(*,*) 'err: bs%npoints'
        stop 1
      endif
      if(bs_3%nsegments /= bs%nsegments) then
        write(*,*) 'err: bs%npoints'
        stop 1
      endif
      if(control_3%dynmat /= control%dynmat) then
        write(*,*) 'err: dynmat'
        stop 1
      endif
      do i = 1, bs_2%nsegments
        diffv(1:3) = bs_3%bvec(1:3,i) - bs%bvec(1:3,i)
        if(vecmag3(diffv) > 1.0d-5) then
          write(*,*) 'warning: bs_3%bvec(1:3,i)=',bs_3%bvec(1:3,i)
          write(*,*) '           bs%bvec(1:3,i)=',bs%bvec(1:3,i)
        endif
      enddo
      do i = 1, bs_2%nsegments
        diffv(1:3) = bs_3%evec(1:3,i) - bs%evec(1:3,i)
        if(vecmag3(diffv) > 1.0d-5) then
          write(*,*) 'warning: bs_3%evec(1:3,i)=',bs_3%evec(1:3,i)
          write(*,*) '           bs%evec(1:3,i)=',bs%evec(1:3,i)
        endif
      enddo
      if(loto_3%loto /= loto%loto) then
        write(*,*) 'err: loto'
        stop 1
      endif
      do i = 1, ns1_3
        if(abs(mass_3(i) - mass(i)) > 1.0d-5) then
          write(*,*) 'err: mass'
          stop 1
        endif
      enddo
    endif

    if(control%gv == 1) then

      absrawdf=trim(control%dir)//'/'//trim(rawdf)
      write(*,'(A)') 'absrawdf is '//trim(absrawdf)
      ns1_2 = ns1
      ns2_2 = ns2
      call readrawd(rawdu,absrawdf,control_2,phyunit_2,ns1_2,ns2_2,m21_2(1),mass_2(1),fconst_2(1,1,1,1),loto_2,s1_2,s2_2,bs_2)
      ns1_3 = ns1
      ns2_3 = ns2
      call readrawd(rawdu,absrawdf,control_3,phyunit_3,ns1_3,ns2_3,m21_3(1),mass_3(1),fconst_3(1,1,1,1),loto_3,s1_3,s2_3,bs_3)
    endif

    if(control%GP == 1 .or. control%gv == 1) then
      nstages = control%dynnstages
      if(nstages < 1) then
        write(*,*) 'err: control%dynnstages= ',control%dynnstages
        stop 1
      endif
      degentol = control%dyndegentol
      if(abs(degentol) < 1.0d-15) then
        write(*,*) 'err: control%dyndegentol= ',control%dyndegentol
        stop 1
      endif

      allocate(pev(nb,0:nstages))
      allocate(aev(nb,0:nstages))
      allocate(aevinvcm(nb,0:nstages))
      allocate(D1(nb,nb))
      allocate(D2(nb,nb))
      allocate(D3(nb,nb))
      allocate(Di1(nb,nb))
      allocate(Di2(nb,nb))
      allocate(DeltaD(nb,nb))
      allocate(tmpv1(nb,nb))
      allocate(tmpv2(nb,nb))
      allocate(tmpv3(nb,nb))
      allocate(T1(nb,nb))
      allocate(T2(nb,nb))
      allocate(w1(nb))
      allocate(w2(nb))
      allocate(w3(nb))
      allocate(indarr(nb))
      allocate(zfreq(nb))
      allocate(pfreq(nb))
      allocate(mfreq(nb))
    endif

    if(control%GP == 1) then
      GPeps = control%GPeps
      if(abs(GPeps) < 1.0d-5) then
        write(*,*) 'err: control%GPeps= ',GPeps
        stop 1
      endif
      GPscalefactor = control%GPscalefactor
      if(abs(GPeps) < 1.0d-5) then
        write(*,*) 'err: control%GPscalefactor= ',GPscalefactor
        stop 1
      endif
      open(unit=freqgpu,file=trim(freqgpf),status='replace')
      write(freqgpu,*) '#np=', npoints*nsegments*nb
      write(freqgpu,*) '#nb=', nb
    endif

    gvdelta = zero
    if(control%gv == 1) then
      gvdelta = control%gvdelta
      if(abs(gvdelta) < 1.0d-5) then
        write(*,*) 'err: control%gvdelta= ',control%gvdelta
        stop 1
      endif
      open(unit=gvbsu,file=trim(gvbsf),status='replace')
      write(gvbsu,'(A)') '# '//trim(N2str(nsegments))//' '//trim(N2str(npoints))//' '//trim(N2str(nb))
    endif

    ind = 0
    do i = 1, nsegments

      bq(1:3) = bs%bvec(1:3,i)
      eq(1:3) = bs%evec(1:3,i)

      do j = 1, npoints
        ind = ind + 1
        write(*,*) 'kpt ind = ',ind, '/ ',totalp

        if(npoints == 1) then
          write(*,*) 'npoints = ',npoints
          write(*,'(A)') 'err: npoints must be greater than 1'
          stop 1
        endif

        qvec(1:3) = bq(1:3) + (j-1d0)*(eq(1:3)-bq(1:3))/(npoints-1d0)
        bs%q(1:3,j,i) = qvec(1:3)

        if(control%GP == 1) then

          call cal_dynmat(control,phyunit,qvec(1),fconst(1,1,1,1),loto,mass(1),dym(1,1),nb,s1,s2,m21(1),ns1,ns2)

          if(control%GPscheme == 1 .or. control%GPscheme == 2) then

            qvec_2(1:3) = bs_2%bvec(1:3,i) + (j-1d0)*(bs_2%evec(1:3,i)-bs_2%bvec(1:3,i))/(npoints-1d0)
            call cal_dynmat(control,phyunit,qvec_2(1),fconst_2(1,1,1,1),loto_2,mass_2(1),dym_2(1,1),nb,s1_2,s2_2,m21_2(1),ns1_2,ns2_2)
          endif
          if(control%GPscheme == 1 .or. control%GPscheme == 3) then

            qvec_3(1:3) = bs_3%bvec(1:3,i) + (j-1d0)*(bs_3%evec(1:3,i)-bs_3%bvec(1:3,i))/(npoints-1d0)
            call cal_dynmat(control,phyunit,qvec_3(1),fconst_3(1,1,1,1),loto_3,mass_3(1),dym_3(1,1),nb,s1_3,s2_3,m21_3(1),ns1_3,ns2_3)
          endif

          if(control%dynmat == 3) then
            write(*,*) 'err: GP is active. We must use either Upper triangular or Hermitian average'
            stop 1
          endif

          if(control%dynmat == 2) then
            do ii = 1, nb
              do jj = 1, ii

                cmplv = (dym(ii,jj) + conjg(dym(jj,ii)))/2d0
                dym(ii,jj) = cmplv
                dym(jj,ii) = conjg(cmplv)

                cmplv = (dym_2(ii,jj) + conjg(dym_2(jj,ii)))/2d0
                dym_2(ii,jj) = cmplv
                dym_2(jj,ii) = conjg(cmplv)

                cmplv = (dym_3(ii,jj) + conjg(dym_3(jj,ii)))/2d0
                dym_3(ii,jj) = cmplv
                dym_3(jj,ii) = conjg(cmplv)
              enddo
            enddo
          endif

          do ii = 1, nb
            do jj = 1, nb
              dym(jj,ii) = dym(jj,ii)*phyunit%omega2nuRydinvh**2
              dym_2(jj,ii) = dym_2(jj,ii)*phyunit%omega2nuRydinvh**2
              dym_3(jj,ii) = dym_3(jj,ii)*phyunit%omega2nuRydinvh**2
            enddo
          enddo

          if(control%GPscheme == 1) then

            call pt_shiftedfreq(aev(1,0),pev(1,0),mfreq(1),pfreq(1),nstages,minstage,degentol,nb,dym_3(1,1),dym_2(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w3(1),w2(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))
            call pt_shiftedfreq(aev(1,0),pev(1,0),mfreq(1),zfreq(1),nstages,minstage,degentol,nb,dym_3(1,1),dym(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w3(1),w1(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))

          else if(control%GPscheme == 2) then
            call pt_shiftedfreq(aev(1,0),pev(1,0),zfreq(1),pfreq(1),nstages,minstage,degentol,nb,dym(1,1),dym_2(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w1(1),w2(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))

          else if(control%GPscheme == 3) then
            call pt_shiftedfreq(aev(1,0),pev(1,0),zfreq(1),mfreq(1),nstages,minstage,degentol,nb,dym(1,1),dym_3(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w1(1),w2(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))
          endif

          do k = 1, nb
            if(abs(zfreq(k)) < 1.0d-3) then
              zfreq(k) = 1.0d-3
            endif
            if(control%GPscheme == 1) then
              gru_bs%omega(k,j,i) = (-GPscalefactor/(four*GPeps))*(pfreq(k)**2-mfreq(k)**2)/zfreq(k)**2
            else if(control%GPscheme == 2) then
              gru_bs%omega(k,j,i) = (-GPscalefactor/(two*GPeps))*(pfreq(k)**2-zfreq(k)**2)/zfreq(k)**2
            else if(control%GPscheme == 3) then
              gru_bs%omega(k,j,i) = (GPscalefactor/(two*GPeps))*(mfreq(k)**2-zfreq(k)**2)/zfreq(k)**2
            endif
            write(freqgpu,*) zfreq(k),gru_bs%omega(k,j,i)
          enddo
        endif

        if(control%gv == 1) then
          do cartesian_comp = 1, 3

            call cal_dynmat(control,phyunit,qvec(1),fconst(1,1,1,1),loto,mass(1),dym(1,1),nb,s1,s2,m21(1),ns1,ns2)

            Cartdq(1:3) = zero
            Cartdq(cartesian_comp) = control%gvdelta
            call matvec3(Binv(1,1),Cartdq(1),dq(1))
            qvec_2(1:3) = qvec(1:3) + dq(1:3)
            call cal_dynmat(control,phyunit,qvec_2(1),fconst_2(1,1,1,1),loto_2,mass_2(1),dym_2(1,1),nb,s1_2,s2_2,m21_2(1),ns1_2,ns2_2)

            Cartdq(1:3) = zero
            Cartdq(cartesian_comp) = -control%gvdelta
            call matvec3(Binv(1,1),Cartdq(1),dq(1))
            qvec_3(1:3) = qvec(1:3) + dq(1:3)
            call cal_dynmat(control,phyunit,qvec_3(1),fconst_3(1,1,1,1),loto_3,mass_3(1),dym_3(1,1),nb,s1_3,s2_3,m21_3(1),ns1_3,ns2_3)

            if(control%dynmat == 2) then
              do ii = 1, nb
                do jj = 1, ii

                  cmplv = (dym(ii,jj) + conjg(dym(jj,ii)))/2.0d0
                  dym(ii,jj) = cmplv
                  dym(jj,ii) = conjg(cmplv)

                  cmplv = (dym_2(ii,jj) + conjg(dym_2(jj,ii)))/2.0d0
                  dym_2(ii,jj) = cmplv
                  dym_2(jj,ii) = conjg(cmplv)

                  cmplv = (dym_3(ii,jj) + conjg(dym_3(jj,ii)))/2.0d0
                  dym_3(ii,jj) = cmplv
                  dym_3(jj,ii) = conjg(cmplv)
                enddo
              enddo
            endif

            do ii = 1, nb
              do jj = 1, nb
                dym(jj,ii) = dym(jj,ii)*phyunit%omega2nuRydinvh**2
                dym_2(jj,ii) = dym_2(jj,ii)*phyunit%omega2nuRydinvh**2
                dym_3(jj,ii) = dym_3(jj,ii)*phyunit%omega2nuRydinvh**2
              enddo
            enddo

            if(control%gvscheme == 1) then
              call pt_shiftedfreq(aev(1,0),pev(1,0),mfreq(1),pfreq(1),nstages,minstage,degentol,nb,dym_3(1,1),dym_2(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w3(1),w2(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))
              call pt_shiftedfreq(aev(1,0),pev(1,0),mfreq(1),zfreq(1),nstages,minstage,degentol,nb,dym_3(1,1),dym(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w3(1),w1(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))

            else if(control%gvscheme == 2) then
              call pt_shiftedfreq(aev(1,0),pev(1,0),zfreq(1),pfreq(1),nstages,minstage,degentol,nb,dym(1,1),dym_2(1,1),deltaD(1,1),tmpv1(1,1),tmpv2(1,1),tmpv3(1,1),T1(1,1),T2(1,1),w1(1),w2(1),indarr(1),CWORK(1),LCWORK,RWORK(1),LRWORK,Di1(1,1),Di2(1,1))
            else
              write(*,*) 'bdiff not implemented.'
              write(*,*) 'err: gv scheme should be either cdiff or fdiff'
              stop 1
            endif

            do k = 1, nb
              if(abs(zfreq(k)) < 1.0d-3) then
                zfreq(k) = 1.0d-3
              endif
              if(control%gvscheme == 1) then

                gv = one/(four*gvdelta)*(pfreq(k)**2-mfreq(k)**2)/zfreq(k) * invcm * phyunit%length

              else if(control%gvscheme == 2) then

                gv = one/(two*gvdelta)*(pfreq(k)**2-zfreq(k)**2)/zfreq(k) * invcm * phyunit%length

              else
                write(*,*) 'Bad control%gvscheme= ',control%gvscheme
                stop 1
              endif
              groupv%velocity(cartesian_comp,k,j,i) = gv
            enddo
          enddo
          do k = 1, nb
            gparr(k) = vecmag3(  groupv%velocity(1,k,j,i) )
          enddo
          write(gvbsu,'(500F20.8)') bs%dis(j,i), gparr(1:nb)

          do k = 1, nb

          enddo
        endif

        call cal_dynmat(control,phyunit,qvec(1),fconst(1,1,1,1),loto,mass(1),dym(1,1),nb,s1,s2,m21(1),ns1,ns2)

        if(control%dynmat == 1 .or. control%dynmat == 2) then
          if(control%dynmat == 2) then
            do ii = 1, nb
              do jj = 1, ii
                cmplv = (dym(ii,jj) + conjg(dym(jj,ii)))/2d0
                dym(ii,jj) = cmplv
                dym(jj,ii) = conjg(cmplv)
              enddo
            enddo
          endif

          if(control%storedyn == 1) then
            write(allqu,'(/A/)') 'Dynamical  Matrix at q point'
            write(allqu,'(A,3F20.10)') 'q = ',qvec(1:3)
            write(allqu,*)
            do ii = 1, nb/3
              do jj = 1, nb/3
                write(allqu,'(I4,I4)') ii, jj
                write(allqu,'(6E20.10)') dym(  3*(ii-1)+1, 3*(jj-1)+1: 3*(jj-1)+3)*phyunit%omega2nuRydinvh**2
                write(allqu,'(6E20.10)') dym(  3*(ii-1)+2, 3*(jj-1)+1: 3*(jj-1)+3)*phyunit%omega2nuRydinvh**2
                write(allqu,'(6E20.10)') dym(  3*(ii-1)+3, 3*(jj-1)+1: 3*(jj-1)+3)*phyunit%omega2nuRydinvh**2
              enddo
            enddo
          endif
          call zheev('V','U',nb,dym(1,1),nb,bs%r(1),CWORK(1),LCWORK,RWORK(1),INFO)
          if(INFO /= 0) then
            write(*,*) 'error found in calling zheev, info is ',info
            write(*,'(A)') 'err: check zheev'
            stop 1
          endif
        else if(control%dynmat == 3) then

          call zgeev('n','v',nb,dym(1,1),nb,CW(1),VL(1,1),LDVL,VR(1,1),LDVR,CWORK(1),LCWORK,RWORK(1),INFO)
          if(INFO /= 0) then
            write(*,*) 'problem encountered in zgeev.'
            write(*,*) 'INFO =',INFO
            write(*,'(A)') 'err: zgeev failed'
            stop 1
          endif
        endif
        do k = 1, nb
          if(control%dynmat == 3) then
            omega(k) = sqrt(CW(k))*phyunit%freqfac
            bs%freq(k,j,i) = omega(k)
          else if(control%dynmat == 1 .or. control%dynmat == 2) then

            romegasq = bs%r(k)
            freq = sqrt(abs(romegasq))*phyunit%freqfac
            if(romegasq >= zero) then
              bs%freq(k,j,i) = freq
            else
              bs%freq(k,j,i) = -freq
            endif
          endif
        enddo

        if(control%dynmat == 1 .or. control%dynmat == 2) then

          do kk = 1, nb
            totalarea = zero
            do mm = 1, nb/3
              zv(1:3) = dym(3*(mm-1)+1: 3*(mm-1) + 3,kk)

              call find_ellipse_area2(zv(1),area,r1,r2)
              if(area > 0.6d0) then
                write(*,*) 'large area: zv(1:3) = ',zv(1:3)
                write(*,'(A,3F20.10)') 'fracvec is ',fracvec(1:3)
                write(*,'(A,I5)') 'band number is ',kk
                write(*,'(A,3F20.10)') 'r1,r2,r2/r1 =',r1,r2,r2/r1
                write(*,*)
              endif
              totalarea = totalarea + area
            enddo

          enddo
        else
          write(*,'(A)') 'Note: eigenvectors file is empty since we do not write the eigenvectors when use_full_mat is true.'
        endif
      enddo
    enddo
    if(control%GP == 1) then
      close(freqgpu)
    endif
    if(control%GP == 1 .or. control%gv == 1) then
      deallocate(fconst_2)
      deallocate(bs_2%bvec)
      deallocate(bs_2%evec)
      deallocate(bs_2%q)
      deallocate(bs_2%freq)
      deallocate(bs_2%dis)
      deallocate(bs_2%r)
      deallocate(m21_2)
      deallocate(mass_2)
      deallocate(loto_2%corfc)
      deallocate(loto_2%zeu)
      deallocate(dym_2)

      deallocate(fconst_3)
      deallocate(bs_3%bvec)
      deallocate(bs_3%evec)
      deallocate(bs_3%q)
      deallocate(bs_3%freq)
      deallocate(bs_3%dis)
      deallocate(bs_3%r)
      deallocate(m21_3)
      deallocate(mass_3)
      deallocate(loto_3%corfc)
      deallocate(loto_3%zeu)
      deallocate(dym_3)

      deallocate(D1)
      deallocate(D2)
      deallocate(D3)

      close(gvbsu)
    endif

    deallocate(gru_bs%omega)
    deallocate(groupv%velocity)
    deallocate(gparr)
    deallocate(omega)
    deallocate(CW)
    deallocate(VL,VR,CWORK)
    deallocate(RWORK)
  end subroutine cal_band_struct

  subroutine assign_mass(mass,ns1,s1)
    type(supercell) :: s1
    integer :: i,ns1
    real(double) :: mass(ns1)
    do i = 1, ns1
      mass(i) = massofa(s1%at(i)%z)
    enddo
  end subroutine assign_mass

  subroutine output_fconst(calfcu,fname,fconst,ns2,ns1)
    character(len=*) :: fname
    integer :: ns2,ns1,calfcu
    real(double) :: fconst(3,ns2,3,ns1)
    integer :: i,j,k,m

    open(unit=calfcu,file=trim(fname),status='replace')
    write(calfcu,*) ns1
    write(calfcu,*) ns2
    do i = 1, ns1
      do j = 1, 3
        do k = 1, ns2
          do m = 1, 3

            write(calfcu,'(A)') 'm='//trim(N2str(m))//'/3, k='//trim(N2str(k))//'/(ns2='//trim(N2str(ns2))//'), j='//trim(N2str(j))//'/3, i='//trim(N2str(i))//'/(ns1='//trim(N2str(ns1))//')'
            write(calfcu,*) fconst(m,k,j,i)
          enddo
        enddo
      enddo
    enddo
    close(unit=calfcu)
    write(*,'(A)') trim(fname)//' is generated'
  end subroutine output_fconst

end module qmod
