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

module extra
  use qmod
  implicit none
contains
  subroutine abc
  end subroutine abc
end module extra

program sample
  use extra
  implicit none

  type(supercell) :: s1,s2
  integer,allocatable :: m12(:),m21(:)
  integer :: ns1,ns2,nb
  complex(double),allocatable :: dym(:,:)
  real(double),allocatable :: zf(:,:),nzf(:,:,:,:,:),fconst(:,:,:,:)
  real(double),allocatable :: mass(:)
  type(bst) :: bs
  integer :: i,nsegments,npoints

  type(unitt) :: phyunit
  character(len=DSL) :: tmpstr
  type(controlt) :: control
  character(len=DSL) :: absparfile
  character(len=DSL) :: absqpathfile
  character(len=DSL) :: abspposcar
  character(len=DSL) :: abssposcar
  character(len=DSL) :: absvasp_qpoints
  character(len=DSL) :: absellipse_areaf
  character(len=DSL) :: absallqf
  character(len=DSL) :: zstarfile
  character(len=DSL) :: filestem,inputformat
  type(lotot) :: loto
  integer,allocatable :: Infinite_NumGr_At(:),Infinite_MemNum_GrAt(:,:)
  integer,allocatable :: Finite_NumGr_At(:),Finite_MemNum_GrAt(:,:)
  real(double),allocatable :: Infinite_Dis_GrAt(:,:)
  real(double),allocatable :: Finite_Dis_GrAt(:,:)
  integer,allocatable :: CompleteShellInd(:)
  type(eigencheckert) :: eigenchecker
  character(len=DSL) :: lengthstr,energystr,massstr,absposcardir,absevecdir
  type(dtype) :: variable,defaultv
  real(double) :: tvec(3)

  write(*,*)
  write(*,*)

  call fargv(1,control%dir)
  if(trim(control%dir) == 'help') then
    write(*,'(A)') 'Usage: qphonon (1) working directory'
    write(*,'(A)') 'Prepare qphonon.par, p.vasp, supercell.vasp, vis.vasp'
  endif
  call fargn(1)

  absposcardir = trim(control%dir)//'/'//trim(poscardir)
  absevecdir = trim(control%dir)//'/'//trim(evecdir)
  call system("rm -rf "//trim(absposcardir))
  call system("rm -rf "//trim(absevecdir))
  call system("mkdir "//trim(absposcardir))
  call system("mkdir "//trim(absevecdir))

  abspposcar=trim(control%dir)//'/'//trim(pposcar)
  write(*,'(A)') 'opening file = '//trim(abspposcar)//', the primitive POSCAR'

  call read_struc(vaspu,trim(abspposcar),filestem,inputformat,s1)
  write(*,*) 's1 reading done.'

  abssposcar=trim(control%dir)//'/'//trim(sposcar)
  write(*,'(A)') 'opening file = '//trim(abssposcar)//', the wrapping supercell'

  call read_struc(vaspu,trim(abssposcar),filestem,inputformat,s2)
  write(*,*) 's2 reading done.'

  ns1 = s1%n
  ns2 = s2%n
  nb = 3*ns1
  write(*,*) 'dimension of the dynamical matrix is ',nb

  absparfile=trim(control%dir)//'/'//trim(parfile)
  open(unit=inu,file=trim(absparfile),action='read')
  write(*,'(A)') trim(absparfile)//' is opened.'
  close(inu)
  write(*,'(A)') trim(absparfile)//' is immediately closed.'

  variable%s = 'generic'
  call read_gendata(stype,inu,trim(absparfile),'title',variable,DEFAULTOPT,defaultv)
  control%title = trim(variable%s)
  write(*,'(A)') 'title is '//trim(control%title)

  write(*,'(A)') 'input units: '

  call read_gendata(stype,inu,trim(absparfile),'unit',variable,NONDEFAULTOPT,defaultv,vec=3)
  lengthstr = trim(variable%svec(1))
  energystr = trim(variable%svec(2))
  massstr = trim(variable%svec(3))

  call get_phys_unit(lengthstr,energystr,massstr,phyunit)

  call read_gendata(itype,inu,trim(absparfile),'frequnit',variable,NONDEFAULTOPT,defaultv)
  control%frequnit = variable%i
  if(control%frequnit == 1) then
    write(*,'(A)') 'Output frequency unit is invcm'
  else if(variable%i == 2) then
    write(*,'(A)') 'Output frequency unit is THz'
  else
    write(*,*) 'err: frequnit'
    stop 1
  endif

  call read_gendata(itype,inu,trim(absparfile),'diffscheme',variable,NONDEFAULTOPT,defaultv)
  control%diffscheme = variable%i
  if(control%diffscheme == 1) then
    write(*,'(A)') 'Central difference'
  else if(control%diffscheme == 2) then
    write(*,'(A)') 'Forward difference'
  else if(control%diffscheme == 3) then
    write(*,'(A)') 'Backward difference'
  else
    write(*,*) 'err: diffscheme'
    stop 1
  endif

  call read_gendata(itype,inu,trim(absparfile),'ASR',variable,NONDEFAULTOPT,defaultv)
  control%ASR = variable%i

  if(control%ASR == 0) then
    write(*,'(A)') 'Newton 3rd Law is not imposed on the force constants'
  else if(control%ASR == 1) then
    write(*,'(A)') 'diagonal-element Newton 3rd Law is imposed on the force constants'
  else if(control%ASR == 2) then
    write(*,'(A)') '1x1y1z-element Newton 3rd Law is imposed on the force constants'
  else
    write(*,'(A)') 'control%ASR = ',control%ASR
    write(*,'(A)') 'Invalid switch for Newton 3rd Law'
    write(*,'(A)') 'err: Check the switch, most likely you need to modify PHONON.par'
    stop 1
  endif

  call read_gendata(itype,inu,trim(absparfile),'iradcut',variable,NONDEFAULTOPT,defaultv)
  control%iradcut = variable%i

  write(*,*) 'control%iradcut = ',control%iradcut
  if(control%iradcut == 1) then
    write(*,*) 'Force constant switch: We use all interactions'
  else if(control%iradcut == 2) then
    write(*,*) 'Force constant switch: We use interactions with complete shells'
  else if(control%iradcut == 3) then
    write(*,*) 'Force constant switch: Advanced user-defined option. Prepare cutoff-radii.dat'
  else
    write(*,'(A)') 'err: control%iradcut is 1, 2, or 3'
    stop 1
  endif

  call read_gendata(itype,inu,trim(absparfile),'dynmat',variable,NONDEFAULTOPT,defaultv)
  control%dynmat = variable%i
  if(control%dynmat == 1) then
    write(*,'(A)') 'dynmat='//trim(N2str(1))//', use the upper half of the matrix'
  else if(control%dynmat == 2) then
    write(*,'(A)') 'dynmat='//trim(N2str(2))//', do a simple Hermitian average'
  else if(control%dynmat == 3) then
    write(*,'(A)') 'dynmat='//trim(N2str(3))//', full dynamical matrix without doing any symmetrization'
  else
    write(*,*) 'err: dynmat'
    stop 1
  endif

  call read_gendata(rtype,inu,trim(absparfile),'delta',variable,DEFAULTOPT,defaultv)
  control%delta = variable%r
  write(*,*) 'delta for atomic displacement is ',control%delta

  call read_gendata(itype,inu,trim(absparfile),'LOTO',variable,NONDEFAULTOPT,defaultv)
  if(variable%i == 0) then
    loto%loto = 0
    write(*,*) 'no LOTO'
  else if(variable%i == 1) then
    loto%loto = 1
    write(*,*) 'with LOTO'
  else
    write(*,*) 'err: loto'
    stop 1
  endif

  absqpathfile=trim(control%dir)//'/'//trim(qpathfile)
  open(unit=kinu,file=trim(absqpathfile),action='read')
  read(kinu,*) tmpstr, bs%nsegments, bs%npoints

  if(trim(tmpstr) == 'noBS') then
    control%BS = 0
    write(*,'(A)') 'We are NOT going to calculate the band structure.'
  else if(trim(tmpstr) == 'BS') then
    control%BS = 1
    write(*,'(A)') 'We are going to calculate the band structure.'
  else
    write(*,'(A)') 'err: Unrecognizable band structure switch, should be BS or noBS.'
    stop 1
  endif
  write(*,*) 'number of segments is ',bs%nsegments
  write(*,*) 'npoints along each segment is ',bs%npoints
  nsegments = bs%nsegments
  npoints = bs%npoints
  allocate(bs%bvec(3,nsegments))
  allocate(bs%evec(3,nsegments))

  allocate(bs%q(3,npoints,nsegments))
  allocate(bs%freq(nb,npoints,nsegments))
  allocate(bs%dis(npoints,nsegments))
  allocate(bs%r(nb))
  do i = 1, bs%nsegments
    read(kinu,*) bs%bvec(1:3,i),bs%evec(1:3,i)
    write(*,'(A,3F5.2,A,3F5.2)') 'read : bs%bvec(1:3,i),bs%evec(1:3,i) =',bs%bvec(1:3,i),', ',bs%evec(1:3,i)
  enddo
  close(kinu)

  call read_gendata(itype,inu,trim(absparfile),'evec',variable,NONDEFAULTOPT,defaultv)
  control%evec = variable%i
  if(control%evec == 0) then
    write(*,'(A)') 'NOT to analyse the eigen solutions for the atomic movements.'
  else if(control%evec == 1) then
    write(*,'(A)') 'To analyse the eigen solutions for the atomic movements.'

    open(unit=visu,file=trim(control%dir)//'/'//trim(visposcar),action='read')
    close(visu)
  else
    write(*,'(A)') 'err: evec'
    stop 1
  endif

  call read_gendata(rtype,inu,trim(absparfile),'qpoint',variable,NONDEFAULTOPT,defaultv,vec=3)
  control%special_kpt(1:3) = (/ variable%rvec(1), variable%rvec(2), variable%rvec(3) /)

  write(*,'(A,3F10.4)') 'read : special_kpt is ',control%special_kpt(1:3)

  defaultv%r = 4.0d0
  call read_gendata(rtype,inu,trim(absparfile),'evecamp',variable,DEFAULTOPT,defaultv)
  eigenchecker%amp_evec = variable%r
  write(*,*) 'eigenchecker%amp_evec = ',eigenchecker%amp_evec

  defaultv%r = 0.015d0
  call read_gendata(rtype,inu,trim(absparfile),'paraamp',variable,DEFAULTOPT,defaultv)
  eigenchecker%amp_deltaU = variable%r
  write(*,*) 'eigenchecker%amp_deltaU = ',eigenchecker%amp_deltaU

  defaultv%i = 11
  call read_gendata(itype,inu,trim(absparfile),'paranum',variable,DEFAULTOPT,defaultv)
  eigenchecker%n_deltaU = variable%i
  write(*,*) 'eigenchecker%n_deltaU = ',eigenchecker%n_deltaU

  defaultv%r = 3.0d0
  call read_gendata(rtype,inu,trim(absparfile),'visamp',variable,DEFAULTOPT,defaultv)
  eigenchecker%amp_vis = variable%r
  write(*,*) 'eigenchecker%amp_vis = ',eigenchecker%amp_vis

  defaultv%i = 15
  call read_gendata(itype,inu,trim(absparfile),'visnum',variable,DEFAULTOPT,defaultv)
  eigenchecker%n_vis = variable%i
  write(*,*) 'eigenchecker%n_vis = ',eigenchecker%n_vis

  if(eigenchecker%n_deltaU < 1 .or. eigenchecker%n_vis < 1) then
    write(*,'(A)') 'err: eigenchecker%n_deltaU,eigenchecker%n_vis= ',eigenchecker%n_deltaU,eigenchecker%n_vis
    stop 1
  endif

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'fcmulti',variable,DEFAULTOPT,defaultv)
  control%fcmulti = variable%i
  write(*,*) 'control%fcmulti = ',control%fcmulti

  defaultv%i=0
  call read_gendata(itype,inu,trim(absparfile),'storedyn',variable,DEFAULTOPT,defaultv)
  control%storedyn = variable%i
  write(*,*) 'control%storedyn = ',control%storedyn

  if(control%dynmat == 3 .and. control%storedyn == 1) then
    write(*,*) 'Warning: Since control%dynmat == 3, we do not store the dynamical matrix'

  endif

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'gv',variable,DEFAULTOPT,defaultv)
  control%gv = variable%i
  write(*,*) 'control%gv = ',control%gv

  defaultv%r = 0
  call read_gendata(rtype,inu,trim(absparfile),'gvdelta',variable,DEFAULTOPT,defaultv)
  control%gvdelta = variable%r
  write(*,*) 'control%gvdelta = ',control%gvdelta

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'gvscheme',variable,DEFAULTOPT,defaultv)
  control%gvscheme = variable%i
  write(*,*) 'control%gvscheme = ',control%gvscheme

  if(control%gvscheme == 0) then

  else if(control%gvscheme == 1) then
    write(*,'(A)') 'gv cdiff'
  else if(control%gvscheme == 2) then
    write(*,'(A)') 'gv fdiff'
  else if(control%gvscheme == 3) then
    write(*,'(A)') 'gv bdiff'
  else
    write(*,*) 'control%gvscheme = ',control%gvscheme
    write(*,'(A)') 'err: diff gvscheme not recognized'
    stop 1
  endif

  defaultv%i = -1
  call read_gendata(itype,inu,trim(absparfile),'GP',variable,DEFAULTOPT,defaultv)
  control%GP = variable%i
  write(*,*) 'control%GP = ',control%GP

  defaultv%i = -1
  call read_gendata(itype,inu,trim(absparfile),'GPscheme',variable,DEFAULTOPT,defaultv)
  control%GPscheme = variable%i
  write(*,*) 'control%GPscheme = ',control%GPscheme

  defaultv%s = 'xx'
  call read_gendata(stype,inu,trim(absparfile),'positivedir',variable,DEFAULTOPT,defaultv)
  control%positivedir = trim(variable%s)
  write(*,*) 'control%positivedir = '//trim(control%positivedir)

  defaultv%s = 'xx'
  call read_gendata(stype,inu,trim(absparfile),'negativedir',variable,DEFAULTOPT,defaultv)
  control%negativedir = trim(variable%s)
  write(*,*) 'control%negativedir = '//trim(control%negativedir)

  defaultv%r = 0
  call read_gendata(rtype,inu,trim(absparfile),'GPeps',variable,DEFAULTOPT,defaultv)
  control%GPeps = variable%r
  write(*,*) 'control%GPeps = ',control%GPeps

  defaultv%r = 0
  call read_gendata(rtype,inu,trim(absparfile),'GPscalefactor',variable,DEFAULTOPT,defaultv)
  control%GPscalefactor = variable%r
  write(*,*) 'control%GPscalefactor = ',control%GPscalefactor

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'dynnstages',variable,DEFAULTOPT,defaultv)
  control%dynnstages = variable%i
  write(*,*) 'control%dynnstages = ',control%dynnstages

  defaultv%r = 0
  call read_gendata(rtype,inu,trim(absparfile),'dyndegentol',variable,DEFAULTOPT,defaultv)
  control%dyndegentol = variable%r
  write(*,*) 'control%dyndegentol = ',control%dyndegentol

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'bc',variable,DEFAULTOPT,defaultv)
  control%bc = variable%i
  write(*,*) 'control%bc = ',control%bc

  defaultv%r = 0
  call read_gendata(rtype,inu,trim(absparfile),'bcdelta',variable,DEFAULTOPT,defaultv)
  control%bcdelta = variable%r
  write(*,*) 'control%bcdelta = ',control%bcdelta

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'bcsubn',variable,DEFAULTOPT,defaultv)
  control%bcsubn = variable%i
  write(*,*) 'control%bcsubn = ',control%bcsubn

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'bcndata',variable,DEFAULTOPT,defaultv)
  control%bcndata = variable%i
  write(*,*) 'control%bcndata = ',control%bcndata

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'bcndeg',variable,DEFAULTOPT,defaultv)
  control%bcndeg = variable%i
  write(*,*) 'control%bcndeg = ',control%bcndeg

  defaultv%i = 0
  call read_gendata(itype,inu,trim(absparfile),'bcminstage',variable,DEFAULTOPT,defaultv)
  control%bcminstage = variable%i
  write(*,*) 'control%bcminstage = ',control%bcminstage

  defaultv%r = 0
  call read_gendata(rtype,inu,trim(absparfile),'bcdegentol',variable,DEFAULTOPT,defaultv)
  control%bcdegentol = variable%r
  write(*,*) 'In invcm, control%bcdegentol = ',control%bcdegentol

  if(control%frequnit == 1) then
    phyunit%freqfac = sqrt(phyunit%energy/(phyunit%length**2*phyunit%mass))*(one/twopi)*Hz2invcm
  else if(control%frequnit == 2) then
    phyunit%freqfac = sqrt(phyunit%energy/(phyunit%length**2*phyunit%mass))*(one/twopi)*Hz2THz
  endif
  write(*,*) 'phyunit%freqfac = ',phyunit%freqfac
  phyunit%omega2nuRydinvh = sqrt(phyunit%energy/(phyunit%length**2*phyunit%mass))*(one/twopi)/Rydoverh

  write(*,*) 'ns1=',ns1
  allocate(mass(ns1))
  call assign_mass(mass(1),ns1,s1)

  allocate(m12(ns1))
  tvec(1:3) = (/zero,zero,zero/)
  call assign_m12_tvec(tvec(1),m12(1),ns1,s1,s2,distance_thres)
  write(*,'(A,200I4)') 'm12(1:ns1) = (',m12(1:ns1)
  write(*,'(A)') ')'

  allocate(m21(ns2))
  call assign_m21(s2,s1,m21(1),ns2)

  allocate(zf(3,ns2),nzf(3,ns2,3,ns1,2))
  write(*,'(A)') 'To read all forces'
  call read_allforces(ns1,ns2,zf(1,1),nzf(1,1,1,1,1),control%dir)

  allocate(Infinite_NumGr_At(ns1))
  allocate(Infinite_MemNum_GrAt(ngmax,ns1))
  allocate(Infinite_Dis_GrAt(ngmax,ns1))
  allocate(Finite_NumGr_At(ns1))
  allocate(Finite_Dis_GrAt(ns2,ns1))
  allocate(Finite_MemNum_GrAt(ns2,ns1))
  allocate(CompleteShellInd(ns1))

  call InfiniteCrystalStat(s1,ns1,Infinite_NumGr_At(1),Infinite_MemNum_GrAt(1,1),Infinite_Dis_GrAt(1,1))
  call FiniteCrystalStat(s1,ns1,Infinite_NumGr_At(1),Infinite_MemNum_GrAt(1,1),Infinite_Dis_GrAt(1,1),s2,ns2,Finite_NumGr_At(1),Finite_MemNum_GrAt(1,1),Finite_Dis_GrAt(1,1),CompleteShellInd(1))

  allocate(fconst(3,ns2,3,ns1))
  call cal_fconst(s1,s2,control,Finite_NumGr_At(1),Finite_Dis_GrAt(1,1),CompleteShellInd(1),m12(1),ns1,ns2,fconst(1,1,1,1),zf(1,1),nzf(1,1,1,1,1))

  call output_fconst(calfcu,trim(control%dir)//'/'//trim(calfconstsfile),fconst(1,1,1,1),ns2,ns1)

  loto%ns1 = ns1
  loto%ns2 = ns2

  allocate(loto%corfc(3,ns2,3,ns1))
  if(loto%loto == 1) then
    allocate(loto%zeu(3,3,ns1))
    zstarfile=trim(control%dir)//'/epsil.dat'
    write(*,'(A)') 'file= '//trim(zstarfile)//' is open'
    call read_epsil_and_Born_effective_charges(zstarfile,loto)
  endif

  allocate(dym(nb,nb))

  absvasp_qpoints=trim(control%dir)//'/'//trim(vasp_qpoints)
  write(*,*)
  write(*,'(A)') 'To open file = '//trim(absvasp_qpoints)
  open(unit=qtravu,file=trim(absvasp_qpoints),status='replace')
  write(qtravu,'(A)') 'q-points along high symmetry lines'
  write(qtravu,'(I6)') bs%npoints
  write(qtravu,'(A)') 'Line-mode'
  write(qtravu,'(A)') 'rec'
  do i = 1, bs%nsegments
    write(qtravu,'(3F20.10)') bs%bvec(1:3,i)
    write(qtravu,'(3F20.10)') bs%evec(1:3,i)
    write(qtravu,*)
  enddo
  close(qtravu)
  write(*,'(A)') 'file = '//trim(absvasp_qpoints)//' is now closed'

  if(control%BS == 1) then

    if(control%bc == 1) then
      call bandconnect(control,phyunit,m21(1),mass(1),fconst(1,1,1,1),loto,ns1,ns2,nb,s1,s2,bs)

    endif

    if(control%storedyn == 1) then
      absallqf=trim(control%dir)//'/'//trim(allqf)
      open(unit=allqu,file=trim(absallqf),status='replace')
    endif

    absellipse_areaf=trim(control%dir)//'/'//trim(evecdir)//'/'//trim(ellipse_areaf)
    open(unit=ellipse_areau,file=trim(absellipse_areaf),status='replace')

    call cal_band_struct(control,phyunit,m21(1),dym(1,1),mass(1),fconst(1,1,1,1),loto,ns1,ns2,nb,s1,s2,bs)

    close(ellipse_areau)
    call output_bs(bs,nb,control%dir)
    write(*,*) 'To plot the phonon dispersion relation, first execute the following command:'
    write(*,*) 'qph2bs or /opt1/bs-babel qphonon . redundant 0.0 0'
    if(control%storedyn == 1) then
      close(unit=allqu)
    endif
  endif

  if(control%bc == 1) then
    write(*,*) 'To plot the phonon connectivity: '
    write(*,*) 'xmgrace -nxy bc-bs.dat 1-anchor.dat'
  endif

end program sample
