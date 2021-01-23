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
module sgsym
  use types_and_constants
  implicit none

  character(len=*),parameter :: sgtable='table-sg'
  character(len=*),parameter :: dsgtable='table-dsg'
  character(len=*),parameter :: pgtable='table-pg'
  character(len=*),parameter :: pgchartab='table-pg-char'

  real(double),parameter :: dis_tol_tryput = 1.0d-6

  real(double),parameter :: dev_tol_get_supercell_atom_map = 1.0d-8

  real(double),parameter :: FCCmat(1:3,1:3) = reshape( (/ zero,half,half, &
                                             half,zero,half,&
                                             half,half,zero/), (/3,3/) )

  real(double),parameter :: BCCmat(1:3,1:3) = reshape( (/ -half,half,half, &
                                             half,-half,half,&
                                             half,half,-half/), (/3,3/) )

  real(double),parameter :: Mat2Trigonal(1:3,1:3) = reshape(  (/ third,-third,third, &
                                                      third, two*third, third,&
                                                   -two*third,-third,third    /), (/3,3/)  )

  real(double),parameter :: Mat6Trigonal(1:3,1:3) = reshape( (/two*third,third,third,&
                                                 -third,third,third,&
                                                 -third,-two*third,third/),(/3,3/))

  real(double),parameter :: BCTmat(1:3,1:3) = reshape( (/ half,-half,half, &
                                             half,half,half,&
                                             -half,-half,half/), (/3,3/) )

  real(double),parameter :: OBaseCmat1(1:3,1:3) = reshape( (/ half,half,zero, &
                                             -half,half,zero,&
                                             zero,zero,one/), (/3,3/) )

  real(double),parameter :: OBaseCmat2(1:3,1:3) = reshape( (/ half,-half,zero, &
                                             half,half,zero,&
                                             zero,zero,one/), (/3,3/) )

  real(double),parameter :: OBaseAmat(1:3,1:3) = reshape( (/ one,zero,zero, &
                                             zero,half,-half,&
                                             zero,half,half/), (/3,3/) )

  real(double),parameter :: OFCmat(1:3,1:3) = reshape( (/ half,zero,half, &
                                             half,half,zero,&
                                             zero,half,half/), (/3,3/) )

  real(double),parameter :: OBCmat(1:3,1:3) = reshape( (/ half,half,half, &
                                             -half,half,half,&
                                             -half,-half,half/), (/3,3/) )

  real(double),parameter :: MonocBaxisBaseCentredCmat(1:3,1:3) = reshape( (/ half,-half,zero, &
                                             half,half,zero,&
                                             zero,zero,one/), (/3,3/) )

contains

  subroutine init_static_sg_database(qtablepath)
    character(len=*) :: qtablepath

    qtablepath='abc'

  end subroutine init_static_sg_database
end module sgsym
