!----------------------------------
  subroutine setdetails(mpi_comm)
!----------------------------------
!
! Read input file and set variables
!
      use modparallel, only : master
      use modecp, only : flagecp
      implicit none
      integer,intent(in) :: mpi_comm
!
! Set defaults before reading input file
!
      call setdefault1
!
! Read input data and open checkpoint file if necessary
!
      if(master) call opendatfile
      call readinput(mpi_comm)
!
! Set basis functions
!
      call setbasis(mpi_comm)
!
! Set ECP functions
!
      if(flagecp) call setecp(mpi_comm)
!
! Set maximum memory size
!
      call maxmemset
!
! Set number of electrons
!
      call setelectron
!
! Reset defaults after reading input file
!
      call setdefault2
!
! Set functional information and adjust the number of DFT grids
!
      call setdft
!
! Set functional information and adjust the number of DFT grids
!
      call setmp2
!
! Write input data
!
      call writecondition
      call writegeom
      call writebasis
      if(flagecp) call writeecp
!
! Set atom charge including dummy atom
!
      call setcharge(mpi_comm)
!
      return
end


