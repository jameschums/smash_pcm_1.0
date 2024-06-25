!---------------------------------
  subroutine setcharge(mpi_comm)
!---------------------------------
!
! Set atom charge
!
      use modparallel, only : master
      use modiofile, only : input, maxline
      use modmolecule, only : znuc, natom
      use modparam, only : mxatom
      implicit none
      integer,intent(in) :: mpi_comm
      integer :: ii, jj, iatom
      real(8) :: znew
      character(len=254) :: line
!
      if(master) then
        rewind(input)
        do ii= 1,maxline
          read(input,'(a)',end=200)line
          if(line(1:6) == 'CHARGE') then
            write(*,'(/," -----------------")')
            write(*,'(  "   Atomic charge")')
            write(*,'(  " -----------------")')
            write(*,'(  "   Atomic charges are set manually.")')
            do jj= 1,mxatom
              read(input,'(a)',end=100) line
              read(line,*,end=100) iatom, znew
              znuc(iatom)= znew
              write(*,'("   Charge of Atom ",i5,"     ",f7.3)')iatom, znew
            enddo
          endif
        enddo
 100    write(*,*)
      endif
 200  continue
!
      call para_bcastr(znuc,natom,0,mpi_comm)
!
      return
end

