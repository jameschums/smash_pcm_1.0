!--------------------
  subroutine setmp2
!--------------------
!
! Set MP2 information
!
      use modmp2, only : ncore
      implicit none
      integer :: ncorecalc
!
      if(ncore == -1) ncore= ncorecalc()
      write(*,'("james setmp2.F90")')
      return
end
