!--------------------------------------------------------------
  subroutine fixdtor(dcoordredun,numbond,numangle,numtorsion)
!--------------------------------------------------------------
!
! Fix dihedral displacements
!
      implicit none
      integer,intent(in) :: numbond, numangle, numtorsion
      integer :: ii
      real(8),parameter :: pi2=6.283185307179586D+00
      real(8),intent(inout) :: dcoordredun(numbond+numangle+numtorsion)
!
!$OMP parallel do
      do ii= numbond+numangle+1,numbond+numangle+numtorsion
        if(abs(dcoordredun(ii)+pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)+pi2
        if(abs(dcoordredun(ii)+pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)+pi2
        if(abs(dcoordredun(ii)-pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)-pi2
        if(abs(dcoordredun(ii)-pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)-pi2
      enddo
!$OMP end parallel do
!
      write(*,'("james mp2 fixdtor.F90 line 26")')
!
      return
end
