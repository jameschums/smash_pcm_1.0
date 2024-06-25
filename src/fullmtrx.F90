!------------------------------------------
 subroutine fullmtrx(trimat,fullmat,ndim)
!------------------------------------------
!
! Copy triangler size matrix to full size matrix
!
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j, ii
      real(8),intent(in) :: trimat((ndim*(ndim+1))/2)
      real(8),intent(out) :: fullmat(ndim,ndim)
!
      write(*,'("james mp2 fullmtrx.F90")')
!
!$OMP parallel do private(ii)
      do i= 1,ndim
        ii= i*(i-1)/2
        do j= 1,i
          fullmat(j,i)= trimat(ii+j)
        enddo
      enddo
!$OMP end parallel do
      return
end

