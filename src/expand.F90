!----------------------------------------
  subroutine expand(array1,array2,ndim)
!----------------------------------------
!
! Expand upper-triangle matrix to full matrix
!
      implicit none
      integer,intent(in) :: ndim
      integer :: ij, i, j
      real(8),intent(in) :: array1(ndim*(ndim+1)/2)
      real(8),intent(out) :: array2(ndim,ndim)
!
      ij= 0
      do i= 1,ndim
        do j= 1,i
          ij= ij+1
          array2(j,i)= array1(ij)
        enddo
      enddo
      return
end


!----------------------------------------
  subroutine expand2(array1,array2,ndim)
!----------------------------------------
!
! Expand upper-triangle matrix to full matrix
!
      implicit none
      integer,intent(in) :: ndim
      integer :: ij, i, j
      real(8),intent(in) :: array1(ndim*(ndim+1)/2)
      real(8),intent(out) :: array2(ndim,ndim)
!
      ij= 0
      do i= 1,ndim
        do j= 1,i
          ij= ij+1
          array2(j,i)= array1(ij)
          array2(i,j)= array1(ij)
        enddo
      enddo
      return
end

