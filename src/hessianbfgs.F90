!----------------------------------------------------------------------
  subroutine hessianbfgs(ehess,coord,coordold,egrad,egradold,vec,ndim)
!----------------------------------------------------------------------
!
! Update hessian matrix using BFGS method
!
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: coord(ndim), coordold(ndim), egrad(ndim), egradold(ndim)
      real(8),intent(inout) :: ehess(ndim*(ndim+1)/2), vec(ndim,3)
      real(8) :: denom, factor, ddot
!
      do i= 1,ndim
        vec(i,1)= coord(i)-coordold(i)
        vec(i,2)= egrad(i)-egradold(i)
        vec(i,3)= zero
      enddo
      denom= ddot(ndim,vec(1,1),1,vec(1,2),1)
      denom= one/denom
!
      do i= 1,ndim
        ii= i*(i-1)/2
        do j= 1,i-1
          vec(i,3)= vec(i,3)+ehess(ii+j)*vec(j,1)
          vec(j,3)= vec(j,3)+ehess(ii+j)*vec(i,1)
        enddo
        vec(i,3)= vec(i,3)+ehess(ii+i)*vec(i,1)
      enddo
      factor= ddot(ndim,vec(1,1),1,vec(1,3),1)
      factor= one/factor
      do i= 1,ndim
        ii= i*(i-1)/2
        do j= 1,i
          ehess(ii+j)= ehess(ii+j)-factor*vec(j,3)*vec(i,3) &
&                                 +denom*vec(j,2)*vec(i,2)
        enddo
      enddo
!
      return
end

