!-----------------------------------------------------------------------------------------
  subroutine calcuewdmtrx(cmoa,cmob,energymoa,energymob,fulldmtrx1,fulldmtrx2,ewdmtrx, &
&                         ndim,neleca,nelecb)
!-----------------------------------------------------------------------------------------
!
! Calculate energy-weighted and normal density matrix for open-shell
!
! In  : cmoa     (Alpha MO coefficient matrix)
!       cmob     (Beta MO coefficient matrix)
!       energymoa(Alpha MO energy)
!       energymob(Beta MO energy)
!       ndim     (Dimension of basis functions)
!       neleca   (Number of alpha electrons)
!       nelecb   (Number of beta electrons)
! Out : fulldmtrx1(Full alpha+beta density matrix)
!       fulldmtrx2(Full alpha-beta density matrix)
!       ewdmtrx   (Upper-triangle energy-weighted density matrix)
!
      implicit none
      integer,intent(in) :: ndim, neleca, nelecb
      integer :: i, j, k, ij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmoa(ndim,ndim), cmob(ndim,ndim), energymoa(ndim), energymob(ndim)
      real(8),intent(out) :: fulldmtrx1(ndim,ndim), fulldmtrx2(ndim,ndim)
      real(8),intent(out) :: ewdmtrx(ndim*(ndim+1)/2)
      real(8) :: dena, denb
!
      fulldmtrx1= transpose(cmoa)
      fulldmtrx2= transpose(cmob)
!$OMP parallel do schedule(guided) private(ij)
      do i= ndim,1,-1
        ij= i*(i-1)/2
        do j= 1,i
          ewdmtrx(ij+j)= zero
          do k= 1,neleca
            ewdmtrx(ij+j)= ewdmtrx(ij+j)-fulldmtrx1(k,i)*fulldmtrx1(k,j)*energymoa(k)
          enddo
          do k= 1,nelecb
            ewdmtrx(ij+j)= ewdmtrx(ij+j)-fulldmtrx2(k,i)*fulldmtrx2(k,j)*energymob(k)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call dgemm('N','T',ndim,ndim,neleca,one,cmoa,ndim,cmoa,ndim,zero,fulldmtrx1,ndim)
      call dgemm('N','T',ndim,ndim,nelecb,one,cmob,ndim,cmob,ndim,zero,fulldmtrx2,ndim)
!
!$OMP parallel do private(dena,denb)
      do i= 1,ndim
        do j= 1,ndim
          dena= fulldmtrx1(j,i)+fulldmtrx2(j,i)
          denb= fulldmtrx1(j,i)-fulldmtrx2(j,i)
          fulldmtrx1(j,i)= dena
          fulldmtrx2(j,i)= denb
        enddo
      enddo
!$OMP end parallel do
!
      return
end


