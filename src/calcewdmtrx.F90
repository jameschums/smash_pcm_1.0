!--------------------------------------------------------------------
  subroutine calcewdmtrx(cmo,energymo,fulldmtrx,ewdmtrx,ndim,nelec)
!--------------------------------------------------------------------
!
! Calculate energy-weighted and normal density matrix
!
! In  : cmo      (MO coefficient matrix)
!       energymo (MO energy)
!       ndim     (Dimension of basis functions)
!       nelec    (Number of electrons)
! Out : fulldmtrx(Full density matrix)
!       ewdmtrx  (Upper-triangle energy-weighted density matrix)
!
      implicit none
      integer,intent(in) :: ndim, nelec
      integer :: i, j, k, ij
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: cmo(ndim,ndim), energymo(ndim)
      real(8),intent(out) :: fulldmtrx(ndim,ndim), ewdmtrx(ndim*(ndim+1)/2)
!
      fulldmtrx= transpose(cmo)
!$OMP parallel do schedule(guided) private(ij)
      do i= ndim,1,-1
        ij= i*(i-1)/2
        do j= 1,i
          ewdmtrx(ij+j)= zero
          do k= 1,nelec
            ewdmtrx(ij+j)= ewdmtrx(ij+j)-fulldmtrx(k,i)*fulldmtrx(k,j)*energymo(k)
          enddo
          ewdmtrx(ij+j)= ewdmtrx(ij+j)*two
        enddo
      enddo
!$OMP end parallel do
!
      call dgemm('N','T',ndim,ndim,nelec,two,cmo,ndim,cmo,ndim,zero,fulldmtrx,ndim)
      return
end
