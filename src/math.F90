! Copyright 2014-2017  Kazuya Ishimura
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!


!-----------------------------------------------------------------------------------------
  subroutine mtrxcanon(ortho,overlap,eigen,ndim,newdim,threshover,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------------------------
!
! Calculate canonicalization matrix
!
! The ortho matrix satisfiles (Ortho)-daggar * S * (Ortho) = I
! where S is the overlap matrix.
!
      implicit none
      integer,intent(in) :: ndim, nproc, myrank, mpi_comm
      integer,intent(out) :: newdim
      integer :: i, j, icount
      real(8),parameter :: one=1.0D+00
      real(8),intent(in) :: threshover
      real(8),intent(inout) :: overlap(ndim,ndim)
      real(8),intent(out) :: ortho(ndim,ndim), eigen(ndim)
      real(8) :: ecanon
!
! Diagonalize ortho matrix
!
      call diag('V','U',ndim,overlap,ndim,eigen,nproc,myrank,mpi_comm)
!
! Eliminate eigenvectors with small eigenvalues
!
      icount= 0
      do i= 1,ndim
        if(eigen(i) < threshover) icount= icount+1
      enddo
      newdim= ndim-icount
!
! Form canonical orthonormal matrix
!
!$OMP parallel do private(ecanon)
      do i= 1,newdim
        ecanon= one/sqrt(eigen(i+icount))
        do j= 1,ndim
          ortho(j,i)=overlap(j,i+icount)*ecanon
        enddo
      enddo
!$OMP end parallel do
      return
end


!---------------------------------------------------
  subroutine canonicalize(hmat,ortho,work,nao,nmo)
!---------------------------------------------------
!
! Canonicalize symmetric Hamiltonian matrix
!   In : hmat(nao*nao), elements are in upper triangle.
!   Out: hmat(nmo*nmo)
!
      implicit none
      integer,intent(in) :: nao, nmo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: ortho(nao*nmo)
      real(8),intent(inout) :: hmat(nao*nao), work(nao*nao)
!
      call dsymm('L','U',nao,nmo,one,hmat,nao,ortho,nao,zero,work,nao)
      call dgemm('T','N',nmo,nmo,nao,one,ortho,nao,work,nao,zero,hmat,nao)
      return
end


!-------------------------------------------------------------------------------------
  subroutine canonicalizep(hmat,ortho,work,work2,nao,nmo,idis,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------------------
!
! Canonicalize symmetric Hamiltonian matrix
!   In : hmat(nao*nao), elements are in upper triangle.
!   Out: hmat(nmo*nmo), matrix size is nao*nao
!
      implicit none
      integer,intent(in) :: nao, nmo, nproc, myrank, mpi_comm, idis(nproc,14)
      integer :: num, istart
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: ortho(nao,nmo)
      real(8),intent(inout) :: hmat(nao*nao), work(nao*nao), work2(*)
!
      num=idis(myrank+1,1)
      istart=idis(myrank+1,2)+1
      if(num > 0)then
        call dsymm('L','U',nao,num,one,hmat,nao,ortho(1,istart),nao,zero,work,nao)
        call dgemm('T','N',nmo,num,nao,one,ortho,nao,work,nao,zero,work2,nao)
      endif
      call para_allgathervr(work2,num*nao,hmat,idis(1,3),idis(1,4),nproc,mpi_comm)
!     call dsymm('L','U',nao,nmo,one,hmat,nao,ortho,nao,zero,work,nao)
!     call dgemm('T','N',nmo,nmo,nao,one,ortho,nao,work,nao,zero,hmat,nao)
      return
end


!--------------------------------------------------------
  subroutine orthonorm(hmat,ndim,nproc,myrank,mpi_comm)
!--------------------------------------------------------
!
! Gram-Schmidt orthonormalization
!
      implicit none
      integer,intent(in) :: ndim, nproc, myrank, mpi_comm
      integer :: i, j, k, ndest
      real(8),intent(inout) :: hmat(ndim,ndim)
      real(8) :: dum
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
!
      do i= 1,ndim
        dum= zero
        do k= 1,ndim
          dum= dum+hmat(k,i)*hmat(k,i)
        enddo
        dum= one/sqrt(dum)
        do k= 1,ndim
          hmat(k,i)= hmat(k,i)*dum
        enddo
        if(i == ndim) cycle
!$OMP parallel do private(dum,k)
        do j= i+1,ndim
          if(mod(j,nproc).ne.myrank)cycle
          dum= zero
          do k= 1,ndim
            dum= dum-hmat(k,i)*hmat(k,j)
          enddo
          do k= 1,ndim
            hmat(k,j)= hmat(k,j)+dum*hmat(k,i)
          enddo
!         dum=-ddot(ndim,hmat(1,i),1,hmat(1,j),1)
!         call daxpy(ndim,dum,hmat(1,i),1,hmat(1,j),1)
        enddo
        ndest= mod(i+1,nproc)
        call para_bcastr(hmat(1,i+1),ndim,ndest,mpi_comm)
      enddo
      return
end


!-------------------------------------------------------------------------
  subroutine diag(jobz,uplo,ndim,vector,lda,eigen,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------
!
! Diagonalize matrix
!
      use modparallel, only : master
      implicit none
      integer,intent(in) :: ndim, lda, nproc, myrank, mpi_comm
      integer :: info
      integer, allocatable :: iwork(:)
      real(8),intent(out) :: eigen(lda)
      real(8),intent(inout) :: vector(*)
      real(8), allocatable :: work(:)
      character(len=1),intent(in) :: jobz, uplo
!
      info= 0
!     call memset(ndim*ndim)
!     allocate(work(ndim*ndim))
!     call dsyev(jobz,uplo,ndim,vector,lda,eigen,work,ndim*ndim,info)
!     deallocate(work)
!     call memunset(ndim*ndim)
      call memset(3*ndim*ndim+45*ndim)
      allocate(iwork(10*ndim),work(3*ndim*ndim+35*ndim))
      call dsyevd(jobz,uplo,ndim,vector,lda,eigen,work,3*ndim*ndim+35*ndim,iwork,10*ndim,info)
      deallocate(iwork,work)
      call memunset(3*ndim*ndim+45*ndim)
!
      if(info /= 0) then
        if(master)write(*,'(" Error! Diagonalization failed, info =",i5)')info
        call iabort
      endif
      return
end


!----------------------------------------------------------------------------------------------
  subroutine mtrxcanoninv(ortho,overinv,overlap,ndim,newdim,threshover,nproc,myrank,mpi_comm)
!----------------------------------------------------------------------------------------------
!
! Calculate canonicalization matrix and inverse matrix of overlap
!
! The ortho matrix satisfiles (Ortho)-daggar * S * (Ortho) = I
! where S is the overlap matrix.
!
      implicit none
      integer,intent(in) :: ndim, nproc, myrank, mpi_comm
      integer,intent(out) :: newdim
      integer :: i, j, icount
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: threshover
      real(8),intent(out) :: ortho(ndim,ndim), overinv(ndim,ndim)
      real(8),intent(inout) :: overlap(ndim,ndim)
      real(8),allocatable :: eigen(:), ecanon(:)
!
! Set arrays
!
      call memset(ndim*2)
      allocate(eigen(ndim),ecanon(ndim))
!
! Diagonalize ortho matrix
!
      call diag('V','U',ndim,overlap,ndim,eigen,nproc,myrank,mpi_comm)
!
! Eliminate eigenvectors with small eigenvalues
!
      icount= 0
      do i= 1,ndim
        if(eigen(i) < threshover) icount= icount+1
      enddo
      newdim= ndim-icount
!
! Form canonicalization and inverse matrix
!
!$OMP parallel
!$OMP do
      do i= 1,ndim
        eigen(i)  = one/eigen(i)
      enddo
!$OMP enddo
!$OMP do
      do i= 1,newdim
        ecanon(i)= sqrt(eigen(i+icount))
      enddo
!$OMP enddo
!$OMP do
      do i= 1,ndim
        do j= 1,ndim
          ortho(j,i) =overlap(j,i)*eigen(i)
        enddo
      enddo
!$OMP enddo
!$OMP end parallel
      call dgemm('N','T',ndim,ndim,ndim,one,overlap,ndim,ortho,ndim,zero,overinv,ndim)
!$OMP parallel
!$OMP do
      do i= 1,newdim
        do j= 1,ndim
          ortho(j,i)=overlap(j,i+icount)*ecanon(i)
        enddo
      enddo
!$OMP enddo
!$OMP end parallel
!
! Unset arrays
!
      call memunset(ndim*2)
      deallocate(eigen,ecanon)
      return
end

