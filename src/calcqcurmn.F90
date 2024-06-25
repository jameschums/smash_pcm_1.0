!------------------------------------------------------------------------------------------
  subroutine calcqcurmn(qcrmna,qcrmnb,qcvec,cmoa,cmob,work,nao,nocca,noccb,nvira,nvirb, &
&                       itdav,maxqcdiagsub)
!------------------------------------------------------------------------------------------
!
! Calculate Rmn for quadratically convergent method of UHF
!
      implicit none
      integer,intent(in) :: nao, nocca, noccb, nvira, nvirb, itdav, maxqcdiagsub
      integer :: ii, jj, ij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: qcvec(nocca*nvira+noccb*nvirb+1,maxqcdiagsub+1,2)
      real(8),intent(in) :: cmoa(nao,nao), cmob(nao,nao)
      real(8),intent(out) :: qcrmna(nao*nao), qcrmnb(nao*nao), work(nao,nao)
!
      call dgemm('N','N',nao,nvira,nocca,one,cmoa,nao,qcvec(2,itdav,1),nocca,zero,qcrmna,nao)
      call dgemm('N','T',nao,nao,nvira,one,qcrmna,nao,cmoa(1,nocca+1),nao,zero,work,nao)
!
      ij= 0
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          qcrmna(ij+jj)= work(jj,ii)+work(ii,jj)
        enddo
      enddo
!$OMP end parallel do
!
      call dgemm('N','N',nao,nvirb,noccb,one,cmob,nao,qcvec(nocca*nvira+2,itdav,1),noccb,zero, &
&                qcrmnb,nao)
      call dgemm('N','T',nao,nao,nvirb,one,qcrmnb,nao,cmob(1,noccb+1),nao,zero,work,nao)
!
      ij= 0
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          qcrmnb(ij+jj)= work(jj,ii)+work(ii,jj)
        enddo
      enddo
!$OMP end parallel do
!
      return
end


