!---------------------------------------------------------------------------------------------
  subroutine uhfqc(focka,fockb,cmoa,cmob,qcrmax,qcgmna,qcgmnb,qcvec, &
&                  qcmat,qcmatsave,qceigen,overlap,xint, &
&                  qcworka,qcworkb,work,hfexchange,nao,nmo,nocca,noccb,nvira,nvirb,nshell, &
&                  maxdim,maxqcdiag,maxqcdiagsub,threshqc, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!---------------------------------------------------------------------------------------------
!
! Driver of Davidson diagonalization for quadratically convergent of UHF
!
      use modparallel, only : master
      implicit none
      integer,intent(in) :: nao, nmo, nocca, noccb, nvira, nvirb, nshell, maxdim, maxqcdiag
      integer,intent(in) :: maxqcdiagsub, nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: itdav, itqcdiag, ii, ij, jj, kk, ia, ib, istart, icount
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: overlap(nao*(nao+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: hfexchange, threshqc
      real(8),intent(inout) :: qcrmax(nshell*(nshell+1)/2)
      real(8),intent(inout) :: qcgmna(nao*(nao+1)/2), qcgmnb(nao*(nao+1)/2)
      real(8),intent(inout) :: qcvec(nocca*nvira+noccb*nvirb+1,maxqcdiagsub+1,2)
      real(8),intent(inout) :: qcmat(maxqcdiagsub,maxqcdiagsub)
      real(8),intent(inout) :: qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2)
      real(8),intent(inout) :: qceigen(maxqcdiagsub)
      real(8),intent(inout) :: qcworka(nao,nao), qcworkb(nao,nao), work(nao,nao)
      real(8),intent(inout) :: focka(nao,nao),fockb(nao,nao), cmoa(nao,nao), cmob(nao,nao)
      real(8) :: tmp, tmpa, tmpb, qcnorm, ddot, rotqc
!
      qcmatsave(:)= zero
!
! Calculate Fock matrix in MO basis
!
      call expand(focka,qcworka,nao)
      call dsymm('L','U',nao,nmo,one,qcworka,nao,cmoa,nao,zero,work,nao)
      call dgemm('T','N',nmo,nmo,nao,one,cmoa,nao,work,nao,zero,focka,nao)
      call expand(fockb,qcworka,nao)
      call dsymm('L','U',nao,nmo,one,qcworka,nao,cmob,nao,zero,work,nao)
      call dgemm('T','N',nmo,nmo,nao,one,cmob,nao,work,nao,zero,fockb,nao)
!
! Calculate initial (first and second) qc vector
!
      qcvec(1,1,1)= one
      qcvec(2:nocca*nvira+noccb*nvirb+1,1,1)= zero
      qcvec(1,1,2)= zero
      qcvec(1,2,1)= zero
!
      qcnorm= zero
!$OMP parallel private(ij) reduction(+:qcnorm)
!$OMP do
      do ii= 1,nvira
        ij=(ii-1)*nocca+1
        do jj= 1,nocca
          qcvec(ij+jj,1,2)= focka(jj,ii+nocca)
          qcvec(ij+jj,2,1)= focka(jj,ii+nocca)
          qcnorm= qcnorm+focka(jj,ii+nocca)*focka(jj,ii+nocca)
        enddo
      enddo
!$OMP end do
!$OMP do
      do ii= 1,nvirb
        ij=(ii-1)*noccb+nocca*nvira+1
        do jj= 1,noccb
          qcvec(ij+jj,1,2)= fockb(jj,ii+noccb)
          qcvec(ij+jj,2,1)= fockb(jj,ii+noccb)
          qcnorm= qcnorm+fockb(jj,ii+noccb)*fockb(jj,ii+noccb)
        enddo
      enddo
!$OMP end do
!$OMP end parallel
      qcnorm= one/sqrt(qcnorm)
!$OMP parallel do
      do ii= 2,nocca*nvira+noccb*nvirb+1
        qcvec(ii,2,1)= qcvec(ii,2,1)*qcnorm
      enddo
!$OMP end parallel do

!
! Start Davidson diagonalization
!
      itdav= 2
      do itqcdiag= 2,maxqcdiag
!
! Calculate Gmn
!
        call calcqcurmn(qcworka,qcworkb,qcvec,cmoa,cmob,work,nao,nocca,noccb,nvira,nvirb, &
&                       itdav,maxqcdiagsub)
        call calcudmax(qcworka,qcworkb,qcrmax,work,nproc2,myrank2,mpi_comm2)
        call calcqcugmn(qcgmna,qcgmnb,work,qcworka,qcworkb,qcrmax,xint,hfexchange,maxdim, &
&                       nao,nshell,nproc1,myrank1,mpi_comm1)
!
! Add two-electron integral contribution
!
        call expand(qcgmna,qcworka,nao)
        call dsymm('L','U',nao,nvira,one,qcworka,nao,cmoa(1,nocca+1),nao,zero,work,nao)
        call dgemm('T','N',nocca,nvira,nao,one,cmoa,nao,work,nao,zero,qcvec(2,itdav,2),nocca)
        call expand(qcgmnb,qcworka,nao)
        call dsymm('L','U',nao,nvirb,one,qcworka,nao,cmob(1,noccb+1),nao,zero,work,nao)
        call dgemm('T','N',noccb,nvirb,nao,one,cmob,nao,work,nao,zero, &
&                  qcvec(nocca*nvira+2,itdav,2),noccb)
!
! Add Fock matrix element contribution
!
        tmp= zero
!$OMP parallel private(kk) reduction(+:tmp)
!$OMP do collapse(2)
        do ia= 1,nvira
          do ii= 1,nocca
            kk= (ia-1)*nocca+ii+1
            tmp= tmp+focka(ii,ia+nocca)*qcvec(kk,itdav,1)
            qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+focka(ii,ia+nocca)*qcvec(1,itdav,1)
            do ib= 1,nvira
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+focka(ib+nocca,ia+nocca) &
&                                                 *qcvec((ib-1)*nocca+ii+1,itdav,1)
            enddo
            do ij= 1,nocca
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)-focka(ij,ii)*qcvec((ia-1)*nocca+ij+1,itdav,1)
            enddo
          enddo
        enddo
!$OMP end do
!$OMP do collapse(2)
        do ia= 1,nvirb
          do ii= 1,noccb
            kk= (ia-1)*noccb+ii+nocca*nvira+1
            tmp= tmp+fockb(ii,ia+noccb)*qcvec(kk,itdav,1)
            qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fockb(ii,ia+noccb)*qcvec(1,itdav,1)
            do ib= 1,nvirb
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fockb(ib+noccb,ia+noccb) &
&                                                 *qcvec((ib-1)*noccb+ii+nocca*nvira+1,itdav,1)
            enddo
            do ij= 1,noccb
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)-fockb(ij,ii) &
&                                                 *qcvec((ia-1)*noccb+ij+nocca*nvira+1,itdav,1)
            enddo
          enddo
        enddo
!$OMP end do
!$OMP end parallel
        qcvec(1,itdav,2)= tmp
!
! Calculate small matrix <b|A|b>
!
        qceigen(:)= zero
!$OMP parallel reduction(+:qceigen)
        do ii= 1,itdav
!$OMP do
          do jj= 1,nocca*nvira+noccb*nvirb+1
            qceigen(ii)= qceigen(ii)+qcvec(jj,ii,1)*qcvec(jj,itdav,2)
          enddo
!$OMP end do
        enddo
!$OMP end parallel
!
        istart= itdav*(itdav-1)/2
        do ii= 1,itdav
          qcmatsave(istart+ii)= qceigen(ii)
        enddo
        icount= 0
        do ii= 1,itdav
          do jj= 1,ii
            icount= icount+1
            qcmat(jj,ii)= qcmatsave(icount)
          enddo
        enddo
!
! Diagonalize small matrix
!
        call diag('V','U',itdav,qcmat,maxqcdiagsub,qceigen,nproc2,myrank2,mpi_comm2)
!
! Form correction vector
!
        qcnorm= zero
!$OMP parallel do reduction(+:qcnorm)
        do jj= 1,nocca*nvira+noccb*nvirb+1
          qcvec(jj,itdav+1,1)= zero
          do ii= 1,itdav
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1) &
&                               +qcmat(ii,1)*(qcvec(jj,ii,2)-qceigen(1)*qcvec(jj,ii,1))
          enddo
          qcnorm= qcnorm+qcvec(jj,itdav+1,1)*qcvec(jj,itdav+1,1)
        enddo
!$OMP end parallel do
        qcnorm= sqrt(qcnorm)
!
        if(master) write(*,'(5x,"QC Cycle ",i2," : Norm = ",1p,d10.3)') itdav-1, qcnorm
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
        qcvec(1,itdav+1,1)=-qcvec(1,itdav+1,1)/qceigen(1)
!$OMP parallel private(ij)
!$OMP do
        do ii= 1,nvira
          ij=(ii-1)*nocca+1
          do jj= 1,nocca
            qcvec(ij+jj,itdav+1,1)= qcvec(ij+jj,itdav+1,1)/ &
&                                  (focka(ii+nocca,ii+nocca)-focka(jj,jj)-qceigen(1))
          enddo
        enddo
!$OMP end do
!$OMP do
        do ii= 1,nvirb
          ij=(ii-1)*noccb+nocca*nvira+1
          do jj= 1,noccb
            qcvec(ij+jj,itdav+1,1)= qcvec(ij+jj,itdav+1,1)/ &
&                                  (fockb(ii+noccb,ii+noccb)-fockb(jj,jj)-qceigen(1))
          enddo
        enddo
!$OMP end do
!$OMP end parallel
!
! Schmidt orthogonalization
!
        do ii= 1,itdav
          tmp= zero
          do jj= 1,nocca*nvira+noccb*nvirb+1
            tmp= tmp-qcvec(jj,ii,1)*qcvec(jj,itdav+1,1)
          enddo
          do jj= 1,nocca*nvira+noccb*nvirb+1
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1)+tmp*qcvec(jj,ii,1)
          enddo
        enddo
!
! Renormalization of new vector
!
        qcnorm= sqrt(ddot(nocca*nvira+noccb*nvirb+1,qcvec(1,itdav+1,1),1,qcvec(1,itdav+1,1),1))
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
        if(itqcdiag ==(maxqcdiag)) then
          if(master) then
            write(*,'(" Error! Number of iteration for Quadratically convergent ",&
&                     "method exceeds maxqcdiag=",i3,".")') maxqcdiag
            write(*,'(" Set larger value for maxqcdiag in scf section.")')
          endif
          call iabort
        endif
!
! Reset Davidson diagonalization
!
        if(itdav == maxqcdiagsub) then
          do jj= 1,nocca*nvira+noccb*nvirb+1
            qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
          enddo
          do ii= 2,itdav
            do jj= 1,nocca*nvira+noccb*nvirb+1
              qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
            enddo
          enddo
!
          itdav=1
          cycle
        endif
!
        qcnorm= one/qcnorm
!$OMP parallel do
        do ii= 1,nocca*nvira+noccb*nvirb+1
          qcvec(ii,itdav+1,1)= qcvec(ii,itdav+1,1)*qcnorm
        enddo
!$OMP end parallel do
        itdav= itdav+1
      enddo
!
! End of Davidson diagonalization
!
! Calculate orbital rotation matrix
!
!$OMP parallel private(rotqc,tmp)
!$OMP do
      do jj= 1,nocca*nvira+noccb*nvirb+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
      enddo
!$OMP end do
      do ii= 2,itdav
!$OMP do
        do jj= 1,nocca*nvira+noccb*nvirb+1
          qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
        enddo
!$OMP end do
      enddo
!
      tmp= one/qcvec(1,1,1)
!$OMP do
      do jj= 2,nocca*nvira+noccb*nvirb+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*tmp
      enddo
!$OMP end do
!
! Rotate occupied MOs
!
!$OMP do
      do ii= 1,nocca
        do jj= 1,nvira
          rotqc= qcvec((jj-1)*nocca+ii+1,1,1)
          do kk= 1,nao
            cmoa(kk,ii)= cmoa(kk,ii)+rotqc*cmoa(kk,jj+nocca)
          enddo
        enddo
      enddo
!$OMP end do
!$OMP do
      do ii= 1,noccb
        do jj= 1,nvirb
          rotqc= qcvec((jj-1)*noccb+ii+nocca*nvira+1,1,1)
          do kk= 1,nao
            cmob(kk,ii)= cmob(kk,ii)+rotqc*cmob(kk,jj+noccb)
          enddo
        enddo
      enddo
!$OMP end do
!$OMP end parallel
!
! Schmidt orthonormalization of new occupied MOs
!
!$OMP parallel private(ij,tmpa,tmpb)
!$OMP do
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          work(jj,ii)= overlap(ij+jj)
          work(ii,jj)= overlap(ij+jj)
        enddo
      enddo
!$OMP end do
!
      do ii= 1,nmo
!$OMP do
        do jj= 1,nao
          qcworka(jj,1)= zero
          qcworkb(jj,1)= zero
          do kk= 1,nao
            qcworka(jj,1)= qcworka(jj,1)+work(kk,jj)*cmoa(kk,ii)
            qcworkb(jj,1)= qcworkb(jj,1)+work(kk,jj)*cmob(kk,ii)
          enddo
        enddo
!$OMP end do
        tmpa= zero
        tmpb= zero
        do kk= 1,nao
          tmpa= tmpa+qcworka(kk,1)*cmoa(kk,ii)
          tmpb= tmpb+qcworkb(kk,1)*cmob(kk,ii)
        enddo
        tmpa= one/sqrt(tmpa)
        tmpb= one/sqrt(tmpb)
!$OMP barrier
!$OMP do
        do kk= 1,nao
          cmoa(kk,ii)= cmoa(kk,ii)*tmpa
          cmob(kk,ii)= cmob(kk,ii)*tmpb
          qcworka(kk,1)= qcworka(kk,1)*tmpa
          qcworkb(kk,1)= qcworkb(kk,1)*tmpb
        enddo
!$OMP end do
        if(ii == nmo) cycle
!
!$OMP do
        do jj= ii+1,nmo
          tmpa= zero
          tmpb= zero
          do kk= 1,nao
            tmpa= tmpa-qcworka(kk,1)*cmoa(kk,jj)
            tmpb= tmpb-qcworkb(kk,1)*cmob(kk,jj)
          enddo
          do kk= 1,nao
            cmoa(kk,jj)= cmoa(kk,jj)+tmpa*cmoa(kk,ii)
            cmob(kk,jj)= cmob(kk,jj)+tmpb*cmob(kk,ii)
          enddo
        enddo
!$OMP end do
      enddo
!$OMP end parallel
!
      return
end

