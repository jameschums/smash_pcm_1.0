!--------------------------------------------------------------------------------------------
  subroutine rhfqc(fock,cmo,qcrmax,qcgmn,qcvec,qcmat,qcmatsave,qceigen,overlap,xint, &
&                  qcwork,work,hfexchange,nao,nmo,nocc,nvir,nshell,maxdim,maxqcdiag, &
&                  maxqcdiagsub,threshqc,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!--------------------------------------------------------------------------------------------
!
! Driver of Davidson diagonalization for quadratically convergent of RHF
!
      use modparallel, only : master
      implicit none
      integer,intent(in) :: nao, nmo, nocc, nvir, nshell, maxdim, maxqcdiag, maxqcdiagsub
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: itdav, itqcdiag, ii, ij, jj, kk, ia, ib, istart, icount
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: overlap(nao*(nao+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: hfexchange, threshqc
      real(8),intent(out) :: qcvec(nocc*nvir+1,maxqcdiagsub+1,2), qcrmax(nshell*(nshell+1)/2)
      real(8),intent(out) :: qcwork(nao,nao), qcmat(maxqcdiagsub,maxqcdiagsub)
      real(8),intent(out) :: qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2), qceigen(maxqcdiagsub)
      real(8),intent(out) :: work(nao,nao)
      real(8),intent(inout) :: fock(nao,nao), cmo(nao,nao)
      real(8),intent(inout) :: qcgmn(nao*(nao+1)/2)
      real(8) :: tmp, qcnorm, ddot, rotqc
!
      qcmatsave(:)= zero
!
! Calculate Fock matrix in MO basis
!
      call expand(fock,qcwork,nao)
      call dsymm('L','U',nao,nocc+nvir,one,qcwork,nao,cmo,nao,zero,work,nao)
      call dgemm('T','N',nocc+nvir,nocc+nvir,nao,one,cmo,nao,work,nao,zero,fock,nao)
!
! Calculate initial (first and second) qc vector
!
      qcvec(1,1,1)= one
      qcvec(2:nocc*nvir+1,1,1)= zero
      qcvec(1,1,2)= zero
      qcvec(1,2,1)= zero
!
      qcnorm= zero
!$OMP parallel do private(ij) reduction(+:qcnorm)
      do ii= 1,nvir
        ij=(ii-1)*nocc+1
        do jj= 1,nocc
          qcvec(ij+jj,1,2)= fock(jj,ii+nocc)
          qcvec(ij+jj,2,1)= fock(jj,ii+nocc)
          qcnorm= qcnorm+fock(jj,ii+nocc)*fock(jj,ii+nocc)
        enddo
      enddo
!$OMP end parallel do
      qcnorm= one/sqrt(qcnorm)
!$OMP parallel do
      do ii= 2,nocc*nvir+1
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
        call calcqcrmn(qcwork,qcvec,cmo,work,nao,nocc,nvir,itdav,maxqcdiagsub)
        call calcrdmax(qcwork,qcrmax,work,nproc2,myrank2,mpi_comm2)
        call formrdftfock(qcgmn,work,qcwork,qcrmax,xint,maxdim,hfexchange, &
&                         nproc1,myrank1,mpi_comm1)
!
! Add two-electron integral contribution
!
        call expand(qcgmn,qcwork,nao)
        call dsymm('L','U',nao,nvir,one,qcwork,nao,cmo(1,nocc+1),nao,zero,work,nao)
        call dgemm('T','N',nocc,nvir,nao,one,cmo,nao,work,nao,zero,qcvec(2,itdav,2),nocc)
!
! Add Fock matrix element contribution
!
        tmp= zero
!$OMP parallel do collapse(2) private(kk) reduction(+:tmp)
        do ia= 1,nvir
          do ii= 1,nocc
            kk= (ia-1)*nocc+ii+1
            tmp= tmp+fock(ii,ia+nocc)*qcvec(kk,itdav,1)
            qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fock(ii,ia+nocc)*qcvec(1,itdav,1)
            do ib= 1,nvir
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fock(ib+nocc,ia+nocc) &
&                                                 *qcvec((ib-1)*nocc+ii+1,itdav,1)
            enddo
            do ij= 1,nocc
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)-fock(ij,ii)*qcvec((ia-1)*nocc+ij+1,itdav,1)
            enddo
          enddo
        enddo
!$OMP end parallel do
        qcvec(1,itdav,2)= tmp
!
! Calculate small matrix <b|A|b>
!
        qceigen(:)= zero
!$OMP parallel reduction(+:qceigen)
        do ii= 1,itdav
!$OMP do
          do jj= 1,nocc*nvir+1
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
        do jj= 1,nocc*nvir+1
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
!$OMP parallel do private(ij)
        do ii= 1,nvir
          ij=(ii-1)*nocc+1
          do jj= 1,nocc
            qcvec(ij+jj,itdav+1,1)= qcvec(ij+jj,itdav+1,1)/ &
&                                  (fock(ii+nocc,ii+nocc)-fock(jj,jj)-qceigen(1))
          enddo
        enddo
!$OMP end parallel do
!
! Schmidt orthogonalization
!
        do ii= 1,itdav
          tmp= zero
          do jj= 1,nocc*nvir+1
            tmp= tmp-qcvec(jj,ii,1)*qcvec(jj,itdav+1,1)
          enddo
          do jj= 1,nocc*nvir+1
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1)+tmp*qcvec(jj,ii,1)
          enddo
        enddo
!
! Renormalization of new vector
!
        qcnorm= sqrt(ddot(nocc*nvir+1,qcvec(1,itdav+1,1),1,qcvec(1,itdav+1,1),1))
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
! Reset Davidson diagonalization
!
        if(itdav == maxqcdiagsub) then
          do jj= 1,nocc*nvir+1
            qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
          enddo
          do ii= 2,itdav
            do jj= 1,nocc*nvir+1
              qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
            enddo
          enddo
!    
          itdav=1
          cycle
        endif
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
        qcnorm= one/qcnorm
!$OMP parallel do
        do ii= 1,nocc*nvir+1
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
      do jj= 1,nocc*nvir+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
      enddo
!$OMP end do
      do ii= 2,itdav
!$OMP do
        do jj= 1,nocc*nvir+1
          qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
        enddo
!$OMP end do
      enddo
!
      tmp= one/qcvec(1,1,1)
!$OMP do
      do jj= 2,nocc*nvir+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*tmp
      enddo
!$OMP end do
!
! Rotate occupied MOs
!
!$OMP do
      do ii= 1,nocc
        do jj= 1,nvir
          rotqc= qcvec((jj-1)*nocc+ii+1,1,1)
          do kk= 1,nao
            cmo(kk,ii)= cmo(kk,ii)+rotqc*cmo(kk,jj+nocc)
          enddo
        enddo
      enddo
!$OMP end do
!$OMP end parallel
!
! Schmidt orthonormalization of new occupied MOs
!
!$OMP parallel private(ij,tmp)
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
          qcwork(jj,1)= zero
          do kk= 1,nao
            qcwork(jj,1)= qcwork(jj,1)+work(kk,jj)*cmo(kk,ii)
          enddo
        enddo
!$OMP end do
        tmp= zero
        do kk= 1,nao
          tmp= tmp+qcwork(kk,1)*cmo(kk,ii)
        enddo
        tmp= one/sqrt(tmp)
!$OMP barrier
!$OMP do
        do kk= 1,nao
          cmo(kk,ii)= cmo(kk,ii)*tmp
          qcwork(kk,1)= qcwork(kk,1)*tmp
        enddo
!$OMP end do
        if(ii == nmo) cycle
!
!$OMP do
        do jj= ii+1,nmo
          tmp= zero
          do kk= 1,nao
            tmp= tmp-qcwork(kk,1)*cmo(kk,jj)
          enddo
          do kk= 1,nao
            cmo(kk,jj)= cmo(kk,jj)+tmp*cmo(kk,ii)
          enddo
        enddo
!$OMP end do
      enddo
!$OMP end parallel
!
      return
end


