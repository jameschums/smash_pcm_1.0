!---------------------------------------------------------------------------------
  subroutine calcqcugmn(gmn1,gmn2,gmn3,rmna,rmnb,rmtrx,xint,hfexchange,maxdim, &
&                       nao,nshell,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------
!
! Driver of Gmn matrix formation from two-electron intgrals
!
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: nshell, nao, maxdim, nproc, myrank, mpi_comm
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: rmna(nao*(nao+1)/2), rmnb(nao*(nao+1)/2)
      real(8),intent(in) :: rmtrx(nshell*(nshell+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: hfexchange
      real(8),intent(out) :: gmn1(nao*(nao+1)/2), gmn2(nao*(nao+1)/2), gmn3(nao*(nao+1)/2)
      real(8) :: xijkl, rmax, twoeri(maxdim**4), rmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      gmn2(:)= zero
      gmn3(:)= zero
!
!      write(*,'("james calcqcugmn.F90 smash23 line 25")')
!
      ncount= 0
      ncount= ncount+(2*nshell**3+3*nshell**2+nshell)/6+myrank
      do ish= 1,nshell
        icount(ish)= ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
      enddo
!
      ish= nshell
      ii= ish*(ish-1)/2
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,rmax1,rmax,twoeri,jj,kk, &
!$OMP kstart,last,ltmp,lnum,ll) firstprivate(ish,ii) reduction(+:gmn2,gmn3)
      do ijsh= nshell*(nshell+1)/2,1,-1
        do ishcheck=1,nshell
          if(ijsh > ii) then
            jsh= ijsh-ii
            exit
          else
            ish= ish-1
            ii= ish*(ish-1)/2
          endif
        enddo
!
        ij= ii+jsh
        jj= jsh*(jsh-1)/2
        kstart=mod(icount(ish)-ish*(jsh-1),nproc)+1
        do ksh= kstart,ish,nproc
          kk= ksh*(ksh-1)/2
          ik= ii+ksh
          jk= jj+ksh
          if(jsh.lt.ksh) jk= kk+jsh
          rmax1=max(two*rmtrx(ij),rmtrx(ik),rmtrx(jk))
          last= ksh
          if(ish == ksh) last= jsh
          ll=min(jsh,ksh)
          lnum=0
!         do lsh= 1,ksh
          do lsh= 1,ll
            kl= kk+lsh
            il= ii+lsh
            jl= jj+lsh
            xijkl=xint(ij)*xint(kl)
            rmax=max(rmax1,two*rmtrx(kl),rmtrx(il),rmtrx(jl))
            if(xijkl*rmax.ge.cutint2) then
              lnum=lnum+1
              ltmp(lnum)=lsh
            endif
          enddo
          do lsh= ll+1,last
            kl= kk+lsh
            il= ii+lsh
            jl= lsh*(lsh-1)/2+jsh
            xijkl=xint(ij)*xint(kl)
            rmax=max(rmax1,two*rmtrx(kl),rmtrx(il),rmtrx(jl))
            if(xijkl*rmax.ge.cutint2) then
              lnum=lnum+1
              ltmp(lnum)=lsh
            endif
          enddo
          do lsh= 1,lnum
            call calc2eri(twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
            call ugmneri(gmn2,gmn3,rmna,rmnb,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim,hfexchange)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      do ii= 1,nao
        gmn2(ii*(ii+1)/2)= gmn2(ii*(ii+1)/2)*two
        gmn3(ii*(ii+1)/2)= gmn3(ii*(ii+1)/2)*two
      enddo
!
      call para_allreducer(gmn2,gmn1,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(gmn3,gmn2,nao*(nao+1)/2,mpi_comm)
      return
end

