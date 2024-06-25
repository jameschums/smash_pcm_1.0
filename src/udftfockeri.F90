!------------------------------------------------------------------------------------
  subroutine udftfockeri(focka,fockb,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,lsh,maxdim, &
&                        hfexchange)
!------------------------------------------------------------------------------------
!
! Form unrestricted DFT Fock matrix from two-electron intgrals
!
      use modbasis, only : nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2, jloc0, kloc0
      real(8),parameter :: half=0.5D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: twoeri(maxdim,maxdim,maxdim,maxdim), hfexchange
      real(8),intent(inout) :: focka(nao*(nao+1)/2), fockb(nao*(nao+1)/2)
      real(8) :: val, val2, val4
      logical :: ieqj, keql, ieqk, jeql, ikandjl, ijorkl
!
      nbfi = mbf(ish)
      nbfj = mbf(jsh)
      nbfk = mbf(ksh)
      nbfl = mbf(lsh)
      locbfi= locbf(ish)
      locbfj= locbf(jsh)
      locbfk= locbf(ksh)
      locbfl= locbf(lsh)
!
      ieqj= ish.eq.jsh
      keql= ksh.eq.lsh
      ieqk= ish.eq.ksh
      jeql= jsh.eq.lsh
      ikandjl= ieqk.and.jeql
      ijorkl= ieqj.or.keql
      jmax= nbfj
      lmax= nbfl
      ij= 0
      jloc0= locbfj*(locbfj-1)/2
      kloc0= locbfk*(locbfk-1)/2
!
      if(.not.ieqk)then
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            nij= iloc2+jloc
            do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
              nik= iloc2+kloc
              njk= jloc2+kloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(keql) lmax= k
              do l= 1,lmax
                val= twoeri(l,k,j,i)
                if(abs(val).lt.cutint2) cycle
                lloc= locbfl+l
                nkl= kloc2+lloc
                nil= iloc2+lloc
                njl= jloc2+lloc
                if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
                if(ijorkl) then
                  if(iloc.eq.jloc) val= val*half
                  if(kloc.eq.lloc) val= val*half
                endif
                val2= val*two*hfexchange
                val4= val*four
!
                focka(nij)= focka(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
                fockb(nij)= fockb(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
                focka(nkl)= focka(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
                fockb(nkl)= fockb(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
!
                focka(nik)= focka(nik)-val2*dmtrxa(njl)
                focka(nil)= focka(nil)-val2*dmtrxa(njk)
                focka(njk)= focka(njk)-val2*dmtrxa(nil)
                focka(njl)= focka(njl)-val2*dmtrxa(nik)
!
                fockb(nik)= fockb(nik)-val2*dmtrxb(njl)
                fockb(nil)= fockb(nil)-val2*dmtrxb(njk)
                fockb(njk)= fockb(njk)-val2*dmtrxb(nil)
                fockb(njl)= fockb(njl)-val2*dmtrxb(nik)
              enddo
            enddo
          enddo
        enddo
      else
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            ij= ij+1
            kl= 0
     kloop: do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
              if(keql) lmax= k
              do l= 1,lmax
                kl= kl+1
                if(ikandjl.and.kl.gt.ij) exit kloop
                val= twoeri(l,k,j,i)
                if(abs(val).lt.cutint2) cycle
                lloc= locbfl+l
                nij= iloc2+jloc
                nkl= kloc2+lloc
                nik= iloc2+kloc
                nil= iloc2+lloc
                njk= jloc2+kloc
                njl= jloc2+lloc
                if(jloc.lt.kloc) njk= kloc2+jloc
                if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
                if(iloc.lt.kloc) then
                  lloc2= lloc*(lloc-1)/2
                  nij= kloc2+lloc
                  nkl= iloc2+jloc
                  nik= kloc2+iloc
                  nil= kloc2+jloc
                  njk= lloc2+iloc
                  njl= lloc2+jloc
                  if(lloc.lt.iloc) njk= iloc2+lloc
                  if(lloc.lt.jloc) njl= jloc2+lloc
                elseif(iloc.eq.kloc.and.jloc.eq.lloc) then
                  val= val*half
                endif
                if(ijorkl) then
                  if(iloc.eq.jloc) val= val*half
                  if(kloc.eq.lloc) val= val*half
                endif
                val2= val*two*hfexchange
                val4= val*four
!
                focka(nij)= focka(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
                fockb(nij)= fockb(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
                focka(nkl)= focka(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
                fockb(nkl)= fockb(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
!
                focka(nik)= focka(nik)-val2*dmtrxa(njl)
                focka(nil)= focka(nil)-val2*dmtrxa(njk)
                focka(njk)= focka(njk)-val2*dmtrxa(nil)
                focka(njl)= focka(njl)-val2*dmtrxa(nik)
!
                fockb(nik)= fockb(nik)-val2*dmtrxb(njl)
                fockb(nil)= fockb(nil)-val2*dmtrxb(njk)
                fockb(njk)= fockb(njk)-val2*dmtrxb(nil)
                fockb(njl)= fockb(njl)-val2*dmtrxb(nik)
              enddo
            enddo kloop
          enddo
        enddo
      endif
!
      return
end


