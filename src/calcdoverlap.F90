!-------------------------------------------------
  subroutine calcdoverlap(egrad,ewdmtrx,ish,jsh)
!-------------------------------------------------
!
! Driver of overlap derivative term
!
! In : ewdmtrx  (Energy-weighted density matrix)
!    : ish, jsh (Shell indices)
! Inout : egrad (Energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      use modmolecule, only : natom, coord
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj
      integer :: nbfi, nbfj, iprim, jprim, ncarti, ncartj, i, j, iang, jang, ii, ij
      integer :: ix, jx, iy, jy, iz, jz, ncart(0:6)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: ewdmtrx(nao*(nao+1)/2)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, exj2, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpij(3,2), cij
      real(8) :: xyzint(3), sx(0:7,0:6,2), sy(0:7,0:6,2), sz(0:7,0:6,2)
      real(8) :: dsint(28,28,3)
      data ncart/1,3,6,10,15,21,28/
!
      iatom = locatom(ish)
      jatom = locatom(jsh)
      if(iatom == jatom) return
!
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      nprimi= mprim(ish)
      nangi = mtype(ish)
      nbfi  = mbf(ish)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
      nprimj= mprim(jsh)
      nangj = mtype(jsh)
      nbfj  = mbf(jsh)
!
      if((nangi > 6).or.(nangj > 6)) then
        write(*,'(" Error! This program supports up to h function in calcdoverlap")')
        call iabort
      endif
      
!       write(*,'("james calcdoverlap.F90 Driver of overlap derivative term",/,&
! &               "In : ewdmtrx  (Energy-weighted density matrix) ish, jsh (Shell indices)",/,&
! &               "Inout : egrad (Energy gradient value)",/)')
!
      do i= 1,3
        xyzij(i)= coord(i,iatom)-coord(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      ncarti= ncart(nangi)
      ncartj= ncart(nangj)
!
      do i= 1,ncarti
        do j= 1,ncartj
          dsint(j,i,1)= zero
          dsint(j,i,2)= zero
          dsint(j,i,3)= zero
        enddo
      enddo
!
! Calculate overlap derivative for each primitive
!
      do iprim= 1,nprimi
        exi= ex(iloc+iprim)
        ci = coeff(iloc+iprim)
        do jprim= 1,nprimj
          exj= ex(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          exj2= exj*two
          ex3= sqrt(ex2)
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff(jloc+jprim)*fac
!
          do iang= 0,nangi
            do jang= 0,nangj+1
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang,1)= xyzint(1)*ex3
              sy(jang,iang,1)= xyzint(2)*ex3
              sz(jang,iang,1)= xyzint(3)*ex3
            enddo
          enddo
          do iang= 0,nangi
            sx(0,iang,2)= sx(1,iang,1)*exj2
            sy(0,iang,2)= sy(1,iang,1)*exj2
            sz(0,iang,2)= sz(1,iang,1)*exj2
          enddo
          do iang= 0,nangi
            do jang= 1,nangj
              sx(jang,iang,2)= sx(jang+1,iang,1)*exj2-sx(jang-1,iang,1)*jang
              sy(jang,iang,2)= sy(jang+1,iang,1)*exj2-sy(jang-1,iang,1)*jang
              sz(jang,iang,2)= sz(jang+1,iang,1)*exj2-sz(jang-1,iang,1)*jang
            enddo
          enddo
          cij= ci*cj
          i= 0
          do ix= nangi,0,-1
            do iy= nangi-ix,0,-1
              iz= nangi-ix-iy
              i= i+1
              j= 0
              do jx= nangj,0,-1
                do jy= nangj-jx,0,-1
                  jz= nangj-jx-jy
                  j= j+1
                  dsint(j,i,1)= dsint(j,i,1)+cij*sx(jx,ix,2)*sy(jy,iy,1)*sz(jz,iz,1)
                  dsint(j,i,2)= dsint(j,i,2)+cij*sx(jx,ix,1)*sy(jy,iy,2)*sz(jz,iz,1)
                  dsint(j,i,3)= dsint(j,i,3)+cij*sx(jx,ix,1)*sy(jy,iy,1)*sz(jz,iz,2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(dsint(1,1,1),nbfi,nbfj,ncarti)
        call nrmlz1(dsint(1,1,2),nbfi,nbfj,ncarti)
        call nrmlz1(dsint(1,1,3),nbfi,nbfj,ncarti)
      endif
!
      do i= 1,nbfi
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,nbfj
          egrad(1,iatom)= egrad(1,iatom)-ewdmtrx(ij+j)*dsint(j,i,1)
          egrad(2,iatom)= egrad(2,iatom)-ewdmtrx(ij+j)*dsint(j,i,2)
          egrad(3,iatom)= egrad(3,iatom)-ewdmtrx(ij+j)*dsint(j,i,3)
          egrad(1,jatom)= egrad(1,jatom)+ewdmtrx(ij+j)*dsint(j,i,1)
          egrad(2,jatom)= egrad(2,jatom)+ewdmtrx(ij+j)*dsint(j,i,2)
          egrad(3,jatom)= egrad(3,jatom)+ewdmtrx(ij+j)*dsint(j,i,3)
        enddo
      enddo
!
      return
end


