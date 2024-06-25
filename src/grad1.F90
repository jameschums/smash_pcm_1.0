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

!-------------------------------------------------------------------------------
  subroutine int1cgmd(egrad,fulldmtrx,exij,cij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
!-------------------------------------------------------------------------------
!
! Driver of first derivative of 1-electron Coulomb integrals (j|(Z/r)'|i) 
! using McMurchie-Davidson method
!
! In : exij     (exponents of basis functions)
!      coij     (coefficients of basis functions)
!      coordij  (x,y,z coordinates of basis functions)
!      coord    (x,y,z coordinates of atoms)
!      znuc     (charges of atoms)
!      natom    (number of atoms)
!      nprimij  (numbers of primitive functions)
!      nangij   (degrees of angular momentum)
!      nbfij    (numbers of basis functions)
!      locbfij  (starting addresses of basis functions)
!      mxprsh   (size of primitive fuction array)
!      threshex (threshold of exponential calculation)
!      iandj    (flag of ish == jsh)
! Out: egrad    (energy gradient value)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), locbfij(2), natom, nao, mxprsh
      integer,parameter :: mxprsh2=30
      integer :: inttyp, nij, iprim, jprim, i, ii, jj, nbfij2(2), locbfij2(2)
      real(8),parameter :: one=1.0D+00, pi2=6.283185307179586D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exij(mxprsh,2), cij(mxprsh,2), coordij(3,2)
      real(8),intent(in) :: coord(3,natom), znuc(natom), threshex
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rij, exi, exj, ci, cj, ex12, ex21, ex2i, ex2j, rij2, pixyz(3)
      real(8) :: exfac(5,mxprsh2*mxprsh2), pijxyz(3,mxprsh2*mxprsh2)
      logical,intent(in) :: iandj
!
      if(mxprsh > mxprsh2) then
        write(*,'(" Error! Parameter mxprsh2 in int1cgmd is small!")')
        call exit
      endif
!
      inttyp=nangij(2)*3+nangij(1)+1
      if(nangij(2).ge.nangij(1)) then
        ii=1
        jj=2
        nbfij2(1)= nbfij(1)
        nbfij2(2)= nbfij(2)
        locbfij2(1)= locbfij(1)
        locbfij2(2)= locbfij(2)
      else
        ii=2
        jj=1
        nbfij2(1)= nbfij(2)
        nbfij2(2)= nbfij(1)
        locbfij2(1)= locbfij(2)
        locbfij2(2)= locbfij(1)
      endif
!
      do i= 1,3
        xyz(i)= coordij(i,ii)-coordij(i,jj)
      enddo
      rij= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
      nij= 0
      do iprim= 1,nprimij(ii)
        exi = exij(iprim,ii)
        ci  = cij(iprim,ii)*pi2
        do i= 1,3
          pixyz(i)= exi*coordij(i,ii)
        enddo
        do jprim= 1,nprimij(jj)
          exj = exij(jprim,jj)
          ex12= exi+exj
          ex21= one/ex12
          ex2i= ex21*exi
          ex2j= ex21*exj
          rij2= rij*ex2i*exj
          if(rij2 > threshex) cycle
          nij= nij+1
          cj = cij(jprim,jj)
          exfac(1,nij)= ex12
          exfac(2,nij)= ex21
          exfac(3,nij)= ex2i
          exfac(4,nij)= ex2j
          exfac(5,nij)= exp(-rij2)*ex21*ci*cj
          do i= 1,3
            pijxyz(i,nij)=(pixyz(i)+exj*coordij(i,jj))*ex21
          enddo
        enddo
      enddo
!
      select case(inttyp)
        case (1)
          call int1gcss(egrad,fulldmtrx,exfac,pijxyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,iandj)
        case (2,4)
          call int1gcps(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2)
        case (5)
          call int1gcpp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,iandj)
        case (3,7)
          call int1gcds(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,nbfij2)
        case (6,8)
          call int1gcdp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,nbfij2)
        case (9)
          call int1gcdd(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,nbfij2,iandj)
      end select
!
      return
end


!-----------------------------------------------------------------------
  subroutine int1gcss(egrad,fulldmtrx,exfac,pijxyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,iandj)
!-----------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <s|V'|s>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2)
      integer :: ijprim, i, k, iatom, igrid
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh)
      real(8),intent(in) :: pijxyz(3,mxprsh*mxprsh), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, c12, pcxyz(3), tinv
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:1), r1(3)
      logical,intent(in) :: iandj
!
      do iatom= 1,natom
        do i= 1,3
          r1(i)= zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= ft(0)*half*tinv
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            ft(0)= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&                 +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&                 +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&                 +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            ft(1)= fgrid(0,1,igrid)      +fgrid( 1,1,igrid)*tval  +fgrid( 2,1,igrid)*tval2 &
&                 +fgrid(3,1,igrid)*tval3+fgrid( 4,1,igrid)*tval4 +fgrid( 5,1,igrid)*tval5 &
&                 +fgrid(6,1,igrid)*tval6+fgrid( 7,1,igrid)*tval7 +fgrid( 8,1,igrid)*tval8 &
&                 +fgrid(9,1,igrid)*tval9+fgrid(10,1,igrid)*tval10
          endif
          ft(1)=-ft(1)*two*ex12*c12*znuc(iatom)
          do i= 1,3
            r1(i)= r1(i)+ft(1)*pcxyz(i)
          enddo
        enddo
        if(.not.iandj) then
          do k= 1,3
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+1,locbfij(1)+1)*r1(k)
          enddo
        else
          do k= 1,3
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+1,locbfij(1)+1)*r1(k)*half
          enddo
        endif
      enddo
!
      return
end





