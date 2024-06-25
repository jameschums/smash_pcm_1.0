!--------------------------------------------------------------------------------
  subroutine int1grys(egrad,fulldmtrx,exij,coij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
!--------------------------------------------------------------------------------
!
! Calculate derivative of 1-electron Coulomb integrals (j|Z/r|i) using Rys quadratures
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), locbfij(2), natom, nao, mxprsh
      integer :: nroots, ncart(0:6), ncarti, ncartj, i, j, k, iprim, jprim, iatom, iroot
      integer :: iang, jang, ix, iy, iz, jx, jy, jz
      real(8),parameter :: sqrtpi2=1.128379167095513D+00 !2.0/sqrt(pi)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, six=6.0D+00, eight=8.0D+00, p24=24.0D+00, eighth=0.125D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064590D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01
      real(8),parameter :: facf1=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facf2=3.87298334620741688D+00 ! sqrt(15)
      real(8),parameter :: facf3=0.61237243569579452D+00 ! sqrt(3/2)/2
      real(8),parameter :: facf4=1.93649167310370844D+00 ! sqrt(15)/2
      real(8),parameter :: facg1=2.95803989154980802D+00 ! sqrt(35)/2
      real(8),parameter :: facg2=2.09165006633518887D+00 ! sqrt(35/2)/2
      real(8),parameter :: facg3=1.11803398874989484D+00 ! sqrt(5)/2
      real(8),parameter :: facg4=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facg5=0.55901699437494742D+00 ! sqrt(5)/4
      real(8),parameter :: facg6=0.73950997288745200D+00 ! sqrt(35)/8
      real(8),intent(in) :: fulldmtrx(nao,nao), exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
      real(8),intent(in) :: coord(3,natom), znuc(natom), threshex
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rij, rij2, exi, exj, ci, cj, ex1, ex2, ex3, ex4, fac, rc, tval
      real(8) :: pixyz(3), pijxyz(3), pcxyz(3), xyzpijk(3,3), trys(13), wrys(13)
      real(8) :: cx(0:6,0:6,7,2), cy(0:6,0:6,7,2), cz(0:6,0:6,7,2), ww, xyzint(3)
      real(8) :: dcint1(28,28,3)
      logical,intent(in) :: iandj
      data ncart /1,3,6,10,15,21,28/
!
      nroots=(nangij(1)+nangij(2)+1)/2+1
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
!
      write(*,'("james Calculate derivative of 1-electron Coulomb integrals (j|Z/r|i) Rys quadratures")')
      write(*,'("james int1grys.F90 ")')
!
      do i= 1,3
        xyz(i)= coordij(i,1)-coordij(i,2)
      enddo
      rij= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
!
      do iatom= 1,natom
        dcint1(1:ncartj,1:ncarti,1:3)= zero
        do iprim= 1,nprimij(1)
          exi= exij(iprim,1)
          ci = coij(iprim,1)*sqrtpi2
          do i= 1,3
            pixyz(i)= exi*coordij(i,1)
          enddo
          do jprim= 1,nprimij(2)
            exj= exij(jprim,2)
            ex1= exi+exj
            ex2= one/ex1
            ex3= exi*exj
            rij2=rij*ex2*ex3
            if(rij2 > threshex) cycle
            cj = coij(jprim,2)
            fac=exp(-rij2)*ex2*ci*cj
!
            do i= 1,3
              pijxyz(i)=(pixyz(i)+exj*coordij(i,2))*ex2
            enddo
            do i= 1,3
              pcxyz(i)= pijxyz(i)-coord(i,iatom)
            enddo
            rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
            tval= ex1*rc
            call rysquad(tval,trys,wrys,nroots)
            do iroot= 1,nroots
              ww=-wrys(iroot)*znuc(iatom)*two*ex1*trys(iroot)/(one-trys(iroot))
              ex4= sqrt(ex2*(one-trys(iroot)))
              do i= 1,3
                xyzpijk(i,1)=(one-trys(iroot))*pijxyz(i)+trys(iroot)*coord(i,iatom)-coordij(i,1)
                xyzpijk(i,2)=(one-trys(iroot))*pijxyz(i)+trys(iroot)*coord(i,iatom)-coordij(i,2)
                xyzpijk(i,3)=(one-trys(iroot))*pijxyz(i)+trys(iroot)*coord(i,iatom)-coord(i,iatom)
              enddo
              do iang= 0,nangij(1)
                do jang= 0,nangij(2)
                  call ghquad(xyzint,ex4,xyzpijk,iang,jang)
                  cx(jang,iang,iroot,1)= xyzint(1)
                  cy(jang,iang,iroot,1)= xyzint(2)
                  cz(jang,iang,iroot,1)= xyzint(3)*ww
                  call dghquad(xyzint,ex4,xyzpijk,iang,jang)
                  cx(jang,iang,iroot,2)= xyzint(1)
                  cy(jang,iang,iroot,2)= xyzint(2)
                  cz(jang,iang,iroot,2)= xyzint(3)*ww
                enddo
              enddo
            enddo
            i= 0
            do ix= nangij(1),0,-1
              do iy= nangij(1)-ix,0,-1
                iz= nangij(1)-ix-iy
                i= i+1
                j= 0
                do jx= nangij(2),0,-1
                  do jy= nangij(2)-jx,0,-1
                    jz= nangij(2)-jx-jy
                    j= j+1
                    do iroot= 1,nroots
                      dcint1(j,i,1)= dcint1(j,i,1) &
&                                   +fac*cx(jx,ix,iroot,2)*cy(jy,iy,iroot,1)*cz(jz,iz,iroot,1)
                      dcint1(j,i,2)= dcint1(j,i,2) &
&                                   +fac*cx(jx,ix,iroot,1)*cy(jy,iy,iroot,2)*cz(jz,iz,iroot,1)
                      dcint1(j,i,3)= dcint1(j,i,3) &
&                                   +fac*cx(jx,ix,iroot,1)*cy(jy,iy,iroot,1)*cz(jz,iz,iroot,2)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!
        if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
          call nrmlz1(dcint1(1,1,1),nbfij(1),nbfij(2),ncarti)
          call nrmlz1(dcint1(1,1,2),nbfij(1),nbfij(2),ncarti)
          call nrmlz1(dcint1(1,1,3),nbfij(1),nbfij(2),ncarti)
        endif
!
        if(.not.iandj) then
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*dcint1(j,i,k)
              enddo
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom) &
&                              +fulldmtrx(locbfij(2)+j,locbfij(1)+i)*dcint1(j,i,k)*half
              enddo
            enddo
          enddo
        endif
!
      enddo
!
      return
end

