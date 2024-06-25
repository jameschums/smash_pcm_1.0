!---------------------------------------------------------------------------
  subroutine int1gcdp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,nbfij)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <d|V'|p>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2), nbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, five=5.0D+00, six=6.0D+00, seven=7.0D+00
      real(8),parameter :: p15=1.5D+00, p25=2.5D+00, p35=3.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:4), ft2two, ft3two, ft2inv
      real(8) :: pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r1(9), r2(18), r3(20), r4(15), xx, yy, zz, xy, xz, yz
      real(8) :: fpc2(6), fpc3(10), cint1(6,3,3), work(6)
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      do iatom= 1,natom
        do i= 1,9
          r1(i) = zero
        enddo
        do i= 1,18
          r2(i) = zero
        enddo
        do i= 1,20
          r3(i) = zero
        enddo
        do i= 1,15
          r4(i) = zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          ex21= exfac(2,ijprim)
          ex2i= exfac(3,ijprim)
          ex2j= exfac(4,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          pxx= pcxyz(1)*pcxyz(1)
          pyy= pcxyz(2)*pcxyz(2)
          pzz= pcxyz(3)*pcxyz(3)
          pxy= pcxyz(1)*pcxyz(2)
          pxz= pcxyz(1)*pcxyz(3)
          pyz= pcxyz(2)*pcxyz(3)
          rc= pxx+pyy+pzz
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
            ft(3)= p25 *tinv*ft(2)
            ft(4)= p35 *tinv*ft(3)
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
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,4
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          ex21 = ex21*half
          ft2two= ft(2)*extwo
          ft3two= ft(3)*extwo
          ft(4) = ft(4)*extwo
          ft2inv= ft(2)*ex21
          do i= 1,3
            r1(i)  = r1(i)  +ft(1)*pcxyz(i)*ex2j
            r1(i+3)= r1(i+3)-ft(1)*pcxyz(i)*ex2i
            r1(i+6)= r1(i+6)+ft(1)*pcxyz(i)*ex2j*ex2i*ex2i*extwo
          enddo
          fpc2(1)= ft2two*pxx-ft(1)
          fpc2(2)= ft2two*pyy-ft(1)
          fpc2(3)= ft2two*pzz-ft(1)
          fpc2(4)= ft2two*pxy
          fpc2(5)= ft2two*pxz
          fpc2(6)= ft2two*pyz
          do i= 1,6
            r2(i   )= r2(i   )-fpc2(i)*ex2j*ex2i
            r2(i+ 6)= r2(i+ 6)+fpc2(i)*ex2i*ex2i
            r2(i+12)= r2(i+12)+fpc2(i)*ex21
          enddo
          fpc3(1) =-ft3two*pxx*pcxyz(1)+ft(2)*pcxyz(1)*three
          fpc3(2) =-ft3two*pyy*pcxyz(2)+ft(2)*pcxyz(2)*three
          fpc3(3) =-ft3two*pzz*pcxyz(3)+ft(2)*pcxyz(3)*three
          fpc3(4) =-ft3two*pxx*pcxyz(2)+ft(2)*pcxyz(2)
          fpc3(5) =-ft3two*pxx*pcxyz(3)+ft(2)*pcxyz(3)
          fpc3(6) =-ft3two*pxy*pcxyz(2)+ft(2)*pcxyz(1)
          fpc3(7) =-ft3two*pyy*pcxyz(3)+ft(2)*pcxyz(3)
          fpc3(8) =-ft3two*pxz*pcxyz(3)+ft(2)*pcxyz(1)
          fpc3(9) =-ft3two*pyz*pcxyz(3)+ft(2)*pcxyz(2)
          fpc3(10)=-ft3two*pxy*pcxyz(3)
          do i= 1,10
            r3(i   )= r3(i   )-fpc3(i)*ex2j
            r3(i+10)= r3(i+10)+fpc3(i)*ex2i
          enddo
          r4( 1)= r4( 1)+ft(4)*pxx*pxx-ft(3)*pxx*six      +ft2inv*three
          r4( 2)= r4( 2)+ft(4)*pyy*pyy-ft(3)*pyy*six      +ft2inv*three
          r4( 3)= r4( 3)+ft(4)*pzz*pzz-ft(3)*pzz*six      +ft2inv*three
          r4( 4)= r4( 4)+ft(4)*pxx*pxy-ft(3)*pxy*three
          r4( 5)= r4( 5)+ft(4)*pxx*pxz-ft(3)*pxz*three
          r4( 6)= r4( 6)+ft(4)*pxy*pyy-ft(3)*pxy*three
          r4( 7)= r4( 7)+ft(4)*pyy*pyz-ft(3)*pyz*three
          r4( 8)= r4( 8)+ft(4)*pxz*pzz-ft(3)*pxz*three
          r4( 9)= r4( 9)+ft(4)*pyz*pzz-ft(3)*pyz*three
          r4(10)= r4(10)+ft(4)*pxx*pyy-ft(3)*pxx-ft(3)*pyy+ft2inv
          r4(11)= r4(11)+ft(4)*pxx*pzz-ft(3)*pxx-ft(3)*pzz+ft2inv
          r4(12)= r4(12)+ft(4)*pyy*pzz-ft(3)*pyy-ft(3)*pzz+ft2inv
          r4(13)= r4(13)+ft(4)*pxx*pyz-ft(3)*pyz
          r4(14)= r4(14)+ft(4)*pxy*pyz-ft(3)*pxz
          r4(15)= r4(15)+ft(4)*pxy*pzz-ft(3)*pxy
        enddo
        cint1(1,1,1)=r4( 1)+r3( 1)*xyz(1)+r3(11)*xyz(1)*two+r2( 1)*xx*two+r2( 7)*xx &
&                   +r2(13)*three+r1(1)*xyz(1)+r1(4)*xyz(1)*two+r1(7)*xx*xyz(1)
        cint1(2,1,1)=r4( 4)+r3( 4)*xyz(1)+r3(14)*xyz(1)+r3(11)*xyz(2)+r2( 4)*xx &
&                   +r2( 1)*xy+r2( 7)*xy+r2(16)+r1(4)*xyz(2)+r1(7)*xy*xyz(1)
        cint1(3,1,1)=r4( 5)+r3( 5)*xyz(1)+r3(15)*xyz(1)+r3(11)*xyz(3)+r2( 5)*xx &
&                   +r2( 1)*xz+r2( 7)*xz+r2(17)+r1(4)*xyz(3)+r1(7)*xz*xyz(1)
        cint1(4,1,1)=r4(10)+r3( 6)*xyz(1)+r3(14)*xyz(2)*two+r2( 4)*xy*two+r2( 7)*yy &
&                   +r2(13)+r1(1)*xyz(1)+r1(7)*yy*xyz(1)
        cint1(5,1,1)=r4(13)+r3(10)*xyz(1)+r3(15)*xyz(2)+r3(14)*xyz(3)+r2( 5)*xy &
&                   +r2( 4)*xz+r2( 7)*yz+r1(7)*yz*xyz(1)
        cint1(6,1,1)=r4(11)+r3( 8)*xyz(1)+r3(15)*xyz(3)*two+r2( 5)*xz*two+r2( 7)*zz &
&                   +r2(13)+r1(1)*xyz(1)+r1(7)*zz*xyz(1)
        cint1(1,2,1)=r4( 4)+r3( 1)*xyz(2)+r3(14)*xyz(1)*two+r2( 1)*xy*two+r2(10)*xx &
&                   +r2(16)+r1(1)*xyz(2)+r1(7)*xx*xyz(2)
        cint1(2,2,1)=r4(10)+r3( 4)*xyz(2)+r3(16)*xyz(1)+r3(14)*xyz(2)+r2( 4)*xy &
&                   +r2( 1)*yy+r2(10)*xy+r2(13)+r1(4)*xyz(1)+r1(7)*xy*xyz(2)
        cint1(3,2,1)=r4(13)+r3( 5)*xyz(2)+r3(20)*xyz(1)+r3(14)*xyz(3)+r2( 5)*xy &
&                   +r2( 1)*yz+r2(10)*xz+r1(7)*xz*xyz(2)
        cint1(4,2,1)=r4( 6)+r3( 6)*xyz(2)+r3(16)*xyz(2)*two+r2( 4)*yy*two+r2(10)*yy &
&                   +r2(16)*three+r1(1)*xyz(2)+r1(4)*xyz(2)*two+r1(7)*yy*xyz(2)
        cint1(5,2,1)=r4(14)+r3(10)*xyz(2)+r3(20)*xyz(2)+r3(16)*xyz(3)+r2( 5)*yy &
&                   +r2( 4)*yz+r2(10)*yz+r2(17)+r1(4)*xyz(3)+r1(7)*yz*xyz(2)
        cint1(6,2,1)=r4(15)+r3( 8)*xyz(2)+r3(20)*xyz(3)*two+r2( 5)*yz*two+r2(10)*zz &
&                   +r2(16)+r1(1)*xyz(2)+r1(7)*zz*xyz(2)
        cint1(1,3,1)=r4( 5)+r3( 1)*xyz(3)+r3(15)*xyz(1)*two+r2( 1)*xz*two+r2(11)*xx &
&                   +r2(17)+r1(1)*xyz(3)+r1(7)*xx*xyz(3)
        cint1(2,3,1)=r4(13)+r3( 4)*xyz(3)+r3(20)*xyz(1)+r3(15)*xyz(2)+r2( 4)*xz &
&                   +r2( 1)*yz+r2(11)*xy+r1(7)*xy*xyz(3)
        cint1(3,3,1)=r4(11)+r3( 5)*xyz(3)+r3(18)*xyz(1)+r3(15)*xyz(3)+r2( 5)*xz &
&                   +r2( 1)*zz+r2(11)*xz+r2(13)+r1(4)*xyz(1)+r1(7)*xz*xyz(3)
        cint1(4,3,1)=r4(14)+r3( 6)*xyz(3)+r3(20)*xyz(2)*two+r2( 4)*yz*two+r2(11)*yy &
&                   +r2(17)+r1(1)*xyz(3)+r1(7)*yy*xyz(3)
        cint1(5,3,1)=r4(15)+r3(10)*xyz(3)+r3(18)*xyz(2)+r3(20)*xyz(3)+r2( 5)*yz &
&                   +r2( 4)*zz+r2(11)*yz+r2(16)+r1(4)*xyz(2)+r1(7)*yz*xyz(3)
        cint1(6,3,1)=r4( 8)+r3( 8)*xyz(3)+r3(18)*xyz(3)*two+r2( 5)*zz*two+r2(11)*zz &
&                   +r2(17)*three+r1(1)*xyz(3)+r1(4)*xyz(3)*two+r1(7)*zz*xyz(3)
        cint1(1,1,2)=r4( 4)+r3( 4)*xyz(1)+r3(14)*xyz(1)*two+r2( 4)*xx*two+r2(10)*xx &
&                   +r2(16)*three+r1(2)*xyz(1)+r1(5)*xyz(1)*two+r1(8)*xx*xyz(1)
        cint1(2,1,2)=r4(10)+r3( 6)*xyz(1)+r3(16)*xyz(1)+r3(14)*xyz(2)+r2( 2)*xx &
&                   +r2( 4)*xy+r2(10)*xy+r2(14)+r1(5)*xyz(2)+r1(8)*xy*xyz(1)
        cint1(3,1,2)=r4(13)+r3(10)*xyz(1)+r3(20)*xyz(1)+r3(14)*xyz(3)+r2( 6)*xx &
&                   +r2( 4)*xz+r2(10)*xz+r2(18)+r1(5)*xyz(3)+r1(8)*xz*xyz(1)
        cint1(4,1,2)=r4( 6)+r3( 2)*xyz(1)+r3(16)*xyz(2)*two+r2( 2)*xy*two+r2(10)*yy &
&                   +r2(16)+r1(2)*xyz(1)+r1(8)*yy*xyz(1)
        cint1(5,1,2)=r4(14)+r3( 7)*xyz(1)+r3(20)*xyz(2)+r3(16)*xyz(3)+r2( 6)*xy &
&                   +r2( 2)*xz+r2(10)*yz+r1(8)*yz*xyz(1)
        cint1(6,1,2)=r4(15)+r3( 9)*xyz(1)+r3(20)*xyz(3)*two+r2( 6)*xz*two+r2(10)*zz &
&                   +r2(16)+r1(2)*xyz(1)+r1(8)*zz*xyz(1)
        cint1(1,2,2)=r4(10)+r3( 4)*xyz(2)+r3(16)*xyz(1)*two+r2( 4)*xy*two+r2( 8)*xx &
&                   +r2(14)+r1(2)*xyz(2)+r1(8)*xx*xyz(2)
        cint1(2,2,2)=r4( 6)+r3( 6)*xyz(2)+r3(12)*xyz(1)+r3(16)*xyz(2)+r2( 2)*xy &
&                   +r2( 4)*yy+r2( 8)*xy+r2(16)+r1(5)*xyz(1)+r1(8)*xy*xyz(2)
        cint1(3,2,2)=r4(14)+r3(10)*xyz(2)+r3(17)*xyz(1)+r3(16)*xyz(3)+r2( 6)*xy &
&                   +r2( 4)*yz+r2( 8)*xz+r1(8)*xz*xyz(2)
        cint1(4,2,2)=r4( 2)+r3( 2)*xyz(2)+r3(12)*xyz(2)*two+r2( 2)*yy*two+r2( 8)*yy &
&                   +r2(14)*three+r1(2)*xyz(2)+r1(5)*xyz(2)*two+r1(8)*yy*xyz(2)
        cint1(5,2,2)=r4( 7)+r3( 7)*xyz(2)+r3(17)*xyz(2)+r3(12)*xyz(3)+r2( 6)*yy &
&                   +r2( 2)*yz+r2( 8)*yz+r2(18)+r1(5)*xyz(3)+r1(8)*yz*xyz(2)
        cint1(6,2,2)=r4(12)+r3( 9)*xyz(2)+r3(17)*xyz(3)*two+r2( 6)*yz*two+r2( 8)*zz &
&                   +r2(14)+r1(2)*xyz(2)+r1(8)*zz*xyz(2)
        cint1(1,3,2)=r4(13)+r3( 4)*xyz(3)+r3(20)*xyz(1)*two+r2( 4)*xz*two+r2(12)*xx &
&                   +r2(18)+r1(2)*xyz(3)+r1(8)*xx*xyz(3)
        cint1(2,3,2)=r4(14)+r3( 6)*xyz(3)+r3(17)*xyz(1)+r3(20)*xyz(2)+r2( 2)*xz &
&                   +r2( 4)*yz+r2(12)*xy+r1(8)*xy*xyz(3)
        cint1(3,3,2)=r4(15)+r3(10)*xyz(3)+r3(19)*xyz(1)+r3(20)*xyz(3)+r2( 6)*xz &
&                   +r2( 4)*zz+r2(12)*xz+r2(16)+r1(5)*xyz(1)+r1(8)*xz*xyz(3)
        cint1(4,3,2)=r4( 7)+r3( 2)*xyz(3)+r3(17)*xyz(2)*two+r2( 2)*yz*two+r2(12)*yy &
&                   +r2(18)+r1(2)*xyz(3)+r1(8)*yy*xyz(3)
        cint1(5,3,2)=r4(12)+r3( 7)*xyz(3)+r3(19)*xyz(2)+r3(17)*xyz(3)+r2( 6)*yz &
&                   +r2( 2)*zz+r2(12)*yz+r2(14)+r1(5)*xyz(2)+r1(8)*yz*xyz(3)
        cint1(6,3,2)=r4( 9)+r3( 9)*xyz(3)+r3(19)*xyz(3)*two+r2( 6)*zz*two+r2(12)*zz &
&                   +r2(18)*three+r1(2)*xyz(3)+r1(5)*xyz(3)*two+r1(8)*zz*xyz(3)
        cint1(1,1,3)=r4( 5)+r3( 5)*xyz(1)+r3(15)*xyz(1)*two+r2( 5)*xx*two+r2(11)*xx &
&                   +r2(17)*three+r1(3)*xyz(1)+r1(6)*xyz(1)*two+r1(9)*xx*xyz(1)
        cint1(2,1,3)=r4(13)+r3(10)*xyz(1)+r3(20)*xyz(1)+r3(15)*xyz(2)+r2( 6)*xx &
&                   +r2( 5)*xy+r2(11)*xy+r2(18)+r1(6)*xyz(2)+r1(9)*xy*xyz(1)
        cint1(3,1,3)=r4(11)+r3( 8)*xyz(1)+r3(18)*xyz(1)+r3(15)*xyz(3)+r2( 3)*xx &
&                   +r2( 5)*xz+r2(11)*xz+r2(15)+r1(6)*xyz(3)+r1(9)*xz*xyz(1)
        cint1(4,1,3)=r4(14)+r3( 7)*xyz(1)+r3(20)*xyz(2)*two+r2( 6)*xy*two+r2(11)*yy &
&                   +r2(17)+r1(3)*xyz(1)+r1(9)*yy*xyz(1)
        cint1(5,1,3)=r4(15)+r3( 9)*xyz(1)+r3(18)*xyz(2)+r3(20)*xyz(3)+r2( 3)*xy &
&                   +r2( 6)*xz+r2(11)*yz+r1(9)*yz*xyz(1)
        cint1(6,1,3)=r4( 8)+r3( 3)*xyz(1)+r3(18)*xyz(3)*two+r2( 3)*xz*two+r2(11)*zz &
&                   +r2(17)+r1(3)*xyz(1)+r1(9)*zz*xyz(1)
        cint1(1,2,3)=r4(13)+r3( 5)*xyz(2)+r3(20)*xyz(1)*two+r2( 5)*xy*two+r2(12)*xx &
&                   +r2(18)+r1(3)*xyz(2)+r1(9)*xx*xyz(2)
        cint1(2,2,3)=r4(14)+r3(10)*xyz(2)+r3(17)*xyz(1)+r3(20)*xyz(2)+r2( 6)*xy &
&                   +r2( 5)*yy+r2(12)*xy+r2(17)+r1(6)*xyz(1)+r1(9)*xy*xyz(2)
        cint1(3,2,3)=r4(15)+r3( 8)*xyz(2)+r3(19)*xyz(1)+r3(20)*xyz(3)+r2( 3)*xy &
&                   +r2( 5)*yz+r2(12)*xz+r1(9)*xz*xyz(2)
        cint1(4,2,3)=r4( 7)+r3( 7)*xyz(2)+r3(17)*xyz(2)*two+r2( 6)*yy*two+r2(12)*yy &
&                   +r2(18)*three+r1(3)*xyz(2)+r1(6)*xyz(2)*two+r1(9)*yy*xyz(2)
        cint1(5,2,3)=r4(12)+r3( 9)*xyz(2)+r3(19)*xyz(2)+r3(17)*xyz(3)+r2( 3)*yy &
&                   +r2( 6)*yz+r2(12)*yz+r2(15)+r1(6)*xyz(3)+r1(9)*yz*xyz(2)
        cint1(6,2,3)=r4( 9)+r3( 3)*xyz(2)+r3(19)*xyz(3)*two+r2( 3)*yz*two+r2(12)*zz &
&                   +r2(18)+r1(3)*xyz(2)+r1(9)*zz*xyz(2)
        cint1(1,3,3)=r4(11)+r3( 5)*xyz(3)+r3(18)*xyz(1)*two+r2( 5)*xz*two+r2( 9)*xx &
&                   +r2(15)+r1(3)*xyz(3)+r1(9)*xx*xyz(3)
        cint1(2,3,3)=r4(15)+r3(10)*xyz(3)+r3(19)*xyz(1)+r3(18)*xyz(2)+r2( 6)*xz &
&                   +r2( 5)*yz+r2( 9)*xy+r1(9)*xy*xyz(3)
        cint1(3,3,3)=r4( 8)+r3( 8)*xyz(3)+r3(13)*xyz(1)+r3(18)*xyz(3)+r2( 3)*xz &
&                   +r2( 5)*zz+r2( 9)*xz+r2(17)+r1(6)*xyz(1)+r1(9)*xz*xyz(3)
        cint1(4,3,3)=r4(12)+r3( 7)*xyz(3)+r3(19)*xyz(2)*two+r2( 6)*yz*two+r2( 9)*yy &
&                   +r2(15)+r1(3)*xyz(3)+r1(9)*yy*xyz(3)
        cint1(5,3,3)=r4( 9)+r3( 9)*xyz(3)+r3(13)*xyz(2)+r3(19)*xyz(3)+r2( 3)*yz &
&                   +r2( 6)*zz+r2( 9)*yz+r2(18)+r1(6)*xyz(2)+r1(9)*yz*xyz(3)
        cint1(6,3,3)=r4( 3)+r3( 3)*xyz(3)+r3(13)*xyz(3)*two+r2( 3)*zz*two+r2( 9)*zz &
&                    +r2(15)*three+r1(3)*xyz(3)+r1(6)*xyz(3)*two+r1(9)*zz*xyz(3)
!
        if(nbfij(2) == 6) then
          do k= 1,3
            do i= 1,3
              cint1(2,i,k)= cint1(2,i,k)*sqrt3
              cint1(3,i,k)= cint1(3,i,k)*sqrt3
              cint1(5,i,k)= cint1(5,i,k)*sqrt3
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,3
              do j= 1,6
                work(j)= cint1(j,i,k)
              enddo
              cint1(1,i,k)= work(2)*sqrt3
              cint1(2,i,k)= work(5)*sqrt3
              cint1(3,i,k)=(work(6)*two-work(1)-work(4))*half
              cint1(4,i,k)= work(3)*sqrt3
              cint1(5,i,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        endif
        do k= 1,3
          do i= 1,3
            do j= 1,nbfij(2)
              egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)
            enddo
          enddo
        enddo
      enddo
!
      return
end


