!---------------------------------------------------------------------------
  subroutine int1gcdd(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,nbfij,iandj)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <d|V'|d>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2), nbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, seven=7.0D+00, nine=9.0D+00
      real(8),parameter :: ten=1.0D+01, fift=1.5D+01
      real(8),parameter :: p15=1.5D+00, p25=2.5D+00, p35=3.5D+00, p45=4.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo, ex2ii, ex2ij, ex2jj
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:5), ft2two, ft3two, ft4two
      real(8) :: ft2inv, ft3inv, pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: pxxx, pyyy, pzzz, pxxy, pxxz, pxyy, pyyz, pxzz, pyzz, pxyz
      real(8) :: tinv, r1(15), r2(24), r3(40), r4(30), r5(21), xx, yy, zz, xy, xz, yz
      real(8) :: xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz2
      real(8) :: fpc2(6), fpc3(10), fpc4(15), cint1(6,6,3), work(6)
      logical,intent(in) :: iandj
!
!      write(*,'("james int1gcdd.F90 line 31")')
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      xxx= xx*xyz(1)
      yyy= yy*xyz(2)
      zzz= zz*xyz(3)
      xxy= xx*xyz(2)
      xxz= xx*xyz(3)
      xyy= yy*xyz(1)
      yyz= yy*xyz(3)
      xzz= zz*xyz(1)
      yzz= zz*xyz(2)
      xyz2= xy*xyz(3)
      do iatom= 1,natom
        do i= 1,15
          r1(i) = zero
        enddo
        do i= 1,24
          r2(i) = zero
        enddo
        do i= 1,40
          r3(i) = zero
        enddo
        do i= 1,30
          r4(i) = zero
        enddo
        do i= 1,21
          r5(i) = zero
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
          pxxx= pxx*pcxyz(1)
          pyyy= pyy*pcxyz(2)
          pzzz= pzz*pcxyz(3)
          pxxy= pxx*pcxyz(2)
          pxxz= pxx*pcxyz(3)
          pxyy= pyy*pcxyz(1)
          pyyz= pyy*pcxyz(3)
          pxzz= pzz*pcxyz(1)
          pyzz= pzz*pcxyz(2)
          pxyz= pxy*pcxyz(3)
          rc= pxx+pyy+pzz
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
            ft(3)= p25 *tinv*ft(2)
            ft(4)= p35 *tinv*ft(3)
            ft(5)= p45 *tinv*ft(4)
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
            do ii= 0,5
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,5
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo = ex12*two
          ex21  = ex21*half
          ft2two= ft(2)*extwo
          ft3two= ft(3)*extwo
          ft4two= ft(4)*extwo
          ft(5) = ft(5)*extwo
          ft2inv= ft(2)*ex21
          ft3inv= ft(3)*ex21
          ex2ii= ex2i*ex2i
          ex2ij= ex2i*ex2j
          ex2jj= ex2j*ex2j
          do i= 1,3
            r1(i   )= r1(i   )-ft(1)*pcxyz(i)*ex2jj*ex2ii*extwo
            r1(i+ 3)= r1(i+ 3)-ft(1)*pcxyz(i)*ex2jj
            r1(i+ 6)= r1(i+ 6)+ft(1)*pcxyz(i)*ex2ij
            r1(i+ 9)= r1(i+ 9)-ft(1)*pcxyz(i)*ex2ii
            r1(i+12)= r1(i+12)-ft(1)*pcxyz(i)*ex21
          enddo
          fpc2(1)= ft2two*pxx-ft(1)
          fpc2(2)= ft2two*pyy-ft(1)
          fpc2(3)= ft2two*pzz-ft(1)
          fpc2(4)= ft2two*pxy
          fpc2(5)= ft2two*pxz
          fpc2(6)= ft2two*pyz
          do i= 1,6
            r2(i   )= r2(i   )+fpc2(i)*ex2jj*ex2i
            r2(i+ 6)= r2(i+ 6)-fpc2(i)*ex2ii*ex2j
            r2(i+12)= r2(i+12)-fpc2(i)*ex2j*ex21
            r2(i+18)= r2(i+18)+fpc2(i)*ex2i*ex21
          enddo
          fpc3(1) =-ft3two*pxxx+ft(2)*pcxyz(1)*three
          fpc3(2) =-ft3two*pyyy+ft(2)*pcxyz(2)*three
          fpc3(3) =-ft3two*pzzz+ft(2)*pcxyz(3)*three
          fpc3(4) =-ft3two*pxxy+ft(2)*pcxyz(2)
          fpc3(5) =-ft3two*pxxz+ft(2)*pcxyz(3)
          fpc3(6) =-ft3two*pxyy+ft(2)*pcxyz(1)
          fpc3(7) =-ft3two*pyyz+ft(2)*pcxyz(3)
          fpc3(8) =-ft3two*pxzz+ft(2)*pcxyz(1)
          fpc3(9) =-ft3two*pyzz+ft(2)*pcxyz(2)
          fpc3(10)=-ft3two*pxyz
          do i= 1,10
            r3(i   )= r3(i   )+fpc3(i)*ex2jj
            r3(i+10)= r3(i+10)-fpc3(i)*ex2ij
            r3(i+20)= r3(i+20)+fpc3(i)*ex2ii
            r3(i+30)= r3(i+30)+fpc3(i)*ex21
          enddo
          fpc4( 1)= ft4two*pxx*pxx-ft(3)*pxx*six      +ft2inv*three
          fpc4( 2)= ft4two*pyy*pyy-ft(3)*pyy*six      +ft2inv*three
          fpc4( 3)= ft4two*pzz*pzz-ft(3)*pzz*six      +ft2inv*three
          fpc4( 4)= ft4two*pxx*pxy-ft(3)*pxy*three
          fpc4( 5)= ft4two*pxx*pxz-ft(3)*pxz*three
          fpc4( 6)= ft4two*pxy*pyy-ft(3)*pxy*three
          fpc4( 7)= ft4two*pyy*pyz-ft(3)*pyz*three
          fpc4( 8)= ft4two*pxz*pzz-ft(3)*pxz*three
          fpc4( 9)= ft4two*pyz*pzz-ft(3)*pyz*three
          fpc4(10)= ft4two*pxx*pyy-ft(3)*pxx-ft(3)*pyy+ft2inv
          fpc4(11)= ft4two*pxx*pzz-ft(3)*pxx-ft(3)*pzz+ft2inv
          fpc4(12)= ft4two*pyy*pzz-ft(3)*pyy-ft(3)*pzz+ft2inv
          fpc4(13)= ft4two*pxx*pyz-ft(3)*pyz
          fpc4(14)= ft4two*pxy*pyz-ft(3)*pxz
          fpc4(15)= ft4two*pxy*pzz-ft(3)*pxy
          do i= 1,15
            r4(i   )= r4(i   )-fpc4(i)*ex2j
            r4(i+15)= r4(i+15)+fpc4(i)*ex2i
          enddo
          r5( 1)= r5( 1)-ft(5)*pxxx*pxx+ft(4)*pxxx*ten             -ft3inv*pcxyz(1)*fift
          r5( 2)= r5( 2)-ft(5)*pyyy*pyy+ft(4)*pyyy*ten             -ft3inv*pcxyz(2)*fift
          r5( 3)= r5( 3)-ft(5)*pzzz*pzz+ft(4)*pzzz*ten             -ft3inv*pcxyz(3)*fift
          r5( 4)= r5( 4)-ft(5)*pxxx*pxy+ft(4)*pxxy*six             -ft3inv*pcxyz(2)*three
          r5( 5)= r5( 5)-ft(5)*pxxx*pxz+ft(4)*pxxz*six             -ft3inv*pcxyz(3)*three
          r5( 6)= r5( 6)-ft(5)*pxyy*pyy+ft(4)*pxyy*six             -ft3inv*pcxyz(1)*three
          r5( 7)= r5( 7)-ft(5)*pyyy*pyz+ft(4)*pyyz*six             -ft3inv*pcxyz(3)*three
          r5( 8)= r5( 8)-ft(5)*pxzz*pzz+ft(4)*pxzz*six             -ft3inv*pcxyz(1)*three
          r5( 9)= r5( 9)-ft(5)*pyzz*pzz+ft(4)*pyzz*six             -ft3inv*pcxyz(2)*three
          r5(10)= r5(10)-ft(5)*pxxx*pyy+ft(4)*pxyy*three+ft(4)*pxxx-ft3inv*pcxyz(1)*three
          r5(11)= r5(11)-ft(5)*pxxx*pyz+ft(4)*pxyz*three
          r5(12)= r5(12)-ft(5)*pxxx*pzz+ft(4)*pxzz*three+ft(4)*pxxx-ft3inv*pcxyz(1)*three
          r5(13)= r5(13)-ft(5)*pxxy*pyy+ft(4)*pxxy*three+ft(4)*pyyy-ft3inv*pcxyz(2)*three
          r5(14)= r5(14)-ft(5)*pxyy*pyz+ft(4)*pxyz*three
          r5(15)= r5(15)-ft(5)*pyyy*pzz+ft(4)*pyzz*three+ft(4)*pyyy-ft3inv*pcxyz(2)*three
          r5(16)= r5(16)-ft(5)*pxxz*pzz+ft(4)*pxxz*three+ft(4)*pzzz-ft3inv*pcxyz(3)*three
          r5(17)= r5(17)-ft(5)*pxyz*pzz+ft(4)*pxyz*three
          r5(18)= r5(18)-ft(5)*pyyz*pzz+ft(4)*pyyz*three+ft(4)*pzzz-ft3inv*pcxyz(3)*three
          r5(19)= r5(19)-ft(5)*pxxy*pyz+ft(4)*pxxz      +ft(4)*pyyz-ft3inv*pcxyz(3)
          r5(20)= r5(20)-ft(5)*pxxy*pzz+ft(4)*pxxy      +ft(4)*pyzz-ft3inv*pcxyz(2)
          r5(21)= r5(21)-ft(5)*pxyy*pzz+ft(4)*pxyy      +ft(4)*pxzz-ft3inv*pcxyz(1)
        enddo
        cint1(1,1,1)=r5( 1)+r4( 1)*xyz(1)*two+r4(16)*xyz(1)*two+r3( 1)*xx+r3(11)*xx*four &
&                   +r3(21)*xx+r3(31)*six+r2( 1)*xxx*two+r2( 7)*xxx*two+r2(13)*xyz(1)*six &
&                   +r2(19)*xyz(1)*six+r1( 1)*xx*xx+r1( 4)*xx+r1( 7)*xx*four+r1(10)*xx &
&                   +r1(13)*three
        cint1(2,1,1)=r5( 4)+r4( 4)*xyz(1)*two+r4(19)*xyz(1)+r4(16)*xyz(2)+r3( 4)*xx &
&                   +r3(14)*xx*two+r3(11)*xy*two+r3(21)*xy+r3(34)*three+r2( 4)*xxx &
&                   +r2( 1)*xxy+r2( 7)*xxy*two+r2(16)*xyz(1)*two+r2(22)*xyz(1) &
&                   +r2(19)*xyz(2)*three+r1( 1)*xx*xy+r1( 7)*xy*two+r1(10)*xy
        cint1(3,1,1)=r5( 5)+r4( 5)*xyz(1)*two+r4(20)*xyz(1)+r4(16)*xyz(3)+r3( 5)*xx &
&                   +r3(15)*xx*two+r3(11)*xz*two+r3(21)*xz+r3(35)*three+r2( 5)*xxx &
&                   +r2( 1)*xxz+r2( 7)*xxz*two+r2(17)*xyz(1)*two+r2(23)*xyz(1) &
&                   +r2(19)*xyz(3)*three+r1( 1)*xx*xz+r1( 7)*xz*two+r1(10)*xz
        cint1(4,1,1)=r5(10)+r4(10)*xyz(1)*two+r4(19)*xyz(2)*two+r3( 6)*xx+r3(14)*xy*four &
&                   +r3(21)*yy+r3(31)+r3(36)+r2( 4)*xxy*two+r2( 7)*xyy*two+r2(13)*xyz(1)*two &
&                   +r2(22)*xyz(2)*two+r1( 1)*xx*yy+r1( 4)*xx+r1(10)*yy+r1(13)
        cint1(5,1,1)=r5(11)+r4(13)*xyz(1)*two+r4(20)*xyz(2)+r4(19)*xyz(3)+r3(10)*xx &
&                   +r3(15)*xy*two+r3(14)*xz*two+r3(21)*yz+r3(40)+r2( 5)*xxy+r2( 4)*xxz &
&                   +r2( 7)*xyz2*two+r2(23)*xyz(2)+r2(22)*xyz(3)+r1( 1)*xx*yz+r1(10)*yz
        cint1(6,1,1)=r5(12)+r4(11)*xyz(1)*two+r4(20)*xyz(3)*two+r3( 8)*xx+r3(15)*xz*four &
&                   +r3(21)*zz+r3(31)+r3(38)+r2( 5)*xxz*two+r2( 7)*xzz*two+r2(13)*xyz(1)*two &
&                   +r2(23)*xyz(3)*two+r1( 1)*xx*zz+r1( 4)*xx+r1(10)*zz+r1(13)
        cint1(1,2,1)=r5( 4)+r4( 4)*xyz(1)+r4( 1)*xyz(2)+r4(19)*xyz(1)*two+r3( 1)*xy &
&                   +r3(14)*xx*two+r3(11)*xy*two+r3(24)*xx+r3(34)*three+r2( 1)*xxy*two &
&                   +r2(10)*xxx+r2( 7)*xxy+r2(16)*xyz(1)+r2(13)*xyz(2)*three &
&                   +r2(22)*xyz(1)*two+r1( 1)*xy*xx+r1( 4)*xy+r1( 7)*xy*two
        cint1(2,2,1)=r5(10)+r4(10)*xyz(1)+r4( 4)*xyz(2)+r4(25)*xyz(1)+r4(19)*xyz(2) &
&                   +r3( 4)*xy+r3(16)*xx+r3(14)*xy+r3(14)*xy+r3(11)*yy+r3(24)*xy+r3(31) &
&                   +r3(36)+r2( 4)*xxy+r2( 1)*xyy+r2(10)*xxy+r2( 7)*xyy+r2(13)*xyz(1) &
&                   +r2(16)*xyz(2)+r2(19)*xyz(1)+r2(22)*xyz(2)+r1( 1)*xy*xy+r1( 7)*xx &
&                   +r1( 7)*yy+r1(13)
        cint1(3,2,1)=r5(11)+r4(13)*xyz(1)+r4( 5)*xyz(2)+r4(28)*xyz(1)+r4(19)*xyz(3) &
&                   +r3( 5)*xy+r3(20)*xx+r3(14)*xz+r3(15)*xy+r3(11)*yz+r3(24)*xz+r3(40) &
&                   +r2( 5)*xxy+r2( 1)*xyz2+r2(10)*xxz+r2( 7)*xyz2+r2(17)*xyz(2) &
&                   +r2(22)*xyz(3)+r1( 1)*xy*xz+r1( 7)*yz
        cint1(4,2,1)=r5(13)+r4( 6)*xyz(1)+r4(10)*xyz(2)+r4(25)*xyz(2)*two+r3( 6)*xy &
&                   +r3(16)*xy*two+r3(14)*yy*two+r3(24)*yy+r3(34)*three+r2( 4)*xyy*two &
&                   +r2(10)*xyy+r2( 7)*yyy+r2(16)*xyz(1)*three+r2(13)*xyz(2) &
&                   +r2(19)*xyz(2)*two+r1( 1)*xy*yy+r1( 4)*xy+r1( 7)*xy*two
        cint1(5,2,1)=r5(19)+r4(14)*xyz(1)+r4(13)*xyz(2)+r4(28)*xyz(2)+r4(25)*xyz(3) &
&                   +r3(10)*xy+r3(20)*xy+r3(16)*xz+r3(15)*yy+r3(14)*yz+r3(24)*yz+r3(35) &
&                   +r2( 5)*xyy+r2( 4)*xyz2+r2(10)*xyz2+r2( 7)*yyz+r2(17)*xyz(1) &
&                   +r2(19)*xyz(3)+r1( 1)*xy*yz+r1( 7)*xz
        cint1(6,2,1)=r5(20)+r4(15)*xyz(1)+r4(11)*xyz(2)+r4(28)*xyz(3)*two+r3( 8)*xy &
&                   +r3(20)*xz*two+r3(15)*yz*two+r3(24)*zz+r3(34)+r2( 5)*xyz2*two+r2(10)*xzz &
&                   +r2( 7)*yzz+r2(16)*xyz(1)+r2(13)*xyz(2)+r1( 1)*xy*zz+r1( 4)*xy
        cint1(1,3,1)=r5( 5)+r4( 5)*xyz(1)+r4( 1)*xyz(3)+r4(20)*xyz(1)*two+r3( 1)*xz &
&                   +r3(15)*xx*two+r3(11)*xz*two+r3(25)*xx+r3(35)*three+r2( 1)*xxz*two &
&                   +r2(11)*xxx+r2( 7)*xxz+r2(17)*xyz(1)+r2(13)*xyz(3)*three &
&                   +r2(23)*xyz(1)*two+r1( 1)*xz*xx+r1( 4)*xz+r1( 7)*xz*two
        cint1(2,3,1)=r5(11)+r4(13)*xyz(1)+r4( 4)*xyz(3)+r4(28)*xyz(1)+r4(20)*xyz(2) &
&                   +r3( 4)*xz+r3(20)*xx+r3(15)*xy+r3(14)*xz+r3(11)*yz+r3(25)*xy+r3(40) &
&                   +r2( 4)*xxz+r2( 1)*xyz2+r2(11)*xxy+r2( 7)*xyz2+r2(16)*xyz(3) &
&                   +r2(23)*xyz(2)+r1( 1)*xz*xy+r1( 7)*yz
        cint1(3,3,1)=r5(12)+r4(11)*xyz(1)+r4( 5)*xyz(3)+r4(26)*xyz(1)+r4(20)*xyz(3) &
&                   +r3( 5)*xz+r3(18)*xx+r3(15)*xz+r3(15)*xz+r3(11)*zz+r3(25)*xz+r3(31) &
&                   +r3(38)+r2( 5)*xxz+r2( 1)*xzz+r2(11)*xxz+r2( 7)*xzz+r2(13)*xyz(1) &
&                   +r2(17)*xyz(3)+r2(19)*xyz(1)+r2(23)*xyz(3)+r1( 1)*xz*xz+r1( 7)*xx &
&                   +r1( 7)*zz+r1(13)
        cint1(4,3,1)=r5(19)+r4(14)*xyz(1)+r4(10)*xyz(3)+r4(28)*xyz(2)*two+r3( 6)*xz &
&                   +r3(20)*xy*two+r3(14)*yz*two+r3(25)*yy+r3(35)+r2( 4)*xyz2*two+r2(11)*xyy &
&                   +r2( 7)*yyz+r2(17)*xyz(1)+r2(13)*xyz(3)+r1( 1)*xz*yy+r1( 4)*xz
        cint1(5,3,1)=r5(20)+r4(15)*xyz(1)+r4(13)*xyz(3)+r4(26)*xyz(2)+r4(28)*xyz(3) &
&                   +r3(10)*xz+r3(18)*xy+r3(20)*xz+r3(15)*yz+r3(14)*zz+r3(25)*yz+r3(34) &
&                   +r2( 5)*xyz2+r2( 4)*xzz+r2(11)*xyz2+r2( 7)*yzz+r2(16)*xyz(1) &
&                   +r2(19)*xyz(2)+r1( 1)*xz*yz+r1( 7)*xy
        cint1(6,3,1)=r5(16)+r4( 8)*xyz(1)+r4(11)*xyz(3)+r4(26)*xyz(3)*two+r3( 8)*xz &
&                   +r3(18)*xz*two+r3(15)*zz*two+r3(25)*zz+r3(35)*three+r2( 5)*xzz*two &
&                   +r2(11)*xzz+r2( 7)*zzz+r2(17)*xyz(1)*three+r2(13)*xyz(3) &
&                   +r2(19)*xyz(3)*two+r1( 1)*xz*zz+r1( 4)*xz+r1( 7)*xz*two
        cint1(1,4,1)=r5(10)+r4( 4)*xyz(2)*two+r4(25)*xyz(1)*two+r3( 1)*yy+r3(14)*xy*four &
&                   +r3(26)*xx+r3(31)+r3(36)+r2( 1)*xyy*two+r2(10)*xxy*two+r2(16)*xyz(2)*two &
&                   +r2(19)*xyz(1)*two+r1( 1)*yy*xx+r1( 4)*yy+r1(10)*xx+r1(13)
        cint1(2,4,1)=r5(13)+r4(10)*xyz(2)*two+r4(21)*xyz(1)+r4(25)*xyz(2)+r3( 4)*yy &
&                   +r3(16)*xy*two+r3(14)*yy*two+r3(26)*xy+r3(34)*three+r2( 4)*xyy &
&                   +r2( 1)*yyy+r2(10)*xyy*two+r2(13)*xyz(2)*two+r2(22)*xyz(1)*three &
&                   +r2(19)*xyz(2)+r1( 1)*yy*xy+r1( 7)*xy*two+r1(10)*xy
        cint1(3,4,1)=r5(19)+r4(13)*xyz(2)*two+r4(29)*xyz(1)+r4(25)*xyz(3)+r3( 5)*yy &
&                   +r3(20)*xy*two+r3(14)*yz*two+r3(26)*xz+r3(35)+r2( 5)*xyy+r2( 1)*yyz &
&                   +r2(10)*xyz2*two+r2(23)*xyz(1)+r2(19)*xyz(3)+r1( 1)*yy*xz+r1(10)*xz
        cint1(4,4,1)=r5( 6)+r4( 6)*xyz(2)*two+r4(21)*xyz(2)*two+r3( 6)*yy+r3(16)*yy*four &
&                   +r3(26)*yy+r3(36)*six+r2( 4)*yyy*two+r2(10)*yyy*two+r2(16)*xyz(2)*six &
&                   +r2(22)*xyz(2)*six+r1( 1)*yy*yy+r1( 4)*yy+r1( 7)*yy*four+r1(10)*yy &
&                   +r1(13)*three
        cint1(5,4,1)=r5(14)+r4(14)*xyz(2)*two+r4(29)*xyz(2)+r4(21)*xyz(3)+r3(10)*yy &
&                   +r3(20)*yy*two+r3(16)*yz*two+r3(26)*yz+r3(40)*three+r2( 5)*yyy &
&                   +r2( 4)*yyz+r2(10)*yyz*two+r2(17)*xyz(2)*two+r2(23)*xyz(2) &
&                   +r2(22)*xyz(3)*three+r1( 1)*yy*yz+r1( 7)*yz*two+r1(10)*yz
        cint1(6,4,1)=r5(21)+r4(15)*xyz(2)*two+r4(29)*xyz(3)*two+r3( 8)*yy+r3(20)*yz*four &
&                   +r3(26)*zz+r3(36)+r3(38)+r2( 5)*yyz*two+r2(10)*yzz*two+r2(16)*xyz(2)*two &
&                   +r2(23)*xyz(3)*two+r1( 1)*yy*zz+r1( 4)*yy+r1(10)*zz+r1(13)
        cint1(1,5,1)=r5(11)+r4( 5)*xyz(2)+r4( 4)*xyz(3)+r4(28)*xyz(1)*two+r3( 1)*yz &
&                   +r3(15)*xy*two+r3(14)*xz*two+r3(30)*xx+r3(40)+r2( 1)*xyz2*two &
&                   +r2(11)*xxy+r2(10)*xxz+r2(17)*xyz(2)+r2(16)*xyz(3)+r1( 1)*yz*xx+r1( 4)*yz
        cint1(2,5,1)=r5(19)+r4(13)*xyz(2)+r4(10)*xyz(3)+r4(29)*xyz(1)+r4(28)*xyz(2) &
&                   +r3( 4)*yz+r3(20)*xy+r3(15)*yy+r3(16)*xz+r3(14)*yz+r3(30)*xy+r3(35) &
&                   +r2( 4)*xyz2+r2( 1)*yyz+r2(11)*xyy+r2(10)*xyz2+r2(13)*xyz(3) &
&                   +r2(23)*xyz(1)+r1( 1)*yz*xy+r1( 7)*xz
        cint1(3,5,1)=r5(20)+r4(11)*xyz(2)+r4(13)*xyz(3)+r4(30)*xyz(1)+r4(28)*xyz(3) &
&                   +r3( 5)*yz+r3(18)*xy+r3(15)*yz+r3(20)*xz+r3(14)*zz+r3(30)*xz+r3(34) &
&                   +r2( 5)*xyz2+r2( 1)*yzz+r2(11)*xyz2+r2(10)*xzz+r2(13)*xyz(2) &
&                   +r2(22)*xyz(1)+r1( 1)*yz*xz+r1( 7)*xy
        cint1(4,5,1)=r5(14)+r4(14)*xyz(2)+r4( 6)*xyz(3)+r4(29)*xyz(2)*two+r3( 6)*yz &
&                   +r3(20)*yy*two+r3(16)*yz*two+r3(30)*yy+r3(40)*three+r2( 4)*yyz*two &
&                   +r2(11)*yyy+r2(10)*yyz+r2(17)*xyz(2)+r2(16)*xyz(3)*three+ &
&                   r2(23)*xyz(2)*two+r1( 1)*yz*yy+r1( 4)*yz+r1( 7)*yz*two
        cint1(5,5,1)=r5(21)+r4(15)*xyz(2)+r4(14)*xyz(3)+r4(30)*xyz(2)+r4(29)*xyz(3) &
&                   +r3(10)*yz+r3(18)*yy+r3(20)*yz+r3(20)*yz+r3(16)*zz+r3(30)*yz+r3(36) &
&                   +r3(38)+r2( 5)*yyz+r2( 4)*yzz+r2(11)*yyz+r2(10)*yzz+r2(16)*xyz(2) &
&                   +r2(17)*xyz(3)+r2(22)*xyz(2)+r2(23)*xyz(3)+r1( 1)*yz*yz+r1( 7)*yy &
&                   +r1( 7)*zz+r1(13)
        cint1(6,5,1)=r5(17)+r4( 8)*xyz(2)+r4(15)*xyz(3)+r4(30)*xyz(3)*two+r3( 8)*yz &
&                   +r3(18)*yz*two+r3(20)*zz*two+r3(30)*zz+r3(40)*three+r2( 5)*yzz*two &
&                   +r2(11)*yzz+r2(10)*zzz+r2(17)*xyz(2)*three+r2(16)*xyz(3) &
&                   +r2(22)*xyz(3)*two+r1( 1)*yz*zz+r1( 4)*yz+r1( 7)*yz*two
        cint1(1,6,1)=r5(12)+r4( 5)*xyz(3)*two+r4(26)*xyz(1)*two+r3( 1)*zz+r3(15)*xz*four &
&                   +r3(28)*xx+r3(31)+r3(38)+r2( 1)*xzz*two+r2(11)*xxz*two+r2(17)*xyz(3)*two &
&                   +r2(19)*xyz(1)*two+r1( 1)*zz*xx+r1( 4)*zz+r1(10)*xx+r1(13)
        cint1(2,6,1)=r5(20)+r4(13)*xyz(3)*two+r4(30)*xyz(1)+r4(26)*xyz(2)+r3( 4)*zz &
&                   +r3(20)*xz*two+r3(15)*yz*two+r3(28)*xy+r3(34)+r2( 4)*xzz+r2( 1)*yzz &
&                   +r2(11)*xyz2*two+r2(22)*xyz(1)+r2(19)*xyz(2)+r1( 1)*zz*xy+r1(10)*xy
        cint1(3,6,1)=r5(16)+r4(11)*xyz(3)*two+r4(23)*xyz(1)+r4(26)*xyz(3)+r3( 5)*zz &
&                   +r3(18)*xz*two+r3(15)*zz*two+r3(28)*xz+r3(35)*three+r2( 5)*xzz &
&                   +r2( 1)*zzz+r2(11)*xzz*two+r2(13)*xyz(3)*two+r2(23)*xyz(1)*three &
&                   +r2(19)*xyz(3)+r1( 1)*zz*xz+r1( 7)*xz*two+r1(10)*xz
        cint1(4,6,1)=r5(21)+r4(14)*xyz(3)*two+r4(30)*xyz(2)*two+r3( 6)*zz+r3(20)*yz*four &
&                   +r3(28)*yy+r3(36)+r3(38)+r2( 4)*yzz*two+r2(11)*yyz*two+r2(17)*xyz(3)*two &
&                   +r2(22)*xyz(2)*two+r1( 1)*zz*yy+r1( 4)*zz+r1(10)*yy+r1(13)
        cint1(5,6,1)=r5(17)+r4(15)*xyz(3)*two+r4(23)*xyz(2)+r4(30)*xyz(3)+r3(10)*zz &
&                   +r3(18)*yz*two+r3(20)*zz*two+r3(28)*yz+r3(40)*three+r2( 5)*yzz &
&                   +r2( 4)*zzz+r2(11)*yzz*two+r2(16)*xyz(3)*two+r2(23)*xyz(2)*three &
&                   +r2(22)*xyz(3)+r1( 1)*zz*yz+r1( 7)*yz*two+r1(10)*yz
        cint1(6,6,1)=r5( 8)+r4( 8)*xyz(3)*two+r4(23)*xyz(3)*two+r3( 8)*zz+r3(18)*zz*four &
&                   +r3(28)*zz+r3(38)*six+r2( 5)*zzz*two+r2(11)*zzz*two+r2(17)*xyz(3)*six &
&                   +r2(23)*xyz(3)*six+r1( 1)*zz*zz+r1( 4)*zz+r1( 7)*zz*four+r1(10)*zz &
&                   +r1(13)*three
        cint1(1,1,2)=r5( 4)+r4( 4)*xyz(1)*two+r4(19)*xyz(1)*two+r3( 4)*xx+r3(14)*xx*four &
&                   +r3(24)*xx+r3(34)*six+r2( 4)*xxx*two+r2(10)*xxx*two+r2(16)*xyz(1)*six &
&                   +r2(22)*xyz(1)*six+r1( 2)*xx*xx+r1( 5)*xx+r1( 8)*xx*four+r1(11)*xx &
&                   +r1(14)*three
        cint1(2,1,2)=r5(10)+r4(10)*xyz(1)*two+r4(25)*xyz(1)+r4(19)*xyz(2)+r3( 6)*xx &
&                   +r3(16)*xx*two+r3(14)*xy*two+r3(24)*xy+r3(36)*three+r2( 2)*xxx &
&                   +r2( 4)*xxy+r2(10)*xxy*two+r2(14)*xyz(1)*two+r2(20)*xyz(1) &
&                   +r2(22)*xyz(2)*three+r1( 2)*xx*xy+r1( 8)*xy*two+r1(11)*xy
        cint1(3,1,2)=r5(11)+r4(13)*xyz(1)*two+r4(28)*xyz(1)+r4(19)*xyz(3)+r3(10)*xx &
&                   +r3(20)*xx*two+r3(14)*xz*two+r3(24)*xz+r3(40)*three+r2( 6)*xxx &
&                   +r2( 4)*xxz+r2(10)*xxz*two+r2(18)*xyz(1)*two+r2(24)*xyz(1) &
&                   +r2(22)*xyz(3)*three+r1( 2)*xx*xz+r1( 8)*xz*two+r1(11)*xz
        cint1(4,1,2)=r5(13)+r4( 6)*xyz(1)*two+r4(25)*xyz(2)*two+r3( 2)*xx+r3(16)*xy*four &
&                   +r3(24)*yy+r3(32)+r3(34)+r2( 2)*xxy*two+r2(10)*xyy*two+r2(16)*xyz(1)*two &
&                   +r2(20)*xyz(2)*two+r1( 2)*xx*yy+r1( 5)*xx+r1(11)*yy+r1(14)
        cint1(5,1,2)=r5(19)+r4(14)*xyz(1)*two+r4(28)*xyz(2)+r4(25)*xyz(3)+r3( 7)*xx &
&                   +r3(20)*xy*two+r3(16)*xz*two+r3(24)*yz+r3(37)+r2( 6)*xxy+r2( 2)*xxz &
&                   +r2(10)*xyz2*two+r2(24)*xyz(2)+r2(20)*xyz(3)+r1( 2)*xx*yz+r1(11)*yz
        cint1(6,1,2)=r5(20)+r4(15)*xyz(1)*two+r4(28)*xyz(3)*two+r3( 9)*xx+r3(20)*xz*four &
&                   +r3(24)*zz+r3(34)+r3(39)+r2( 6)*xxz*two+r2(10)*xzz*two+r2(16)*xyz(1)*two &
&                   +r2(24)*xyz(3)*two+r1( 2)*xx*zz+r1( 5)*xx+r1(11)*zz+r1(14)
        cint1(1,2,2)=r5(10)+r4(10)*xyz(1)+r4( 4)*xyz(2)+r4(25)*xyz(1)*two+r3( 4)*xy &
&                   +r3(16)*xx*two+r3(14)*xy*two+r3(26)*xx+r3(36)*three+r2( 4)*xxy*two &
&                   +r2( 8)*xxx+r2(10)*xxy+r2(14)*xyz(1)+r2(16)*xyz(2)*three &
&                   +r2(20)*xyz(1)*two+r1( 2)*xy*xx+r1( 5)*xy+r1( 8)*xy*two
        cint1(2,2,2)=r5(13)+r4( 6)*xyz(1)+r4(10)*xyz(2)+r4(21)*xyz(1)+r4(25)*xyz(2) &
&                   +r3( 6)*xy+r3(12)*xx+r3(16)*xy+r3(16)*xy+r3(14)*yy+r3(26)*xy+r3(32) &
&                   +r3(34)+r2( 2)*xxy+r2( 4)*xyy+r2( 8)*xxy+r2(10)*xyy+r2(16)*xyz(1) &
&                   +r2(14)*xyz(2)+r2(22)*xyz(1)+r2(20)*xyz(2)+r1( 2)*xy*xy+r1( 8)*xx &
&                   +r1( 8)*yy+r1(14)
        cint1(3,2,2)=r5(19)+r4(14)*xyz(1)+r4(13)*xyz(2)+r4(29)*xyz(1)+r4(25)*xyz(3) &
&                   +r3(10)*xy+r3(17)*xx+r3(16)*xz+r3(20)*xy+r3(14)*yz+r3(26)*xz+r3(37) &
&                   +r2( 6)*xxy+r2( 4)*xyz2+r2( 8)*xxz+r2(10)*xyz2+r2(18)*xyz(2)+r2(20)*xyz(3) &
&                   +r1( 2)*xy*xz+r1( 8)*yz
        cint1(4,2,2)=r5( 6)+r4( 2)*xyz(1)+r4( 6)*xyz(2)+r4(21)*xyz(2)*two+r3( 2)*xy &
&                   +r3(12)*xy*two+r3(16)*yy*two+r3(26)*yy+r3(36)*three+r2( 2)*xyy*two &
&                   +r2( 8)*xyy+r2(10)*yyy+r2(14)*xyz(1)*three+r2(16)*xyz(2) &
&                   +r2(22)*xyz(2)*two+r1( 2)*xy*yy+r1( 5)*xy+r1( 8)*xy*two
        cint1(5,2,2)=r5(14)+r4( 7)*xyz(1)+r4(14)*xyz(2)+r4(29)*xyz(2)+r4(21)*xyz(3) &
&                   +r3( 7)*xy+r3(17)*xy+r3(12)*xz+r3(20)*yy+r3(16)*yz+r3(26)*yz+r3(40) &
&                   +r2( 6)*xyy+r2( 2)*xyz2+r2( 8)*xyz2+r2(10)*yyz+r2(18)*xyz(1)+r2(22)*xyz(3) &
&                   +r1( 2)*xy*yz+r1( 8)*xz
        cint1(6,2,2)=r5(21)+r4(12)*xyz(1)+r4(15)*xyz(2)+r4(29)*xyz(3)*two+r3( 9)*xy &
&                   +r3(17)*xz*two+r3(20)*yz*two+r3(26)*zz+r3(36)+r2( 6)*xyz2*two+r2( 8)*xzz &
&                   +r2(10)*yzz+r2(14)*xyz(1)+r2(16)*xyz(2)+r1( 2)*xy*zz+r1( 5)*xy
        cint1(1,3,2)=r5(11)+r4(13)*xyz(1)+r4( 4)*xyz(3)+r4(28)*xyz(1)*two+r3( 4)*xz &
&                   +r3(20)*xx*two+r3(14)*xz*two+r3(30)*xx+r3(40)*three+r2( 4)*xxz*two &
&                   +r2(12)*xxx+r2(10)*xxz+r2(18)*xyz(1)+r2(16)*xyz(3)*three &
&                   +r2(24)*xyz(1)*two+r1( 2)*xz*xx+r1( 5)*xz+r1( 8)*xz*two
        cint1(2,3,2)=r5(19)+r4(14)*xyz(1)+r4(10)*xyz(3)+r4(29)*xyz(1)+r4(28)*xyz(2) &
&                   +r3( 6)*xz+r3(17)*xx+r3(20)*xy+r3(16)*xz+r3(14)*yz+r3(30)*xy+r3(37) &
&                   +r2( 2)*xxz+r2( 4)*xyz2+r2(12)*xxy+r2(10)*xyz2+r2(14)*xyz(3) &
&                   +r2(24)*xyz(2)+r1( 2)*xz*xy+r1( 8)*yz
        cint1(3,3,2)=r5(20)+r4(15)*xyz(1)+r4(13)*xyz(3)+r4(30)*xyz(1)+r4(28)*xyz(3) &
&                   +r3(10)*xz+r3(19)*xx+r3(20)*xz+r3(20)*xz+r3(14)*zz+r3(30)*xz+r3(34) &
&                   +r3(39)+r2( 6)*xxz+r2( 4)*xzz+r2(12)*xxz+r2(10)*xzz+r2(16)*xyz(1) &
&                   +r2(18)*xyz(3)+r2(22)*xyz(1)+r2(24)*xyz(3)+r1( 2)*xz*xz+r1( 8)*xx &
&                   +r1( 8)*zz+r1(14)
        cint1(4,3,2)=r5(14)+r4( 7)*xyz(1)+r4( 6)*xyz(3)+r4(29)*xyz(2)*two+r3( 2)*xz &
&                   +r3(17)*xy*two+r3(16)*yz*two+r3(30)*yy+r3(40)+r2( 2)*xyz2*two+r2(12)*xyy &
&                   +r2(10)*yyz+r2(18)*xyz(1)+r2(16)*xyz(3)+r1( 2)*xz*yy+r1( 5)*xz
        cint1(5,3,2)=r5(21)+r4(12)*xyz(1)+r4(14)*xyz(3)+r4(30)*xyz(2)+r4(29)*xyz(3) &
&                   +r3( 7)*xz+r3(19)*xy+r3(17)*xz+r3(20)*yz+r3(16)*zz+r3(30)*yz+r3(36) &
&                   +r2( 6)*xyz2+r2( 2)*xzz+r2(12)*xyz2+r2(10)*yzz+r2(14)*xyz(1) &
&                   +r2(22)*xyz(2)+r1( 2)*xz*yz+r1( 8)*xy
        cint1(6,3,2)=r5(17)+r4( 9)*xyz(1)+r4(15)*xyz(3)+r4(30)*xyz(3)*two+r3( 9)*xz &
&                   +r3(19)*xz*two+r3(20)*zz*two+r3(30)*zz+r3(40)*three+r2( 6)*xzz*two &
&                   +r2(12)*xzz+r2(10)*zzz+r2(18)*xyz(1)*three+r2(16)*xyz(3) &
&                   +r2(22)*xyz(3)*two+r1( 2)*xz*zz+r1( 5)*xz+r1( 8)*xz*two
        cint1(1,4,2)=r5(13)+r4(10)*xyz(2)*two+r4(21)*xyz(1)*two+r3( 4)*yy+r3(16)*xy*four &
&                   +r3(22)*xx+r3(32)+r3(34)+r2( 4)*xyy*two+r2( 8)*xxy*two+r2(14)*xyz(2)*two &
&                   +r2(22)*xyz(1)*two+r1( 2)*yy*xx+r1( 5)*yy+r1(11)*xx+r1(14)
        cint1(2,4,2)=r5( 6)+r4( 6)*xyz(2)*two+r4(17)*xyz(1)+r4(21)*xyz(2)+r3( 6)*yy &
&                   +r3(12)*xy*two+r3(16)*yy*two+r3(22)*xy+r3(36)*three+r2( 2)*xyy &
&                   +r2( 4)*yyy+r2( 8)*xyy*two+r2(16)*xyz(2)*two+r2(20)*xyz(1)*three &
&                   +r2(22)*xyz(2)+r1( 2)*yy*xy+r1( 8)*xy*two+r1(11)*xy
        cint1(3,4,2)=r5(14)+r4(14)*xyz(2)*two+r4(22)*xyz(1)+r4(21)*xyz(3)+r3(10)*yy &
&                   +r3(17)*xy*two+r3(16)*yz*two+r3(22)*xz+r3(40)+r2( 6)*xyy+r2( 4)*yyz &
&                   +r2( 8)*xyz2*two+r2(24)*xyz(1)+r2(22)*xyz(3)+r1( 2)*yy*xz+r1(11)*xz
        cint1(4,4,2)=r5( 2)+r4( 2)*xyz(2)*two+r4(17)*xyz(2)*two+r3( 2)*yy+r3(12)*yy*four &
&                   +r3(22)*yy+r3(32)*six+r2( 2)*yyy*two+r2( 8)*yyy*two+r2(14)*xyz(2)*six &
&                   +r2(20)*xyz(2)*six+r1( 2)*yy*yy+r1( 5)*yy+r1( 8)*yy*four+r1(11)*yy &
&                   +r1(14)*three
        cint1(5,4,2)=r5( 7)+r4( 7)*xyz(2)*two+r4(22)*xyz(2)+r4(17)*xyz(3)+r3( 7)*yy &
&                   +r3(17)*yy*two+r3(12)*yz*two+r3(22)*yz+r3(37)*three+r2( 6)*yyy &
&                   +r2( 2)*yyz+r2( 8)*yyz*two+r2(18)*xyz(2)*two+r2(24)*xyz(2) &
&                   +r2(20)*xyz(3)*three+r1( 2)*yy*yz+r1( 8)*yz*two+r1(11)*yz
        cint1(6,4,2)=r5(15)+r4(12)*xyz(2)*two+r4(22)*xyz(3)*two+r3( 9)*yy+r3(17)*yz*four &
&                   +r3(22)*zz+r3(32)+r3(39)+r2( 6)*yyz*two+r2( 8)*yzz*two+r2(14)*xyz(2)*two &
&                   +r2(24)*xyz(3)*two+r1( 2)*yy*zz+r1( 5)*yy+r1(11)*zz+r1(14)
        cint1(1,5,2)=r5(19)+r4(13)*xyz(2)+r4(10)*xyz(3)+r4(29)*xyz(1)*two+r3( 4)*yz &
&                   +r3(20)*xy*two+r3(16)*xz*two+r3(27)*xx+r3(37)+r2( 4)*xyz2*two+r2(12)*xxy &
&                   +r2( 8)*xxz+r2(18)*xyz(2)+r2(14)*xyz(3)+r1( 2)*yz*xx+r1( 5)*yz
        cint1(2,5,2)=r5(14)+r4(14)*xyz(2)+r4( 6)*xyz(3)+r4(22)*xyz(1)+r4(29)*xyz(2) &
&                   +r3( 6)*yz+r3(17)*xy+r3(20)*yy+r3(12)*xz+r3(16)*yz+r3(27)*xy+r3(40) &
&                   +r2( 2)*xyz2+r2( 4)*yyz+r2(12)*xyy+r2( 8)*xyz2+r2(16)*xyz(3) &
&                   +r2(24)*xyz(1)+r1( 2)*yz*xy+r1( 8)*xz
        cint1(3,5,2)=r5(21)+r4(15)*xyz(2)+r4(14)*xyz(3)+r4(27)*xyz(1)+r4(29)*xyz(3) &
&                   +r3(10)*yz+r3(19)*xy+r3(20)*yz+r3(17)*xz+r3(16)*zz+r3(27)*xz+r3(36) &
&                   +r2( 6)*xyz2+r2( 4)*yzz+r2(12)*xyz2+r2( 8)*xzz+r2(16)*xyz(2) &
&                   +r2(20)*xyz(1)+r1( 2)*yz*xz+r1( 8)*xy
        cint1(4,5,2)=r5( 7)+r4( 7)*xyz(2)+r4( 2)*xyz(3)+r4(22)*xyz(2)*two+r3( 2)*yz &
&                   +r3(17)*yy*two+r3(12)*yz*two+r3(27)*yy+r3(37)*three+r2( 2)*yyz*two &
&                   +r2(12)*yyy+r2( 8)*yyz+r2(18)*xyz(2)+r2(14)*xyz(3)*three &
&                   +r2(24)*xyz(2)*two+r1( 2)*yz*yy+r1( 5)*yz+r1( 8)*yz*two
        cint1(5,5,2)=r5(15)+r4(12)*xyz(2)+r4( 7)*xyz(3)+r4(27)*xyz(2)+r4(22)*xyz(3) &
&                   +r3( 7)*yz+r3(19)*yy+r3(17)*yz+r3(17)*yz+r3(12)*zz+r3(27)*yz+r3(32) &
&                   +r3(39)+r2( 6)*yyz+r2( 2)*yzz+r2(12)*yyz+r2( 8)*yzz+r2(14)*xyz(2) &
&                   +r2(18)*xyz(3)+r2(20)*xyz(2)+r2(24)*xyz(3)+r1( 2)*yz*yz+r1( 8)*yy &
&                   +r1( 8)*zz+r1(14)
        cint1(6,5,2)=r5(18)+r4( 9)*xyz(2)+r4(12)*xyz(3)+r4(27)*xyz(3)*two+r3( 9)*yz &
&                   +r3(19)*yz*two+r3(17)*zz*two+r3(27)*zz+r3(37)*three+r2( 6)*yzz*two &
&                   +r2(12)*yzz+r2( 8)*zzz+r2(18)*xyz(2)*three+r2(14)*xyz(3) &
&                   +r2(20)*xyz(3)*two+r1( 2)*yz*zz+r1( 5)*yz+r1( 8)*yz*two
        cint1(1,6,2)=r5(20)+r4(13)*xyz(3)*two+r4(30)*xyz(1)*two+r3( 4)*zz+r3(20)*xz*four &
&                   +r3(29)*xx+r3(34)+r3(39)+r2( 4)*xzz*two+r2(12)*xxz*two+r2(18)*xyz(3)*two &
&                   +r2(22)*xyz(1)*two+r1( 2)*zz*xx+r1( 5)*zz+r1(11)*xx+r1(14)
        cint1(2,6,2)=r5(21)+r4(14)*xyz(3)*two+r4(27)*xyz(1)+r4(30)*xyz(2)+r3( 6)*zz &
&                   +r3(17)*xz*two+r3(20)*yz*two+r3(29)*xy+r3(36)+r2( 2)*xzz+r2( 4)*yzz &
&                   +r2(12)*xyz2*two+r2(20)*xyz(1)+r2(22)*xyz(2)+r1( 2)*zz*xy+r1(11)*xy
        cint1(3,6,2)=r5(17)+r4(15)*xyz(3)*two+r4(24)*xyz(1)+r4(30)*xyz(3)+r3(10)*zz &
&                   +r3(19)*xz*two+r3(20)*zz*two+r3(29)*xz+r3(40)*three+r2( 6)*xzz &
&                   +r2( 4)*zzz+r2(12)*xzz*two+r2(16)*xyz(3)*two+r2(24)*xyz(1)*three &
&                   +r2(22)*xyz(3)+r1( 2)*zz*xz+r1( 8)*xz*two+r1(11)*xz
        cint1(4,6,2)=r5(15)+r4( 7)*xyz(3)*two+r4(27)*xyz(2)*two+r3( 2)*zz+r3(17)*yz*four &
&                   +r3(29)*yy+r3(32)+r3(39)+r2( 2)*yzz*two+r2(12)*yyz*two+r2(18)*xyz(3)*two &
&                   +r2(20)*xyz(2)*two+r1( 2)*zz*yy+r1( 5)*zz+r1(11)*yy+r1(14)
        cint1(5,6,2)=r5(18)+r4(12)*xyz(3)*two+r4(24)*xyz(2)+r4(27)*xyz(3)+r3( 7)*zz &
&                   +r3(19)*yz*two+r3(17)*zz*two+r3(29)*yz+r3(37)*three+r2( 6)*yzz &
&                   +r2( 2)*zzz+r2(12)*yzz*two+r2(14)*xyz(3)*two+r2(24)*xyz(2)*three &
&                   +r2(20)*xyz(3)+r1( 2)*zz*yz+r1( 8)*yz*two+r1(11)*yz
        cint1(6,6,2)=r5( 9)+r4( 9)*xyz(3)*two+r4(24)*xyz(3)*two+r3( 9)*zz+r3(19)*zz*four &
&                   +r3(29)*zz+r3(39)*six+r2( 6)*zzz*two+r2(12)*zzz*two+r2(18)*xyz(3)*six &
&                   +r2(24)*xyz(3)*six+r1( 2)*zz*zz+r1( 5)*zz+r1( 8)*zz*four+r1(11)*zz &
&                   +r1(14)*three
        cint1(1,1,3)=r5( 5)+r4( 5)*xyz(1)*two+r4(20)*xyz(1)*two+r3( 5)*xx+r3(15)*xx*four &
&                   +r3(25)*xx+r3(35)*six+r2( 5)*xxx*two+r2(11)*xxx*two+r2(17)*xyz(1)*six &
&                   +r2(23)*xyz(1)*six+r1( 3)*xx*xx+r1( 6)*xx+r1( 9)*xx*four+r1(12)*xx &
&                   +r1(15)*three
        cint1(2,1,3)=r5(11)+r4(13)*xyz(1)*two+r4(28)*xyz(1)+r4(20)*xyz(2)+r3(10)*xx &
&                   +r3(20)*xx*two+r3(15)*xy*two+r3(25)*xy+r3(40)*three+r2( 6)*xxx &
&                   +r2( 5)*xxy+r2(11)*xxy*two+r2(18)*xyz(1)*two+r2(24)*xyz(1) &
&                   +r2(23)*xyz(2)*three+r1( 3)*xx*xy+r1( 9)*xy*two+r1(12)*xy
        cint1(3,1,3)=r5(12)+r4(11)*xyz(1)*two+r4(26)*xyz(1)+r4(20)*xyz(3)+r3( 8)*xx &
&                   +r3(18)*xx*two+r3(15)*xz*two+r3(25)*xz+r3(38)*three+r2( 3)*xxx &
&                   +r2( 5)*xxz+r2(11)*xxz*two+r2(15)*xyz(1)*two+r2(21)*xyz(1) &
&                   +r2(23)*xyz(3)*three+r1( 3)*xx*xz+r1( 9)*xz*two+r1(12)*xz
        cint1(4,1,3)=r5(19)+r4(14)*xyz(1)*two+r4(28)*xyz(2)*two+r3( 7)*xx+r3(20)*xy*four &
&                   +r3(25)*yy+r3(35)+r3(37)+r2( 6)*xxy*two+r2(11)*xyy*two+r2(17)*xyz(1)*two &
&                   +r2(24)*xyz(2)*two+r1( 3)*xx*yy+r1( 6)*xx+r1(12)*yy+r1(15)
        cint1(5,1,3)=r5(20)+r4(15)*xyz(1)*two+r4(26)*xyz(2)+r4(28)*xyz(3)+r3( 9)*xx &
&                   +r3(18)*xy*two+r3(20)*xz*two+r3(25)*yz+r3(39)+r2( 3)*xxy+r2( 6)*xxz &
&                   +r2(11)*xyz2*two+r2(21)*xyz(2)+r2(24)*xyz(3)+r1( 3)*xx*yz+r1(12)*yz
        cint1(6,1,3)=r5(16)+r4( 8)*xyz(1)*two+r4(26)*xyz(3)*two+r3( 3)*xx+r3(18)*xz*four &
&                   +r3(25)*zz+r3(33)+r3(35)+r2( 3)*xxz*two+r2(11)*xzz*two+r2(17)*xyz(1)*two &
&                   +r2(21)*xyz(3)*two+r1( 3)*xx*zz+r1( 6)*xx+r1(12)*zz+r1(15)
        cint1(1,2,3)=r5(11)+r4(13)*xyz(1)+r4( 5)*xyz(2)+r4(28)*xyz(1)*two+r3( 5)*xy &
&                   +r3(20)*xx*two+r3(15)*xy*two+r3(30)*xx+r3(40)*three+r2( 5)*xxy*two &
&                   +r2(12)*xxx+r2(11)*xxy+r2(18)*xyz(1)+r2(17)*xyz(2)*three &
&                   +r2(24)*xyz(1)*two+r1( 3)*xy*xx+r1( 6)*xy+r1( 9)*xy*two
        cint1(2,2,3)=r5(19)+r4(14)*xyz(1)+r4(13)*xyz(2)+r4(29)*xyz(1)+r4(28)*xyz(2) &
&                   +r3(10)*xy+r3(17)*xx+r3(20)*xy+r3(20)*xy+r3(15)*yy+r3(30)*xy+r3(35) &
&                   +r3(37)+r2( 6)*xxy+r2( 5)*xyy+r2(12)*xxy+r2(11)*xyy+r2(17)*xyz(1) &
&                   +r2(18)*xyz(2)+r2(23)*xyz(1)+r2(24)*xyz(2)+r1( 3)*xy*xy+r1( 9)*xx &
&                   +r1( 9)*yy+r1(15)
        cint1(3,2,3)=r5(20)+r4(15)*xyz(1)+r4(11)*xyz(2)+r4(30)*xyz(1)+r4(28)*xyz(3) &
&                   +r3( 8)*xy+r3(19)*xx+r3(20)*xz+r3(18)*xy+r3(15)*yz+r3(30)*xz+r3(39) &
&                   +r2( 3)*xxy+r2( 5)*xyz2+r2(12)*xxz+r2(11)*xyz2+r2(15)*xyz(2) &
&                   +r2(24)*xyz(3)+r1( 3)*xy*xz+r1( 9)*yz
        cint1(4,2,3)=r5(14)+r4( 7)*xyz(1)+r4(14)*xyz(2)+r4(29)*xyz(2)*two+r3( 7)*xy &
&                   +r3(17)*xy*two+r3(20)*yy*two+r3(30)*yy+r3(40)*three+r2( 6)*xyy*two &
&                   +r2(12)*xyy+r2(11)*yyy+r2(18)*xyz(1)*three+r2(17)*xyz(2) &
&                   +r2(23)*xyz(2)*two+r1( 3)*xy*yy+r1( 6)*xy+r1( 9)*xy*two
        cint1(5,2,3)=r5(21)+r4(12)*xyz(1)+r4(15)*xyz(2)+r4(30)*xyz(2)+r4(29)*xyz(3) &
&                   +r3( 9)*xy+r3(19)*xy+r3(17)*xz+r3(18)*yy+r3(20)*yz+r3(30)*yz+r3(38) &
&                   +r2( 3)*xyy+r2( 6)*xyz2+r2(12)*xyz2+r2(11)*yyz+r2(15)*xyz(1) &
&                   +r2(23)*xyz(3)+r1( 3)*xy*yz+r1( 9)*xz
        cint1(6,2,3)=r5(17)+r4( 9)*xyz(1)+r4( 8)*xyz(2)+r4(30)*xyz(3)*two+r3( 3)*xy &
&                   +r3(19)*xz*two+r3(18)*yz*two+r3(30)*zz+r3(40)+r2( 3)*xyz2*two+r2(12)*xzz &
&                   +r2(11)*yzz+r2(18)*xyz(1)+r2(17)*xyz(2)+r1( 3)*xy*zz+r1( 6)*xy
        cint1(1,3,3)=r5(12)+r4(11)*xyz(1)+r4( 5)*xyz(3)+r4(26)*xyz(1)*two+r3( 5)*xz &
&                   +r3(18)*xx*two+r3(15)*xz*two+r3(28)*xx+r3(38)*three+r2( 5)*xxz*two &
&                   +r2( 9)*xxx+r2(11)*xxz+r2(15)*xyz(1)+r2(17)*xyz(3)*three &
&                   +r2(21)*xyz(1)*two+r1( 3)*xz*xx+r1( 6)*xz+r1( 9)*xz*two
        cint1(2,3,3)=r5(20)+r4(15)*xyz(1)+r4(13)*xyz(3)+r4(30)*xyz(1)+r4(26)*xyz(2) &
&                   +r3(10)*xz+r3(19)*xx+r3(18)*xy+r3(20)*xz+r3(15)*yz+r3(28)*xy+r3(39) &
&                   +r2( 6)*xxz+r2( 5)*xyz2+r2( 9)*xxy+r2(11)*xyz2+r2(18)*xyz(3) &
&                   +r2(21)*xyz(2)+r1( 3)*xz*xy+r1( 9)*yz
        cint1(3,3,3)=r5(16)+r4( 8)*xyz(1)+r4(11)*xyz(3)+r4(23)*xyz(1)+r4(26)*xyz(3) &
&                   +r3( 8)*xz+r3(13)*xx+r3(18)*xz+r3(18)*xz+r3(15)*zz+r3(28)*xz+r3(33) &
&                   +r3(35)+r2( 3)*xxz+r2( 5)*xzz+r2( 9)*xxz+r2(11)*xzz+r2(17)*xyz(1) &
&                   +r2(15)*xyz(3)+r2(23)*xyz(1)+r2(21)*xyz(3)+r1( 3)*xz*xz+r1( 9)*xx &
&                   +r1( 9)*zz+r1(15)
        cint1(4,3,3)=r5(21)+r4(12)*xyz(1)+r4(14)*xyz(3)+r4(30)*xyz(2)*two+r3( 7)*xz &
&                   +r3(19)*xy*two+r3(20)*yz*two+r3(28)*yy+r3(38)+r2( 6)*xyz2*two+r2( 9)*xyy &
&                   +r2(11)*yyz+r2(15)*xyz(1)+r2(17)*xyz(3)+r1( 3)*xz*yy+r1( 6)*xz
        cint1(5,3,3)=r5(17)+r4( 9)*xyz(1)+r4(15)*xyz(3)+r4(23)*xyz(2)+r4(30)*xyz(3) &
&                   +r3( 9)*xz+r3(13)*xy+r3(19)*xz+r3(18)*yz+r3(20)*zz+r3(28)*yz+r3(40) &
&                   +r2( 3)*xyz2+r2( 6)*xzz+r2( 9)*xyz2+r2(11)*yzz+r2(18)*xyz(1)+r2(23)*xyz(2) &
&                   +r1( 3)*xz*yz+r1( 9)*xy
        cint1(6,3,3)=r5( 8)+r4( 3)*xyz(1)+r4( 8)*xyz(3)+r4(23)*xyz(3)*two+r3( 3)*xz &
&                   +r3(13)*xz*two+r3(18)*zz*two+r3(28)*zz+r3(38)*three+r2( 3)*xzz*two &
&                   +r2( 9)*xzz+r2(11)*zzz+r2(15)*xyz(1)*three+r2(17)*xyz(3) &
&                   +r2(23)*xyz(3)*two+r1( 3)*xz*zz+r1( 6)*xz+r1( 9)*xz*two
        cint1(1,4,3)=r5(19)+r4(13)*xyz(2)*two+r4(29)*xyz(1)*two+r3( 5)*yy+r3(20)*xy*four &
&                   +r3(27)*xx+r3(35)+r3(37)+r2( 5)*xyy*two+r2(12)*xxy*two+r2(18)*xyz(2)*two &
&                   +r2(23)*xyz(1)*two+r1( 3)*yy*xx+r1( 6)*yy+r1(12)*xx+r1(15)
        cint1(2,4,3)=r5(14)+r4(14)*xyz(2)*two+r4(22)*xyz(1)+r4(29)*xyz(2)+r3(10)*yy &
&                   +r3(17)*xy*two+r3(20)*yy*two+r3(27)*xy+r3(40)*three+r2( 6)*xyy &
&                   +r2( 5)*yyy+r2(12)*xyy*two+r2(17)*xyz(2)*two+r2(24)*xyz(1)*three &
&                   +r2(23)*xyz(2)+r1( 3)*yy*xy+r1( 9)*xy*two+r1(12)*xy
        cint1(3,4,3)=r5(21)+r4(15)*xyz(2)*two+r4(27)*xyz(1)+r4(29)*xyz(3)+r3( 8)*yy &
&                   +r3(19)*xy*two+r3(20)*yz*two+r3(27)*xz+r3(38)+r2( 3)*xyy+r2( 5)*yyz &
&                   +r2(12)*xyz2*two+r2(21)*xyz(1)+r2(23)*xyz(3)+r1( 3)*yy*xz+r1(12)*xz
        cint1(4,4,3)=r5( 7)+r4( 7)*xyz(2)*two+r4(22)*xyz(2)*two+r3( 7)*yy+r3(17)*yy*four &
&                   +r3(27)*yy+r3(37)*six+r2( 6)*yyy*two+r2(12)*yyy*two+r2(18)*xyz(2)*six &
&                   +r2(24)*xyz(2)*six+r1( 3)*yy*yy+r1( 6)*yy+r1( 9)*yy*four+r1(12)*yy &
&                   +r1(15)*three
        cint1(5,4,3)=r5(15)+r4(12)*xyz(2)*two+r4(27)*xyz(2)+r4(22)*xyz(3)+r3( 9)*yy &
&                   +r3(19)*yy*two+r3(17)*yz*two+r3(27)*yz+r3(39)*three+r2( 3)*yyy &
&                   +r2( 6)*yyz+r2(12)*yyz*two+r2(15)*xyz(2)*two+r2(21)*xyz(2) &
&                   +r2(24)*xyz(3)*three+r1( 3)*yy*yz+r1( 9)*yz*two+r1(12)*yz
        cint1(6,4,3)=r5(18)+r4( 9)*xyz(2)*two+r4(27)*xyz(3)*two+r3( 3)*yy+r3(19)*yz*four &
&                   +r3(27)*zz+r3(33)+r3(37)+r2( 3)*yyz*two+r2(12)*yzz*two+r2(18)*xyz(2)*two &
&                   +r2(21)*xyz(3)*two+r1( 3)*yy*zz+r1( 6)*yy+r1(12)*zz+r1(15)
        cint1(1,5,3)=r5(20)+r4(11)*xyz(2)+r4(13)*xyz(3)+r4(30)*xyz(1)*two+r3( 5)*yz &
&                   +r3(18)*xy*two+r3(20)*xz*two+r3(29)*xx+r3(39)+r2( 5)*xyz2*two+r2( 9)*xxy &
&                   +r2(12)*xxz+r2(15)*xyz(2)+r2(18)*xyz(3)+r1( 3)*yz*xx+r1( 6)*yz
        cint1(2,5,3)=r5(21)+r4(15)*xyz(2)+r4(14)*xyz(3)+r4(27)*xyz(1)+r4(30)*xyz(2) &
&                   +r3(10)*yz+r3(19)*xy+r3(18)*yy+r3(17)*xz+r3(20)*yz+r3(29)*xy+r3(38) &
&                   +r2( 6)*xyz2+r2( 5)*yyz+r2( 9)*xyy+r2(12)*xyz2+r2(17)*xyz(3) &
&                   +r2(21)*xyz(1)+r1( 3)*yz*xy+r1( 9)*xz
        cint1(3,5,3)=r5(17)+r4( 8)*xyz(2)+r4(15)*xyz(3)+r4(24)*xyz(1)+r4(30)*xyz(3) &
&                   +r3( 8)*yz+r3(13)*xy+r3(18)*yz+r3(19)*xz+r3(20)*zz+r3(29)*xz+r3(40) &
&                   +r2( 3)*xyz2+r2( 5)*yzz+r2( 9)*xyz2+r2(12)*xzz+r2(17)*xyz(2) &
&                   +r2(24)*xyz(1)+r1( 3)*yz*xz+r1( 9)*xy
        cint1(4,5,3)=r5(15)+r4(12)*xyz(2)+r4( 7)*xyz(3)+r4(27)*xyz(2)*two+r3( 7)*yz &
&                   +r3(19)*yy*two+r3(17)*yz*two+r3(29)*yy+r3(39)*three+r2( 6)*yyz*two &
&                   +r2( 9)*yyy+r2(12)*yyz+r2(15)*xyz(2)+r2(18)*xyz(3)*three &
&                   +r2(21)*xyz(2)*two+r1( 3)*yz*yy+r1( 6)*yz+r1( 9)*yz*two
        cint1(5,5,3)=r5(18)+r4( 9)*xyz(2)+r4(12)*xyz(3)+r4(24)*xyz(2)+r4(27)*xyz(3) &
&                   +r3( 9)*yz+r3(13)*yy+r3(19)*yz+r3(19)*yz+r3(17)*zz+r3(29)*yz+r3(33) &
&                   +r3(37)+r2( 3)*yyz+r2( 6)*yzz+r2( 9)*yyz+r2(12)*yzz+r2(18)*xyz(2) &
&                   +r2(15)*xyz(3)+r2(24)*xyz(2)+r2(21)*xyz(3)+r1( 3)*yz*yz+r1( 9)*yy &
&                   +r1( 9)*zz+r1(15)
        cint1(6,5,3)=r5( 9)+r4( 3)*xyz(2)+r4( 9)*xyz(3)+r4(24)*xyz(3)*two+r3( 3)*yz &
&                   +r3(13)*yz*two+r3(19)*zz*two+r3(29)*zz+r3(39)*three+r2( 3)*yzz*two &
&                   +r2( 9)*yzz+r2(12)*zzz+r2(15)*xyz(2)*three+r2(18)*xyz(3) &
&                   +r2(24)*xyz(3)*two+r1( 3)*yz*zz+r1( 6)*yz+r1( 9)*yz*two
        cint1(1,6,3)=r5(16)+r4(11)*xyz(3)*two+r4(23)*xyz(1)*two+r3( 5)*zz+r3(18)*xz*four &
&                   +r3(23)*xx+r3(33)+r3(35)+r2( 5)*xzz*two+r2( 9)*xxz*two+r2(15)*xyz(3)*two &
&                   +r2(23)*xyz(1)*two+r1( 3)*zz*xx+r1( 6)*zz+r1(12)*xx+r1(15)
        cint1(2,6,3)=r5(17)+r4(15)*xyz(3)*two+r4(24)*xyz(1)+r4(23)*xyz(2)+r3(10)*zz &
&                   +r3(19)*xz*two+r3(18)*yz*two+r3(23)*xy+r3(40)+r2( 6)*xzz+r2( 5)*yzz &
&                   +r2( 9)*xyz2*two+r2(24)*xyz(1)+r2(23)*xyz(2)+r1( 3)*zz*xy+r1(12)*xy
        cint1(3,6,3)=r5( 8)+r4( 8)*xyz(3)*two+r4(18)*xyz(1)+r4(23)*xyz(3)+r3( 8)*zz &
&                   +r3(13)*xz*two+r3(18)*zz*two+r3(23)*xz+r3(38)*three+r2( 3)*xzz &
&                   +r2( 5)*zzz+r2( 9)*xzz*two+r2(17)*xyz(3)*two+r2(21)*xyz(1)*three &
&                   +r2(23)*xyz(3)+r1( 3)*zz*xz+r1( 9)*xz*two+r1(12)*xz
        cint1(4,6,3)=r5(18)+r4(12)*xyz(3)*two+r4(24)*xyz(2)*two+r3( 7)*zz+r3(19)*yz*four &
&                   +r3(23)*yy+r3(33)+r3(37)+r2( 6)*yzz*two+r2( 9)*yyz*two+r2(15)*xyz(3)*two &
&                   +r2(24)*xyz(2)*two+r1( 3)*zz*yy+r1( 6)*zz+r1(12)*yy+r1(15)
        cint1(5,6,3)=r5( 9)+r4( 9)*xyz(3)*two+r4(18)*xyz(2)+r4(24)*xyz(3)+r3( 9)*zz &
&                   +r3(13)*yz*two+r3(19)*zz*two+r3(23)*yz+r3(39)*three+r2( 3)*yzz &
&                   +r2( 6)*zzz+r2( 9)*yzz*two+r2(18)*xyz(3)*two+r2(21)*xyz(2)*three &
&                   +r2(24)*xyz(3)+r1( 3)*zz*yz+r1( 9)*yz*two+r1(12)*yz
        cint1(6,6,3)=r5( 3)+r4( 3)*xyz(3)*two+r4(18)*xyz(3)*two+r3( 3)*zz+r3(13)*zz*four &
&                   +r3(23)*zz+r3(33)*six+r2( 3)*zzz*two+r2( 9)*zzz*two+r2(15)*xyz(3)*six &
&                   +r2(21)*xyz(3)*six+r1( 3)*zz*zz+r1( 6)*zz+r1( 9)*zz*four+r1(12)*zz &
&                   +r1(15)*three
!
        if(nbfij(1) == 6) then
          do k= 1,3
            do j= 1,6
              cint1(j,2,k)= cint1(j,2,k)*sqrt3
              cint1(j,3,k)= cint1(j,3,k)*sqrt3
              cint1(j,5,k)= cint1(j,5,k)*sqrt3
            enddo
          enddo
        else
          do k= 1,3
            do j= 1,6
              do i= 1,6
                work(i)= cint1(j,i,k)
              enddo
              cint1(j,1,k)= work(2)*sqrt3
              cint1(j,2,k)= work(5)*sqrt3
              cint1(j,3,k)=(work(6)*two-work(1)-work(4))*half
              cint1(j,4,k)= work(3)*sqrt3
              cint1(j,5,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        endif
        if(nbfij(2) == 6) then
          do k= 1,3
            do i= 1,nbfij(1)
              cint1(2,i,k)= cint1(2,i,k)*sqrt3
              cint1(3,i,k)= cint1(3,i,k)*sqrt3
              cint1(5,i,k)= cint1(5,i,k)*sqrt3
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,nbfij(1)
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
!
        if(.not.iandj) then
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)
              enddo
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom) &
&                              +fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)*half
              enddo
            enddo
          enddo
        endif
      enddo
!
      return
end



