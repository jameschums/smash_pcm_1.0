!---------------------------------------------------------------------------
  subroutine int1gcpp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,iandj)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <p|V'|p>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, five=5.0D+00, p15=1.5D+00, p25=2.5D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:3), fttwo, fpc2(6)
      real(8) :: pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r1(6), r2(12), r3(10), xx, yy, zz, xy, xz, yz, cint1(3,3,3)
      logical,intent(in) :: iandj
!
!      write(*,'("james int1gcpp.F90 Calculate Helmann-Feynman gradient term <p|V´|p>")')
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      do iatom= 1,natom
        do i= 1,6
          r1(i) = zero
        enddo
        do i= 1,12
          r2(i) = zero
        enddo
        do i= 1,10
          r3(i) = zero
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
            do ii= 0,3
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,3
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          fttwo= ft(2)*extwo
          ft(3)= ft(3)*extwo
          do i= 1,3
            r1(i)  = r1(i)  +ft(1)*pcxyz(i)*ex2j*ex2i*extwo
            r1(i+3)= r1(i+3)-ft(1)*pcxyz(i)
          enddo
          fpc2(1)= fttwo*pxx-ft(1)
          fpc2(2)= fttwo*pyy-ft(1)
          fpc2(3)= fttwo*pzz-ft(1)
          fpc2(4)= fttwo*pxy
          fpc2(5)= fttwo*pxz
          fpc2(6)= fttwo*pyz
          do i= 1,6
            r2(i)  = r2(i)  -fpc2(i)*ex2j
            r2(i+6)= r2(i+6)+fpc2(i)*ex2i
          enddo
          r3(1) = r3(1) -ft(3)*pxx*pcxyz(1)+ft(2)*pcxyz(1)*three
          r3(2) = r3(2) -ft(3)*pyy*pcxyz(2)+ft(2)*pcxyz(2)*three
          r3(3) = r3(3) -ft(3)*pzz*pcxyz(3)+ft(2)*pcxyz(3)*three
          r3(4) = r3(4) -ft(3)*pxx*pcxyz(2)+ft(2)*pcxyz(2)
          r3(5) = r3(5) -ft(3)*pxx*pcxyz(3)+ft(2)*pcxyz(3)
          r3(6) = r3(6) -ft(3)*pxy*pcxyz(2)+ft(2)*pcxyz(1)
          r3(7) = r3(7) -ft(3)*pyy*pcxyz(3)+ft(2)*pcxyz(3)
          r3(8) = r3(8) -ft(3)*pxz*pcxyz(3)+ft(2)*pcxyz(1)
          r3(9) = r3(9) -ft(3)*pyz*pcxyz(3)+ft(2)*pcxyz(2)
          r3(10)= r3(10)-ft(3)*pxy*pcxyz(3)
        enddo
        cint1(1,1,1)= r3( 1)+r2(1)*xyz(1)+r2( 7)*xyz(1)+r1(1)*xx+r1(4)
        cint1(2,1,1)= r3( 4)+r2(4)*xyz(1)+r2( 7)*xyz(2)+r1(1)*xy
        cint1(3,1,1)= r3( 5)+r2(5)*xyz(1)+r2( 7)*xyz(3)+r1(1)*xz
        cint1(1,2,1)= r3( 4)+r2(1)*xyz(2)+r2(10)*xyz(1)+r1(1)*xy
        cint1(2,2,1)= r3( 6)+r2(4)*xyz(2)+r2(10)*xyz(2)+r1(1)*yy+r1(4)
        cint1(3,2,1)= r3(10)+r2(5)*xyz(2)+r2(10)*xyz(3)+r1(1)*yz
        cint1(1,3,1)= r3( 5)+r2(1)*xyz(3)+r2(11)*xyz(1)+r1(1)*xz
        cint1(2,3,1)= r3(10)+r2(4)*xyz(3)+r2(11)*xyz(2)+r1(1)*yz
        cint1(3,3,1)= r3( 8)+r2(5)*xyz(3)+r2(11)*xyz(3)+r1(1)*zz+r1(4)
        cint1(1,1,2)= r3( 4)+r2(4)*xyz(1)+r2(10)*xyz(1)+r1(2)*xx+r1(5)
        cint1(2,1,2)= r3( 6)+r2(2)*xyz(1)+r2(10)*xyz(2)+r1(2)*xy
        cint1(3,1,2)= r3(10)+r2(6)*xyz(1)+r2(10)*xyz(3)+r1(2)*xz
        cint1(1,2,2)= r3( 6)+r2(4)*xyz(2)+r2( 8)*xyz(1)+r1(2)*xy
        cint1(2,2,2)= r3( 2)+r2(2)*xyz(2)+r2( 8)*xyz(2)+r1(2)*yy+r1(5)
        cint1(3,2,2)= r3( 7)+r2(6)*xyz(2)+r2( 8)*xyz(3)+r1(2)*yz
        cint1(1,3,2)= r3(10)+r2(4)*xyz(3)+r2(12)*xyz(1)+r1(2)*xz
        cint1(2,3,2)= r3( 7)+r2(2)*xyz(3)+r2(12)*xyz(2)+r1(2)*yz
        cint1(3,3,2)= r3( 9)+r2(6)*xyz(3)+r2(12)*xyz(3)+r1(2)*zz+r1(5)
        cint1(1,1,3)= r3( 5)+r2(5)*xyz(1)+r2(11)*xyz(1)+r1(3)*xx+r1(6)
        cint1(2,1,3)= r3(10)+r2(6)*xyz(1)+r2(11)*xyz(2)+r1(3)*xy
        cint1(3,1,3)= r3( 8)+r2(3)*xyz(1)+r2(11)*xyz(3)+r1(3)*xz
        cint1(1,2,3)= r3(10)+r2(5)*xyz(2)+r2(12)*xyz(1)+r1(3)*xy
        cint1(2,2,3)= r3( 7)+r2(6)*xyz(2)+r2(12)*xyz(2)+r1(3)*yy+r1(6)
        cint1(3,2,3)= r3( 9)+r2(3)*xyz(2)+r2(12)*xyz(3)+r1(3)*yz
        cint1(1,3,3)= r3( 8)+r2(5)*xyz(3)+r2( 9)*xyz(1)+r1(3)*xz
        cint1(2,3,3)= r3( 9)+r2(6)*xyz(3)+r2( 9)*xyz(2)+r1(3)*yz
        cint1(3,3,3)= r3( 3)+r2(3)*xyz(3)+r2( 9)*xyz(3)+r1(3)*zz+r1(6)
        if(.not.iandj) then
          do k= 1,3
            do i= 1,3
              do j= 1,3
                egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)
              enddo
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,3
              do j= 1,3
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

