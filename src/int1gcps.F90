!---------------------------------------------------------------------------
  subroutine int1gcps(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <p|V'|s>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, p15=1.5D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:2), tinv, r1(3), r2(6), cint1(3,3)
!
      do iatom= 1,natom
        do i= 1,3
          r1(i)= zero
        enddo
        do i= 1,6
          r2(i)= zero
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
          rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
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
            do ii= 0,2
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
          do i= 1,2
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          do i= 1,3
            r1(i)= r1(i)-ft(1)*pcxyz(i)*ex2i*extwo
          enddo
          r2(1)= r2(1)+ft(2)*pcxyz(1)*pcxyz(1)*extwo-ft(1)
          r2(2)= r2(2)+ft(2)*pcxyz(2)*pcxyz(2)*extwo-ft(1)
          r2(3)= r2(3)+ft(2)*pcxyz(3)*pcxyz(3)*extwo-ft(1)
          r2(4)= r2(4)+ft(2)*pcxyz(1)*pcxyz(2)*extwo
          r2(5)= r2(5)+ft(2)*pcxyz(1)*pcxyz(3)*extwo
          r2(6)= r2(6)+ft(2)*pcxyz(2)*pcxyz(3)*extwo
        enddo
        cint1(1,1)= r2(1)+r1(1)*xyz(1)
        cint1(2,1)= r2(4)+r1(1)*xyz(2)
        cint1(3,1)= r2(5)+r1(1)*xyz(3)
        cint1(1,2)= r2(4)+r1(2)*xyz(1)
        cint1(2,2)= r2(2)+r1(2)*xyz(2)
        cint1(3,2)= r2(6)+r1(2)*xyz(3)
        cint1(1,3)= r2(5)+r1(3)*xyz(1)
        cint1(2,3)= r2(6)+r1(3)*xyz(2)
        cint1(3,3)= r2(3)+r1(3)*xyz(3)
        do k= 1,3
          do j= 1,3
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+1)*cint1(j,k)
          enddo
        enddo
      enddo
!
      return
end

