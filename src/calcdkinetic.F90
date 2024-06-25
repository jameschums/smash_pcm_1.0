!---------------------------------------------------
  subroutine calcdkinetic(egrad,fulldmtrx,ish,jsh)
!---------------------------------------------------
!
! Driver of kinetic derivative term
!
! In : fulldmtrx (Density matrix)
!    : ish, jsh  (Shell indices)
! Inout : egrad  (Energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      use modmolecule, only : natom, coord
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj
      integer :: nbfi, nbfj, iprim, jprim, ncarti, ncartj, i, j, iang, jang, ii
      integer :: ix, jx, iy, jy, iz, jz, ncart(0:6)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exi2, exj, exj2, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpij(3,2), cij
      real(8) :: xyzint(3), sx(0:7,0:8,4), sy(0:7,0:8,4), sz(0:7,0:8,4)
      real(8) :: dtint(28,28,3)
      data ncart/1,3,6,10,15,21,28/
!
      iatom = locatom(ish)
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      nprimi= mprim(ish)
      nangi = mtype(ish)
      nbfi  = mbf(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
      nprimj= mprim(jsh)
      nangj = mtype(jsh)
      nbfj  = mbf(jsh)
!
      if((nangi > 6).or.(nangj > 6))then
        write(*,'(" Error! This program supports up to h function in calcdkinetic.")')
        call iabort
      endif
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
          dtint(j,i,1)= zero
          dtint(j,i,2)= zero
          dtint(j,i,3)= zero
        enddo
      enddo
!
! Calculate overlap derivative for each primitive
!
      do iprim= 1,nprimi
        exi= ex(iloc+iprim)
        ci = coeff(iloc+iprim)
        exi2=-two*exi*exi
        do jprim= 1,nprimj
          exj= ex(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          exj2= exj*two
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff(jloc+jprim)*fac
!
! overlap term
!
          do iang= 0,nangi+2
            do jang= 0,nangj+1
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang,1)= xyzint(1)*ex3
              sy(jang,iang,1)= xyzint(2)*ex3
              sz(jang,iang,1)= xyzint(3)*ex3
            enddo
          enddo
!
! kinetic term
!
          do iang= 0,nangi
            do jang= 0,nangj+1
              sx(jang,iang,2)= sx(jang,iang+2,1)*exi2+sx(jang,iang,1)*(two*iang+one)*exi
              sy(jang,iang,2)= sy(jang,iang+2,1)*exi2+sy(jang,iang,1)*(two*iang+one)*exi
              sz(jang,iang,2)= sz(jang,iang+2,1)*exi2+sz(jang,iang,1)*(two*iang+one)*exi
              if(iang >= 2) then
                sx(jang,iang,2)=sx(jang,iang,2)-sx(jang,iang-2,1)*half*iang*(iang-1)
                sy(jang,iang,2)=sy(jang,iang,2)-sy(jang,iang-2,1)*half*iang*(iang-1)
                sz(jang,iang,2)=sz(jang,iang,2)-sz(jang,iang-2,1)*half*iang*(iang-1)
              endif
            enddo
          enddo
!
! overlap derivative term
!
          do iang= 0,nangi
            sx(0,iang,3)= sx(1,iang,1)*exj2
            sy(0,iang,3)= sy(1,iang,1)*exj2
            sz(0,iang,3)= sz(1,iang,1)*exj2
          enddo
          do iang= 0,nangi
            do jang= 1,nangj
              sx(jang,iang,3)= sx(jang+1,iang,1)*exj2-sx(jang-1,iang,1)*jang
              sy(jang,iang,3)= sy(jang+1,iang,1)*exj2-sy(jang-1,iang,1)*jang
              sz(jang,iang,3)= sz(jang+1,iang,1)*exj2-sz(jang-1,iang,1)*jang
            enddo
          enddo
!
! kinetic derivative term
!
          do iang= 0,nangi
            sx(0,iang,4)= sx(1,iang,2)*exj2
            sy(0,iang,4)= sy(1,iang,2)*exj2
            sz(0,iang,4)= sz(1,iang,2)*exj2
          enddo
          do iang= 0,nangi
            do jang= 1,nangj
              sx(jang,iang,4)= sx(jang+1,iang,2)*exj2-sx(jang-1,iang,2)*jang
              sy(jang,iang,4)= sy(jang+1,iang,2)*exj2-sy(jang-1,iang,2)*jang
              sz(jang,iang,4)= sz(jang+1,iang,2)*exj2-sz(jang-1,iang,2)*jang
            enddo
          enddo
!
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
                  dtint(j,i,1)= dtint(j,i,1)+cij*(sx(jx,ix,4)*sy(jy,iy,1)*sz(jz,iz,1) &
&                                                +sx(jx,ix,3)*sy(jy,iy,2)*sz(jz,iz,1) &
&                                                +sx(jx,ix,3)*sy(jy,iy,1)*sz(jz,iz,2)) 
                  dtint(j,i,2)= dtint(j,i,2)+cij*(sx(jx,ix,2)*sy(jy,iy,3)*sz(jz,iz,1) &
&                                                +sx(jx,ix,1)*sy(jy,iy,4)*sz(jz,iz,1) &
&                                                +sx(jx,ix,1)*sy(jy,iy,3)*sz(jz,iz,2)) 
                  dtint(j,i,3)= dtint(j,i,3)+cij*(sx(jx,ix,2)*sy(jy,iy,1)*sz(jz,iz,3) &
&                                                +sx(jx,ix,1)*sy(jy,iy,2)*sz(jz,iz,3) &
&                                                +sx(jx,ix,1)*sy(jy,iy,1)*sz(jz,iz,4)) 
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(dtint(1,1,1),nbfi,nbfj,ncarti)
        call nrmlz1(dtint(1,1,2),nbfi,nbfj,ncarti)
        call nrmlz1(dtint(1,1,3),nbfi,nbfj,ncarti)
      endif
!
      do i= 1,nbfi
        ii= ilocbf+i
        do j= 1,nbfj
          egrad(1,jatom)= egrad(1,jatom)+fulldmtrx(jlocbf+j,ii)*dtint(j,i,1)
          egrad(2,jatom)= egrad(2,jatom)+fulldmtrx(jlocbf+j,ii)*dtint(j,i,2)
          egrad(3,jatom)= egrad(3,jatom)+fulldmtrx(jlocbf+j,ii)*dtint(j,i,3)
        enddo
      enddo
!
      return
end


!--------------------------------------------------------
  subroutine calcdcoulomb(egrad,fulldmtrx,ish,jsh,len1)
!--------------------------------------------------------
!
! Driver of Coulomb derivative term
!
! In : fulldmtrx (density matrix)
!    : ish, jsh  (shell indices)
! Inout : egrad  (energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom, iloc, jloc, ilocbf, jlocbf
      integer :: iprim, jprim, i, j, k, ii, ncart(0:6)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, eight=8.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p24=24.0D+00, p30=30.0D+00
      real(8),parameter :: p40=40.0D+00
      real(8),parameter :: third=3.333333333333333D-01, eighth=0.125D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrtthird=0.5773502691896258D+00, sqrtfifth=0.4472135954999579D+00
      real(8),parameter :: sqrt3fifth=0.7745966692414834D+00, sqrtseventh=0.3779644730092272D+00
      real(8),parameter :: sqrtinv35=0.1690308509457033D+00, sqrt3inv35=0.2927700218845599D+00
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064591D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrtinv15=2.581988897471611D-01, sqrtinv21=2.182178902359924D-01
      real(8),parameter :: sqrtinv63=1.259881576697424D-01, sqrtinv105=9.759000729485332D-02
      real(8),parameter :: sqrtinv11=3.015113445777636D-01, sqrtinv33=1.740776559556978D-01
      real(8),parameter :: sqrtinv99=1.005037815259212D-01, sqrtinv231=6.579516949597690D-02
      real(8),parameter :: sqrtinv385=5.096471914376255D-02, sqrt5inv231=1.471224715841249D-01
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
      real(8),parameter :: fach1=0.70156076002011400D+00 ! sqrt(63/2)/8
      real(8),parameter :: fach2=2.21852991866235601D+00 ! sqrt(315)/8
      real(8),parameter :: fach3=0.52291251658379721D+00 ! sqrt(35/2)/8
      real(8),parameter :: fach4=2.56173769148989959D+00 ! sqrt(105)/4
      real(8),parameter :: fach5=0.48412291827592711D+00 ! sqrt(15)/8
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: exij(mxprsh,2), cij(mxprsh,2), coordij(3,2)
      real(8) :: cint1(len1,len1), cint2(len1,len1), dcint(28,28,3), work(21)
      data ncart/1,3,6,10,15,21,28/
!
      nangij(1)= mtype(ish)
      nangij(2)= mtype(jsh)
      nprimij(1)= mprim(ish)
      nprimij(2)= mprim(jsh)
      nbfij(1)  = mbf(ish)
      nbfij(2)  = mbf(jsh)
      iatom = locatom(ish)
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
!
      do i= 1,3
        coordij(i,1)= coord(i,iatom)
        coordij(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        cij(iprim,1) = coeff(iloc+iprim)
      enddo
!
      nangij(2)= nangij(2)+1
      nbfij(2)= ncart(nangij(2))
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        cij(jprim,2) = coeff(jloc+jprim)*two*ex(jloc+jprim)
      enddo
!
      if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
        call int1cmd(cint1,exij,cij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      else
        if((nangij(1) > 6).or.(nangij(2) > 7))then
          write(*,'(" Error! This program supports up to h function in int1c")')
          call iabort
        endif
        call int1rys(cint1,exij,cij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      endif
!
      if(mtype(jsh) >= 1) then
        nangij(2)= mtype(jsh)-1
        nbfij(2)= ncart(nangij(2))
        do jprim= 1,nprimij(2)
          cij(jprim,2) = coeff(jloc+jprim)
        enddo
        if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
          call int1cmd(cint2,exij,cij,coordij,coord,znuc,natom, &
&                      nprimij,nangij,nbfij,len1,mxprsh,threshex)
        else
          if((nangij(1) > 6).or.(nangij(2) > 5))then
            write(*,'(" Error! This program supports up to h function in int1c")')
            call iabort
          endif
!
          call int1rys(cint2,exij,cij,coordij,coord,znuc,natom, &
&                      nprimij,nangij,nbfij,len1,mxprsh,threshex)
!
        endif
      else
        cint2(1:nbfij(2),1:nbfij(1))= zero
      endif
!
      select case(mtype(jsh))
        case (0)
          do i= 1,nbfij(1)
            dcint(1,i,1)= cint1(1,i)
            dcint(1,i,2)= cint1(2,i)
            dcint(1,i,3)= cint1(3,i)
          enddo
        case (1)
          do i= 1,nbfij(1)
            dcint(1,i,1)= cint1(1,i)          -cint2(1,i)
            dcint(2,i,1)= cint1(2,i)*sqrtthird
            dcint(3,i,1)= cint1(3,i)*sqrtthird
            dcint(1,i,2)= cint1(2,i)*sqrtthird
            dcint(2,i,2)= cint1(4,i)          -cint2(1,i)
            dcint(3,i,2)= cint1(5,i)*sqrtthird
            dcint(1,i,3)= cint1(3,i)*sqrtthird
            dcint(2,i,3)= cint1(5,i)*sqrtthird
            dcint(3,i,3)= cint1(6,i)          -cint2(1,i)
          enddo
        case (2)
          do i= 1,nbfij(1)
            dcint(1,i,1)= cint1( 1,i)           -cint2(1,i)*two
            dcint(2,i,1)= cint1( 2,i)*sqrt3fifth-cint2(2,i)*sqrt3
            dcint(3,i,1)= cint1( 3,i)*sqrt3fifth-cint2(3,i)*sqrt3
            dcint(4,i,1)= cint1( 4,i)*sqrtfifth
            dcint(5,i,1)= cint1( 5,i)*sqrtfifth
            dcint(6,i,1)= cint1( 6,i)*sqrtfifth
            dcint(1,i,2)= cint1( 2,i)*sqrtfifth
            dcint(2,i,2)= cint1( 4,i)*sqrt3fifth-cint2(1,i)*sqrt3
            dcint(3,i,2)= cint1( 5,i)*sqrtfifth
            dcint(4,i,2)= cint1( 7,i)           -cint2(2,i)*two
            dcint(5,i,2)= cint1( 8,i)*sqrt3fifth-cint2(3,i)*sqrt3
            dcint(6,i,2)= cint1( 9,i)*sqrtfifth
            dcint(1,i,3)= cint1( 3,i)*sqrtfifth
            dcint(2,i,3)= cint1( 5,i)*sqrtfifth
            dcint(3,i,3)= cint1( 6,i)*sqrt3fifth-cint2(1,i)*sqrt3
            dcint(4,i,3)= cint1( 8,i)*sqrtfifth
            dcint(5,i,3)= cint1( 9,i)*sqrt3fifth-cint2(2,i)*sqrt3
            dcint(6,i,3)= cint1(10,i)           -cint2(3,i)*two
          enddo
        case (3)
          do i= 1,nbfij(1)
            dcint( 1,i,1)= cint1( 1,i)            -cint2(1,i)*three
            dcint( 2,i,1)= cint1( 2,i)*sqrtseventh-cint2(2,i)*two*sqrtthird
            dcint( 3,i,1)= cint1( 3,i)*sqrtseventh-cint2(3,i)*two*sqrtthird
            dcint( 4,i,1)= cint1( 4,i)*sqrt3inv35 -cint2(4,i)
            dcint( 5,i,1)= cint1( 5,i)*sqrtinv35  -cint2(5,i)*sqrtthird
            dcint( 6,i,1)= cint1( 6,i)*sqrt3inv35 -cint2(6,i)
            dcint( 7,i,1)= cint1( 7,i)*sqrtseventh
            dcint( 8,i,1)= cint1( 8,i)*sqrtinv35
            dcint( 9,i,1)= cint1( 9,i)*sqrtinv35
            dcint(10,i,1)= cint1(10,i)*sqrtseventh
            dcint( 1,i,2)= cint1( 2,i)*sqrtseventh
            dcint( 2,i,2)= cint1( 4,i)*sqrt3inv35 -cint2(1,i)
            dcint( 3,i,2)= cint1( 5,i)*sqrtinv35
            dcint( 4,i,2)= cint1( 7,i)*sqrtseventh-cint2(2,i)*two*sqrtthird
            dcint( 5,i,2)= cint1( 8,i)*sqrtinv35  -cint2(3,i)*sqrtthird
            dcint( 6,i,2)= cint1( 9,i)*sqrtinv35
            dcint( 7,i,2)= cint1(11,i)            -cint2(4,i)*three
            dcint( 8,i,2)= cint1(12,i)*sqrtseventh-cint2(5,i)*two*sqrtthird
            dcint( 9,i,2)= cint1(13,i)*sqrt3inv35 -cint2(6,i)
            dcint(10,i,2)= cint1(14,i)*sqrtseventh
            dcint( 1,i,3)= cint1( 3,i)*sqrtseventh
            dcint( 2,i,3)= cint1( 5,i)*sqrtinv35
            dcint( 3,i,3)= cint1( 6,i)*sqrt3inv35 -cint2(1,i)
            dcint( 4,i,3)= cint1( 8,i)*sqrtinv35
            dcint( 5,i,3)= cint1( 9,i)*sqrtinv35  -cint2(2,i)*sqrtthird
            dcint( 6,i,3)= cint1(10,i)*sqrtseventh-cint2(3,i)*two*sqrtthird
            dcint( 7,i,3)= cint1(12,i)*sqrtseventh
            dcint( 8,i,3)= cint1(13,i)*sqrt3inv35 -cint2(4,i)
            dcint( 9,i,3)= cint1(14,i)*sqrtseventh-cint2(5,i)*two*sqrtthird
            dcint(10,i,3)= cint1(15,i)            -cint2(6,i)*three
          enddo
        case (4)
          do i= 1,nbfij(1)
            dcint( 1,i,1)= cint1( 1,i)           -cint2( 1,i)*four
            dcint( 2,i,1)= cint1( 2,i)*third     -cint2( 2,i)*three*sqrtfifth
            dcint( 3,i,1)= cint1( 3,i)*third     -cint2( 3,i)*three*sqrtfifth
            dcint( 4,i,1)= cint1( 4,i)*sqrtinv21 -cint2( 4,i)*two*sqrtfifth
            dcint( 5,i,1)= cint1( 5,i)*sqrtinv63 -cint2( 5,i)*two*sqrtinv15
            dcint( 6,i,1)= cint1( 6,i)*sqrtinv21 -cint2( 6,i)*two*sqrtfifth
            dcint( 7,i,1)= cint1( 7,i)*sqrtinv21 -cint2( 7,i)
            dcint( 8,i,1)= cint1( 8,i)*sqrtinv105-cint2( 8,i)*sqrtfifth
            dcint( 9,i,1)= cint1( 9,i)*sqrtinv105-cint2( 9,i)*sqrtfifth
            dcint(10,i,1)= cint1(10,i)*sqrtinv21 -cint2(10,i)
            dcint(11,i,1)= cint1(11,i)*third
            dcint(12,i,1)= cint1(12,i)*sqrtinv63
            dcint(13,i,1)= cint1(13,i)*sqrtinv105
            dcint(14,i,1)= cint1(14,i)*sqrtinv63
            dcint(15,i,1)= cint1(15,i)*third
            dcint( 1,i,2)= cint1( 2,i)*third
            dcint( 2,i,2)= cint1( 4,i)*sqrtinv21 -cint2( 1,i)
            dcint( 3,i,2)= cint1( 5,i)*sqrtinv63
            dcint( 4,i,2)= cint1( 7,i)*sqrtinv21 -cint2( 2,i)*two*sqrtfifth
            dcint( 5,i,2)= cint1( 8,i)*sqrtinv105-cint2( 3,i)*sqrtfifth
            dcint( 6,i,2)= cint1( 9,i)*sqrtinv105
            dcint( 7,i,2)= cint1(11,i)*third     -cint2( 4,i)*three*sqrtfifth
            dcint( 8,i,2)= cint1(12,i)*sqrtinv63 -cint2( 5,i)*two*sqrtinv15
            dcint( 9,i,2)= cint1(13,i)*sqrtinv105-cint2( 6,i)*sqrtfifth
            dcint(10,i,2)= cint1(14,i)*sqrtinv63
            dcint(11,i,2)= cint1(16,i)           -cint2( 7,i)*four
            dcint(12,i,2)= cint1(17,i)*third     -cint2( 8,i)*three*sqrtfifth
            dcint(13,i,2)= cint1(18,i)*sqrtinv21 -cint2( 9,i)*two*sqrtfifth
            dcint(14,i,2)= cint1(19,i)*sqrtinv21 -cint2(10,i)
            dcint(15,i,2)= cint1(20,i)*third
            dcint( 1,i,3)= cint1( 3,i)*third
            dcint( 2,i,3)= cint1( 5,i)*sqrtinv63
            dcint( 3,i,3)= cint1( 6,i)*sqrtinv21 -cint2( 1,i)
            dcint( 4,i,3)= cint1( 8,i)*sqrtinv105
            dcint( 5,i,3)= cint1( 9,i)*sqrtinv105-cint2( 2,i)*sqrtfifth
            dcint( 6,i,3)= cint1(10,i)*sqrtinv21 -cint2( 3,i)*two*sqrtfifth
            dcint( 7,i,3)= cint1(12,i)*sqrtinv63
            dcint( 8,i,3)= cint1(13,i)*sqrtinv105-cint2( 4,i)*sqrtfifth
            dcint( 9,i,3)= cint1(14,i)*sqrtinv63 -cint2( 5,i)*two*sqrtinv15
            dcint(10,i,3)= cint1(15,i)*third     -cint2( 6,i)*three*sqrtfifth
            dcint(11,i,3)= cint1(17,i)*third
            dcint(12,i,3)= cint1(18,i)*sqrtinv21 -cint2( 7,i)
            dcint(13,i,3)= cint1(19,i)*sqrtinv21 -cint2( 8,i)*two*sqrtfifth
            dcint(14,i,3)= cint1(20,i)*third     -cint2( 9,i)*three*sqrtfifth
            dcint(15,i,3)= cint1(21,i)           -cint2(10,i)*four
          enddo
        case (5)
          do i= 1,nbfij(1)
            dcint( 1,i,1)= cint1( 1,i)            -cint2( 1,i)*five
            dcint( 2,i,1)= cint1( 2,i)*sqrtinv11  -cint2( 2,i)*four*sqrtseventh
            dcint( 3,i,1)= cint1( 3,i)*sqrtinv11  -cint2( 3,i)*four*sqrtseventh
            dcint( 4,i,1)= cint1( 4,i)*sqrtinv33  -cint2( 4,i)*three*sqrt3inv35
            dcint( 5,i,1)= cint1( 5,i)*sqrtinv99  -cint2( 5,i)*three*sqrtinv35
            dcint( 6,i,1)= cint1( 6,i)*sqrtinv33  -cint2( 6,i)*three*sqrt3inv35
            dcint( 7,i,1)= cint1( 7,i)*sqrt5inv231-cint2( 7,i)*two*sqrtseventh
            dcint( 8,i,1)= cint1( 8,i)*sqrtinv231 -cint2( 8,i)*two*sqrtinv35
            dcint( 9,i,1)= cint1( 9,i)*sqrtinv231 -cint2( 9,i)*two*sqrtinv35
            dcint(10,i,1)= cint1(10,i)*sqrt5inv231-cint2(10,i)*two*sqrtseventh
            dcint(11,i,1)= cint1(11,i)*sqrtinv33  -cint2(11,i)
            dcint(12,i,1)= cint1(12,i)*sqrtinv231 -cint2(12,i)*sqrtseventh
            dcint(13,i,1)= cint1(13,i)*sqrtinv385 -cint2(13,i)*sqrt3inv35
            dcint(14,i,1)= cint1(14,i)*sqrtinv231 -cint2(14,i)*sqrtseventh
            dcint(15,i,1)= cint1(15,i)*sqrtinv33  -cint2(15,i)
            dcint(16,i,1)= cint1(16,i)*sqrtinv11
            dcint(17,i,1)= cint1(17,i)*sqrtinv99
            dcint(18,i,1)= cint1(18,i)*sqrtinv231
            dcint(19,i,1)= cint1(19,i)*sqrtinv231
            dcint(20,i,1)= cint1(20,i)*sqrtinv99
            dcint(21,i,1)= cint1(21,i)*sqrtinv11
            dcint( 1,i,2)= cint1( 2,i)*sqrtinv11
            dcint( 2,i,2)= cint1( 4,i)*sqrtinv33  -cint2( 1,i)
            dcint( 3,i,2)= cint1( 5,i)*sqrtinv99
            dcint( 4,i,2)= cint1( 7,i)*sqrt5inv231-cint2( 2,i)*two*sqrtseventh
            dcint( 5,i,2)= cint1( 8,i)*sqrtinv231 -cint2( 3,i)*sqrtseventh
            dcint( 6,i,2)= cint1( 9,i)*sqrtinv231
            dcint( 7,i,2)= cint1(11,i)*sqrtinv33  -cint2( 4,i)*three*sqrt3inv35
            dcint( 8,i,2)= cint1(12,i)*sqrtinv231 -cint2( 5,i)*two*sqrtinv35
            dcint( 9,i,2)= cint1(13,i)*sqrtinv385 -cint2( 6,i)*sqrt3inv35
            dcint(10,i,2)= cint1(14,i)*sqrtinv231
            dcint(11,i,2)= cint1(16,i)*sqrtinv11  -cint2( 7,i)*four*sqrtseventh
            dcint(12,i,2)= cint1(17,i)*sqrtinv99  -cint2( 8,i)*three*sqrtinv35
            dcint(13,i,2)= cint1(18,i)*sqrtinv231 -cint2( 9,i)*two*sqrtinv35
            dcint(14,i,2)= cint1(19,i)*sqrtinv231 -cint2(10,i)*sqrtseventh
            dcint(15,i,2)= cint1(20,i)*sqrtinv99
            dcint(16,i,2)= cint1(22,i)            -cint2(11,i)*five
            dcint(17,i,2)= cint1(23,i)*sqrtinv11  -cint2(12,i)*four*sqrtseventh
            dcint(18,i,2)= cint1(24,i)*sqrtinv33  -cint2(13,i)*three*sqrt3inv35
            dcint(19,i,2)= cint1(25,i)*sqrt5inv231-cint2(14,i)*two*sqrtseventh
            dcint(20,i,2)= cint1(26,i)*sqrtinv33  -cint2(15,i)
            dcint(21,i,2)= cint1(27,i)*sqrtinv11
            dcint( 1,i,3)= cint1( 3,i)*sqrtinv11
            dcint( 2,i,3)= cint1( 5,i)*sqrtinv99
            dcint( 3,i,3)= cint1( 6,i)*sqrtinv33  -cint2( 1,i)
            dcint( 4,i,3)= cint1( 8,i)*sqrtinv231
            dcint( 5,i,3)= cint1( 9,i)*sqrtinv231 -cint2( 2,i)*sqrtseventh
            dcint( 6,i,3)= cint1(10,i)*sqrt5inv231-cint2( 3,i)*two*sqrtseventh
            dcint( 7,i,3)= cint1(12,i)*sqrtinv231
            dcint( 8,i,3)= cint1(13,i)*sqrtinv385 -cint2( 4,i)*sqrt3inv35
            dcint( 9,i,3)= cint1(14,i)*sqrtinv231 -cint2( 5,i)*two*sqrtinv35
            dcint(10,i,3)= cint1(15,i)*sqrtinv33  -cint2( 6,i)*three*sqrt3inv35
            dcint(11,i,3)= cint1(17,i)*sqrtinv99
            dcint(12,i,3)= cint1(18,i)*sqrtinv231 -cint2( 7,i)*sqrtseventh
            dcint(13,i,3)= cint1(19,i)*sqrtinv231 -cint2( 8,i)*two*sqrtinv35
            dcint(14,i,3)= cint1(20,i)*sqrtinv99  -cint2( 9,i)*three*sqrtinv35
            dcint(15,i,3)= cint1(21,i)*sqrtinv11  -cint2(10,i)*four*sqrtseventh
            dcint(16,i,3)= cint1(23,i)*sqrtinv11
            dcint(17,i,3)= cint1(24,i)*sqrtinv33  -cint2(11,i)
            dcint(18,i,3)= cint1(25,i)*sqrt5inv231-cint2(12,i)*two*sqrtseventh
            dcint(19,i,3)= cint1(26,i)*sqrtinv33  -cint2(13,i)*three*sqrt3inv35
            dcint(20,i,3)= cint1(27,i)*sqrtinv11  -cint2(14,i)*four*sqrtseventh
            dcint(21,i,3)= cint1(28,i)            -cint2(15,i)*five
          enddo
      end select
!
      nbfij(2)  = mbf(jsh)
      select case(nbfij(2))
        case(5)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,6
                work(j)= dcint(j,i,k)
              enddo
              dcint(1,i,k)= work(2)
              dcint(2,i,k)= work(5)
              dcint(3,i,k)=(work(6)*two-work(1)-work(4))*half
              dcint(4,i,k)= work(3)
              dcint(5,i,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        case(7)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,10
                work(j)= dcint(j,i,k)
              enddo
              dcint(1,i,k)=(-work(7)+work(2)*three                   )*facf1
              dcint(2,i,k)=  work(5)                                  *facf2
              dcint(3,i,k)=(-work(7)-work(2)+work(9)*four            )*facf3
              dcint(4,i,k)=( work(10)*two-work(3)*three-work(8)*three)*half
              dcint(5,i,k)=(-work(1)-work(4)+work(6)*four            )*facf3
              dcint(6,i,k)=( work(3)-work(8)                         )*facf4
              dcint(7,i,k)=( work(1)-work(4)*three                   )*facf1
            enddo
          enddo
        case(10)
          do k= 1,3
            do i= 1,nbfij(1)
              dcint(2,i,k)= dcint(2,i,k)*sqrt5
              dcint(3,i,k)= dcint(3,i,k)*sqrt5
              dcint(4,i,k)= dcint(4,i,k)*sqrt5
              dcint(5,i,k)= dcint(5,i,k)*sqrt15
              dcint(6,i,k)= dcint(6,i,k)*sqrt5
              dcint(8,i,k)= dcint(8,i,k)*sqrt5
              dcint(9,i,k)= dcint(9,i,k)*sqrt5
            enddo
          enddo
        case(9)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,15
                work(j)= dcint(j,i,k)
              enddo
              dcint(1,i,k)=(work(2)-work(7))*facg1
              dcint(2,i,k)=(-work(12)+work(5)*three)*facg2
              dcint(3,i,k)=(-work(2)-work(7)+work(9)*six)*facg3
              dcint(4,i,k)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              dcint(5,i,k)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                          -work(6)*p24-work(13)*p24)*eighth
              dcint(6,i,k)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              dcint(7,i,k)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
              dcint(8,i,k)=(work(3)-work(8)*three)*facg2
              dcint(9,i,k)=(work(1)+work(11)-work(4)*six)*facg6
            enddo
          enddo
        case(15)
          do k= 1,3
            do i= 1,nbfij(1)
              dcint( 2,i,k)= dcint( 2,i,k)*sqrt7
              dcint( 3,i,k)= dcint( 3,i,k)*sqrt7
              dcint( 4,i,k)= dcint( 4,i,k)*sqrt35third
              dcint( 5,i,k)= dcint( 5,i,k)*sqrt35
              dcint( 6,i,k)= dcint( 6,i,k)*sqrt35third
              dcint( 7,i,k)= dcint( 7,i,k)*sqrt7
              dcint( 8,i,k)= dcint( 8,i,k)*sqrt35
              dcint( 9,i,k)= dcint( 9,i,k)*sqrt35
              dcint(10,i,k)= dcint(10,i,k)*sqrt7
              dcint(12,i,k)= dcint(12,i,k)*sqrt7
              dcint(13,i,k)= dcint(13,i,k)*sqrt35third
              dcint(14,i,k)= dcint(14,i,k)*sqrt7
            enddo
          enddo
        case(11)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,21
                work(j)= dcint(j,i,k)
              enddo
              dcint( 1,i,k)=(work(2)*five-work(7)*ten+work(16))*fach1
              dcint( 2,i,k)=(work(5)*four-work(12)*four)*fach2
              dcint( 3,i,k)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                           -work(18)*eight)*fach3
              dcint( 4,i,k)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
              dcint( 5,i,k)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                           +work(20)*eight)*fach5
              dcint( 6,i,k)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15 &
&                           -work(19)*p40+work(21)*eight)*eighth
              dcint( 7,i,k)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                           +work(15)*eight)*fach5
              dcint( 8,i,k)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
              dcint( 9,i,k)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                           -work(13)*p24)*fach3
              dcint(10,i,k)=(work(3)-work(8)*six+work(17))*fach2
              dcint(11,i,k)=(work(1)-work(4)*ten+work(11)*five)*fach1
            enddo
          enddo
        case(21)
          do k= 1,3
            do i= 1,nbfij(1)
              dcint( 2,i,k)= dcint( 2,i,k)*three
              dcint( 3,i,k)= dcint( 3,i,k)*three
              dcint( 4,i,k)= dcint( 4,i,k)*sqrt21
              dcint( 5,i,k)= dcint( 5,i,k)*sqrt63
              dcint( 6,i,k)= dcint( 6,i,k)*sqrt21
              dcint( 7,i,k)= dcint( 7,i,k)*sqrt21
              dcint( 8,i,k)= dcint( 8,i,k)*sqrt105
              dcint( 9,i,k)= dcint( 9,i,k)*sqrt105
              dcint(10,i,k)= dcint(10,i,k)*sqrt21
              dcint(11,i,k)= dcint(11,i,k)*three
              dcint(12,i,k)= dcint(12,i,k)*sqrt63
              dcint(13,i,k)= dcint(13,i,k)*sqrt105
              dcint(14,i,k)= dcint(14,i,k)*sqrt63
              dcint(15,i,k)= dcint(15,i,k)*three
              dcint(17,i,k)= dcint(17,i,k)*three
              dcint(18,i,k)= dcint(18,i,k)*sqrt21
              dcint(19,i,k)= dcint(19,i,k)*sqrt21
              dcint(20,i,k)= dcint(20,i,k)*three
            enddo
          enddo
      end select
!
      do i= 1,nbfij(1)
        ii= ilocbf+i
        do j= 1,nbfij(2)
          egrad(1,jatom)= egrad(1,jatom)+fulldmtrx(jlocbf+j,ii)*dcint(j,i,1)
          egrad(2,jatom)= egrad(2,jatom)+fulldmtrx(jlocbf+j,ii)*dcint(j,i,2)
          egrad(3,jatom)= egrad(3,jatom)+fulldmtrx(jlocbf+j,ii)*dcint(j,i,3)
        enddo
      enddo
      return
end


