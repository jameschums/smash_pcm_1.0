!---------------------------------------------------------------------------------------
  subroutine calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work1, &
&                            work2,work3,work4,workv,iopt,iredun,isizered, &
&                            maxredun,numbond,numangle,numtorsion,numredun, &
&                            nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Calculate new Cartesian coordinate with gradient, hessian and redundant coordinate
! using Rational Function Optimization (RFO) method
!
! Parameters for force constants (bond strech (up to 3rd period), angle bend, torsion):
!                          H. B. Schlegel, Theoret. Chim. Acta, 333 (1984) 66.
! Parameters for force constants (bond strech (4th - 6th period)):
!                          J. M. Wittbrodt, H. B. Schlegel, J. Mol. Strut., 398 (1997) 55.
! Covalent raddi (H - Kr): H. B. Schlegel, Theoret. Chim. Acta, 333 (1984) 66.
! Covalent radii (Rb- Cn): P. Pyykko, M. Atsumi, Chem. Eur. J., 186 (2009) 15.
!
      use modparallel, only : master
      use modprint, only : iprint
      use modunit, only : toang, tobohr
      use modmolecule, only : numatomic, natom
      implicit none
      integer,parameter :: maxiterdx=100, maxiterrfo=10000
      integer,intent(in) :: iopt, isizered, maxredun, iredun(4,isizered/4)
      integer,intent(in) :: numbond, numangle, numtorsion, numredun, nproc, myrank, mpi_comm
      integer :: irow(112)
      integer :: natom3, ii, jj, ij, kk, iatom, jatom, katom, iterrfo, iterdx
      integer :: numdim
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, convl=1.0D-08, convrms=1.0D-06
      real(8),parameter :: rad2deg=5.729577951308232D+01
      real(8),intent(inout) :: coord(natom*3), coordold(natom*3), coordredun(numredun,2)
      real(8),intent(inout) :: egrad(natom*3), egradredun(numredun,2)
      real(8),intent(inout) :: ehess(numredun*(numredun+1)/2), work1(maxredun,maxredun)
      real(8),intent(inout) :: work2(maxredun,maxredun), work3(maxredun,maxredun)
      real(8),intent(inout) :: work4(maxredun,maxredun), workv(maxredun,3)
      real(8) :: parambond(7,7), rij, paramb
      real(8) :: rcov(112), rjk, suml, rlambda, rmsdx, rmsqx, ddot
      character(len=33) :: paramred
      character(len=5) :: chartmp(5)
      character(len=3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data irow/2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 26*7/
      data parambond/ &
&     -0.2440D+00, 0.3520D+00, 0.6600D+00, 0.7126D+00, 0.8335D+00, 0.9491D+00, 1.0D+00, &
&      0.3520D+00, 1.0850D+00, 1.5220D+00, 1.4725D+00, 1.6549D+00, 1.7190D+00, 1.8D+00, &
&      0.6600D+00, 1.5220D+00, 2.0680D+00, 1.8238D+00, 2.1164D+00, 2.3185D+00, 2.5D+00, &
&      0.7126D+00, 1.4725D+00, 1.8238D+00, 2.0203D+00, 2.2137D+00, 2.5206D+00, 2.8D+00, &
&      0.8335D+00, 1.6549D+00, 2.1164D+00, 2.2137D+00, 2.3718D+00, 2.5110D+00, 2.7D+00, &
&      0.9491D+00, 1.7190D+00, 2.3185D+00, 2.5206D+00, 2.5110D+00, 2.5D+00,    2.5D+00, &
&      1.0D+00,    1.8D+00,    2.5D+00,    2.8D+00,    2.7D+00,    2.5D+00,    2.5D+00/
      data rcov/ &
&     0.32D+00, 0.60D+00, 1.20D+00, 1.05D+00, 0.81D+00, 0.77D+00, 0.74D+00, 0.74D+00, 0.72D+00, &
&     0.72D+00, 1.50D+00, 1.40D+00, 1.30D+00, 1.17D+00, 1.10D+00, 1.04D+00, 0.99D+00, 0.99D+00, &
&     1.80D+00, 1.60D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, &
&     1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.30D+00, 1.20D+00, 1.20D+00, 1.10D+00, 1.10D+00, &
&     2.10D+00, 1.85D+00, 1.63D+00, 1.54D+00, 1.47D+00, 1.38D+00, 1.28D+00, 1.25D+00, 1.25D+00, &
&     1.20D+00, 1.28D+00, 1.36D+00, 1.42D+00, 1.40D+00, 1.40D+00, 1.36D+00, 1.33D+00, 1.31D+00, &
&     2.32D+00, 1.96D+00, 1.80D+00, 1.63D+00, 1.76D+00, 1.74D+00, 1.73D+00, 1.72D+00, 1.68D+00, &
&     1.69D+00, 1.68D+00, 1.67D+00, 1.66D+00, 1.65D+00, 1.64D+00, 1.70D+00, 1.62D+00, 1.52D+00, &
&     1.46D+00, 1.37D+00, 1.31D+00, 1.29D+00, 1.22D+00, 1.23D+00, 1.24D+00, 1.33D+00, 1.44D+00, &
&     1.44D+00, 1.51D+00, 1.45D+00, 1.47D+00, 1.42D+00, 2.23D+00, 2.01D+00, 1.86D+00, 1.75D+00, &
&     1.69D+00, 1.70D+00, 1.71D+00, 1.72D+00, 1.66D+00, 1.66D+00, 1.68D+00, 1.68D+00, 1.65D+00, &
&     1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00, 1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, &
&     1.29D+00, 1.28D+00, 1.21D+00, 1.22D+00/
!
      natom3= natom*3
!
! Calculate B-matrix
!
      call calcbmatrix(coordredun,coord,work1,iredun,numbond,numangle,numtorsion)
!
! Calculate G=B*Bt
!
      call dgemm('N','T',numredun,numredun,natom3,one,work1,numredun,work1,numredun, &
&                zero,work2,maxredun)
!
! Calculate G-inverse and K-matrix
!
      call diag('V','U',numredun,work2,maxredun,workv(1,3),nproc,myrank,mpi_comm)
      numdim=0
      do ii= 1,numredun
        if(abs(workv(ii,3)) >= 1.0D-5) then
          numdim= numdim+1
          workv(ii,3)= one/workv(ii,3)
        endif
      enddo
!$OMP parallel do
      do ii= numredun-numdim+1,numredun
        do jj= 1,numredun
          work3(jj,ii)= work2(jj,ii)*workv(ii,3)
        enddo
      enddo
!$OMP end parallel do
      call dgemm('N','T',numredun,numredun,numdim,one,work2(1,numredun-numdim+1),maxredun, &
&                work3(1,numredun-numdim+1),maxredun,zero,work4,maxredun)
!
! Calculate (G-inverse)*B
!
      call dgemm('N','N',numredun,natom3,numredun,one,work4,maxredun,work1,numredun, &
&                zero,work3,maxredun)
!
! Calculate Fred=(G-inverse)*B*Fcart
!
      call dgemv('N',numredun,natom3,one,work3,maxredun,egrad,1,zero,egradredun,1)
!
      if(iopt >= 2) then
!
! Update Hessian
!
        call hessianbfgsred(ehess,coordredun,coordredun(1,2),egradredun,egradredun(1,2), &
&                           workv,numredun,numbond,numangle,numtorsion)
!
        do ii= 1,numredun
          ij= ii*(ii-1)/2
          do jj= 1,ii
            work4(jj,ii)=ehess(ij+jj)
          enddo
        enddo
      else
!
! Set initial Hessian
!
        do ii= 1,numbond
          iatom= numatomic(iredun(1,ii))
          jatom= numatomic(iredun(2,ii))
          rij= coordredun(ii,1)
          paramb= parambond(irow(iatom),irow(jatom))
          work4(ii,ii)= 1.734D+00/((rij-paramb)*(rij-paramb)*(rij-paramb))
        enddo
        do ii= numbond+1,numbond+numangle
          iatom= numatomic(iredun(1,ii))
          katom= numatomic(iredun(3,ii))
          if((iatom == 1).or.(katom == 1)) then
            work4(ii,ii)= 0.16D+00
          else
            work4(ii,ii)= 0.25D+00
          endif
        enddo
        do ii= numbond+numangle+1,numredun
          jatom= numatomic(iredun(2,ii))
          katom= numatomic(iredun(3,ii))
          jj=(iredun(2,ii)-1)*3
          kk=(iredun(3,ii)-1)*3
          rjk=sqrt((coord(jj+1)-coord(kk+1))**2+(coord(jj+2)-coord(kk+2))**2 &
&                 +(coord(jj+3)-coord(kk+3))**2)-(rcov(jatom)+rcov(katom))*tobohr
          if(rjk > zero) rjk= zero
          work4(ii,ii)= 0.0023D+00-0.07D+00*rjk
        enddo
!
        do ii= 1,numredun
          ij= ii*(ii-1)/2
          do jj= 1,ii-1
            work4(jj,ii)= zero
            ehess(ij+jj)= zero
          enddo
          ehess(ij+ii)= work4(ii,ii)
        enddo
      endif
!
! Copy coordinates and gradients
!
      do ii= 1,natom3
        coordold(ii)= coord(ii)
      enddo
      do ii= 1,numredun
        coordredun(ii,2)= coordredun(ii,1)
        egradredun(ii,2)= egradredun(ii,1)
      enddo
!
! Calculate Kt*FC*K and diagonalize, and then obtain eigenvalues and orthogonal basis
!
      call dsymm('L','U',numredun,numdim,one,work4,maxredun,work2(1,numredun-numdim+1),maxredun, &
&                zero,work1,maxredun)
      call dgemm('T','N',numdim,numdim,numredun,one,work2(1,numredun-numdim+1),maxredun, &
&                work1,maxredun,zero,work4,maxredun)
      call diag('V','U',numdim,work4,maxredun,workv(1,3),nproc,myrank,mpi_comm)
!
! Calculate K*(orthogonal basis)
!
      call dgemm('N','N',numredun,numdim,numdim,one,work2(1,numredun-numdim+1),maxredun, &
&                work4,maxredun,zero,work1,maxredun)
!
! Calculate (K*(orthogonal basis))t*Fred
!
      call dgemv('T',numredun,numdim,one,work1,maxredun,egradredun,1,zero,workv,1)
!
! Rational Function Optimization (RFO) step
!
      rlambda= zero
      do iterrfo= 1,maxiterrfo
        suml= zero
        do ii= 1,numdim
          suml= suml+workv(ii,1)*workv(ii,1)/(rlambda-workv(ii,3))
        enddo
        if(master.and.(iprint >= 3)) then
          write(*,'(" Lambda iteration of RFO",i5,3x,"Lambda=",1p,d15.8,4x,"Sum=",1p,d15.8)') &
&               iterrfo,rlambda,suml
        endif
        if(abs(rlambda-suml) <= convl) exit
        rlambda= suml
        if(iterrfo == maxiterrfo) then
          if(master) write(*,'(" Error! RFOnot converge. ",i5," and conv ",d15.8," pants")') iterrfo,convl
          call iabort
        endif
      enddo
!
! Calculate displacement work(*,3) in coordinate where Hessian is diagonal
!
      do ii=1,numdim
        workv(ii,1)= workv(ii,1)/(rlambda-workv(ii,3))
      enddo
!
! Calculate displacement work(*,2) in redundant coordinate
!
      call dgemv('N',numredun,numdim,one,work1,maxredun,workv,1,zero,workv(1,2),1)
!
! Calculate displacement work(*,1) in Cartesian coordinate
!
      call dgemv('T',numredun,natom3,one,work3,maxredun,workv(1,2),1,zero,workv,1)
!
! Calculate next Cartesian coordinate
!
      do iterdx= 1,maxiterdx
!
! Update Catesian coordinate
!
        do ii= 1,natom3
          coord(ii)=coord(ii)+workv(ii,1)
        enddo
!
! Check convergence of displacement in Cartesian coordinate
!
        rmsdx= sqrt(ddot(natom3,workv,1,workv,1)/natom3)
        rmsqx= sqrt(ddot(numredun,workv(1,2),1,workv(1,2),1)/numredun)
        if(master.and.(iprint >= 3)) then
          write(*,'(" Displacement Iteration",i3,2x,"RMS(Cart)=",1p,d10.3,4x, &
&                   "RMS(Red)=",1p,d10.3)') iterdx, rmsdx, rmsqx
        endif
        if(rmsdx < convrms) exit
!
! Calculate B-matrix and new displacement in redundant coordinate
!
! delta-q workv(*,3)
!
        call calcbmatrix(workv,coord,work1,iredun,numbond,numangle,numtorsion)
        do ii= 1,numredun
          workv(ii,3)= workv(ii,1)-coordredun(ii,1)
        enddo
        call fixdtor(workv(1,3),numbond,numangle,numtorsion)
!
! delta(delta-q) workv(*,2)
!
        do ii= 1,numredun
          workv(ii,2)=workv(ii,2)-workv(ii,3)
        enddo
        call fixdtor(workv(1,2),numbond,numangle,numtorsion)
!
        do ii= 1,numredun
          coordredun(ii,1)= workv(ii,1)
        enddo
!
! Calculate G=B*Bt
!
        call dgemm('N','T',numredun,numredun,natom3,one,work1,numredun,work1,numredun, &
&                  zero,work2,maxredun)
!
!   Calculate G-inverse
!
        call diag('V','U',numredun,work2,maxredun,workv(1,1),nproc,myrank,mpi_comm)
        numdim=0
        do ii= 1,numredun
          if(abs(workv(ii,1)) >= 1.0D-5) then
            numdim= numdim+1
            workv(ii,1)= one/workv(ii,1)
          endif
        enddo
        do ii= numredun-numdim+1,numredun
          do jj= 1,numredun
            work3(jj,ii)= work2(jj,ii)*workv(ii,1)
          enddo
        enddo
        call dgemm('N','T',numredun,numredun,numdim,one,work2(1,numredun-numdim+1),maxredun, &
&                  work3(1,numredun-numdim+1),maxredun,zero,work4,maxredun)
!
! Calculate Bt*(G-inverse)
!
        call dgemm('T','N',natom3,numredun,numredun,one,work1,numredun,work4,maxredun, &
&                  zero,work3,maxredun)
!
! Calculate displacement work(*,1) in Cartesian coordinate
!
        call dgemv('N',natom3,numredun,one,work3,maxredun,workv(1,2),1,zero,workv,1)
!
        if(iterdx == maxiterdx) then
          if(master) then
            write(*,'(" Error! Transformation from redundant to Cartesian did not converge.")')
          endif
          call iabort
        endif
      enddo
      if(master) then
        if(iprint >= 3) write(*,*)
        write(*,'(" ---------------------------------------------------------------")')
        write(*,'("   Redundant coordinate parameters (Angstrom and Degree)")')
        write(*,'("                                        New           Old")')
        write(*,'(" ---------------------------------------------------------------")')
        do ii= 1,numbond
          write(chartmp(1:3),'(i5)')ii,iredun(1:2,ii)
          paramred= trim(trim("Bond"//adjustl(chartmp(1)) //"   ("//adjustl(chartmp(2)))//"," &
&                                  //adjustl(chartmp(3)))//")"
          write(*,'(3x,a33,f9.4,5x,f9.4)')paramred,coordredun(ii,1)*toang, &
&                                                  coordredun(ii,2)*toang
        enddo
        do ii= numbond+1,numbond+numangle
          write(chartmp(1:4),'(i5)')ii-numbond,iredun(1:3,ii)
          paramred= trim(trim(trim("Angle"//adjustl(chartmp(1)) //"  ("//adjustl(chartmp(2))) &
&                                    //","//adjustl(chartmp(3)))//","//adjustl(chartmp(4)))//")"
          write(*,'(3x,a33,f9.4,5x,f9.4)')paramred,coordredun(ii,1)*rad2deg, &
&                                                  coordredun(ii,2)*rad2deg
        enddo
        do ii= numbond+numangle+1,numbond+numangle+numtorsion
          write(chartmp(1:5),'(i5)')ii-numbond-numangle,iredun(1:4,ii)
          paramred= trim(trim(trim(trim("Torsion"//adjustl(chartmp(1)) //"("// &
&                   adjustl(chartmp(2)))//","//adjustl(chartmp(3)))//","// &
&                   adjustl(chartmp(4)))//","//adjustl(chartmp(5)))//")"
          write(*,'(3x,a33,f9.4,5x,f9.4)')paramred,coordredun(ii,1)*rad2deg, &
&                                                  coordredun(ii,2)*rad2deg
        enddo
        write(*,'(" ---------------------------------------------------------------")')
      endif

!
! Print delta xyz
!
      if(master.and.(iprint >= 2)) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Delta xyz (Angstrom)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do ii= 1,natom3/3
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(ii)), &
&              ((coord((ii-1)*3+jj)-coordold((ii-1)*3+jj))*toang,jj=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
      return
end

