!-----------------------
  subroutine readbasis
!-----------------------
!
! Read basis set
!
      use modparallel, only : master
      use modparam, only : mxprim, mxshell
      use modiofile, only : input, maxline
      use modbasis, only : exgen, coeffgen, locgenprim, mgenprim, mgentype, locgenshell, &
&                          ngenshell, atombasis
      implicit none
      integer :: ii, jj, iprim, ishell, ll, ielem(-9:112), nelem, kprim, numprim, natomshell
      character(len=3) :: element(-9:112)
      character(len=100) :: line
      character(len=16) :: symbol
      character(len=3) :: table(-9:112)= &
&     (/'BQ9','BQ8','BQ7','BQ6','BQ5','BQ4','BQ3','BQ2','BQ ','X  ',&
&       'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
&       'S  ','CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',&
&       'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ','RH ',&
&       'PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ',&
&       'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ',&
&       'OS ','IR ','PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',&
&       'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ','MD ','NO ','LR ','RF ','DB ',&
&       'SG ','BH ','HS ','MT ','UUN','UUU','UUB'/)
!
      atombasis(:)= ''
      ngenshell(:)= 0
      iprim= 0
      ishell= 0
!
      if(master) then
        rewind(input)
        do ii= 1,maxline
          read(input,'(a)',end=9999)line
          if(line(1:5) == "BASIS") exit
          if(ii == maxline) then
            write(*,'(" Error! Keyword BASIS is not found.")')
            call iabort
          endif
        enddo
!
        do ll= -9,112
          line=''
          read(input,'(a)',end=300)line
          if(len_trim(line) == 0) exit
          element(:)=''
!
! Read elements
!
          read(line,*,end=100)(element(ii),ii=-9,112)
 100      continue
          nelem= 0
!
! Check elements
!
          do ii= -9,112
            if((element(ii) == '0').or.(element(ii) == '')) exit
            if(element(ii) == 'BQ1') element(ii)= 'BQ'
            do jj= -9,112
              if(element(ii) == table(jj)) then
                ielem(ii)= jj
                nelem= nelem+1
                exit
              endif
              if(jj == 112) then
                write(*,'(" Error! This program does not support Atom",a3,".")')element(ii)
                call iabort
              endif
            enddo
          enddo
!
! Read basis functions
!
          natomshell= 0
          do ii= -9,nelem-10
            locgenshell(ielem(ii))= ishell
          enddo
          do jj= 1,maxline
            symbol= ''
            read(input,'(a)',err=200,end=200) line
            read(line,*,end=200,err=9998) symbol, numprim
            ishell= ishell+1
            natomshell= natomshell+1
            locgenprim(ishell)= iprim
            mgenprim(ishell)= numprim
            select case(symbol)
              case('S')
                mgentype(ishell)= 0
              case('P')
                mgentype(ishell)= 1
              case('D')
                mgentype(ishell)= 2
              case('F')
                mgentype(ishell)= 3
              case('G')
                mgentype(ishell)= 4
              case('H')
                mgentype(ishell)= 5
              case('I')
                mgentype(ishell)= 6
              case('SP')
                mgentype(ishell)  = 0
                mgentype(ishell+1)= 1
              case default
                write(*,'(" Error! The angular momentum ",a2," is not supported.")') symbol
                call iabort
            end select
            if(symbol /= 'SP') then
              do kprim= 1,numprim 
                iprim= iprim+1
                read(input,*,end=9998,err=9998) exgen(iprim), coeffgen(iprim)
              enddo
            else
              do kprim= 1,numprim 
                iprim= iprim+1
                read(input,*,end=9998,err=9998) exgen(iprim), coeffgen(iprim), &
&                                                             coeffgen(iprim+numprim)
                exgen(iprim+numprim)= exgen(iprim)
              enddo
              ishell= ishell+1
              locgenprim(ishell)= iprim
              mgenprim(ishell)= numprim
              iprim= iprim+numprim
              natomshell= natomshell+1
            endif
            cycle
!
 200        if(symbol(1:2) == '**') then
              do ii= -9,nelem-10
                ngenshell(ielem(ii))= natomshell
              enddo
              exit
            elseif(symbol == '') then
              write(*,'(" Error! End of basis functions is not found.")')
              call iabort
            endif
            do ii= -9,nelem-10
              atombasis(ielem(ii))= symbol
            enddo
          enddo
!
        enddo
 300    continue
      endif
      return
!
9999  write(*,'(" Error! Keyword BASIS is not found.")')
      call iabort
9998  write(*,'(" Error! Format of basis functions is incorrect.")')
      call iabort
end
