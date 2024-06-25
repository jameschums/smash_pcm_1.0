!-------------------------
  subroutine setdefault2
!-------------------------
!
! Reset defaults after reading input file
!
      use modparallel, only : master
      use modthresh, only : precision, cutint2, threshweight, threshrho, threshdfock, threshdftao
      use moddft, only : nrad, nleb
      use modscf, only : dconv
      implicit none
      real(8),parameter :: zero= 0.0D+00

      select case(precision)
        case('HIGH')
          if(cutint2 < zero) cutint2= 1.0D-12
          if(dconv   < zero) dconv  = 5.0D-06
          if(threshweight < zero) threshweight=1.0D-08
          if(threshrho    < zero) threshrho   =1.0D-06
          if(threshdfock  < zero) threshdfock =1.0D-05
          if(threshdftao  < zero) threshdftao =1.0D-04
          if(nrad == 0) nrad= 150
          if(nleb == 0) nleb= 590
        case('MEDIUM')
          if(cutint2 < zero) cutint2= 1.0D-11
          if(dconv   < zero) dconv  = 5.0D-06
          if(threshweight < zero) threshweight=1.0D-08
          if(threshrho    < zero) threshrho   =1.0D-05
          if(threshdfock  < zero) threshdfock =1.0D-04
          if(threshdftao  < zero) threshdftao =1.0D-03
          if(nrad == 0) nrad= 96
          if(nleb == 0) nleb= 302
        case('LOW')
          if(cutint2 < zero) cutint2= 1.0D-10
          if(dconv   < zero) dconv  = 1.0D-05
          if(threshweight < zero) threshweight=1.0D-08
          if(threshrho    < zero) threshrho   =1.0D-04
          if(threshdfock  < zero) threshdfock =1.0D-04
          if(threshdftao  < zero) threshdftao =1.0D-02
          if(nrad == 0) nrad= 72
          if(nleb == 0) nleb= 302
        case default
          if(master) write(*,'(" Error! This program does not support precision= ", &
&                          a16,".")') precision
          call iabort
      end select
!
      return
end

