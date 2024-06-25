!-------------------------
  subroutine setdefault1
!-------------------------
!
! Set defaults before reading input file
!
      use modiofile, only : check
      use modwarn, only : nwarn
      use modguess, only : spher_g, guess
      use modmemory, only : memmax, memused, memusedmax, memory
      use modprint, only : iprint
      use modunit, only : bohr
      use modbasis, only : spher, basis
      use modscf, only : maxiter, dconv, fdiff, scfconv, maxdiis, maxsoscf, maxqc, &
&                        maxqcdiag, maxqcdiagsub, extrap
      use modthresh, only : precision, cutint2, threshsoscf, threshqc, threshover, threshatom, &
&                           threshdiis, threshweight, threshrho, threshdfock, threshdftao, &
&                           threshmp2cphf
      use moddft, only : idftex, idftcor, nrad, nleb, bqrad
      use modopt, only : nopt, optconv, cartesian
      use modecp, only : ecp, flagecp
      use modjob, only : scftype, runtype, method
      use modmolecule, only : multi, charge
      use modmp2, only : ncore, nvfz, maxmp2diis, maxmp2iter
      use modprop, only : octupole
      implicit none
!
      nwarn  = 0
      memmax = 1000000000
      memused= 0
      memusedmax= 0
      memory = ''
      maxiter= 150
      maxdiis= 20
      maxsoscf= 20
      maxqc   = 15
      maxqcdiag= 100
      maxqcdiagsub= 10
      fdiff  =.true.
      scfconv='DIIS'
      extrap =.false.
      threshsoscf= 0.25D+00
      threshqc   = 1.0D-05
      threshover = 1.0D-06
      threshatom = 2.0D-01
      threshdiis = 6.0D-01
      threshmp2cphf=1.0D-10
      idftex = 0
      idftcor= 0
      iprint = 2
      bohr   =.false.
      spher  =.true.
      spher_g=.true.
      nopt   = 100
      optconv= 1.0D-04
      cartesian=.false.
      multi  = 1
      charge = 0.0D+00
      bqrad(:)=1.0D+00
      nvfz= 0
      maxmp2diis= 20
      maxmp2iter= 100
!
      cutint2=-1.0d+00
      nrad = 0
      nleb = 0
      ncore= -1
      dconv=-1.0D+00
      threshweight=-1.0D+00
      threshrho=-1.0D+00
      threshdfock=-1.0D+00
      threshdftao=-1.0D+00
!
      precision='MEDIUM'
      flagecp= .false.
      scftype='RHF'
      method='HARTREE-FOCK'
      runtype='ENERGY'
      basis='STO-3G'
      guess='HUCKEL'
      ecp=''
      check=''
      octupole=.false.
!
      return
end

