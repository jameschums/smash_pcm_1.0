#'C:/Program Files (x86)/Intel/oneAPI/setvars.bat' 21may2024
##F90 = gfortran -DnoMPI -lm -llapack64 -fallow-argument-mismatch -fPIC -fallow-invalid-boz -frecursive -DILP64 -m64 -g
#F90 = nvfortran -DnoMPI -lm -llapack64 -fPIC -stdpar=gpu 
F90 = nvfortran -DnoMPI -lm -llapack64 -fPIC -fast -DILP64 -m64 -g -Wall #-stdpar=gpu -Minfo=accel
##gfortran  OPT = -fopenmp -fdefault-integer-8 # -fallow-argument-mismatch is needed for gfortran 10
#nvfortran F90  = nvfortran -DnoMPI #nvfortran
LIB  = -lm -llapack64 -lblas -lopenblas #nvfortran
OPT  = -mp -i8 -fast -fast -DILP64 #nvfortran -O3 -stdpar=gpu -Minfo=accel
#LIB =-lm  -llapack64 -lblas -lopenblas -fallow-invalid-boz -frecursive # high performance libraries are recommended #gfortran 
#OPT = -fopenmp -fdefault-integer-8 #gfortran
SRC =   basis.F90 geom.F90 guess.F90 fileio.F90 int1.F90 machine.F90 main.F90 math.F90 \
	memory.F90 scf.F90 scflib.F90 int2.F90 int2elec.F90 int2sp.F90 int2spd1.F90 \
	int2spd2.F90 int2spd3.F90 int2spd4.F90 gradient.F90 rysquad.F90 \
	dft.F90 dftfunc.F90 lebedev.F90 mp2.F90 mp2grad.F90 grad1.F90 grad2.F90 ecp.F90 \
	ecpfunc.F90 ecpder.F90 parallel.F90 setdefault.F90 setelectron.F90 \
	setparallel.F90 setdetails.F90 calcrenergy.F90 calcuenergy.F90 setdefault1.F90 \
	calcrgradient.F90 calcugradient.F90 calcrgeometry.F90 calcugeometry.F90 \
	setnextopt.F90 setdft.F90 setmp2.F90 tstamp.F90 parallelinfo.F90 para_waitall.F90 \
	para_irecvr.F90 para_isendr.F90 para_init.F90 para_sendrecvr.F90 \
	para_allgathervr.F90 para_allreducei.F90 para_allreducer.F90 calcbmatrix.F90 \
	readbasis.F90 calcnewcoordred.F90 calcnewcoord.F90 fixdtor.F90 writegeom.F90 \
	writebasis.F90 writeecp.F90 fullmtrx.F90 hessianbfgsred.F90 udftfockeri.F90 \
	hessianbfgs.F90 distarray.F90 expand.F90 ugmneri.F90 calcqcugmn.F90 \
	uhfqc.F90 setcharge.F90 rhfqc.F90 calcqcurmn.F90 gradoneei.F90 \
	calcuewdmtrx.F90 calcdkinetic.F90 calcewdmtrx.F90 calchelfey.F90 \
	calcdoverlap.F90 dghquad.F90 int1grys.F90 int1gcdd.F90 int1gcdp.F90 \
	int1gcds.F90 int1gcpp.F90 int1gcps.F90 mkrhs.F90 llgnew.F90 \
	forces_dd.F90 efld.F90 matvec.F90 cosmo.F90 jacobi_diis.F90 cosmomain.F90
MOD = module.F90 modulefmt.F90 modulerys.F90 ddcosmo.F90 
OBJ = $(addprefix obj/,$(SRC:.F90=.o))
OBJM= $(addprefix obj/,$(MOD:.F90=.o))
OBJS= $(OBJM) $(OBJ)
OBJDIR= obj
.SUFFIXES:
.SUFFIXES: .F .F90 .o

smashpcm : $(OBJS)
#$(F90) -o bin/smashpcm -I $(OBJDIR) $(OPT) $(OBJS) $(LIB)
	$(F90) -o bin/smashpcm $(OPT) $(OBJS) $(LIB)


smashpcm.a : $(OBJS)
	ar rv obj/smashpcm.a $(OBJS)

$(OBJDIR)/%.o : src/%.F90
#	$(F90) $(OPT) -J$(OBJDIR) -o $@ -c $<
	$(F90) $(OPT) -module $(OBJDIR) -o $@ -c $<
$(OBJ): $(addprefix src/,$(MOD))

clean :
	rm -f $(OBJDIR)/*.mod $(OBJDIR)/*.o bin/smashpcm molden.xyz input.dat* *.out Input.txt
