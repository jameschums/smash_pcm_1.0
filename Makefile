#'C:/Program Files (x86)/Intel/oneAPI/setvars.bat' 21may2024
F90 = gfortran -DnoMPI -lm -llapack64 -fallow-argument-mismatch -fPIC -fallow-invalid-boz -frecursive
#F90 = nvfortran -DnoMPI -lm -llapack64 -fPIC -stdpar=gpu
#gfortran  OPT = -fopenmp -fdefault-integer-8 # -fallow-argument-mismatch is needed for gfortran 10
#nvfortran F90  = nvfortran -DnoMPI
#nvfortran LIB  = -llapack -lblas
#OPT  = -mp -i8 -O3 #nvfortran
LIB =-lm  -llapack64 -lblas -lopenblas -fallow-invalid-boz -frecursive # high performance libraries are recommended
OPT = -fopenmp -fdefault-integer-8 #gfortran
SRC = basis.F90 geom.F90 guess.F90 fileio.F90 int1.F90 machine.F90 main.F90 math.F90 \
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
	int1gcds.F90 int1gcpp.F90 int1gcps.F90
MOD = module.F90 modulefmt.F90 modulerys.F90
OBJ = $(addprefix obj/,$(SRC:.F90=.o))
OBJM= $(addprefix obj/,$(MOD:.F90=.o))
OBJS= $(OBJM) $(OBJ)
OBJDIR= obj
.SUFFIXES:
.SUFFIXES: .F .F90 .o

smash : $(OBJS)
	$(F90) -o bin/smash -I $(OBJDIR) $(OPT) $(OBJS) $(LIB)
#nvfortran	$(F90) -o bin/smash $(OPT) $(OBJS) $(LIB)


smash.a : $(OBJS)
	ar rv obj/smash.a $(OBJS)

$(OBJDIR)/%.o : src/%.F90
	$(F90) $(OPT) -J$(OBJDIR) -o $@ -c $<
#nvfortran	$(F90) $(OPT) -module $(OBJDIR) -o $@ -c $<
$(OBJ): $(addprefix src/,$(MOD))

clean :
	rm -f $(OBJDIR)/*.mod $(OBJDIR)/*.o bin/smash *.xyz *.chk *.log input.dat* *.out
