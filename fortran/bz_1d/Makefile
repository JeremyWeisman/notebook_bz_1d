# Makefile

FC = gfortran
LD = gfortran

FCFLAGS = -J ./obj -I/usr/local/miniconda3/envs/mumps_seq/include
LDFLAGS = 
LIBS = -L/usr/local/miniconda3/envs/mumps_seq/lib -ldmumps -lmumps_common -lpord -lmpiseq -lmetis -lblas -lesmumps -lscotch

EXE = run/bz_1d

OBJDIR = ./obj
RUNDIR = ./run

OBJS = obj/mod_precision.o \
       obj/mod_cartesian_grid.o \
       obj/mod_coo_matrix.o \
       obj/mod_utils.o \
       obj/mod_bz_2eq_1d.o \
       obj/mod_bz_3eq_1d.o \
       obj/radau5.o \
       obj/decsol.o \
       obj/dc_decsol.o \
       obj/mod_radau.o \
       obj/rock4.o \
       obj/rho.o \
       obj/mod_rock.o \
       obj/mod_strang.o \
       obj/mod_imex.o \
       obj/mod_integration.o \
       obj/bz_1d_main.o

default: $(OBJDIR) $(RUNDIR) $(EXE)

$(OBJDIR):
	@echo "Creation du repertoire $@"
	@mkdir -p $@
	@echo

$(RUNDIR):
	@echo "Creation du repertoire $@"
	@mkdir -p $@
	@echo

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

obj/%.o : src/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

obj/%.o : src/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

#dependencies

obj/mod_utils.o: obj/mod_precision.o obj/mod_cartesian_grid.o

obj/mo1_bz_2eq_1d.o: obj/mod_precision.o obj/mod_cartesian_grid.o obj/mod_csr_matrix.o

obj/mod_integration : obj/mod_precision.o obj/mod_radau.o obj/mod_strang.o obj/mod_bz_2eq_1d.o

obj/bz_2eq_1d_main.o: obj/mod_precision.o obj/mod_cartesian_grid.o obj/mod_utils.o obj/mod_bz_2eq_1d.o obj/mod_integration.o

clean:
	rm -f $(OBJS) obj/*.mod $(EXE) 

cleanall:
	rm -rf $(RUNDIR) $(OBJDIR) 
