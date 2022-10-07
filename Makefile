#### define source file
SUFFIX = _test_eo
MAIN              = main
INIT              = init
INITSRC = init_art_radtest
MULTI             = multi_timestep
TIMESTEP          = timestep
PRIM              = primordial
RADSOURCE         = radiation_source_daisu
RANDOM            = mt19937
STOCHASTICSTELLAR = stochastic_stellar_evolv
ART               = art
CHEMISTRY         = chemistry
EXTERNF           = externalForce_art
GRID              = grid
FG2CG             = fg2cg
IO                = io
IOUTIL            = io_util
OUTPUTD = outputdata_art_radtest
BOUNDARY          = boundary_fixed
GRID_BOUNDARY     = grid_boundary
REFLUX            = reflux
STRING            = string
REFINE            = refine
REFINECOND = refineCond_nest_jeans
UNIT              = unit
KINZOKU           = kinzoku
KINZOKU2          = kinzoku2
RADTR             = radtr_anisotropy_eo
#RADTR_CHEM = radtr_chem
RADTR_CHEM        = radtr_chem_dustsubl_Tsubl
PARAM             = parameter
MODPAR = modelParameter_art_radtest
EOS               = eos_adiabatic_Kim_art
MPILIB            = mpilib
UTIL              = util
SYSTEMCALL        = systemcall_gnu
SORT              = sort
PACKARR           = packarr
DIFFUSE           = diffuse
OB                = overBlockCoordinates
OBINT             = ob_interp
UG                = uniformgrid
WS                = writeSnap
SC                = softCluster
SP      = sinkParticle_art
NANINF            = naninf
RESCUE            = rescue_art
ANALYSIS          = analysis_art
MEMSTAT           = memstat
##### FMG
FMGDATA	          = fmg_data
FMGDATAIF         = fmg_data_if
RSTRCT            = restriction
INTERP            = interpolation
FMGCONV           = fmg_converge
FMGBND            = fmg_boundary
FMGBNDP           = fmg_boundary_phys_periodic
FMGGHCL           = fmg_ghostcell
FMGREFLUX         = fmg_reflux
VMGINT            = vmg_interpol
FMGINTC           = fmg_interpol_cubic
MGDATA	          = mg_data
MGINTC            = mg_interpol_cubic
MG                = mg
VMG               = vmg
FMG               = fmg

##### auto-created header file
RECLH             = recl.h

##### HEALPix
HEALPIX_DIR       = ./Healpix
VPATH             =   .
VPATH             += $(HEALPIX_DIR)

HP_BIT            = bit_manipulation
HP_TYPES          = healpix_types
HP_INDMED         = indmed
HP_NUMREC         = num_rec
HP_STAT           = statistics
HP_EXT            = extension
HP_LONG           = long_intrinsic
HP_MISC           = misc_utils
HP_PIXTOOLS       = pix_tools
HP_CGETENV        = cgetEnvironment
OBJECT_HP         = $(HP_TYPES).o \
$(HP_CGETENV).o \
$(HP_EXT).o \
$(HP_MISC).o \
$(HP_INDMED).o \
$(HP_NUMREC).o \
$(HP_LONG).o \
$(HP_STAT).o \
$(HP_BIT).o \
$(HP_PIXTOOLS).o


##### intel Fortran + MPICH2
# FC	            = mpif90
# FFLAGS          = -traceback -g -warn all -check all -debug all
# FFLAGS          = -traceback -g -warn all -warn nointerfaces -check all -debug all -mcmodel=large -i-dynamic
# FFLAGS          = -u -O3 -shared-intel -mcmodel=large -fno-alias -fno-fnalias
# FFLAGS          = -u -O3 -shared-intel -mcmodel=large -fno-alias -fno-fnalias -opt-report-file$*.opt
# FFLAGS          = -u -CB -g
# CPPFLAGS        = -nostdinc
LDFLAGS	          =$(FFLAGS)

##### PGI Fortran + MPICH2
# FC	            = ftn
# FFLAGS          = -O3 -fastsse -Mlist -Mipa
# #FFLAGS         = -O3 -Mlist -Mipa
# CPPFLAGS        =
# LDFLAGS	        =$(FFLAGS)

##### HA8000
# FC	            = mpif90
# FFLAGS          = -Oss -noparallel -lf90c -w
# CPPFLAGS        =
# LDFLAGS	        =$(FFLAGS)

#### XC30 intel Fortran
CC                = cc
FC                = ftn
FFLAGS            = -O3 -mieee-fp
#FFLAGS = -traceback -g -warn all -warn nointerfaces -check all -debug all -check uninit
CPPFLAGS          = -nostdinc
LDFLAGS	          =$(FFLAGS)

ARCH              = ./Configs/objects.defs
include $(ARCH)

.SUFFIXES:
.SUFFIXES: .o .f90 .F90
.f90.o:; $(FC) $(FFLAGS) -c $<
.F90.o:; cat $< | ./precpp.pl | cpp -P -C $(CPPFLAGS) 2> /dev/null | ./postcpp.pl | ./column.pl > $*_cpp.f90; $(FC) -c $(FFLAGS) -o $*.o $*_cpp.f90 # ; rm $*_cpp.f90

####HEALPix compilation
%.o: $(HEALPIX_DIR)/%.F90
	$(FC) $(FFLAGS) -DNOCFITSIO -c $<
%.o: $(HEALPIX_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

all: $(MAIN) $(INIT)

OBJECT_MAIN       = $(OBJECT_HP) $(OBJECT) $(MAIN).o

BINDIR            = sfumato/ART_sublimation/bin/

$(MAIN): $(OBJECT_MAIN)
	$(FC) $(LDFLAGS) -o $(MAIN)$(SUFFIX) $(OBJECT_MAIN)
	mkdir -p /work/$$USER/$(BINDIR) ; \
	cp $(MAIN)$(SUFFIX) /work/$$USER/$(BINDIR)/$(MAIN)$(SUFFIX) ; \
	cp set_restartfiles.sh /work/$$USER/$(BINDIR) ; \
	cp set_restartfiles_step.sh /work/$$USER/$(BINDIR)


OBJECT_INIT       = $(OBJECT_HP) $(OBJECT) $(INITSRC).o

$(INIT): $(OBJECT_INIT)
	$(FC) $(LDFLAGS) -o $(INIT)$(SUFFIX) $(OBJECT_INIT)
	mkdir -p /work/$$USER/$(BINDIR) ; \
	cp $(INIT)$(SUFFIX) /work/$$USER/$(BINDIR)/$(INIT)$(SUFFIX)

$(MAIN).o: $< config.h

$(INITSRC).o: $< config.h

$(MULTI).o: $< config.h

$(TIMESTEP).o: $< config.h

$(BOUNDARY).o: $< config.h

$(GRID_BOUNDARY).o: $< config.h

$(REFINE).o: $< config.h

$(REFINECOND).o: $< config.h

$(GRID).o: $< config.h

$(IO).o: $< config.h

$(IOUTIL).o: $< config.h

$(OUTPUTD).o: $< config.h

$(EOS).o: $< config.h

$(PARAM).o: $< config.h

$(MODPAR).o: $< config.h

$(UNIT).o: $< config.h

$(KINZOKU).o: $< config.h

$(KINZOKU2).o: $< config.h

$(RADTR).o: $< config.h

$(RADTR_CHEM).o: $< config.h

$(MPILIB).o: $< config.h

$(UTIL).o: $< config.h

$(SYSTEMCALL).o: $< config.h

$(STRING).o: $< config.h

$(FG2CG).o: $< config.h

$(TIMESLICE).o: $< config.h

$(REFLUX).o: $< config.h

$(PACKARR).o: $< config.h

$(FMGDATA).o: $< config.h

$(RSTRCT).o: $< config.h

$(INTERP).o: $< config.h

$(FMGDIFF).o: $< config.h

$(FMGDATAIF).o: $< config.h

$(FMGCONV).o: $< config.h

$(FMGBND).o: $< config.h

$(FMGBNDP).o: $< config.h

$(FMGGHCL).o: $< config.h

$(FMGREFLUX).o: $< config.h

$(VMGINT).o: $< config.h

$(FMGINTC).o: $< config.h

$(MGINTC).o: $< config.h

$(MGDATA).o: $< config.h

$(MG).o: $< config.h mg_od.F90 mg_poisson.F90 mg_diffusion.F90 mg_ad.F90 debug_fmg.h

$(VMG).o: $< config.h vmg_od.F90 vmg_poisson.F90 vmg_diffusion.F90 vmg_ad.F90 debug_fmg.h

$(FMG).o: $< config.h fmg_test.F90 fmg_poisson.F90 fmg_diffusion.F90 fmg_od.F90 fmg_ad.F90 fmg_data_if.F90 debug_fmg.h

$(NANINF).o: $< config.h
	cat $*.F90 | ./precpp.pl | cpp -P -C $(CPPFLAGS) 2> /dev/null | ./postcpp.pl | ./column.pl > $*_cpp.f90; \
	$(FC) -c $(FFLAGS) -o $*.o $*_cpp.f90

$(RESCUE).o: $< config.h

$(DIFFUSE).o: $< config.h

$(ANALYSIS).o: $< config.h

$(MEMSTAT).o: $< config.h

$(OB).o : $< config.h

$(OBINT).o : $< config.h

$(UG).o : $< config.h # $(RECLH)

# $(RECLH): recl/recltest.h
# 	(cd recl; make)
# 	recl/make_headerfile.sh

$(WS).o : $< config.h

$(SP).o : $< config.h

$(SC).o : $< config.h

$(PRIM).o : $< config.h

$(STOCHASTICSTELLAR).o : $< config.h

$(RADSOURCE).o : $< config.h $(PRIM).o

$(ART).o : $< config.h

$(CHEMISTRY).o : $< config.h $(PRIM).o

healpix: $(OBJECT_HP)

clean:
	-rm $(MAIN)$(SUFFIX) $(INIT)$(SUFFIX) $(OBJECT) $(MAIN).o $(INITSRC).o *.mod $(OBJECT_HP) *_cpp.f90 # $(RECLH)
#	(cd recl; make clean)
