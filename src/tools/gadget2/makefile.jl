makefilebody = """
OPTIONS =  \$(OPTIMIZE) \$(OPT) \n
EXEC   = Gadget2                \n
OBJS   = main.o  run.o  predict.o begrun.o endrun.o global.o  \\
	     timestep.o  init.o restart.o  io.o    \\
	     accel.o   read_ic.o  ngb.o  \\
	     system.o  allocate.o  density.o  \\
	     gravtree.o hydra.o  driftfac.o  \\
	     domain.o  allvars.o potential.o  \\
         forcetree.o   peano.o gravtree_forcetest.o \\
	     pm_periodic.o pm_nonperiodic.o longrange.o           \n
INCL   = allvars.h  proto.h  tags.h  Makefile                 \n
CFLAGS = \$(OPTIONS) \$(GSL_INCL) \$(FFTW_INCL) \$(HDF5INCL)  \n
ifeq (NOTYPEPREFIX_FFTW,\$(findstring NOTYPEPREFIX_FFTW,\$(OPT)))
  FFTW_LIB = \$(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,\$(findstring DOUBLEPRECISION_FFTW,\$(OPT)))
  FFTW_LIB = \$(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = \$(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif \n
LIBS   =   \$(HDF5LIB) -g  \$(MPICHLIB)  \$(GSL_LIBS) -lgsl -lgslcblas -lm \$(FFTW_LIB) \n 
\$(EXEC): \$(OBJS) 
\t\t\$(CC) \$(OBJS) \$(LIBS)   -o  \$(EXEC)  \n
\$(OBJS): \$(INCL) \n
clean:
\t\trm -f \$(OBJS) \$(EXEC)
"""

function write_gadget2_makefile(filename::String, opts...;
    defaultopts = ["PEANOHILBERT", "WALLCLOCK", "SYNCHRONIZATION", "NOSTOP_WHEN_BELOW_MINTIMESTEP", "DOUBLEPRECISION", "DOUBLEPRECISION_FFTW"],
    compiler = "mpicc",
    optimize = "-O3 -Wall",
    gsl_include = "-I/usr/local/include",
    gsl_lib = "-L/usr/local/lib",
    fftw_include = "-I/usr/local/include",
    fftw_lib = "-L/usr/local/lib",
    mpi_lib = "-L/usr/lib/x86_64-linux-gnu/openmpi/lib/ -lmpi",
    hdf5_include = "",
    hdf5_lib = "",
    makefilebody = makefilebody,

)
    f = open(filename, "w")

    makefileopt = join("OPT += -D" .* defaultopts .* "\n")
    makefileopt *= join("OPT += -D" .* opts .* "\n")

    makefilelib = """
    CC       = $(compiler)
    OPTIMIZE = $(optimize)
    GSL_INCL = $(gsl_include)
    GSL_LIBS = $(gsl_lib)
    FFTW_INCL= $(fftw_include)
    FFTW_LIBS= $(fftw_lib)
    MPICHLIB = $(mpi_lib)
    HDF5INCL = $(hdf5_include)
    HDF5LIB  = $(hdf5_lib)
    """

    printstyled("Compile Options:", color = :cyan)
    printstyled(makefileopt, color = :cyan)

    printstyled("Compiler Settings:", color = :light_green)
    printstyled(makefilelib, color = :light_green)

    write(f, makefileopt)
    write(f, makefilelib)
    write(f, makefilebody)
    close(f)
end