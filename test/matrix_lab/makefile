# Build against the fork's CUDA + Hypre PETSc:
#   module purge; module load openmpi/4.1.3 cuda/12.3.2; make
# Paths are explicit: the install's own petscvariables point at /g/data1b, which no longer exists,
# so at run time LD_LIBRARY_PATH must also include $(P)/lib (libpetsc's rpath is stale too).
P := /g/data/au88/jh7060/ISSM/externalpackages/petsc/install-gpu-hypre
CUDA := /apps/cuda/12.3.2/lib64

solve_dumped: solve_dumped.c
	mpicc -O2 -Wall -I$(P)/include -o solve_dumped solve_dumped.c \
	  -L$(P)/lib -Wl,-rpath,$(P)/lib -Wl,-rpath-link,$(P)/lib -Wl,-rpath-link,$(CUDA) -lpetsc -lm
