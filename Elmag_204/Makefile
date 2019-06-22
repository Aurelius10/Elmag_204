###############################################################################
# makefile for ELMAG 2.03
#
mpi: modules204.f90 bturb.f90 user204.f90 main_mpi204.f90 init203.f90 elmag204.f90 aux202.f90 dgl_numrec.f
	mpif90 -JModules -O -w modules204.f90 bturb.f90 user204.f90 main_mpi204.f90 init203.f90 elmag204.f90 aux202.f90 dgl_numrec.f
#
single: modules204.f90 bturb.f90 user204.f90 main_sp204.f90 init203.f90 elmag204.f90 aux202.f90 dgl_numrec.f
	gfortran -JModules -O modules204.f90 bturb.f90 user204.f90 main_sp204.f90 init203.f90 elmag204.f90 aux202.f90 dgl_numrec.f
#
test: modules204.f90 bturb.f90 user204.f90 main_sp204.f90 init203.f90 elmag204.f90 aux202.f90 dgl_numrec.f
	mpif90 -JModules  -Wall -C -g  -fbacktrace -ffpe-trap={underflow,overflow,invalid,denormal} modules203.f90 user203.f90 init203.f90 elmag203.f90 aux202.f90
#
###############################################################################
