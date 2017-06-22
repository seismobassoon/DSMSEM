program = ../../bin/TraPSV 
objs =   modules.o trial.o rk3.o glu2.o dcsymbdl.o dcsymbdl3.o others.o calmat.o TraPSV-MPI.o
FC = ifort
#option = -Oss -noparallel -dcfns
#option =  -CB -traceback -g -assume byterecl -check all 
#option = -auto -ftrapuv -check all -warn all -fpe0 -traceback -parallel
#option = -O3 -assume byterecl  -I/usr/local/papi/include -L/usr/local/papi/lib -lpapi -check all -g -warn all 
#option = -O3 -assume byterecl -I/usr/local/papi/include -check all  
#option = -O3 -assume byterecl -check all  -g -warn all 
#option = -O2 -ftz -fpe0 -traceback -ftrapuv -vec-report0 -implicitnone -check all -assume byterecl -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage 
#option = -check nobounds -O3 -ftz -xSSE4.1 -fpe0 -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage
option = -check nobounds -O3 -ftz -xSSE4.1 -fpe0 -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage
.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(option)
.f90.o:
	$(FC) -o $@ -c $< $(option)

.PHONY: clean
clean:
	rm $(program) $(objs)  
