program =../../TraSH 
#objs =    glu2.o  dclib.o dclib3.o others.o trial.o calmat.o SGTMaker-SH-0.8.1-rev.o
objs =    glu2.o  dclib.o dclib3.o others.o trial.o calmat.o TraSH-0.1.0.o
FC = ifort
#option = -Oss -noparallel -dcfns
#option =  -O -assum byterecl 
#option = -O2 -ftz -fpe0 -traceback -ftrapuv -vec-report0 -implicitnone -check all -assume byterecl -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage 
option = -check nobounds -O3 -ftz -xSSE4.1 -fpe0 -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage


.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(option)
.f90.o:
	$(FC) -o $@ -c $< $(option)

.PHONY: clean
clean:
	rm $(program) $(objs)  
