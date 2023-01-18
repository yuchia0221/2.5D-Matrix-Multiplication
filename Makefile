EXEC=multiply
MATRIX_SIZE=7500
W = 32
OBJ = $(EXEC) $(EXEC)-debug $(EXEC)-trace

VIEWER=jumpshot

# flags
OPT=-O2 -g
DEBUG=-O0 -g

all: $(OBJ)

# build the debug parallel version of the program
$(EXEC)-debug: $(EXEC).cpp
	mpicxx $(DEBUG) $(OMP) -o $(EXEC)-debug $(EXEC).cpp -lrt


# build the version of the program that records MPE traces for jumpshot
$(EXEC)-trace: $(EXEC).cpp
	mpecxx -mpilog $(OPT) -o $(EXEC)-trace $(EXEC).cpp -lrt 

# build the optimized parallel version of the program
$(EXEC): $(EXEC).cpp
	mpicxx $(OPT) $(OMP) -o $(EXEC) $(EXEC).cpp -lrt

runp:
	srun -n $(W) $(EXEC) $(MATRIX_SIZE)

jumpshot:
	/bin/rm -f $(EXEC).clog
	srun -n 8 $(EXEC)-trace 600 2
	/bin/cp Unknown.clog2 $(EXEC).clog2
	/bin/rm Unknown.clog2

run-hpc:
	/bin/rm -rf $(EXEC).m $(EXEC).d
	srun -n 8 hpcrun -e REALTIME@1000 -t -o $(EXEC).m ./$(EXEC) 7500 2
	hpcstruct $(EXEC)
	hpcprof -S $(EXEC).hpcstruct -o $(EXEC).d $(EXEC).m

#view a trace with jumpshot 
view:
	$(VIEWER) $(EXEC).clog2

clean:
	/bin/rm -rf $(OBJ) 

clean-hpc:
	/bin/rm -r $(EXEC).d
	/bin/rm -r $(EXEC).m
	/bin/rm $(EXEC).hpcstruct
