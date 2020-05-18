
NVCC        = nvcc
NVCC_FLAGS  = -O3
OBJ         = main.o graph.o
EXE         = truss


default: $(EXE)

%.o: %.cu
	$(NVCC) $(NVCC_FLAGS) -c -o $@ $<

$(EXE): $(OBJ)
	$(NVCC) $(NVCC_FLAGS) $(OBJ) -o $(EXE)

clean:
	del $(OBJ) $(EXE) truss.*

