# make file to include stuff 
# in MY_CPP 
Inc= -I/home/walter/Desktop/4Simulation/My_CPP
CC=g++
OBJs=DumpSurfaceNode.o  DumpPts.o

default: DumpMesh
DumpSurfaceNode.o: DumpSurfaceNode.C
	$(CC) -c $< -o $@ $(Inc)
DumpPts.o: DumpPts.C
	$(CC) -c $< -o $@ $(Inc)
DumpMesh: $(OBJs)
	$(CC)  $^ -o $@
test: test.c
	$(CC)  $^ -o $@ $(Inc)
.PHONEY:clean

clean:
	rm -f DumpMesh *.o
	
