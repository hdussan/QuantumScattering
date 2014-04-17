#
#  makefile 
#  Solve Schroedinger Equation using the numerov algorithm
#
OBJ = main.o  scattering.o numerov.o righthandside.o\
      kMesh.o readBound_wf.o
CFLAGS = -c -O3
C++ = g++
GSLFLAGS=-lgsl -lgslcblas
#Eigen Path (other machines)
EigenPath =/usr/local/include/Eigen/
InclEigen=-I$(EigenPath) 

#List of headers for dpendencies
mainH  = DifferentialEq.hpp
scatterH = scattering.hpp
numerovH = numerov.hpp
rhsH   = righthandside.hpp
kmeshH = kMesh.hpp
readWFH= readBound_wf.hpp
#List of Objects
mainFile    = main.o
scatterFile = scattering.o
numerovFile = numerov.o
rhsFile     = righthandside.o
kmeshFile   = kMesh.o
readWFFile  = readBound_wf.o
# Compilation
solvRadialSchroedinger:$(OBJ)
	$(C++) $(GSLFLAGS) -o solvRadialSchroedinger $(OBJ)
$(OBJ):%.o: %.cpp
	$(C++) $(CFLAGS) $< -o $@

# Dependencies
$(mainFile):    $(mainH) 
$(scatterFile): $(scatterH) 
$(numerovFile): $(numerovH)
$(rhsFile):	$(rhsH)
$(kmeshFile):   $(kmeshH)
$(readWFFile):  $(readWFH)
#
# Cleaning
#
clean:   
	rm *.o
cleanAll:
	rm solvDiffEq *.o