cpp = g++

cppFlags = -std=c++11 -fopenmp $(dFlags)
cppLibs = -lnetcdf -lnetcdf_c++
oFlags = 
dFlags = -g -Wall -Wextra

srcs = main.cpp Distribution.cpp FileInterface.cpp System.cpp PriorityQueue.cpp Reaction.cpp
makeObjectFile = $(cpp) $(cppFlags) -c #-o $@

.phony: clean
clean:
	rm -rf *.o RxnSysSim debug html latex

doc: *.cpp
	doxygen doxygenConfig

all: RxnSysSim

debug: $(srcs)

RxnSysSim: $(srcs)
	$(cpp) $(cppFlags) -o $@ $(srcs) $(cppLibs)

main.o: main.h main.cpp
	$(makeObjectFile) main.cpp

Distribution.o: Distribution.h Distribution.cpp
	$(makeObjectFile) Distribution.cpp 

FileInterface.o: FileInterface.h FileInterface.cpp
	$(makeObjectFile) FileInterface.cpp

System.o: System.h System.cpp
	$(makeObjectFile) System.cpp

PriorityQueue.o: PriorityQueue.h PriorityQueue.cpp
	$(makeObjectFile) PriorityQueue.cpp	

Reaction.o: Reaction.h Reaction.cpp
	$(makeObjectFile) Reaction.cpp
