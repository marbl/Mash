CXXFLAGS += -I/usr/local/include -std=c++11 -L/usr/local/lib

all: *.cpp MinHash.capnp.h
	$(CXX) $(CXXFLAGS) -g -O2 minimap.cpp MurmurHash3.cpp CommandList.cpp Command.cpp CommandFind.cpp CommandIndex.cpp Index.cpp MinHash.capnp.c++ -o minimap -lstdc++ -lz -lcapnp -lkj

MinHash.capnp.h : MinHash.capnp
	capnp compile -oc++ MinHash.capnp

clean:
	rm -f *.o minimap
