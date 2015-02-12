all: minimap.cpp MurmurHash3.cpp
	$(CXX) -g -O2 minimap.cpp MurmurHash3.cpp -o minimap -lz

clean:
	rm -f *.o minimap
