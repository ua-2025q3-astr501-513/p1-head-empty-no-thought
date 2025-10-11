INC=../include
CXXFLAGS= -c -g -Wall -I$(INC) -std=c++11

all: body node tree bh
	g++ body.o node.o tree.o barnes-hut.o -o runME

bh: body node tree
	g++ $(CXXFLAGS) barnes-hut.cpp

tree: body node
	g++ $(CXXFLAGS) tree.cpp 

node: body
	g++ $(CXXFLAGS) node.cpp 

body: 
	g++ $(CXXFLAGS) body.cpp 

clean:
	rm -rf *.o *.mod runME