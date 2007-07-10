createInitFromDAT+MEA : main.cc vertex.o cell.o wall.o
	c++ main.cc vertex.o cell.o wall.o -o createInitFromDAT+MEA -Wall -pedantic -ansi -g

vertex.o : vertex.cc vertex.h
	c++ -c vertex.cc -o vertex.o -Wall -pedantic -ansi -g

cell.o : cell.cc cell.h
	c++ -c cell.cc -o cell.o -Wall -pedantic -ansi -g

wall.o : wall.cc wall.h
	c++ -c wall.cc -o wall.o -Wall -pedantic -ansi -g
