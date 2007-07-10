createInitFromDAT+MEA : main.cc vertex.o cell.o wall.o direction.o
	c++ main.cc vertex.o cell.o wall.o direction.o -o createInitFromDAT+MEA -Wall -pedantic -ansi -g

vertex.o : vertex.cc vertex.h
	c++ -c vertex.cc -o vertex.o -Wall -pedantic -ansi -g

cell.o : cell.cc cell.h
	c++ -c cell.cc -o cell.o -Wall -pedantic -ansi -g

wall.o : wall.cc wall.h
	c++ -c wall.cc -o wall.o -Wall -pedantic -ansi -g

direction.o : direction.cc direction.h
	c++ -c direction.cc -o direction.o -Wall -pedantic -ansi -g
