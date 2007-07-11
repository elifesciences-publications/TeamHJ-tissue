#ifndef _WALL_H_
#define _WALL_H_

#include <vector>

class Cell;
class Vertex;

class Wall
{
public:
	Wall(Vertex *v1, Vertex *v2);
	~Wall();

	static std::vector<Wall *> getWalls(void);
	void addCell(Cell *cell);
	std::vector<Cell *> getCells(void);
	void addVertex(Vertex *vertex);
	void removeVertex(Vertex *vertex);
	std::vector<Vertex *> getVertices(void);
	//	void setIndex(size_t index);
	size_t getIndex(void);
	void swapVertices(void);
	
private:
	static std::vector<Wall *> walls;
	std::vector<Cell *> cells;
	std::vector<Vertex *> vertices;
	size_t index;
};

#endif /* _WALL_H_ */ 
