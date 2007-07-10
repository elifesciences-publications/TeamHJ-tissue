#ifndef _CELL_H_
#define _CELL_H_

#include <vector>

class Direction;
class Vertex;
class Wall;

class Cell
{
public:
	Cell();
	~Cell();

	static std::vector<Cell *> getCells(void);
	static Cell * getCellWithIndex(size_t index);

	void addWall(Wall *wall);
	std::vector<Wall *> getWalls(void);
	void addVertex(Vertex *vertex);
	std::vector<Vertex *> getVertices(void);
	void setDirection(Direction *direction);
	Direction *getDirection(void);
	void setIndex(size_t index);
	size_t getIndex(void);
	double area(void);
private:
	static std::vector<Cell *> cells;
	std::vector<Wall *> walls;
	std::vector<Vertex *> vertices;
	Direction *direction;
	size_t index;
};

#endif /* _CELL_H_ */
