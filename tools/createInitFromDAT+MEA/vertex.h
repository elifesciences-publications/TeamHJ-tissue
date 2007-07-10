#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <vector>

class Cell;
class Wall;

class Vertex
{
public:
	Vertex();
	~Vertex();

	static std::vector<Vertex *> getVertices(void);

	void addCell(Cell *cell);
	void removeCell(Cell *cell);
	std::vector<Cell *> getCells(void);
	void addWall(Wall *wall);
	void removeWall(Wall *wall);
	std::vector<Wall *> getWalls(void);
	void setPosition(double x, double y);
	void setX(double x);
	void setY(double y);
	double getX(void);
	double getY(void);
	void setIndex(size_t index);
	size_t getIndex(void);

private:
	static std::vector<Vertex *> vertices;
	std::vector<Cell *> cells;
	std::vector<Wall *> walls;
	double x;
	double y;
	size_t index;
};

#endif /* _VERTEX_H_ */
