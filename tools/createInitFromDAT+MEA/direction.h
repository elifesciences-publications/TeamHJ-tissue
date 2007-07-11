#ifndef _DIRECTION_H_
#define _DIRECTION_H_

#include <vector>
class Cell;

class Direction
{
 public:
	Direction();
	~Direction();

	static std::vector<Direction *> getDirections(void);
	static Direction * getDirectionWithIndex(size_t index);

	void setIndex(size_t index);
	size_t getIndex(void);
	void setX(double x);
	double getX(void);
	void setY(double y);
	double getY(void);
	void setCell(Cell *cell);
	Cell *getCell(void);
		
 private:
	static std::vector<Direction *> directions;

	size_t index;
	double x;
	double y;
	Cell *cell;
};


#endif /* _DIRECTION_H_ */
