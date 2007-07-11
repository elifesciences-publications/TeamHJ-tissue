#include <algorithm>
#include <iostream>
#include "direction.h"

std::vector<Direction *> Direction::directions;

Direction::Direction()
{
	directions.push_back(this);
}

Direction::~Direction()
{
	std::vector<Direction *>::iterator iterator;
	iterator = find(directions.begin(), directions.end(), this);
	if (iterator == directions.end()) {
		std::cerr << "Internal error: Unable to find cell in static vector." << std::endl;
		exit(EXIT_FAILURE);
	} else {
		directions.erase(iterator);
	}
}

std::vector<Direction *> Direction::getDirections(void)
{
	return directions;
}

Direction *Direction::getDirectionWithIndex(size_t index)
{
	std::vector<Direction *>::iterator iterator;
	for (iterator = directions.begin(); iterator != directions.end(); ++iterator) {
		if ((*iterator)->getIndex() == index) {
			return *iterator;
		}
	}
	return NULL;
}

void Direction::setIndex(size_t index)
{
	this->index = index;
}

size_t Direction::getIndex(void)
{
	return index;
}

void Direction::setX(double x)
{
	this->x = x;
}

double Direction::getX(void)
{
	return x;
}

void Direction::setY(double y)
{
	this->y = y;
}

double Direction::getY(void)
{
	return y;
}

void Direction::setCell(Cell *cell)
{
	this->cell = cell;
}

Cell *Direction::getCell(void)
{
	return cell;
}

