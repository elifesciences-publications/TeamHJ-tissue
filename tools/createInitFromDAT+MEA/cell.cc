#include <algorithm>
#include <iostream>
#include "cell.h"

std::vector<Cell *> Cell::cells;

Cell::Cell()
{
	direction = NULL;
	cells.push_back(this);
}

Cell::~Cell()
{
	std::vector<Cell *>::iterator iterator;
	iterator = find(cells.begin(), cells.end(), this);
	if (iterator == cells.end()) {
		std::cerr << "Internal error: Unable to find cell in static vector." << std::endl;
		exit(EXIT_FAILURE);
	} else {
		cells.erase(iterator);
	}
}

std::vector<Cell *> Cell::getCells(void)
{
	return cells;
}

Cell *Cell::getCellWithIndex(size_t index)
{
	std::vector<Cell *>::iterator iterator;
	for (iterator = cells.begin(); iterator != cells.end(); ++iterator) {
		if ((*iterator)->getIndex() == index) {
			return *iterator;
		}
	}
	return NULL;
}


void Cell::addWall(Wall *wall)
{
	std::vector<Wall *>::iterator iterator;
	iterator = find(walls.begin(), walls.end(), wall);
	if (iterator != walls.end()) {
		std::cout << "Warning: Trying to add the same wall twice to cell " << index << std::endl;
	} else {
		walls.push_back(wall);
	}
}

std::vector<Wall *> Cell::getWalls(void)
{
	return walls;
}

void Cell::addVertex(Vertex *vertex)
{
	std::vector<Vertex *>::iterator iterator;
	iterator = find(vertices.begin(), vertices.end(), vertex);
	if (iterator != vertices.end()) {
		std::cout << "Warning: Trying to add the same vertex twice to cell " << index << std::endl;
	} else {
		vertices.push_back(vertex);
	}
}

std::vector<Vertex *> Cell::getVertices(void)
{
	return vertices;
}

void Cell::setDirection(Direction *direction)
{
	this->direction = direction;
}

Direction *Cell::getDirection(void)
{
	return direction;
}

void Cell::setIndex(size_t index)
{
	this->index = index;
}

size_t Cell::getIndex(void)
{
	return index;
}

double Cell::area(void)
{
	// TODO: Calculate area.
	return 1.0;
}
