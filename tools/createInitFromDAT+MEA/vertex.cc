#include <algorithm>
#include <iostream>
#include "cell.h"
#include "vertex.h"
#include "wall.h"

std::vector<Vertex *> Vertex::vertices;

Vertex::Vertex()
{
	vertices.push_back(this);
}

Vertex::~Vertex()
{
	std::vector<Vertex *>::iterator iterator;
	iterator = find(vertices.begin(), vertices.end(), this);
	if (iterator == vertices.end()) {
		std::cerr << "Internal error: Unable to find vertex in static vector." << std::endl;
		exit(EXIT_FAILURE);
	} else {
		vertices.erase(iterator);
	}
}

std::vector<Vertex *> Vertex::getVertices(void)
{
	return vertices;
}

void Vertex::addCell(Cell *cell)
{
	std::vector<Cell *>::iterator iterator;
	iterator = find(cells.begin(), cells.end(), cell);
	if (iterator != cells.end()) {
		std::cout << "Warning: Trying to add the same cell twice to vertex " << index << std::endl;
	} else {
		cells.push_back(cell);
	}
}

void Vertex::removeCell(Cell *cell)
{
 	std::vector<Cell *>::iterator iterator;
	iterator = find(cells.begin(), cells.end(), cell);
	if (iterator == cells.end()) {
		std::cout << "Warning: Cell " << cell->getIndex() << " is not know for vertex " << index << ". Unable to remove." << std::endl;
	} else {
		cells.erase(iterator);
	}
}

std::vector<Cell *> Vertex::getCells(void)
{
	return cells;
}


void Vertex::addWall(Wall *wall)
{
	std::vector<Wall *>::iterator iterator;
	iterator = find(walls.begin(), walls.end(), wall);
	if (iterator != walls.end()) {
		std::cout << "Warning: Trying to add the same wall twice to vertex " << index << std::endl;
	} else {
		walls.push_back(wall);
	}
}

void Vertex::removeWall(Wall *wall)
{
 	std::vector<Wall *>::iterator iterator;
	iterator = find(walls.begin(), walls.end(), wall);
	if (iterator == walls.end()) {
		std::cout << "Warning: Wall " << wall->getIndex() << " is not know for vertex " << index << ". Unable to remove." << std::endl;
	} else {
		walls.erase(iterator);
	}
}

std::vector<Wall *> Vertex::getWalls(void)
{
	return walls;
}

void Vertex::setPosition(double x, double y)
{
	this->x = x;
	this->y = y;
}
void Vertex::setX(double x)
{
	this->x = x;
}

void Vertex::setY(double y)
{
	this->y = y;   
}

double Vertex::getX(void)
{
	return x;
}

double Vertex::getY(void)
{
	return y;
}

void Vertex::setIndex(size_t index)
{
	this->index = index;
}

size_t Vertex::getIndex(void)
{
	return index;
}


