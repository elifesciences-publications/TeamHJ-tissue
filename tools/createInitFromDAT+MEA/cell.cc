#include <algorithm>
#include <iostream>
#include <cmath>
#include "cell.h"
#include "vertex.h"
#include "wall.h"

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
	std::vector<Wall *> sortedWalls;
	sortedWalls.push_back(walls[0]);
	while (sortedWalls.size() != walls.size()) {
		size_t sw = sortedWalls[sortedWalls.size() - 1]->getVertices()[1]->getIndex();

		bool foundWall = false;
		for (size_t i = 0; i < walls.size(); ++i) {
			Wall *wall = walls[i];

			std::vector<Vertex *> vertices = wall->getVertices();
			std::vector<Wall *>::iterator iterator = find(sortedWalls.begin(), sortedWalls.end(), wall);
			
			if (iterator != sortedWalls.end())
				continue;

			if (vertices.size() != 2) {
				std::cerr << "Error: Wall " << wall->getIndex() << " does not have exactly two vertices. Unable to calculate area." << std::endl;
				exit(EXIT_FAILURE);
			}
			
			size_t v1 = wall->getVertices()[0]->getIndex();
			size_t v2 = wall->getVertices()[1]->getIndex();
			
			if (v1 == sw) {
				sortedWalls.push_back(wall);
				foundWall = true;
			}
			if (v2 == sw) {
				wall->swapVertices();
				sortedWalls.push_back(wall);
				foundWall = true;
			}
		}
		if (foundWall == true)
			foundWall = false;
		else 
			break;
	}

	// Check
	if (sortedWalls[0]->getVertices()[0] != sortedWalls[sortedWalls.size() - 1]->getVertices()[1]) {
		std::cerr << "Error: Unexpected error while sorting walls for area calculation." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Calculate area.
	double area = 0.0;
	for (size_t i = 0; i < sortedWalls.size(); ++i) {
		double x1x = sortedWalls[i]->getVertices()[0]->getX();
		double x1y = sortedWalls[i]->getVertices()[0]->getY();
		
		double x2x = sortedWalls[(i + 1) % sortedWalls.size()]->getVertices()[0]->getX();
		double x2y = sortedWalls[(i + 1) % sortedWalls.size()]->getVertices()[0]->getY();

		area += x1x * x2y - x2x * x1y;
	}
	area /= 2.0;

	return std::abs(area);
}
