#include <algorithm>
#include <iostream>
#include "vertex.h"
#include "wall.h"

std::vector<Wall *> Wall::walls;

Wall::Wall(Vertex *v1, Vertex *v2)
{
	walls.push_back(this);
	vertices.push_back(v1);
	vertices.push_back(v2);
	index = walls.size();
}

Wall::~Wall()
{
	std::vector<Wall *>::iterator iterator;
	iterator = find(walls.begin(), walls.end(), this);
	if (iterator == walls.end()) {
		std::cerr << "Internal error: Unable to find cell in static vector." << std::endl;
		exit(EXIT_FAILURE);
	} else {
		walls.erase(iterator);
	}
}

std::vector<Wall *> Wall::getWalls(void)
{
	return walls;
}

void Wall::addCell(Cell *cell)
{
	std::vector<Cell *>::iterator iterator;
	iterator = find(cells.begin(), cells.end(), cell);
	if (iterator != cells.end()) {
		std::cout << "Warning: Trying to add the same cell twice to wall " << index << std::endl;
	} else if (cells.size() == 2) {
		std::cerr << "Error: Trying to add a third wall to cell " << index << std::endl;
		exit(EXIT_FAILURE);
	} else {
		cells.push_back(cell);
	}
}

std::vector<Cell *> Wall::getCells(void)
{
	return cells;
}

void Wall::addVertex(Vertex *vertex)
{
	std::vector<Vertex *>::iterator iterator;
	iterator = find(vertices.begin(), vertices.end(), vertex);
	if (iterator != vertices.end()) {
		std::cout << "Warning: Trying to add the same vertex twice to wall " << index << std::endl;
	} else {
		vertices.push_back(vertex);
	}
}

void Wall::removeVertex(Vertex *vertex)
{
  	std::vector<Vertex *>::iterator iterator;
 	iterator = find(vertices.begin(), vertices.end(), vertex);
 	if (iterator == vertices.end()) {
 		std::cout << "Warning: Vertex " << vertex->getIndex() << " is not know for wall " << index << ". Unable to remove." << std::endl;
	} else {
 		vertices.erase(iterator);
 	}
}

std::vector<Vertex *> Wall::getVertices(void)
{
	return vertices;
}

void Wall::swapVertices(void)
{
	if (vertices.size() != 2) {
		std::cerr << "Error: Wall " << index << " does not have exactly two vertices. Unable to swap vertices." << std::endl;
		exit(EXIT_FAILURE);
	} else {
		Vertex *v1 = vertices[0];
		Vertex *v2 = vertices[1];
		vertices[0] = v2;
		vertices[1] = v1;
	}
}


// void Wall::setIndex(size_t index)
// {
// 	this->index = index;
// }

size_t Wall::getIndex(void)
{
	return index;
}
