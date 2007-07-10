#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>
#include <cmath>
#include "cell.h"
#include "direction.h"
#include "vertex.h"
#include "wall.h"

template <class T>
T stringTo(std::string string);

template <class T>
size_t getIndexForElement(std::vector<T> &vector, T element);

void readDATFile(std::istream &input);
void removeNonValidCells(void);
void createWalls(void);
void splitFourVertices(void);
void readMEAFile(std::istream &input);
void checkEntities(void);
void printInit(std::ostream &output);

int main(int argc, char *argv[])
{
	if (argc != 4) {
		std::cerr << "Usage: " << argv[0] << " input.dat input.mea output.init" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ifstream datfile;
	datfile.open(argv[1]);
	if (!datfile) {
		std::cerr << "Error: Unable to open file '" << argv[1] << "' for reading." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ifstream meafile;
	meafile.open(argv[2]);
	if (!meafile) {
		std::cerr << "Error: Unable to open file '" << argv[2] << "' for reading." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ofstream initfile;
	initfile.open(argv[3]);
	if (!initfile) {
		std::cerr << "Error: Unable to open file '" << argv[3] << "' for writing." << std::endl;
		exit(EXIT_FAILURE);
	}

	// TODO:
	//  Flip
	//  Rescale

	readDATFile(datfile);

	removeNonValidCells();

	createWalls();

	splitFourVertices();

	readMEAFile(meafile);

	checkEntities();

	printInit(initfile);


	initfile.close();
	meafile.close();
	datfile.close();
	exit(EXIT_SUCCESS);
}

template <class T>
T stringTo(std::string string)
{
	std::stringstream stream(string);
	T value;
	stream >> value;
	return value;
}

template <class T>
size_t getIndexForElement(std::vector<T> &vector, T element)
{
	for (size_t i = 0; i < vector.size(); ++i) {
		if (vector[i] == element)
			return i;
	}
	std::cerr << "Error: Unable to find element in vector." << std::endl;
	exit(EXIT_FAILURE);
}

void readDATFile(std::istream &input)
{
	std::cout << "Reading .dat file." << std::endl;

	std::string trash;
	size_t numberOfVertices;
	input >> trash 
		 >> trash
		 >> numberOfVertices;
	
	for (size_t n = 0; n < numberOfVertices; ++n) {
		Vertex *vertex = new Vertex;
		std::string tmp;

		input >> tmp; // Vertex index;
		vertex->setIndex(stringTo<size_t>(tmp));

		input >> tmp; // Vertex position;

		std::string::iterator start, end;

		start = find(tmp.begin(), tmp.end(), '<');
		end = find(tmp.begin(), tmp.end(), ',');

 		if (start == tmp.end() || end == tmp.end()) {
 			std::cerr << "Error: Corrupted .dat file." << std::endl;
 			exit(EXIT_FAILURE);
 		}
		std::string x(++start, end);

 		start = find(tmp.begin(), tmp.end(), ',');
 		end = find(tmp.begin(), tmp.end(), '>');

  		if (start == tmp.end() || end == tmp.end()) {
  			std::cerr << "Error: Corrupted .dat file." << std::endl;
  			exit(EXIT_FAILURE);
  		}
  		std::string y(++start, end);

		vertex->setPosition(stringTo<double>(x), stringTo<double>(y));

		size_t numberOfCells;
		input >> numberOfCells;
		for (size_t i = 0; i < numberOfCells; ++i) {
			size_t cellIndex;
			input >> cellIndex;

			Cell *cell = Cell::getCellWithIndex(cellIndex);
			if (cell == NULL) {
				cell = new Cell();
				cell->setIndex(cellIndex);
			}

			cell->addVertex(vertex);
			vertex->addCell(cell);
		}
	}
}

void removeNonValidCells(void)
{
	bool hasRemovedCell = false;
	while (true) {
		std::vector<Cell *> cells = Cell::getCells();
		std::vector<Cell *>::iterator cellIterator;
		for (cellIterator = cells.begin(); cellIterator != cells.end(); ++cellIterator) {
			Cell *cell = *cellIterator;
			if (cell->getVertices().size() < 3) {
				std::cout << "Warning: Found cell with less than three vertices." << std::endl;

				std::vector<Vertex *> vertices = cell->getVertices();
				std::cout << " Cell: " << cell->getIndex() << std::endl;
				std::cout << " Vertices: ";

				for (size_t i = 0; i < vertices.size(); ++i)
					std::cout << vertices[i]->getIndex() << " ";
				std::cout << std::endl;


				for (size_t i = 0; i < vertices.size(); ++i) {
					std::cout << " Removing reference from vertex: " << vertices[i]->getIndex() << std::endl;
					vertices[i]->removeCell(cell);
				}

				std::cout << " Removing cell " << cell->getIndex() << std::endl << std::endl;
				delete cell;

				hasRemovedCell = true;
				continue;
			}
		}
		
		if (hasRemovedCell == true) {
			hasRemovedCell = false;
			continue;
		} else {
			break;
		}
	}
}

void createWalls(void)
{
	// 1. Loop over all cells.
	// 2. Find center of cell.
	// 3. Find polar angle coordinate to every vertex.
	// 4. Sort vertices according to angles.
	// 5. Tie together vertices with walls.
	// 6. Check for already existing walls.

	// 1. Loop over all cells.
	std::vector<Cell *> cells = Cell::getCells();
	std::vector<Cell *>::iterator cellIterator;
	for (cellIterator = cells.begin(); cellIterator != cells.end(); ++cellIterator) {
		Cell *cell = *cellIterator;
		std::vector<Vertex *> vertices = cell->getVertices();
		std::vector<Vertex *>::iterator vertexIterator;

		// 3. Find center of cell
		double x = 0.0;
		double y = 0.0;
		for (vertexIterator = vertices.begin(); vertexIterator != vertices.end(); ++vertexIterator) {
			Vertex *vertex = *vertexIterator;
			x += vertex->getX();
			y += vertex->getY();
		}
		x /= vertices.size();
		y /= vertices.size();


		// 3. Find polar angle coordinate to every vertex.
		std::map<double, Vertex *> angles;

		for (vertexIterator = vertices.begin(); vertexIterator != vertices.end(); ++vertexIterator) {
			Vertex *vertex = *vertexIterator;
			
			double vx = vertex->getX() - x;
			double vy = - (vertex->getY() - y);
			double angle;
			if (vx > 0 && vy == 0) {
				angle = 0;
			} else if (vx > 0 && vy > 0) {
				angle = std::atan(vy / vx);
			} else if (vx == 0 && vy > 0) {
				angle = M_PI / 2.0;
			} else if (vx < 0 && vy > 0) {
				angle = M_PI + std::atan(vy / vx);
			} else if (vx < 0 && vy == 0) {
				angle = M_PI;
			} else if (vx < 0 && vy < 0) {
				angle = M_PI + std::atan(vy / vx);
			} else if (vx == 0 && vy < 0) {
				angle = 1.5 * M_PI;
			} else if (vx > 0 && vy < 0) {
				angle = 2 * M_PI + std::atan(vy / vx);
			} else {
				std::cerr << "Internal error: Unable to determine angle." << std::endl;
				exit(EXIT_FAILURE);
			}
			angles[angle] = vertex;
		}


		// 4. Sort vertices according to angles.
 		std::vector<Vertex *> sortedVertices;
 		std::map<double, Vertex *>::iterator angleIterator = angles.lower_bound(0.0);
 		while (angleIterator != angles.end()) {
 			sortedVertices.push_back(angleIterator->second);
 			angles.erase(angleIterator);
 			angleIterator = angles.lower_bound(0.0);
 		}


		// Loop over sorted vertices.
		for (size_t n = 0; n < sortedVertices.size(); ++n) {
			Vertex *vertex1 = sortedVertices[n];
			Vertex *vertex2 = sortedVertices[(n + 1) % sortedVertices.size()];

			std::vector<Wall *> walls = vertex1->getWalls();
			bool foundWallBetweenVertices = false;

			// 6. Check for already existing walls.
			for (size_t i = 0; i < walls.size(); ++i) {
				std::vector<Vertex *> wallVertices = walls[i]->getVertices();
				for (size_t j = 0; j < wallVertices.size(); ++j) {
					if (wallVertices[j] == vertex2) {
						// Found a wall between vertex1 and vertex2.
						// Then add current cell to that wall.
						walls[i]->addCell(cell);
						cell->addWall(walls[i]);

						if (foundWallBetweenVertices == true) {
							std::cerr << "Error: (At least) two walls between vertex " 
									<< vertex1->getIndex() << " and vertex " << vertex2->getIndex() << " were found." 
									<< std::endl;
							exit(EXIT_FAILURE);
						} else {					   
							foundWallBetweenVertices = true;
						}
					}
				}
			}

			// 5. Tie together vertices with walls.
			if (foundWallBetweenVertices == false) {
				// No wall was found between vertex1 and vertex2. 
				// Create wall and add it to current cell and vertex1 and vertex2.
				Wall *wall = new Wall(vertex1, vertex2);
				wall->addCell(cell);
				cell->addWall(wall);
				vertex1->addWall(wall);
				vertex2->addWall(wall);
			}
		}
	}
}

void splitFourVertices(void)
{
	std::vector<Vertex *>::iterator vertexIterator;
	bool foundFourVertex = false;
	while (true) {
		std::vector<Vertex *> vertices = Vertex::getVertices();
		for (vertexIterator = vertices.begin(); vertexIterator != vertices.end(); ++vertexIterator) {
			Vertex *vertex = *vertexIterator;
			std::vector<Wall *> walls = vertex->getWalls();

			if (walls.size() == 4) {
				std::cout << "Warning: Found four-vertices at vertex " << vertex->getIndex() << std::endl;
				
				Wall *wall1 = walls[0];
				Wall *wall2 = NULL;
				Cell *cell1;
				std::vector<Wall *> walls = vertex->getWalls();
				bool foundWall = false;
				for (size_t n = 1; n < walls.size(); ++n) {
					for (size_t i = 0; i < wall1->getCells().size(); ++i) {
						for (size_t j = 0; j < walls[n]->getCells().size(); ++j) {
							if (wall1->getCells()[i] == walls[n]->getCells()[j]) {
								cell1 = wall1->getCells()[i];
								wall2 = walls[n];
								foundWall = true;
							}
							if (foundWall)
								break;
						}
						if (foundWall)
							break;
					}
					if (foundWall)
						break;
				}
				
				if (wall2 == NULL) {
					std::cerr << "Internal error: Could not find a pair of neighboring walls for four vertices split." << std::endl;
					exit(EXIT_FAILURE);
				}
				
				// Remove vertex from wall
				wall1->removeVertex(vertex);
				vertex->removeWall(wall1);
				wall2->removeVertex(vertex);
				vertex->removeWall(wall2);

				// Create new vertex.
				Vertex *newVertex = new Vertex();
				newVertex->setIndex(Vertex::getVertices().size());
				double dx;
				double dy;
				if (wall2->getVertices()[0] == vertex) {
					dx = wall2->getVertices()[1]->getX() - vertex->getX();
					dy = wall2->getVertices()[1]->getY() - vertex->getY();
				} else {
					dx = wall2->getVertices()[0]->getX() - vertex->getX();
					dy = wall2->getVertices()[0]->getY() - vertex->getY();
				}
  				newVertex->setX(vertex->getX() + 0.1 * dx);
 				newVertex->setY(vertex->getY() + 0.1 * dy);

				wall1->addVertex(newVertex);
				wall2->addVertex(newVertex);

				// Create new wall.
				Wall *newWall = new Wall(vertex, newVertex);
				vertex->addWall(newWall);
				newVertex->addWall(newWall);
				newVertex->addWall(wall1);
				newVertex->addWall(wall2);

				// Add cell to new vertex and new wall.
				for (size_t i = 0; i < wall1->getCells().size(); ++i) {
					if (wall1->getCells()[i] != cell1) {
						newWall->addCell(wall1->getCells()[i]);
						wall1->getCells()[i]->addWall(newWall);
						newVertex->addCell(wall1->getCells()[i]);
						wall1->getCells()[i]->addVertex(newVertex);
					}
				}
				for (size_t i = 0; i < wall2->getCells().size(); ++i) {
					if (wall2->getCells()[i] != cell1) {
						newWall->addCell(wall2->getCells()[i]);
						wall2->getCells()[i]->addWall(newWall);
						newVertex->addCell(wall2->getCells()[i]);
						wall2->getCells()[i]->addVertex(newVertex);
					}
				}

				foundFourVertex = true;
			}
		}
		if (foundFourVertex == true) {
			foundFourVertex = false;
		} else {
			break;
		}		
	}
}

// Temporary container class.
class Data {
public:
	double fromx;
	double fromy;
	double tox;
	double toy;
};

void readMEAFile(std::istream &input)
{
	std::vector<Data *> data;

	std::string trash;
	input >> trash >> trash;
	
	size_t numberOfDirections;
	input >> numberOfDirections;

	for (size_t n = 0; n < numberOfDirections; ++n) {
		Data *tmp = new Data;
		input >> trash >> trash;
		input >> tmp->fromx >> tmp->fromy >> trash >> tmp->tox >> tmp->toy >> trash;
		data.push_back(tmp);
	}

	bool foundMatch = false;
	while (true) {
		for (size_t n = 0; n < numberOfDirections; ++n) {
			if (Direction::getDirectionWithIndex(n) != NULL)
				continue;

			double fromx = data[n]->fromx;
			double fromy = data[n]->fromy;
			double tox = data[n]->tox;
			double toy = data[n]->toy;
			
			std::vector<Cell *> candidateCells;
			std::vector<Cell *> cells = Cell::getCells();
			for (size_t i = 0; i < cells.size(); ++i) {
				double xmin = std::numeric_limits<double>::max();
				double xmax = std::numeric_limits<double>::min();
				double ymin = std::numeric_limits<double>::max();
				double ymax = std::numeric_limits<double>::min();
				
				Cell *cell = cells[i];
				if (cell->getDirection() != NULL)
					continue;
			
				std::vector<Vertex *> vertices = cell->getVertices();
				for (size_t j = 0; j < vertices.size(); ++j) {
					Vertex *vertex = vertices[j];
					
					if (vertex->getX() > xmax)
						xmax = vertex->getX();
					if (vertex->getX() < xmin)
						xmin = vertex->getX();
					if (vertex->getY() > ymax)
						ymax = vertex->getY();
					if (vertex->getY() < ymin)
						ymin = vertex->getY();
				}
				
				if ((fromx >= xmin && tox >= xmin) && (fromx <= xmax && tox <= xmax) &&
				    (fromy >= ymin && toy >= ymin) && (fromy <= ymax && toy <= ymax)) {
					candidateCells.push_back(cell);
					foundMatch = true;
				}
			}
			
			if (candidateCells.size() == 1) {
				Cell *cell = candidateCells[0];
				Direction *direction = Direction::getDirectionWithIndex(n);
				if (direction != NULL) {
					std::cerr << "Error: Direction " << n << " not uniquely assigned to one cell. "
							<< "Candidates are " << direction->getCell()->getIndex() 
							<< " and " << cell->getIndex() << std::endl;
					exit(EXIT_FAILURE);
				} else {
					direction = new Direction();
				}
								
				double dx = tox - fromx;
				double dy = toy - fromy;
				double A = std::sqrt(dx * dx + dy * dy);
				dx /= A;
				dy /= A;
				
				direction->setIndex(n);
				direction->setX(dx);
				direction->setY(dy);
				direction->setCell(cell);
				cell->setDirection(direction);

				std::cout << "Matching direction " << direction->getIndex() << " with cell " << cell->getIndex() << std::endl;
			}
		}
		if (foundMatch == true)
			foundMatch = false;
		else 
			break;
	}

	if (Direction::getDirections().size() != numberOfDirections) {
		std::cerr << "Error: Total number of directions does not match number of directions in .mea-file." << std::endl;
		exit(EXIT_FAILURE);
	}
}
	
void checkEntities(void)
{
	std::vector<Cell *> cells = Cell::getCells();
	std::vector<Vertex *> vertices = Vertex::getVertices();
 	std::vector<Wall *> walls = Wall::getWalls();
	
	for (size_t i = 0; i < cells.size(); ++i) {
		Cell *cell = cells[i];
		if (cell->getWalls().size() < 3) {
			std::cerr << "Error: Cell " << cell->getIndex() << " has less than three walls." << std::endl;
			exit(EXIT_FAILURE);
		}
		if (cell->getVertices().size() != cell->getWalls().size()) {
			std::cerr << "Error: The number of vertices (" << cell->getVertices().size() 
					<< ") does not match the number of walls (" << cell->getWalls().size()
					<< ") for cell " << cell->getIndex() << std::endl;

			exit(EXIT_FAILURE);
		}
	}

	for (size_t i = 0; i < vertices.size(); ++i) {
		Vertex *vertex = vertices[i];
		if (vertex->getWalls().size() < 2) {
			std::cerr << "Error: Vertex " << vertex->getIndex() << " has less than two walls." << std::endl;
			exit(EXIT_FAILURE);
		}
		if (vertex->getCells().size() < 1) {
			std::cerr << "Error: Vertex " << vertex->getIndex() << " has less no cells." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	for (size_t i = 0; i < walls.size(); ++i) {
		Wall *wall = walls[i];
		if (wall->getVertices().size() != 2) {
			std::cerr << "Error: Wall " << wall->getIndex() << " has " << wall->getVertices().size() << " vertices." << std::endl;
			exit(EXIT_FAILURE);
		}
		if (wall->getCells().size() != 2 && wall->getCells().size() != 1) {
			std::cerr << "Error: Wall " << wall->getIndex() << " has " << wall->getCells().size() << " cells." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}


void printInit(std::ostream &output)
{
 	std::vector<Cell *> cells = Cell::getCells();
 	std::vector<Vertex *> vertices = Vertex::getVertices();
 	std::vector<Wall *> walls = Wall::getWalls();

 	output << cells.size() << " " << walls.size() << " " << vertices.size() << std::endl;
	
 	for (size_t i = 0; i < walls.size(); ++i) {
 		Wall *wall = walls[i];
		std::vector<Cell *> wallCells = wall->getCells();
 		output << i << " ";
		for (size_t j = 0; j < wallCells.size(); ++j) {
			output << getIndexForElement<Cell *>(cells, wallCells[j]) << " ";
		}
		if (wallCells.size() == 1)
			output << "-1" << " ";

		std::vector<Vertex *> wallVertices = wall->getVertices();
		for (size_t j = 0; j < wallVertices.size(); ++j) {
			output << getIndexForElement<Vertex *>(vertices, wallVertices[j]) << " ";
		}

		output << std::endl;
 	}

	output << std::endl;

	output << vertices.size() << " " << 2 << std::endl;
	for (size_t i = 0; i < vertices.size(); ++i) {
		Vertex *vertex = vertices[i];
		output << vertex->getX() << " " << vertex->getY() << std::endl;
	}

	output << std::endl;

	output << walls.size() << " " << 1 << " " << 2 << std::endl;
	for (size_t i = 0; i < walls.size(); ++i) {
		Wall *wall = walls[i];
		Vertex *v1 = wall->getVertices()[0];
		Vertex *v2 = wall->getVertices()[1];
		double dx = v2->getX() - v1->getX();
		double dy = v2->getY() - v1->getY();
		double distance = std::sqrt(dx * dx + dy * dy);
		output << distance << " " << 99 << " " << 99 << std::endl;
	}

	output << std::endl;
	output << cells.size() << " " << 4 << std::endl;
	for (size_t i = 0; i < cells.size(); ++i) {
		Cell *cell = cells[i];
		Direction *direction = cell->getDirection();
		if (direction != NULL) {
			output << direction->getX() << " " << direction->getY() << " " << 1 << " " << cell->area() << std::endl;
		} else {
			output << "0 0 0 " << cell->area() << std::endl;
		}
	}
}
