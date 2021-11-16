#include "Grid.h"

void grid::printNodes() const{
	std::cout << "\nNode values:\n";
	for (int i = 0; i < this->nN; i++) {
		std::cout << "Node " << i + 1
			<< ":\t" << this->nodes[i].x
			<< "\t" << this->nodes[i].y
			<< "\t" << this->nodes[i].BC
			<< "\n";
	}
}

void grid::printElements() const{
	std::cout << "\nElement ID values (id1, id2, id3, id4):\n";
	for (int i = 0; i < this->nE; i++) {
		std::cout << "Element " << i + 1
			<< ":\t" << this->elements[i].ID[0]
			<< "\t" << this->elements[i].ID[1]
			<< "\t" << this->elements[i].ID[2]
			<< "\t" << this->elements[i].ID[3]
			<< "\n";
	}
}


