#include "Grid.h"

void Grid::printNodes() const{
	std::cout << "\nNode values:\n";
	for (int i = 0; i < this->nN; i++) {
		std::cout << "Node " << i + 1 <<": "
			<< std::setw(12) << this->nodes[i].x
			<< std::setw(12) << this->nodes[i].y
			<< std::setw(12) << this->nodes[i].BC
			<< "\n";
	}
}

void Grid::printElements() const{
	std::cout << "\nElement ID values (id1, id2, id3, id4):\n";
	for (int i = 0; i < this->nE; i++) {
		std::cout << "Element " << i + 1 << ": "
			<< std::setw(12) << this->elements[i].ID[0]
			<< std::setw(12) << this->elements[i].ID[1]
			<< std::setw(12) << this->elements[i].ID[2]
			<< std::setw(12) << this->elements[i].ID[3]
			<< "\n";
	}
}


void Grid::printGlobalH(double** globalH, Grid G) {
	std::cout << std::endl << "::::::::::Global H matrix::::::::::" << std::endl;
	for (int i = 0; i < G.nN; i++) {
		for (int j = 0; j < G.nN; j++) {
			std::cout << std::setw(8) << std::setprecision(4) << globalH[i][j];
		}
		std::cout << std::endl;
	}
}

void Grid::printGlobalC(double** globalC, Grid G) {
	std::cout << std::endl << "::::::::::Global C matrix::::::::::" << std::endl;
	for (int i = 0; i < G.nN; i++) {
		for (int j = 0; j < G.nN; j++) {
			std::cout << std::setw(8) << std::setprecision(4) << globalC[i][j];
		}
		std::cout << std::endl;
	}
}

void Grid::printGlobalP(double* globalP, Grid G) {
	std::cout << std::endl << "::::::::::Global P vector::::::::::" << std::endl;
	for (int j = 0; j < G.nN; j++) {
		std::cout << std::setprecision(12) << globalP[j] << std::endl;
	}
}

