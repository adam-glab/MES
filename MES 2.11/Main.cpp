#include <iostream>
#include "Grid.h"
#include "Element4_2D.h"
#include "Node.h"
#include "Jacobian.h"

#define LAB 4

void calc_jacobian_test(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D E, grid G) {

	double sumX = 0.;
	double sumY = 0.;
	for (int i = 0; i < nE; i++) {
		for (int j = 0; j < nIP; j++) {
			for (int k = 0; k < nIP; k++) {
				sumX += E.dNdE[j][k] * G.nodes[G.elements[i].ID[k] - 1].x;
				sumY += E.dNdn[j][k] * G.nodes[G.elements[i].ID[k] - 1].y;
			}
			J->j_matrix[0][0] = sumX;
			J->j_matrix[0][1] = 0.;
			J->j_matrix[1][0] = 0.;
			J->j_matrix[1][1] = sumY;

			J_inv->j_matrix[0][0] = 1 / sumX;
			J_inv->j_matrix[0][1] = 0.;
			J_inv->j_matrix[1][0] = 0.;
			J_inv->j_matrix[1][1] = 1 / sumY;

			std::cout << "=========================" << std::endl;
			std::cout << "::::::::::j_inv::::::::::" << std::endl;
			std::cout << J_inv->j_matrix[0][0] << '\t' << J_inv->j_matrix[0][1] << '\n' << J_inv->j_matrix[1][0] << '\t' << J_inv->j_matrix[1][1] << std::endl << std::endl;
			std::cout << "::::::::::::j::::::::::::" << std::endl;
			std::cout << J->j_matrix[0][0] << '\t' << J->j_matrix[0][1] << '\n' << J->j_matrix[1][0] << '\t' << J->j_matrix[1][1] << std::endl << std::endl;
			sumX = 0.;
			sumY = 0.;
		}
	}
}

int main() {

	jacobian J;
	jacobian J_inv;
	Element4_2D E(4);
	int nIP = 4;
	// Create grid
	grid G(0.2,0.1,5,4);
	std::cout << "\nNode values:\n";
	for (int i = 0; i < G.nN; i++) {
		std::cout << "Node " << i + 1
			<< ":\t" << G.nodes[i].x
			<< "\t" << G.nodes[i].y
			<< "\n";
	}
	std::cout << "\nElement ID values (id1, id2, id3, id4):\n";
	for (int i = 0; i < G.nE; i++) {
		std::cout << "Element " << i
			<< ":\t" << G.elements[i].ID[0]
			<< "\t" << G.elements[i].ID[1]
			<< "\t" << G.elements[i].ID[2]
			<< "\t" << G.elements[i].ID[3]
			<< "\n";
	}
	calc_jacobian_test(G.nE, nIP, &J, &J_inv, E, G);
	return 0;
}