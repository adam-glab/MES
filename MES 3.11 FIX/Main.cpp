#include <iostream>
#include <iomanip>
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

void calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D E, grid G) {

	double sumXX = 0.;
	double sumXY = 0.;
	double sumYX = 0.;
	double sumYY = 0.;

			for (int k = 0; k < 4; k++) {
				sumXX += E.dNdE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
				sumXY += E.dNdE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
				sumYX += E.dNdn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
				sumYY += E.dNdn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
			}

			J->j_matrix[0][0] = sumXX;
			J->j_matrix[0][1] = sumXY;
			J->j_matrix[1][0] = sumYX;
			J->j_matrix[1][1] = sumYY;

			double detJ = J->j_matrix[0][0] * J->j_matrix[1][1] - J->j_matrix[0][1] * J->j_matrix[1][0];
			J_inv->j_matrix[0][0] = 1 / detJ * sumXX;
			J_inv->j_matrix[0][1] = 1 / detJ * sumXY;
			J_inv->j_matrix[1][0] = 1 / detJ * sumYX;
			J_inv->j_matrix[1][1] = 1 / detJ * sumYY;

			std::cout << "=========================" << std::endl;
			std::cout << "::::::::::j_inv::::::::::" << std::endl;
			std::cout << J_inv->j_matrix[0][0] << '\t' << J_inv->j_matrix[0][1] << '\n' << J_inv->j_matrix[1][0] << '\t' << J_inv->j_matrix[1][1] << std::endl << std::endl;
			std::cout << "::::::::::::j::::::::::::" << std::endl;
			std::cout << J->j_matrix[0][0] << '\t' << J->j_matrix[0][1] << '\n' << J->j_matrix[1][0] << '\t' << J->j_matrix[1][1] << std::endl << std::endl;
			sumXX = 0.;
			sumXY = 0.;
			sumYX = 0.;
			sumYY = 0.;
}


void calcH(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D E, grid G, int k, double detJ, double **sumArray) {
	double dNdX[4][4] = { 0. };
	double dNdY[4][4] = { 0. };
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			dNdX[z][x] = J_inv->j_matrix[0][0] * E.dNdE[z][x] + J_inv->j_matrix[1][0] * E.dNdn[z][x];
			dNdY[z][x] = J_inv->j_matrix[0][1] * E.dNdE[z][x] + J_inv->j_matrix[1][1] * E.dNdn[z][x];
		}
	}

	/*
	* debug -> print dN/dx, dN/dy
	*/

	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			std::cout << dNdX[z][x] << "  ";
		}
		std::cout << " || ";
		for (int x = 0; x < 4; x++) {
			std::cout << dNdY[z][x] << "  ";
		}
		std::cout << std::endl;
	}

	std::cout << "=========================" << std::endl;
	
	double N_X[4][4] = { 0. };
	double N_Y[4][4] = { 0. };
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			N_X[z][x] = dNdX[nIP][z] * dNdX[nIP][x];
			N_Y[z][x] = dNdY[nIP][z] * dNdY[nIP][x];
		}
	}

	/*
	* debug -> print N_X, N_Y
	*/
	
	/*
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			std::cout << N_X[z][x] << "  ";
		}
		std::cout << " || ";
		for (int x = 0; x < 4; x++) {
			std::cout << N_Y[z][x] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << "=========================" << std::endl;
	*/
	double resArray[4][4] = { 0. };
	std::cout << "Hpc" << nIP + 1 << std::endl;
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			resArray[z][x] = k * (N_X[z][x] + N_Y[z][x]) * detJ;
			std::cout << resArray[z][x] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			double tmp = resArray[z][x];
			sumArray[z][x] += tmp;
		}
	}
}

int main() {
	jacobian J;
	jacobian J_inv;
	Element4_2D E(4);
	int nIP = 4;
	// Create grid
	grid G(0.2, 0.1, 5, 4);
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
	//calc_jacobian_test(G.nE, nIP, &J, &J_inv, E, G);

	double k = 30.;

	double** sumOfH = new double*[nIP];
	for (int i = 0; i < nIP; i++) {
		sumOfH[i] = new double[nIP];
	}
	for (int i = 0; i < nIP; i++) {
		for (int j = 0; j < nIP; j++) {
			sumOfH[i][j] = 0.;
		}
	}

	for (int i = 0; i < G.nE; i++) {
		for (int j = 0; j < nIP; j++) {

			calcJacobian(i, j, &J, &J_inv, E, G);
			double detJ = J.j_matrix[0][0] * J.j_matrix[1][1] - J.j_matrix[0][1] * J.j_matrix[1][0];
			std::cout << "detJ: " << detJ << std::endl;
			calcH(i, j, &J, &J_inv, E, G, k, detJ,sumOfH);
			/*J.j_matrix[0][0] = 0.;
			J.j_matrix[1][0] = 0.;
			J.j_matrix[1][1] = 0.;
			J.j_matrix[0][1] = 0.;*/

		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << sumOfH[i][j] << "  ";
		}
		std::cout << std::endl;
	}

	return 0;
}