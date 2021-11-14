#include <iostream>
#include <iomanip>
#include "Grid.h"
#include "Element4_2D.h"
#include "Node.h"
#include "Jacobian.h"

#define LAB 4

void calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D E, grid G) {

	if (E.nIP == 4) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				J->j_matrix[i][j] = 0.;
				J_inv->j_matrix[i][j] = 0.;
			}
		}

		for (int k = 0; k < 4; k++) {
			J->j_matrix[0][0] += E.dNdE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
			J->j_matrix[0][1] += E.dNdE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
			J->j_matrix[1][0] += E.dNdn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
			J->j_matrix[1][1] += E.dNdn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
		}
	}

	if (E.nIP == 9) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				J->j_matrix[i][j] = 0.;
				J_inv->j_matrix[i][j] = 0.;
			}
		}

		for (int k = 0; k < 4; k++) {
			J->j_matrix[0][0] += E.dN_dE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
			J->j_matrix[0][1] += E.dN_dE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
			J->j_matrix[1][0] += E.dN_dn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
			J->j_matrix[1][1] += E.dN_dn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
		}
	}

	double detJ = J->j_matrix[0][0] * J->j_matrix[1][1] - J->j_matrix[0][1] * J->j_matrix[1][0];
	J_inv->j_matrix[0][0] = (1 / detJ) * J->j_matrix[1][1];
	J_inv->j_matrix[0][1] = -1.0 * ((1 / detJ) * J->j_matrix[1][0]);
	J_inv->j_matrix[1][0] = -1.0 * ((1 / detJ) * J->j_matrix[0][1]);
	J_inv->j_matrix[1][1] = (1 / detJ) * J->j_matrix[0][0];

	std::cout << "Point: " << nIP + 1 << std::endl;
	std::cout << "=========================" << std::endl;
	std::cout << "::::::::::j_inv::::::::::" << std::endl;
	std::cout << J_inv->j_matrix[0][0] << '\t' << J_inv->j_matrix[0][1] << '\n' << J_inv->j_matrix[1][0] << '\t' << J_inv->j_matrix[1][1] << std::endl << std::endl;
	std::cout << "::::::::::::j::::::::::::" << std::endl;
	std::cout << J->j_matrix[0][0] << '\t' << J->j_matrix[0][1] << '\n' << J->j_matrix[1][0] << '\t' << J->j_matrix[1][1] << std::endl << std::endl;
}

void calcHTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D E, grid G, int k, double detJ, double** sumArray) {

	if (E.nIP == 4) {
		double dNdX[4] = { 0. };
		double dNdY[4] = { 0. };
		for (int i = 0; i < 4; i++) {
			dNdX[i] = J_inv->j_matrix[0][0] * E.dNdE[nIP][i] + J_inv->j_matrix[1][0] * E.dNdn[nIP][i];
			dNdY[i] = J_inv->j_matrix[0][1] * E.dNdE[nIP][i] + J_inv->j_matrix[1][1] * E.dNdn[nIP][i];
		}
		double N_X[4][4] = { 0. };
		double N_Y[4][4] = { 0. };
		for (int z = 0; z < 4; z++) {
			for (int x = 0; x < 4; x++) {
				N_X[z][x] = dNdX[z] * dNdX[x];
				N_Y[z][x] = dNdY[z] * dNdY[x];
			}
		}
		double resArray[4][4] = { 0. };
		for (int z = 0; z < 4; z++) {
			for (int x = 0; x < 4; x++) {
				resArray[z][x] = k * (N_X[z][x] + N_Y[z][x]) * detJ;
			}
		}
		for (int z = 0; z < 4; z++) {
			for (int x = 0; x < 4; x++) {
				double tmp = resArray[z][x];
				sumArray[z][x] += tmp;
			}
		}
	}

	if (E.nIP == 9) {
		double dNdX[4] = { 0. };
		double dNdY[4] = { 0. };
		for (int i = 0; i < 4; i++) {
			dNdX[i] = J_inv->j_matrix[0][0] * E.dN_dE[nIP][i] + J_inv->j_matrix[1][0] * E.dN_dn[nIP][i];
			dNdY[i] = J_inv->j_matrix[0][1] * E.dN_dE[nIP][i] + J_inv->j_matrix[1][1] * E.dN_dn[nIP][i];
		}
		double N_X[4][4] = { 0. };
		double N_Y[4][4] = { 0. };
		for (int z = 0; z < 4; z++) {
			for (int x = 0; x < 4; x++) {
				N_X[z][x] = dNdX[z] * dNdX[x];
				N_Y[z][x] = dNdY[z] * dNdY[x];
			}
		}
		// include weight values for integral points
		// w1w3 w2w3 w3w3
		// w1w2 w2w2 w3w2
		// w1w1 w2w1 w3w1
		int col = -1;
		int row = 0;
		double weightCol = E.g->wP[0];
		double weightRow = E.g->wP[0];
		double resArray[4][4] = { 0. };
		for (int z = 0; z < 4; z++) {
			if (nIP < 3) {
				weightRow = E.g->wP[0];
				weightCol = E.g->wP[nIP];
			}
			if (nIP > 2 && nIP < 6) {
				weightRow = E.g->wP[1];
				weightCol = E.g->wP[nIP-3];
			}
			if (nIP > 5) {
				weightRow = E.g->wP[2];
				weightCol = E.g->wP[nIP - 6];
			}
			for (int x = 0; x < 4; x++) {
				resArray[z][x] = weightCol * weightRow * k * (N_X[z][x] + N_Y[z][x]) * detJ;
			}
		}
		for (int z = 0; z < 4; z++) {
			for (int x = 0; x < 4; x++) {
				double tmp = resArray[z][x];
				sumArray[z][x] += tmp;
			}
		}
	}
}

int main() {
	//std::cout << std::fixed << std::setprecision(6);
	jacobian J;
	jacobian J_inv;
	int nIP = 9;
	Element4_2D E(nIP);
	
	// Create grid
	grid G(0.025, 0.025, 2, 2);

	std::cout << "\nNode values:\n";
	for (int i = 0; i < G.nN; i++) {
		std::cout << "Node " << i + 1
			<< ":\t" << G.nodes[i].x
			<< "\t" << G.nodes[i].y
			<< "\t" << G.nodes[i].BC
			<< "\n";
	}
	std::cout << "\nElement ID values (id1, id2, id3, id4):\n";
	for (int i = 0; i < G.nE; i++) {
		std::cout << "Element " << i+1
			<< ":\t" << G.elements[i].ID[0]
			<< "\t" << G.elements[i].ID[1]
			<< "\t" << G.elements[i].ID[2]
			<< "\t" << G.elements[i].ID[3]
			<< "\n";
	}
	
	double k = 30.;

	double** sumOfH = new double*[4];
	for (int i = 0; i < 4; i++) {
		sumOfH[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfH[i][j] = 0.;
		}
	}

	for (int i = 0; i < G.nE; i++) {
		for (int j = 0; j < nIP; j++) {
			calcJacobian(i, j, &J, &J_inv, E, G);
			double detJ = J.j_matrix[0][0] * J.j_matrix[1][1] - J.j_matrix[0][1] * J.j_matrix[1][0];
			std::cout << "detJ: " << detJ << std::endl;
			calcHTest(i, j, &J, &J_inv, E, G, k, detJ,sumOfH);
		}
		std::cout << "H matrix for element E" << i + 1 << std::endl;
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				std::cout << sumOfH[x][z] << "  ";
			}
			std::cout << std::endl;
		}

		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				sumOfH[x][z] = 0.;
			}
			std::cout << std::endl;
		}	
	}
	return 0;
}