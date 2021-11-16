#include <iostream>
#include <iomanip>
#include "Grid.h"
#include "Element4_2D.h"
#include "Node.h"
#include "Jacobian.h"

void calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G);
void calcHTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double k, double detJ, double** sumArray);
void calcHbcTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double k);

const int nIP = 2;
const double k = 30.;

int main() {

	jacobian J;
	jacobian J_inv;
	Element4_2D E(nIP);
	
	// Create grid
	//const grid G(0.025, 0.025, 2, 2);
	const grid G(0.2, 0.1, 5, 4);

	E.printGauss();
	G.printNodes();
	G.printElements();
		
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
		std::cout << "::::::::ELEMENT " << i + 1 << "::::::::" << std::endl;
		for (int j = 0; j < nIP*nIP; j++) {
			calcJacobian(i, j, &J, &J_inv, &E, G);
			double detJ = J.j_matrix[0][0] * J.j_matrix[1][1] - J.j_matrix[0][1] * J.j_matrix[1][0];
			//std::cout << "detJ: " << detJ << std::endl <<std::endl;
			calcHTest(i, j, &J, &J_inv, &E, G, k, detJ,sumOfH);
			calcHbcTest(i, j, &J, &J_inv, &E, G, 25);
		}
		std::cout << "H matrix for element E" << i + 1 << std::endl;
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				std::cout << sumOfH[x][z] << "  ";
			}
			std::cout << std::endl;
		}
		// Reset array for next element
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				sumOfH[x][z] = 0.;
			}
		}	
	}
	for (int i = 0; i < 4; i++) {
		delete[] sumOfH[i];
	}
	delete[] sumOfH;

	return 0;
}

void calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G) {

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			J->j_matrix[i][j] = 0.;
			J_inv->j_matrix[i][j] = 0.;
		}
	}
	for (int k = 0; k < 4; k++) {
		J->j_matrix[0][0] += E->dN_dE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
		J->j_matrix[0][1] += E->dN_dE[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
		J->j_matrix[1][0] += E->dN_dn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].x;
		J->j_matrix[1][1] += E->dN_dn[nIP][k] * G.nodes[G.elements[nE].ID[k] - 1].y;
	}
	double detJ = J->j_matrix[0][0] * J->j_matrix[1][1] - J->j_matrix[0][1] * J->j_matrix[1][0];
	J_inv->j_matrix[0][0] = (1 / detJ) * J->j_matrix[1][1];
	J_inv->j_matrix[0][1] = -1.0 * ((1 / detJ) * J->j_matrix[1][0]);
	J_inv->j_matrix[1][0] = -1.0 * ((1 / detJ) * J->j_matrix[0][1]);
	J_inv->j_matrix[1][1] = (1 / detJ) * J->j_matrix[0][0];

	//std::cout << "Integral Point Number: " << nIP + 1 << std::endl;
	//std::cout << "::::::::::j_inv::::::::::" << std::endl;
	//J_inv->printJacobian();
	//std::cout << "::::::::::::j::::::::::::" << std::endl;
	//J->printJacobian();
}

void calcHTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double k, double detJ, double** sumArray) {

	double dNdX[4] = { 0. };
	double dNdY[4] = { 0. };
	for (int i = 0; i < 4; i++) {
		dNdX[i] = J_inv->j_matrix[0][0] * E->dN_dE[nIP][i] + J_inv->j_matrix[1][0] * E->dN_dn[nIP][i];
		dNdY[i] = J_inv->j_matrix[0][1] * E->dN_dE[nIP][i] + J_inv->j_matrix[1][1] * E->dN_dn[nIP][i];
	}
	double N_X[4][4] = { 0. };
	double N_Y[4][4] = { 0. };
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			N_X[z][x] = dNdX[z] * dNdX[x];
			N_Y[z][x] = dNdY[z] * dNdY[x];
		}
	}
	/*
	* include weight values for integral points
	* ==============
	* w1w2 w2w2
	* w1w1 w2w1
	* ==============
	* w1w3 w2w3 w3w3
	* w1w2 w2w2 w3w2
	* w1w1 w2w1 w3w1
	* ==============
	* using mod E.nIP: index = 0,..E.nIP-1
	*/
	double resArray[4][4] = { 0. };
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			resArray[z][x] = E->g->wP[nIP % E->nIP] * E->g->wP[nIP / E->nIP] * k * (N_X[z][x] + N_Y[z][x]) * detJ;
		}
	}
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			double tmp = resArray[z][x];
			sumArray[z][x] += tmp;
		}
	}
}

void calcHbcTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double k) {

	std::cout << "Hbc for IP" << nIP+1 << std::endl;
	double my_detJ = 0.;
	double NNT[4][4] = { 0. };
	double array_detJ[4] = { 0. };
	for (int i = 0; i < 4; i++) {
		if (G.nodes[G.elements[nE].ID[(i+3)%4] - 1].BC == 0 || G.nodes[G.elements[nE].ID[i % 4] - 1].BC == 0) {
			array_detJ[i] = 0.;
		}
		else {
			array_detJ[i] = 0.5 * sqrt(pow((G.nodes[G.elements[nE].ID[(i+3)%4] - 1].x) - (G.nodes[G.elements[nE].ID[i%4] - 1].x), 2) + pow((G.nodes[G.elements[nE].ID[(i+3)%4] - 1].y - (G.nodes[G.elements[nE].ID[i%4] - 1].y)), 2));
		}
	}
	
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			NNT[z][x] = k * array_detJ[nIP/E->nIP] * ((E->N[nIP][0][z] * E->N[nIP][0][x]) + (E->N[nIP][1][z] * E->N[nIP][1][x]));
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << NNT[i][j] << " ";
		}
		std::cout << std::endl;
	}

}