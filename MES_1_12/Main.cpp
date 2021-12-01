#include <iostream>
#include <iomanip>
#include "Solver.h"

// pass to grid
const int nIP = 2;
const double H = 0.1, B = 0.1;
const int nH = 4, nB = 4;
// factor values
const double k = 25., alpha = 300., t_env = 1200, c = 700, ro = 7800, dTau = 50;

// helper print funcs
/*
*
* 
* 
*/

int main() {

	jacobian J;
	jacobian J_inv;
	Element4_2D E(nIP);

	// Grid management
	const grid G(H,B,nH,nB);

	//E.printGauss();
	//G.printNodes();
	//G.printElements();

	// Global array of H-vals in nodes
	// nN x nN
	double** globalH = new double* [G.nH*G.nB];
	for (int i = 0; i < G.nH*G.nB; i++) {
		globalH[i] = new double[G.nH*G.nB];
	}
	for (int i = 0; i < G.nH * G.nB; i++) {
		for (int j = 0; j < G.nH * G.nB; j++) {
			globalH[i][j] = 0.;
		}
	}
	// Global array of P-vals in nodes
	// 1 x nN
	double* globalP = new double [G.nH * G.nB];
	for (int i = 0; i < G.nH * G.nB; i++) {
		globalP[i] = 0.;
	}

	// Global C-matrix
	// nN x nN
	double** globalC = new double* [G.nH * G.nB];
	for (int i = 0; i < G.nH * G.nB; i++) {
		globalC[i] = new double[G.nH * G.nB];
	}
	for (int i = 0; i < G.nH * G.nB; i++) {
		for (int j = 0; j < G.nH * G.nB; j++) {
			globalC[i][j] = 0.;
		}
	}

	// sumOfH -> 2d array to hold H matrixes of each element
	// H matrix in each integral point is calculated inside calcHTest function as a local matrix and then added to the sum
	double** sumOfH = new double*[4];
	for (int i = 0; i < 4; i++) {
		sumOfH[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfH[i][j] = 0.;
		}
	}

	// sumOfP -> 2d array to hold P matrixes of each element
	double* sumOfP = new double[4];
	for (int i = 0; i < 4; i++) {
		sumOfP[i] = 0;
	}

	// sumOfHbc -> 2d array to hold Hbc matrixes of each element
	double** sumOfHbc = new double* [4];
	for (int i = 0; i < 4; i++) {
		sumOfHbc[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfHbc[i][j] = 0.;
		}
	}

	// sumOfC -> 2d array to hold C matrixes of each element
	double** sumOfC = new double* [4];
	for (int i = 0; i < 4; i++) {
		sumOfC[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfC[i][j] = 0.;
		}
	}

	for (int i = 0; i < G.nE; i++) {
		std::cout << "::::::::ELEMENT " << i + 1 << "::::::::" << std::endl;
		for (int j = 0; j < nIP*nIP; j++) {
			Solver::calcJacobian(i, j, &J, &J_inv, &E, G);
			double detJ = J.j_matrix[0][0] * J.j_matrix[1][1] - J.j_matrix[0][1] * J.j_matrix[1][0];
			Solver::calcHTest(i, j, &J, &J_inv, &E, G, k, detJ, sumOfH, globalH);
			Solver::calcHbcTest(i, j, &J, &J_inv, &E, G, alpha, sumOfHbc);
			Solver::calcPTest(i, j, &J, &J_inv, &E, G, alpha, t_env, sumOfP);
			Solver::calcCTest(i, j, &J, &J_inv, &E, G, c, ro, detJ, sumOfC, globalC);
		}

		//std::cout << "H matrix for element E" << i + 1 << std::endl;
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				// add to global H matrix
				globalH[G.elements[i].ID[x] - 1][G.elements[i].ID[z] - 1] += sumOfH[x][z] /*+ sumOfHbc[x][z]*/;
				globalC[G.elements[i].ID[x] - 1][G.elements[i].ID[z] - 1] += sumOfC[x][z];
				//std::cout <<std::setw(12)<< sumOfH[x][z];
			}
			//std::cout << std::endl;
		}

		/*std::cout << "Hbc matrix for element E" << i + 1 << std::endl;
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				std::cout << std::setw(12) << sumOfHbc[x][z];
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;*/
		/*std::cout << "C matrix for element E" << i + 1 << std::endl;
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				std::cout << std::setw(12) << sumOfC[x][z];
			}
			std::cout << std::endl;
		}*/

		std::cout << std::endl;
		std::cout << "P vector for element E" << i + 1 << std::endl;
		for (int z = 0; z < 4; z++) {
			globalP[G.elements[i].ID[z] - 1] += sumOfP[z];
		std::cout << std::setw(12) << sumOfP[z];
		}
		std::cout << std::endl;

		// Reset array for next element
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				sumOfH[x][z] = 0.;
				sumOfHbc[x][z] = 0.;
				sumOfC[x][z] = 0.;
			}
		}
		for (int z = 0; z < 4; z++) {
				sumOfP[z] = 0.;
		}
	}

	std::cout << "::::::::::Global H matrix::::::::::" << std::endl;
	for (int i = 0; i < G.nH * G.nB; i++) {
		for (int j = 0; j < G.nH * G.nB; j++) {
			std::cout << std::setw(8) << std::setprecision(4) << globalH[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << "::::::::::Global P vector::::::::::" << std::endl;
	for (int j = 0; j < G.nH * G.nB; j++) {
		std::cout << std::setprecision(12) << globalP[j] << std::endl;
	}
	std::cout << std::endl;
	std::cout << "::::::::::Global C matrix::::::::::" << std::endl;
	for (int i = 0; i < G.nH * G.nB; i++) {
		for (int j = 0; j < G.nH * G.nB; j++) {
			std::cout << std::setw(8) << std::setprecision(4) << globalC[i][j];
		}
		std::cout << std::endl;
	}

	Solver::includeTimeH(G,globalH,globalC,globalP,dTau,25);
	std::cout << "::::::::::[H] = [H]+[C]/dT::::::::::" << std::endl;
	for (int i = 0; i < G.nH * G.nB; i++) {
		for (int j = 0; j < G.nH * G.nB; j++) {
			std::cout << std::setw(8) << std::setprecision(4) << globalH[i][j];
		}
		std::cout << std::endl;
	}
	
	for (int i = 0; i < 4; i++) {
		delete[] sumOfH[i];
	}
	delete[] sumOfH;

	for (int i = 0; i < 4; i++) {
		delete[] sumOfC[i];
	}
	delete[] sumOfC;

	for (int i = 0; i < 4; i++) {
		delete[] sumOfHbc[i];
	}
	delete[] sumOfHbc;

	for (int i = 0; i < G.nH * G.nB; i++) {
		delete[] globalH[i];
	}
	delete[] globalH;

	for (int i = 0; i < G.nH * G.nB; i++) {
		delete[] globalC[i];
	}
	delete[] globalC;

	delete[] globalP;
	delete[] sumOfP;
	return 0;
}
