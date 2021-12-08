#include "Solver.h"

// pass to grid
const int nIP = 3;
const double H = 0.1, B = 0.1;
const int nH = 4, nB = 4;
// factor values
const double k = 25., alpha = 300., t_env = 1200, c = 700, ro = 7800, dTau = 50, simTime = 500, T0 = 100.;

int main() {

	jacobian J;
	jacobian J_inv;
	Element4_2D E(nIP);
	const Grid G(H, B, nH, nB, T0);
//****************************************************************
// Grid management
//****************************************************************
	//G.printNodes();
	//G.printElements();
	//E.printElementData();

//****************************************************************
// MAIN LOOP
//****************************************************************
	for (int i = 0; i < G.nE; i++) {
		//std::cout << "::::::::ELEMENT " << i + 1 << "::::::::" << std::endl;
		for (int j = 0; j < nIP*nIP; j++) {
			Solver::calcJacobian(i, j, &J, &J_inv, &E, G);
			double detJ = J.j_matrix[0][0] * J.j_matrix[1][1] - J.j_matrix[0][1] * J.j_matrix[1][0];
			Solver::calcHTest(i, j, &J, &J_inv, &E, G, k, detJ, E.sumOfH, G.globalH);
			Solver::calcCTest(i, j, &J, &J_inv, &E, G, c, ro, detJ, E.sumOfC, G.globalC);
		}
	//****************************************************************
	// For each wall of element calculate Hbc and P vector
	//****************************************************************
		for (int x = 0; x < 4; x++) {
			Solver::calcHbcTest(i, x, &J, &J_inv, &E, G, alpha, E.sumOfHbc);
			Solver::calcPTest(i, x, &J, &J_inv, &E, G, alpha, t_env, E.sumOfP);
		}
	//****************************************************************
	// Add to global H matrix and C matrix
	//****************************************************************
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				G.globalH[G.elements[i].ID[x] - 1][G.elements[i].ID[z] - 1] += E.sumOfH[x][z] + E.sumOfHbc[x][z];
				G.globalC[G.elements[i].ID[x] - 1][G.elements[i].ID[z] - 1] += E.sumOfC[x][z];
			}
		}

		Element4_2D::printH(E.sumOfH, i);
		Element4_2D::printHbc(E.sumOfHbc, i);
		Element4_2D::printC(E.sumOfC, i);

	//****************************************************************
	// Add to global P vector
	//****************************************************************
		for (int z = 0; z < 4; z++) {
			G.globalP[G.elements[i].ID[z] - 1] += E.sumOfP[z];
		}

		Element4_2D::printP(E.sumOfP, i);

	//****************************************************************
	// Reset arrays for next element
	//****************************************************************
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				E.sumOfH[x][z] = 0.;
				E.sumOfHbc[x][z] = 0.;
				E.sumOfC[x][z] = 0.;
			}
		}
		for (int z = 0; z < 4; z++) {
			E.sumOfP[z] = 0.;
		}
	}

	//Grid::printGlobalH(G.globalH, G);
	//Grid::printGlobalP(G.globalP, G);
	//Grid::printGlobalC(G.globalC, G);

//****************************************************************
// Operations on completed H,C matrixes and P, T0, T1 vectors
//****************************************************************
	Solver::includeTimeH(G, G.globalH, G.globalC, G.globalP, dTau);
	Solver::calcNodeTemp(G.globalH, G.globalC, G.globalP, G, simTime, dTau);

	return 0;
}