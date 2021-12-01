#include "Solver.h"

void Solver::calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G) {

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

void Solver::calcHTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double k, double detJ, double** sumArray, double** globalArray) {

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

void Solver::calcHbcTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double alpha, double** sumArray) {

	//std::cout << ":::::Hbc for IP" << nIP + 1 << ":::::" << std::endl;
	double my_detJ = 0.;
	double NNT[4][4] = { 0. }; // hold sum(N_transposed * N) * alpha * detJ
	double array_detJ[4] = { 0. };
	for (int i = 0; i < 4; i++) {
		// is boundary condition met
		// (i+1)%4, i%4 -> start from downmost side and go counter-clockwise
		if (G.nodes[G.elements[nE].ID[(i + 1) % 4] - 1].BC == false || G.nodes[G.elements[nE].ID[(i) % 4] - 1].BC == false) {
			array_detJ[i] = 0.;
		}
		else {
			// length of side
			array_detJ[i] = 0.5 * sqrt(pow((G.nodes[G.elements[nE].ID[(i + 1) % 4] - 1].x - (G.nodes[G.elements[nE].ID[(i) % 4] - 1].x)), 2)
				+ pow((G.nodes[G.elements[nE].ID[(i + 1) % 4] - 1].y - (G.nodes[G.elements[nE].ID[(i) % 4] - 1].y)), 2)
			);
		}
		// std::cout << array_detJ[i] << " ";
	}
	//std::cout << std::endl;

	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			NNT[z][x] = alpha * array_detJ[nIP % (E->nIP + nIP)] * ((E->N[nIP][0][z] * E->N[nIP][0][x]) + (E->N[nIP][1][z] * E->N[nIP][1][x]));
		}
	}

	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			double tmp = NNT[z][x];
			sumArray[z][x] += tmp;
		}
	}
	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << std::setw(12) << NNT[i][j];
		}
		std::cout << std::endl;
	}*/
}

void Solver::calcPTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double alpha, double t_env, double* sumArray) {

	std::cout << ":::::P for IP" << nIP + 1 << ":::::" << std::endl;
	double my_detJ = 0.;
	double NT[4] = { 0. }; // hold sum(N_transposed * t_env) * alpha * detJ
	double array_detJ[4] = { 0. };
	for (int i = 0; i < 4; i++) {
		// is boundary condition met
		if (G.nodes[G.elements[nE].ID[(i + 1) % 4] - 1].BC == false || G.nodes[G.elements[nE].ID[(i) % 4] - 1].BC == false) {
			array_detJ[i] = 0.;
		}
		else {
			// length of side
			array_detJ[i] = 0.5 * sqrt(pow((G.nodes[G.elements[nE].ID[(i + 1) % 4] - 1].x - (G.nodes[G.elements[nE].ID[(i) % 4] - 1].x)), 2)
				+ pow((G.nodes[G.elements[nE].ID[(i + 1) % 4] - 1].y - (G.nodes[G.elements[nE].ID[(i) % 4] - 1].y)), 2)
			);
		}
		//std::cout << array_detJ[i] << " ";
	}
	std::cout << std::endl;

	for (int x = 0; x < 4; x++) {
		NT[x] = alpha * array_detJ[nIP % (E->nIP + nIP)] * ((E->N[nIP][0][x] * t_env) + (E->N[nIP][1][x] * t_env));
	}

	for (int x = 0; x < 4; x++) {
		double tmp = NT[x];
		sumArray[x] += tmp;
	}
	for (int j = 0; j < 4; j++) {
		std::cout << std::setw(12) << NT[j];
	}
	std::cout << std::endl;
}

void Solver::calcCTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double c, double ro, double detJ, double** sumArray, double** globalArray){
	

	double localC[4][4] = { 0. };
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			localC[z][x] = c * ro * detJ * E->N_ofIP[nIP][z] * E->N_ofIP[nIP][x];
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
	for (int z = 0; z < 4; z++) {
		for (int x = 0; x < 4; x++) {
			double tmp = localC[z][x];
			sumArray[z][x] += tmp;
		}
	}
}

void Solver::includeTimeH(grid G, double** matrixH, double** matrixC, double* vectorP, double dTau, double T0){
	for (int i = 0; i < G.nH * G.nB; i++) {
		for (int j = 0; j < G.nH * G.nB; j++) {
			//std::cout << std::setw(8) << std::setprecision(4) << matrixH[i][j];
			matrixH[i][j] = matrixH[i][j] + (matrixC[i][j] / dTau);
		}
		//std::cout << std::endl;
	}
	/*std::cout << "::::::::::{P} = {P}+{[C]/dT}*{T0}::::::::::" << std::endl;
	for (int j = 0; j < G.nH * G.nB; j++) {
		std::cout << std::setprecision(12) << vectorP[j] << std::endl;
	}*/
}
