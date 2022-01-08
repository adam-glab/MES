//#include "Solver.h"
#include "Data.h"

//****************************************************************
// In order to change data read mode (be it file stream or user values)
// uncomment/comment marked blocks in main() below
//****************************************************************


//****************************************************************
// SET THE AMOUNT OF INTEGRAL POINTS
//****************************************************************
const int nIP = 3;

// Pass dimensions to grid
const double 
H = 0.1, B = 0.1;
const int 
nH = 4, nB = 4;
// Pass to functions: factor values
double 
k = 25., 
alpha = 300., 
t_env = 1200, 
c = 700, 
ro = 7800, 
dTau = 50, 
simTime = 500, 
T0 = 100.;

Grid G;
jacobian J;
jacobian J_inv;
Element4_2D E(nIP);

int main() {

//****************************************************************
// READ VALUES FROM FILE
//****************************************************************
	Data d;
	ifstream MyData("Test1_4_4_my.txt");

	//****************************************************************
	// comment G = { H,B,nH,nB,T0 }; to avoid data overwrite
	// uncomment code below this comment to enable file stream
	//****************************************************************

	/*d.getParametersTest(MyData);
	d.setParameters(&G);
	simTime = d.simTime;
	dTau = d.dTau;
	k = d.k;
	alpha = d.alpha;
	t_env = d.t_env;
	T0 = d.T0;
	ro = d.ro;
	c = d.c;*/

//****************************************************************
// USER CHANGES VALUES OF CONSTS
//****************************************************************

	//****************************************************************
	// uncomment G = { H,B,nH,nB,T0 }; to apply global data to grid
	// comment code in previous section to avoid factor overwrite from file stream
	//****************************************************************

	G = { H,B,nH,nB,T0 };

//****************************************************************
// Grid management
//****************************************************************
	//G.printNodes();
    //G.printElements();
	//E.printElementData();
	
//****************************************************************
// Global array of H-vals in nodes nN x nN
//****************************************************************
	double** globalH = new double* [G.nN];
	for (int i = 0; i < G.nN; i++) {
		globalH[i] = new double[G.nN];
	}
	for (int i = 0; i < G.nN; i++) {
		for (int j = 0; j < G.nN; j++) {
			globalH[i][j] = 0.;
		}
	}
//****************************************************************
// Global array of P-vals in nodes {1 x nN}
//**************************************************************** 
	double* globalP = new double [G.nN];
	for (int i = 0; i < G.nN; i++) {
		globalP[i] = 0.;
	}
//****************************************************************
// Global C-matrix [nN x nN]
//****************************************************************
	double** globalC = new double* [G.nN];
	for (int i = 0; i < G.nN; i++) {
		globalC[i] = new double[G.nN];
	}
	for (int i = 0; i < G.nN; i++) {
		for (int j = 0; j < G.nN; j++) {
			globalC[i][j] = 0.;
		}
	}
//****************************************************************
// sumOfH -> 2d array to hold H matrixes of each element
// H matrix in each integral point is calculated inside calcHTest function as a local matrix and then added to the sum
//****************************************************************
	double** sumOfH = new double*[4];
	for (int i = 0; i < 4; i++) {
		sumOfH[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfH[i][j] = 0.;
		}
	}
//****************************************************************
// sumOfP -> 2d array to hold P matrixes of each element
//****************************************************************
	double* sumOfP = new double[4];
	for (int i = 0; i < 4; i++) {
		sumOfP[i] = 0;
	}
//****************************************************************
// sumOfHbc -> 2d array to hold Hbc matrixes of each element
//****************************************************************
	double** sumOfHbc = new double* [4];
	for (int i = 0; i < 4; i++) {
		sumOfHbc[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfHbc[i][j] = 0.;
		}
	}
//****************************************************************
// sumOfC -> 2d array to hold C matrixes of each element
//****************************************************************
	double** sumOfC = new double* [4];
	for (int i = 0; i < 4; i++) {
		sumOfC[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sumOfC[i][j] = 0.;
		}
	}
//****************************************************************
// MAIN LOOP
//****************************************************************
	for (int i = 0; i < G.nE; i++) {
		//std::cout << "::::::::ELEMENT " << i + 1 << "::::::::" << std::endl;
		for (int j = 0; j < nIP*nIP; j++) {
			Solver::calcJacobian(i, j, &J, &J_inv, &E, G);
			double detJ = J.j_matrix[0][0] * J.j_matrix[1][1] - J.j_matrix[0][1] * J.j_matrix[1][0];
			Solver::calcHTest(i, j, &J, &J_inv, &E, G, k, detJ, sumOfH, globalH);
			Solver::calcCTest(i, j, &J, &J_inv, &E, G, c, ro, detJ, sumOfC, globalC);
		}
	//****************************************************************
	// For each wall of element calculate Hbc and P vector
	//****************************************************************
		for (int x = 0; x < 4; x++) {
			Solver::calcHbcTest(i, x, &J, &J_inv, &E, G, alpha, sumOfHbc);
			Solver::calcPTest(i, x, &J, &J_inv, &E, G, alpha, t_env, sumOfP);
		}
	//****************************************************************
	// Add to global H matrix and C matrix
	//****************************************************************
		for (int x = 0; x < 4; x++) {
			for (int z = 0; z < 4; z++) {
				globalH[G.elements[i].ID[x] - 1][G.elements[i].ID[z] - 1] += sumOfH[x][z] + sumOfHbc[x][z];
				globalC[G.elements[i].ID[x] - 1][G.elements[i].ID[z] - 1] += sumOfC[x][z];
			}
		}

		//Element4_2D::printH(sumOfH, i);
		//Element4_2D::printHbc(sumOfHbc, i);
		//Element4_2D::printC(sumOfC, i);

	//****************************************************************
	// Add to global P vector
	//****************************************************************
		for (int z = 0; z < 4; z++) {
			globalP[G.elements[i].ID[z] - 1] += sumOfP[z];
		}

		//Element4_2D::printP(sumOfP, i);

	//****************************************************************
	// Reset arrays for next element
	//****************************************************************
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

	//Grid::printGlobalH(globalH, G);
	//Grid::printGlobalP(globalP, G);
	//Grid::printGlobalC(globalC, G);

//****************************************************************
// Operations on completed H,C matrixes and P, T0, T1 vectors
//****************************************************************
	Solver::includeTimeH(G, globalH, globalC, globalP, dTau);
	Solver::calcNodeTemp(globalH, globalC, globalP, G, simTime, dTau);

//****************************************************************
// Free memory
//****************************************************************
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
	for (int i = 0; i < G.nN; i++) {
		delete[] globalH[i];
	}
	delete[] globalH;
	for (int i = 0; i < G.nN; i++) {
		delete[] globalC[i];
	}
	delete[] globalC;
	delete[] globalP;
	delete[] sumOfP;

	MyData.close();

	return 0;
}