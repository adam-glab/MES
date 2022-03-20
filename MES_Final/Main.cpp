#include "Data.h"

int nIP = 3;

// Pass dimensions to grid
double
	H = 0.1, B = 0.1;
int
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


jacobian J;
jacobian J_inv;
Element4_2D E(nIP);


int main() {
//****************************************************************
// READ VALUES FROM FILE
//****************************************************************
	//****************************************************************
	// comment code in USER CHANGES VALUES OF CONSTS section to avoid data overwrite
	// uncomment code below this comment to enable file stream
	//****************************************************************

	/*Data dataModel;
	ifstream dataStream("Test1_4_4_my.txt");

	dataModel.getParametersTest(dataStream);
	Grid G(dataModel.nN, dataModel.nE);
	dataModel.setParameters(G);

	simTime = dataModel.simTime;
	dTau = dataModel.dTau;
	k = dataModel.k;
	alpha = dataModel.alpha;
	t_env = dataModel.t_env;
	T0 = dataModel.T0;
	ro = dataModel.ro;
	c = dataModel.c;*/

//****************************************************************
// USER CHANGES VALUES OF CONSTS
//****************************************************************
	//****************************************************************
	// uncomment code below to apply global data to grid
	// comment code in previous section to avoid factor overwrite from file stream
	//****************************************************************

	Grid G( H, B, nH, nB, nH*nB, T0 );

//****************************************************************
// Grid management
//****************************************************************
	//G.printNodes();
	//G.printElements();
	//E.printElementData();

//****************************************************************
// SOLVE FEM
//****************************************************************

	Solver::solveFEM(nIP, k, alpha, t_env, c, ro, dTau, simTime, T0, G, J, J_inv, E);

//****************************************************************
// Operations on completed H,C matrixes and P, T0, T1 vectors
//****************************************************************
	Solver::includeTimeH(G, G.globalH, G.globalC, G.globalP, dTau);
	Solver::calcNodeTemp(G.globalH, G.globalC, G.globalP, G, simTime, dTau);
	return 0;
}