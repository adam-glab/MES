#pragma once
#include "Solver.h"

using namespace std;

struct Data {
	// pass to Grid
	int 
		nIP = 2;
	double 
		H = 0.1,
		B = 0.1;
	int 
		nH = 4, 
		nB = 4;
	// factor values
	double 
		k = 25., 
		alpha = 25., 
		t_env = 1200,
		c = 700, 
		ro = 7800, 
		simTime = 500,
		T0 = 100.,
		dTau = 50;

	int params[10] = { 0 };
	double** d_nodes; // params[8] x 2
	int** d_elems; // params[9] x 4
	int* bcs; // params
	int bc_counter = 0;

	void getParametersTest(ifstream& MyReadFile);
	void setParameters(Grid* G);
	
};
