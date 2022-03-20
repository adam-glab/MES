#pragma once
#include "Solver.h"

using namespace std;

struct Data {
	int
		nN,
		nE;
	// factor values
	double 
		k, 
		alpha, 
		t_env,
		c, 
		ro, 
		simTime,
		T0,
		dTau;

	int params[10] = { 0 };
	double** d_nodes; // params[8] x 2
	int** d_elems; // params[9] x 4
	int* bcs; // params
	int bc_counter = 0;

	void getParametersTest(ifstream& MyReadFile);
	void setParameters(Grid &G);
};
