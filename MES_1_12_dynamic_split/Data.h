#pragma once
#include "Solver.h"

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
		dTau = 50;
};
