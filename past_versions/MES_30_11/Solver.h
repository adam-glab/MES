#pragma once
#include "Grid.h"
#include "Element4_2D.h"
#include "Node.h"
#include "Jacobian.h"

struct Solver {
	static void calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G);
	static void calcHTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double k, double detJ, double** sumArray, double** globalArray);
	static void calcHbcTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double alpha, double** sumArray);
	static void calcPTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, grid G, double alpha, double t_env, double* sumArray);
};
