#pragma once
#include "Grid.h"
#include "Element4_2D.h"
#include "Node.h"
#include "Jacobian.h"

struct Solver {
	static void calcJacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, Grid &G);
	static void calcHTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, Grid &G, double k, double detJ, double** sumArray, double** globalArray);
	static void calcHbcTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, Grid &G, double alpha, double** sumArray);
	static void calcPTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, Grid &G, double alpha, double t_env, double* sumArray);
	static void calcCTest(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D* E, Grid &G, double c, double ro, double detJ, double** sumArray, double** globalArray);
	static void includeTimeH(Grid &G, double** matrixH, double** matrixC, double* vectorP, double dTau);
	static double* gaussScheme(double** matrix, double* vector, int size);
	static void calcNodeTemp(double** Hmatrix, double** Cmatrix, double* Pvector, Grid &G, double simTime, double timeStep);
	static void getMinMax(double* arr, int size, int step);
	static void solveFEM(int nIP, double k, double alpha, double t_env, double c, double ro, double dTau, double simTime, double T0, Grid& G, jacobian& J, jacobian& J_inv, Element4_2D& E);
};
