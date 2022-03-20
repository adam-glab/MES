#pragma once
#include "Node.h"
#include "Element.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

struct Grid {
	double H, B; // heigth, length
	int nH, nB; // nodes /heigth, nodes /length
	int nN;	// node count
	int nE;	// element count
	double T0; // T0 in node

// Global array of H-vals in nodes [nN x nN]
	double** globalH = new double* [nN];
// Global array of P-vals in nodes {1 x nN}
	double* globalP = new double[nN];
// Global C-matrix [nN x nN]
	double** globalC = new double* [nN];

	Node* nodes;
	Element* elements;

// Constructor to be used in the case of file input
	Grid(int noN, int noE) :  nN(noN), nE(noE) {
		H = 0.;
		B = 0.;
		nH = nB = 0;
		T0 = 0.;

		for (int i = 0; i < nN; i++) {
			globalH[i] = new double[nN];
		}
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				globalH[i][j] = 0.;
			}
		}

		for (int i = 0; i < nN; i++) {
			globalP[i] = 0.;
		}

		for (int i = 0; i < nN; i++) {
			globalC[i] = new double[nN];
		}
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				globalC[i][j] = 0.;
			}
		}
	};

// Constructor to be used in the case of user input
	Grid(double H0, double B0, int noH, int noB, int noN, double temp0) : H(H0), B(B0), nH(noH), nB(noB), nN(noN), T0(temp0) {

		nE = (noH - 1) * (noB - 1);
		nodes = new Node[nN];
		elements = new Element[nE];

		for (int i = 0; i < nN; i++) {
			globalH[i] = new double[nN];
		}
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				globalH[i][j] = 0.;
			}
		}

		for (int i = 0; i < nN; i++) {
			globalP[i] = 0.;
		}

		for (int i = 0; i < nN; i++) {
			globalC[i] = new double[nN];
		}
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				globalC[i][j] = 0.;
			}
		}

		// increment dx = B / (nB - 1)
		// increment dy = H / (nH - 1)

		int increment = 0; // step to preserve proper node value going upwards in the grid

		// #1 - Fill in node values
		for (int i = 0; i < nN; i++) {
			if (i == 0) {
				increment = 0;
			}
			else if (i % nH == 0) {
				increment++;
			}
			nodes[i].x = increment * B / (nB - 1.0); // 0*dx, 1*dx, 2*dx ..
			nodes[i].y = (i % nH) * H / (nH - 1.0);

			if (i % nH == 0 || i % nH == nH - 1 || i < nH || i > nN - nH) {
				nodes[i].BC = true;
			}
			else {
				nodes[i].BC = false;
			}
			nodes[i].t0 = T0;
		}

		// #2 - Fill in element ids
		increment = 0; // mirror step in element array
		for (int i = 0; i < nE; i++) {
			if (i == 0) {
				increment = 0;
			}
			else if (i % (nH - 1) == 0) { // nH element goes to another row
				increment++;
			}
			elements[i].ID[0] = i + 1 + increment;
			elements[i].ID[1] = i + 1 + nH + increment;
			elements[i].ID[2] = i + 2 + nH + increment;
			elements[i].ID[3] = i + 2 + increment;

		}
	};

	~Grid() {
		std::cout << "Grid destructor initialized. End of programme.\n";
		for (int i = 0; i < nN; i++) {
			delete[] globalH[i];
		}
		delete[] globalH;
		for (int i = 0; i < nN; i++) {
			delete[] globalC[i];
		}
		delete[] globalC;
		delete[] globalP;
	}

	void printNodes() const;
	void printElements() const;
	static void printGlobalH(double** globalH, Grid G);
	static void printGlobalC(double** globalC, Grid G);
	static void printGlobalP(double* globalP, Grid G);
};