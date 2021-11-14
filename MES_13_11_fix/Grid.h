#pragma once
#include "Node.h"
#include "Element.h"

struct grid {
	double H, B; // heigth, length
	int nH, nB; // nodes /heigth, nodes /length
	int nN;	// node count
	int nE;	// element count
	node* nodes;
	element* elements;

	grid(double H0, double B0, int noH, int noB) : H(H0), B(B0), nH(noH), nB(noB) {

		nN = noH * noB;
		nE = (noH - 1) * (noB - 1);
		nodes = new node[nN];
		elements = new element[nE];

		// dx = B / (nB - 1)
		// dy = H / (nH - 1)

		int increment = 0; // step to preserve proper node value going upwards in the grid

		// #1 - Fill in node values
		for (int i = 0; i < nN; i++) {
			if (i == 0) {
				increment = 0;
			}
			else if (i % nH == 0) {
				increment++;
			}
			nodes[i].x = increment * B / (nB - 1); // 0*dx, 1*dx, 2*dx ..
			nodes[i].y = (i % nH) * H / (nH - 1);

			if (i % nH == 0 || i % nH == nH - 1 || i < nH || i > nN - nH) {
				nodes[i].BC = true;
			}
			else {
				nodes[i].BC = false;
			}
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
};