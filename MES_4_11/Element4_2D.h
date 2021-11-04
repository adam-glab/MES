#pragma once
#include "Gauss.h"

struct Element4_2D {
	int nIP;	// no. of integral pts
	gauss* g;	// vals of integral pts
	double ksi, eta;
	double dNdE[4][4] = { 0. };
	double dNdn[4][4] = { 0. };
	Element4_2D(int n0) : nIP(n0) {

		if (n0 == 4) {
			g = new gauss(n0 / 2); // init vals
			std::cout << "dN1dE\tdN2dE\tdN3dE\tdN4dE\n";
			int tmp = 0;
			for (int i = 0; i < n0; i++) {
				if (i == 0 || i == 1) {
					tmp = 0;
				}
				else {
					tmp = 1;
				}
				dNdE[i][0] = -0.25 * (1 - g->xC[tmp]);
				dNdE[i][1] = 0.25 * (1 - g->xC[tmp]);
				dNdE[i][2] = 0.25 * (1 + g->xC[tmp]);
				dNdE[i][3] = -0.25 * (1 + g->xC[tmp]);
			}
			for (int i = 0; i < n0; i++) {
				for (int j = 0; j < n0; j++) {
					std::cout << dNdE[i][j] << '\t';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
			std::cout << "dN1dn\tdN2dn\tdN3dn\tdN4dn\n";
			tmp = 0;
			for (int i = 0; i < n0; i++) {
				if (i == 0 || i == 3) {
					tmp = 0;
				}
				else {
					tmp = 1;
				}
				dNdn[i][0] = -0.25 * (1 - g->xC[tmp]);
				dNdn[i][1] = -0.25 * (1 + g->xC[tmp]);
				dNdn[i][2] = 0.25 * (1 + g->xC[tmp]);
				dNdn[i][3] = 0.25 * (1 - g->xC[tmp]);
			}
			
			for (int i = 0; i < n0; i++) {
				for (int j = 0; j < n0; j++) {
					std::cout << dNdn[i][j] << '\t';
				}
				std::cout << std::endl;
			}
		}
	};

	~Element4_2D() {
	}
};