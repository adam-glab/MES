#pragma once
#include "Gauss.h"

struct Element4_2D {
	int nIP;	// no. of integral pts
	gauss* g;	// vals of integral pts
	gauss* g_hbc; // vals of integral pts for hbc matrix
	double ksi, eta;
	double dNdE[4][4] = { 0. };
	double dNdn[4][4] = { 0. };
	double N[2][4] = { 0. };
	double** dN_dE = new double* [nIP];
	double** dN_dn = new double* [nIP];
	
	Element4_2D(int n0) : nIP(n0) {

		if (n0 == 4) {
			g = new gauss(n0 / 2);
			
			std::cout << "\tdNdE\n";
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
				std::cout << "P" << i + 1 << '\t';
				for (int j = 0; j < n0; j++) {
					std::cout << dNdE[i][j] << '\t';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
			std::cout << "\tdNdn\n";
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
				std::cout << "P" << i + 1 << '\t';
				for (int j = 0; j < n0; j++) {
					std::cout << dNdn[i][j] << '\t';
				}
				std::cout << std::endl;
			}

		}

		if (n0 == 9) {
			g = new gauss(n0 / 3); // init vals
			for (int i = 0; i < nIP; i++) {
				dN_dE[i] = new double[nIP];
				dN_dn[i] = new double[nIP];
			}
			for (int i = 0; i < nIP; i++) {
				for (int j = 0; j < 4; j++) { //N1,N2,N3,N4
					dN_dE[i][j] = 0.;
					dN_dn[i][j] = 0.;
				}
			}
			std::cout << "\tdNdE\n";
			int tmp = 0;
			for (int i = 0; i < n0; i++) {
				if (i == 0 || i == 3 || i == 6) {
					tmp = 0; // - sq3/5
				}
				if(i == 2 || i == 5 || i == 8) {
					tmp = 2; // 0
				}
				if(i== 1 || i == 4 || i == 7) {
					tmp = 1; // sq3/5
				}
				dN_dE[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dE[i][1] = 0.25 * (1 - g->xC[tmp]);
				dN_dE[i][2] = 0.25 * (1 + g->xC[tmp]);
				dN_dE[i][3] = -0.25 * (1 + g->xC[tmp]);
			}
			for (int i = 0; i < n0; i++) {
				std::cout <<"P" << i + 1 << '\t';
				for (int j = 0; j < 4; j++) {
					std::cout << dN_dE[i][j] << '\t';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
			std::cout << "\tdNdn\n";
			tmp = 0;
			for (int i = 0; i < n0; i++) {
				if (i % 9 < 3) {
					tmp = 0; // - sq3/5
				}
				if (i % 9 > 2 && i % 9 < 6) {
					tmp = 1; // 0
				}
				if (i % 9 > 5) {
					tmp = 2; // sq3/5
				}
				dN_dn[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dn[i][1] = -0.25 * (1 + g->xC[tmp]);
				dN_dn[i][2] = 0.25 * (1 + g->xC[tmp]);
				dN_dn[i][3] = 0.25 * (1 - g->xC[tmp]);
			}
			
			for (int i = 0; i < n0; i++) {
				std::cout << "P" << i + 1 << '\t';
				for (int j = 0; j < 4; j++) {
					std::cout << dN_dn[i][j] << '\t';
				}
				std::cout << std::endl;
			}

			//g_hbc = new gauss(n0 / 2); // init vals
			//for (int i = 0; i < n0/2; i++) {
			//	if(i == 0)
			//	N[i][0] = 0.25 * (1 - g->xC[0]) * (1 - (-1));
			//	N[i][1] = 0.25 * (1 + g->xC[0]) * (1 - (-1));
			//	N[i][2] = 0.25 * (1 + g->xC[0]) * (1 + (-1));
			//	N[i][3] = 0.25 * (1 - g->xC[0]) * (1 + (-1));
			//}
		}
		
	};
};