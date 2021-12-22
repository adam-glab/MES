#pragma once
#include "Gauss.h"

struct Element4_2D {
	int nIP;	// no. of integral pts
	gauss* g;	// vals of integral pts
	gauss* g_hbc; // vals of integral pts for hbc matrix
	double N[4][2][4] = { 0. };
	double** dN_dE = new double* [nIP * nIP];
	double** dN_dn = new double* [nIP * nIP];

	void printGauss() const;
	
	Element4_2D(int n0) : nIP(n0) {

		if (n0 == 2) {
			g = new gauss(n0);
			for (int i = 0; i < nIP*nIP; i++) {
				dN_dE[i] = new double[nIP * nIP];
				dN_dn[i] = new double[nIP * nIP];
			}
			for (int i = 0; i < nIP * nIP; i++) {
				for (int j = 0; j < 4; j++) { //N1,N2,N3,N4
					dN_dE[i][j] = 0.;
					dN_dn[i][j] = 0.;
				}
			}
			
			int tmp = 0;
			for (int i = 0; i < nIP * nIP; i++) {
				tmp = i / nIP;
				dN_dE[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dE[i][1] = 0.25 * (1 - g->xC[tmp]);
				dN_dE[i][2] = 0.25 * (1 + g->xC[tmp]);
				dN_dE[i][3] = -0.25 * (1 + g->xC[tmp]);
			}
			tmp = 0;
			for (int i = 0; i < nIP * nIP; i++) {
				if (i == 0 || i == 3) {
					tmp = 0;
				}
				else {
					tmp = 1;
				}
				dN_dn[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dn[i][1] = -0.25 * (1 + g->xC[tmp]);
				dN_dn[i][2] = 0.25 * (1 + g->xC[tmp]);
				dN_dn[i][3] = 0.25 * (1 - g->xC[tmp]);
			}


			double tmp1 = 0.;
			double tmp2 = 0.;
			g_hbc = new gauss(n0);
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 2; j++) {
					if (i % n0 == 0) {
						tmp1 = g_hbc->xC[j % nIP];
						tmp2 = -1. + i*1.0;
						if (i > 1) {
							tmp1 = g_hbc->xC[(j + 1) % nIP];
						}
					}
					if (i % n0 != 0) {
						tmp1 = pow(-1, i);
						if (i < 2) {
							tmp1 = pow(-1, i + 1);
						}
						tmp2 = g_hbc->xC[(j+1) % nIP];
						if (i > 2) {
							tmp2 = g_hbc->xC[(j + 1) % nIP];
						}
					}
					N[i][j][0] = 0.25 * ((1 - tmp1) * (1 - tmp2));
					N[i][j][1] = 0.25 * ((1 + tmp1) * (1 - tmp2));
					N[i][j][2] = 0.25 * ((1 + tmp1) * (1 + tmp2));
					N[i][j][3] = 0.25 * ((1 - tmp1) * (1 + tmp2));
				}
			}
		}

		if (n0 == 3) {
			g = new gauss(n0); // init vals
			for (int i = 0; i < nIP*nIP; i++) {
				dN_dE[i] = new double[nIP * nIP];
				dN_dn[i] = new double[nIP * nIP];
			}
			for (int i = 0; i < nIP * nIP; i++) {
				for (int j = 0; j < 4; j++) { //N1,N2,N3,N4
					dN_dE[i][j] = 0.;
					dN_dn[i][j] = 0.;
				}
			}
			int tmp = 0;
			for (int i = 0; i < nIP * nIP; i++) {
				if (i % nIP == 0) {
					tmp = 0;
				}
				if (i % nIP == 2) {
					tmp = 2;
				}
				if (i % nIP == 1) {
					tmp = 1;
				}
				dN_dE[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dE[i][1] = 0.25 * (1 - g->xC[tmp]);
				dN_dE[i][2] = 0.25 * (1 + g->xC[tmp]);
				dN_dE[i][3] = -0.25 * (1 + g->xC[tmp]);
			}
			tmp = 0;
			for (int i = 0; i < nIP * nIP; i++) {
				tmp = i / nIP;
				dN_dn[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dn[i][1] = -0.25 * (1 + g->xC[tmp]);
				dN_dn[i][2] = 0.25 * (1 + g->xC[tmp]);
				dN_dn[i][3] = 0.25 * (1 - g->xC[tmp]);
			}
		}
		
	};
	~Element4_2D() {
		for (int i = 0; i < nIP * nIP; i++) {
			delete[] dN_dE[i];
			delete[] dN_dn[i];
		}
		delete[] dN_dE;
		delete[] dN_dn;
	}
};