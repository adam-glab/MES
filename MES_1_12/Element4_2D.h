#pragma once
#include "Gauss.h"

struct Element4_2D {
	int nIP;	// no. of integral pts
	gauss* g;	// vals of integral pts
	double*** N_shape = new double**[4]; // [4(walls)] [nIP] [4(N-func value)]
	double** N_ofIP = new double* [nIP*nIP]; // [nIP*nIP][4(N-func value)]
	double** dN_dE = new double* [nIP * nIP];// [nIP*nIP][4(N-func value)]
	double** dN_dn = new double* [nIP * nIP];// [nIP*nIP][4(N-func value)]

	void printGauss() const;
	
	Element4_2D(int n0) : nIP(n0) {

		// ====================== //
		//		nIP = 2			  //
		// ====================== //

		if (n0 == 2) {
			// initialize 2d arrays
			g = new gauss(n0);
			for (int i = 0; i < nIP*nIP; i++) {
				dN_dE[i] = new double[4];
				dN_dn[i] = new double[4];
				N_ofIP[i] = new double[4];
			}
			for (int i = 0; i < nIP * nIP; i++) {
				for (int j = 0; j < 4; j++) {
					//dN1,dN2,dN3,dN4
					dN_dE[i][j] = 0.;
					dN_dn[i][j] = 0.;
					//N1,N2,N3,N4
					N_ofIP[i][j] = 0.;
				}
			}
			// initialize 3d array
			for (int i = 0; i < 4; i++) {
				N_shape[i] = new double* [nIP];
				for (int j = 0; j < nIP; j++) {
					N_shape[i][j] = new double[4];
				}
			}
			for (int x = 0; x < nIP * nIP; x++) {
				for (int y = 0; y < nIP; y++) {
					for (int z = 0; z < 4; z++) {
						N_shape[x][y][z] = 0.;
					}
				}
			}
			// dN-vals of integral points for Jacobian and H-matrix
			/*
			* 01 11
			* 00 10
			*/
			int tmp = 0;
			for (int i = 0; i < nIP * nIP; i++) {
				tmp = i / nIP;
				dN_dE[i][0] = -0.25 * (1 - g->xC[tmp]);
				dN_dE[i][1] =  0.25 * (1 - g->xC[tmp]);
				dN_dE[i][2] =  0.25 * (1 + g->xC[tmp]);
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
				dN_dn[i][2] =  0.25 * (1 + g->xC[tmp]);
				dN_dn[i][3] =  0.25 * (1 - g->xC[tmp]);
			}

			// N-vals of integral points for Hbc-matrix
			/*
			* start from downmost surface and go counter clockwise (0:down -> 1:right -> 2:up -> 3:left)
			* counter clockwise rotation applies to gauss points,
			* code below includes conditions for shifts of xC indexes when going through surfaces 2:up and 3:left
			*/
			double tmp1 = 0.;
			double tmp2 = 0.;
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 2; j++) {
					if (i % 2 == 0) {
						tmp1 = g->xC[j % nIP];
						tmp2 = -1. + i*1.0;
						if (i > 1) {
							tmp1 = g->xC[(j + 1) % nIP];
						}
					}
					if (i % 2 != 0) {
						tmp1 = pow(-1, i);
						if (i < 2) {
							tmp1 = pow(-1, i + 1);
						}
						tmp2 = g->xC[(j+1) % nIP];
						if (i > 2) {
							tmp2 = g->xC[(j+1) % nIP];
						}
					}
					N_shape[i][j][0] = 0.25 * ((1 - tmp1) * (1 - tmp2));
					N_shape[i][j][1] = 0.25 * ((1 + tmp1) * (1 - tmp2));
					N_shape[i][j][2] = 0.25 * ((1 + tmp1) * (1 + tmp2));
					N_shape[i][j][3] = 0.25 * ((1 - tmp1) * (1 + tmp2));
				}
			}

			// N-vals of integral points for C-matrix
			for (int i = 0; i < nIP * nIP; i++) {
				if (i < 2) {
					tmp1 = g->xC[(2 - i) % 2]; // 0 1
				}
				if (i > 1) {
					tmp1 = g->xC[(1 + i) % 2]; // 1 0
				}
				tmp2 = g->xC[i / 2];
				N_ofIP[i][0] = 0.25 * ((1 - tmp1) * (1 - tmp2));
				N_ofIP[i][1] = 0.25 * ((1 + tmp1) * (1 - tmp2));
				N_ofIP[i][2] = 0.25 * ((1 + tmp1) * (1 + tmp2));
				N_ofIP[i][3] = 0.25 * ((1 - tmp1) * (1 + tmp2));
			}
		}

		// ====================== //
		//		nIP = 3			  //
		// ====================== //

		if (n0 == 3) {
			g = new gauss(n0); // init vals
			// initialize 2d arrays
			for (int i = 0; i < nIP * nIP; i++) {
				dN_dE[i] = new double[4];
				dN_dn[i] = new double[4];
				N_ofIP[i] = new double[4];
			}
			for (int i = 0; i < nIP * nIP; i++) {
				for (int j = 0; j < 4; j++) {
					//dN1,dN2,dN3,dN4
					dN_dE[i][j] = 0.;
					dN_dn[i][j] = 0.;
					//N1,N2,N3,N4
					N_ofIP[i][j] = 0.;
				}
			}
			// initialize 3d array
			for (int i = 0; i < 4; i++) {
				N_shape[i] = new double* [nIP];
				for (int j = 0; j < nIP; j++) {
					N_shape[i][j] = new double[4];
				}
			}
			for (int x = 0; x < 4; x++) {
				for (int y = 0; y < nIP; y++) {
					for (int z = 0; z < 4; z++) {
						N_shape[x][y][z] = 0.;
					}
				}
			}
			// dN-vals of integral points for Jacobian and H-matrix
			for (int i = 0; i < nIP * nIP; i++) {
				dN_dE[i][0] = -0.25 * (1 - g->xC[i % nIP]);
				dN_dE[i][1] =  0.25 * (1 - g->xC[i % nIP]);
				dN_dE[i][2] =  0.25 * (1 + g->xC[i % nIP]);
				dN_dE[i][3] = -0.25 * (1 + g->xC[i % nIP]);
			}
			for (int i = 0; i < nIP * nIP; i++) {
				dN_dn[i][0] = -0.25 * (1 - g->xC[i / nIP]);
				dN_dn[i][1] = -0.25 * (1 + g->xC[i / nIP]);
				dN_dn[i][2] = 0.25 * (1 + g->xC[i / nIP]);
				dN_dn[i][3] = 0.25 * (1 - g->xC[i / nIP]);
			}

			// N-vals of integral points for Hbc-matrix
			double tmp1 = 0.;
			double tmp2 = 0.;
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 3; j++) {
					if (i % 2 == 0) {
						tmp1 = g->xC[j % nIP];
						tmp2 = -1. + i * 1.0;
						if (i > 1) {
							tmp1 = g->xC[(j) % nIP];
						}
					}
					if (i % 2 != 0) {
						tmp1 = pow(-1, i);
						if (i < 2) {
							tmp1 = pow(-1, i + 1);
						}
						tmp2 = g->xC[(j) % nIP];
						if (i > 2) {
							tmp2 = g->xC[(j) % nIP];
						}
					}
					N_shape[i][j][0] = 0.25 * ((1 - tmp1) * (1 - tmp2));
					N_shape[i][j][1] = 0.25 * ((1 + tmp1) * (1 - tmp2));
					N_shape[i][j][2] = 0.25 * ((1 + tmp1) * (1 + tmp2));
					N_shape[i][j][3] = 0.25 * ((1 - tmp1) * (1 + tmp2));
				}
			}

			// N-vals of integral points for C-matrix
			/**
			* 02 12 22
			* 01 11 21
			* 00 10 20
			*/
			for (int i = 0; i < nIP * nIP; i++) {
				tmp1 = g->xC[i % 3];
				tmp2 = g->xC[i / 3];
				N_ofIP[i][0] = 0.25 * ((1 - tmp1) * (1 - tmp2));
				N_ofIP[i][1] = 0.25 * ((1 + tmp1) * (1 - tmp2));
				N_ofIP[i][2] = 0.25 * ((1 + tmp1) * (1 + tmp2));
				N_ofIP[i][3] = 0.25 * ((1 - tmp1) * (1 + tmp2));
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
	
		for (int x = 0; x < 4; x++) {
			for (int y = 0; y < nIP; y++) {
				delete[] N_shape[x][y];
			}
			delete[] N_shape[x];
		}
		delete[] N_shape;
	}
};