#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>

struct gauss {
	// gauss x const
	double* xC;
	// weight points
	double* wP;
	int N;
	gauss(int N0) : N(N0) {
		if (N == 2) {
			xC = new double[N];
			xC[0] = -sqrt(1.0 / 3.0);
			xC[1] = sqrt(1.0 / 3.0);

			wP = new double[N];
			wP[0] = 1.0;
			wP[1] = 1.0;
		}
		else if (N == 3) {
			xC = new double[N];
			xC[0] = -sqrt(3.0 / 5.0);
			xC[1] = 0.0;
			xC[2] = sqrt(3.0 / 5.0);

			wP = new double[N];
			wP[0] = 5./9.;
			wP[1] = 8./9.;
			wP[2] = 5./9.;
		}
		else {
			throw std::out_of_range("Illegal integral point count passed to constructor\n");
			// further wP,xC
		}
	}

	~gauss() {
		delete[] xC;
		delete[] wP;
	}
};