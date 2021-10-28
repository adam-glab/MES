#pragma once

struct gauss {

	// gauss x const
	double* xC;

	// weight points
	double* wP;

	int N;
	gauss(int N0) : N(N0) {
		if (N == 2) {
			xC = new double[N];
			xC[0] = -0.577350;
			xC[1] = 0.577350;

			wP = new double[N];
			wP[0] = 1.0;
			wP[1] = 1.0;
		}
		else if (N == 3) {
			xC = new double[N];
			xC[0] = -0.774597;
			xC[1] = 0.0;
			xC[2] = 0.774597;

			wP = new double[N];
			wP[0] = 0.555555;
			wP[1] = 0.888888;
			wP[2] = 0.555555;
		}
		else {
			// further wP,xC
		}
	}

	~gauss() {
		std::cout << "Gauss - destructor...\n";
		delete[] xC;
		delete[] wP;
	}
};