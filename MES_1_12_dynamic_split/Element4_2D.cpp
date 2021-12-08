#include "Element4_2D.h"

void Element4_2D::printElementData() const{

	for (int i = 0; i < 4; i++) {
		std::cout << std::setw(16) << "dNdE" << i + 1;
	}
	std::cout << std::endl;
	for (int i = 0; i < nIP * nIP; i++) {
		std::cout << "P" << i + 1;
		for (int j = 0; j < 4; j++) {
			std::cout << std::setw(16) << dN_dE[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;
	for (int i = 0; i < 4; i++) {
		std::cout << std::setw(16) << "dNdn" << i + 1;
	}
	std::cout << std::endl;
	for (int i = 0; i < nIP * nIP; i++) {
		std::cout << "P" << i + 1;
		for (int j = 0; j < 4; j++) {
			std::cout << std::setw(16) << dN_dn[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << "N_shape function vals for Hbc\n";
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < nIP; j++) {
			for (int x = 0; x < 4; x++) {
				std::cout << std::setw(16) << N_shape[i][j][x];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << "N_shape function vals for C\n";
	for (int i = 0; i < nIP*nIP; i++) {
		for (int j = 0; j < 4; j++) {
				std::cout << std::setw(16) << N_ofIP[i][j];
			}
		std::cout << std::endl;
	}
}

void Element4_2D::printH(double** Hmatrix, int el) {
	std::cout << "H matrix for element E" << el + 1 << std::endl;
	for (int x = 0; x < 4; x++) {
		for (int z = 0; z < 4; z++) {
			std::cout << std::setw(12) << Hmatrix[x][z];
		}
		std::cout << std::endl;
	}
}

void Element4_2D::printHbc(double** Hbcmatrix, int el) {
	std::cout << "Hbc matrix for element E" << el + 1 << std::endl;
	for (int x = 0; x < 4; x++) {
		for (int z = 0; z < 4; z++) {
			std::cout << std::setw(12) << Hbcmatrix[x][z];
		}
		std::cout << std::endl;
	}
}

void Element4_2D::printP(double* Pvector, int el) {
	std::cout << "P vector for element E" << el + 1 << std::endl;
	for (int z = 0; z < 4; z++) {
		std::cout << std::setw(12) << Pvector[z];
	}
	std::cout << std::endl;
}

void Element4_2D::printC(double** Cmatrix, int el) {
	std::cout << "C matrix for element E" << el + 1 << std::endl;
	for (int x = 0; x < 4; x++) {
		for (int z = 0; z < 4; z++) {
			std::cout << std::setw(12) << Cmatrix[x][z];
		}
		std::cout << std::endl;
	}
}