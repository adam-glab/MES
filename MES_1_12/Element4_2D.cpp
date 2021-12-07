#include "Element4_2D.h"

void Element4_2D::printGauss() const{

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
