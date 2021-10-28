#include <iostream>
#include "grid.h"
#include "Element4_2D.h"
#include "node.h"
#include "jacobian.h"

#define LAB 4

double funX(double x) {
	return 5 * x * x + 3 * x + 6;
}

double funXY(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}


double Gauss1D(gauss val, int n) { //w argumencie - stopien
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += val.wP[i] * funX(val.xC[i]);
	}
	return sum;
}

double Gauss2D(gauss val, int n) {
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			sum += val.wP[i] * val.wP[j] * funXY(val.xC[i], val.xC[j]);
		}
	}
	return sum;
}

void calc_jacobian(int nE, int nIP, jacobian* J, jacobian* J_inv, Element4_2D E, double x[], double y[]) {
	double sumX = 0.;
	double sumY = 0.;
	for (int i = 0; i < nE; i++) {
		for (int j = 0; j < nIP; j++) {
			sumX += E.dNdE[i][j] * x[j];
			sumY += E.dNdn[i][j] * y[j];
		}
		J->j_matrix[0][0] = sumX;
		J->j_matrix[0][1] = 0.;
		J->j_matrix[1][0] = 0.;
		J->j_matrix[1][1] = sumY;
		std::cout << "::::::::::RESULT START::::::::::" << std::endl;
		std::cout << "::::::::::1/J::::::::::" << std::endl;
		std::cout << 1/J->j_matrix[0][0] << '\t' << J->j_matrix[0][1] << '\n' << J->j_matrix[1][0] << '\t' << 1/J->j_matrix[1][1] << std::endl<<std::endl;
		std::cout << "::::::::::J::::::::::" << std::endl;
		std::cout << J->j_matrix[0][0] << '\t' << J->j_matrix[0][1] << '\n' << J->j_matrix[1][0] << '\t' << J->j_matrix[1][1] << std::endl << std::endl;
		std::cout << "::::::::::RESULT END::::::::::" << std::endl;
		sumX = 0.;
		sumY = 0.;
	}
}

int main() {

#if LAB == 1
	double hei, len;
	int n_h, n_l;

	std::cout << "Enter height, length (greater than 0.0):\n";
	std::cin >> hei >> len;
	std::cout << "Enter nodes on height, on length (greater than 0):\n";
	std::cin >> n_h >> n_l;

	// Create grid
	grid G(hei, len, n_h, n_l);

	std::cout << "\nNode values:\n";
	for (int i = 0; i < G.nN; i++) {
		std::cout << "Node " << i + 1
			<< ":\t" << G.nodes[i].x
			<< "\t" << G.nodes[i].y
			<< "\n";
	}

	std::cout << "\nElement ID values (id1, id2, id3, id4):\n";
	for (int i = 0; i < G.nE; i++) {
		std::cout << "Element " << i + 1
			<< ":\t" << G.elements[i].ID[0]
			<< "\t" << G.elements[i].ID[1]
			<< "\t" << G.elements[i].ID[2]
			<< "\t" << G.elements[i].ID[3]
			<< "\n";
	}
	return 0;
#endif

#if LAB == 2
	int n = 2;
	gauss* values = new gauss(n);
	std::cout << "1D, 2 points:\t" << Gauss1D(*values, n) << "\n";
	values = new gauss(n);
	std::cout << "2D, 2 points:\t" << Gauss2D(*values, n) << "\n";
	n = 3;
	values = new gauss(n);
	std::cout << "1D, 3 points:\t" << Gauss1D(*values, n) << "\n";
	values = new gauss(n);
	std::cout << "2D, 3 points:\t" << Gauss2D(*values, n) << "\n";
	return 0;
#endif

//#if LAB == 3
	//Element4_2D test(4);
	//return 0;
//#endif

#if LAB == 4

	jacobian J;
	jacobian J_inv;
	Element4_2D E(4);

	int nE = 1;
	int nIP = 4;

	double x[4] = { 0., 0.025, 0.025, 0. };
	double y[4] = { 0., 0., 0.025, 0.025 };

	/*
	for (int i = 0; i < nE; i++) {
		for (int j = 0; j < nIP; j++) {
			calc_jacobian(i, j, &J, &J_inv, E, x, y);
		}
	}
	*/
	calc_jacobian(4, 4, &J, &J_inv, E, x, y);
	return 0;
#endif
}