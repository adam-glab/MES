#include <iostream>

double funX(double x) {
	return 5 * x * x + 3 * x + 6;
}

double funXY(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}

// z tablic
//x_k

double m2[2] = { -0.577350,0.577350 };
double m3[3] = { -0.774597,0,0.774597 };
//double m4[4] = { -0.861136,-0.339981,0.339981,0.861136 };
//double m5[5] = { -0.906180,-0.538469,0,0.538469,0.906180 };

//A_k

double A2[2] = { 1.0,1.0 };
double A3[3] = { 0.555555,0.888888,0.555555 };
//double A4[4] = { 0.347855,0.652145,0.652145,0.347855 };
//double A5[5] = { 0.236927,0.478629,0.568889,0.478629,0.236927 };

double Gauss1D(int n) //w argumencie - stopien
{
	double sum = 0.0;
	double alfa = 1;
	double beta = 0;

	// nowe granice
	/*if (a != -1 || b != 1)
	{
		alfa = (b - a) / 2;
		beta = (a + b) / 2;
	}*/
	if (n == 2) {
		for (int i = 0; i < n; i++)
		{
			sum += alfa * A2[i] * funX(alfa * m2[i] - beta);
		}
	}
	if (n == 3) {
		for (int i = 0; i < n; i++)
		{
			sum += alfa * A3[i] * funX(alfa * m3[i] - beta);
		}
	}
	return sum;
}

double Gauss2D(int n) {
	double sum = 0.0;
	//double alfa = 1;
	//double beta = 0;

	// nowe granice
	/*if (a != -1 || b != 1)
	{
		alfa = (b - a) / 2;
		beta = (a + b) / 2;
	}*/

	//23 26
	if (n == 2) {
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				sum += A2[i] * A2[j] * funXY(m2[i], m2[j]);
			}

		}
	}
	if (n == 3) {
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				sum += A3[i] * A3[j] * funXY(m3[i], m3[j]);
			}

		}
	}
	return sum;
}

int main() {

	double res = Gauss1D(3);
	std::cout << res<< "\n";
	res = Gauss2D(3);
	std::cout << res << "\n";
	return 0;
}

#include <iostream>




double funX(double x) {
	return 5 * x * x + 3 * x + 6;
}

double funXY(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}

// z tablic
//x_k

double m2[2] = { -0.577350,0.577350 };
double m3[3] = { -0.774597,0,0.774597 };
//double m4[4] = { -0.861136,-0.339981,0.339981,0.861136 };
//double m5[5] = { -0.906180,-0.538469,0,0.538469,0.906180 };

//A_k

double A2[2] = { 1.0,1.0 };
double A3[3] = { 0.555555,0.888888,0.555555 };
//double A4[4] = { 0.347855,0.652145,0.652145,0.347855 };
//double A5[5] = { 0.236927,0.478629,0.568889,0.478629,0.236927 };

double Gauss1D(int n) //w argumencie - stopien
{
	double sum = 0.0;
	double alfa = 1;
	double beta = 0;

	// nowe granice
	/*if (a != -1 || b != 1)
	{
		alfa = (b - a) / 2;
		beta = (a + b) / 2;
	}*/
	if (n == 2) {
		for (int i = 0; i < n; i++)
		{
			sum += alfa * A2[i] * funX(alfa * m2[i] - beta);
		}
	}
	if (n == 3) {
		for (int i = 0; i < n; i++)
		{
			sum += alfa * A3[i] * funX(alfa * m3[i] - beta);
		}
	}
	return sum;
}

double Gauss2D(int n) {
	double sum = 0.0;
	//double alfa = 1;
	//double beta = 0;

	// nowe granice
	/*if (a != -1 || b != 1)
	{
		alfa = (b - a) / 2;
		beta = (a + b) / 2;
	}*/

	//23 26
	if (n == 2) {
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				sum += A2[i] * A2[j] * funXY(m2[i], m2[j]);
			}

		}
	}
	if (n == 3) {
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				sum += A3[i] * A3[j] * funXY(m3[i], m3[j]);
			}

		}
	}
	return sum;
}

int main() {

	double res = Gauss1D(3);
	std::cout << res << "\n";
	res = Gauss2D(3);
	std::cout << res << "\n";
	return 0;
}
