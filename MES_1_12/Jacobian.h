#pragma once

struct jacobian {
	double j_matrix[2][2] = { 0.0 };
	void printJacobian();
};