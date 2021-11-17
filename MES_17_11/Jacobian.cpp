#include "Jacobian.h"

void jacobian::printJacobian(){
	std::cout << "=========================" << std::endl;
	std::cout << this->j_matrix[0][0] << '\t' << this->j_matrix[0][1] << '\n' << this->j_matrix[1][0] << '\t' << this->j_matrix[1][1] << std::endl << std::endl;
	std::cout << "=========================" << std::endl;
}
