#include <math.h>
#include <iostream>
#include "matrix.h"
#include "vldb.h"

Vldb::Vldb(const std::vector<bool>& labelList, int p, double delta, int maxIterations)
{
	this->labelList = labelList;
	this->p = p;
	this->delta = delta;
	this->maxIterations = maxIterations;
}

matrix Vldb::initOLSParams(const_matrix xMatrix, const_matrix yMatrix)
{
	const_matrix transX = ublas::trans(xMatrix);

	const_matrix middleMatrix = ublas::prec_prod(transX, xMatrix);

	const_matrix inversed = inverse(middleMatrix);
	//std::cout << inversed << std::endl << std::endl;

	const_matrix prec = ublas::prec_prod(inversed, transX);
	return ublas::prec_prod(prec, yMatrix);
}

matrix Vldb::initICParams(const_matrix A, const_matrix B)
{
	return ublas::prec_prod(inverse(A), B);
}

matrix Vldb::combine(const_matrix phi, const_matrix xMatrix)
{
	return ublas::prec_prod(xMatrix, phi);
}

int Vldb::repair(const_matrix yhatMatrix, const_matrix yMatrix)
{
	int rowNum = yhatMatrix.size1();
	matrix residualMatrix = yhatMatrix - yMatrix;

	double aMin = MINVAL;
	int targetIndex = -1;

	for (int i = 0; i < rowNum; ++i) {
		if (labelList.at(i + p)) {
			continue;
		}
		if (fabs(residualMatrix(i, 0)) < delta) {
			continue;
		}

		double yhat = yhatMatrix(i, 0);
		double yhatabs = fabs(yhat);

		if (yhatabs < aMin) { // no need to > 0
			aMin = yhatabs;
			targetIndex = i;
		}
	}

	return targetIndex;
}


