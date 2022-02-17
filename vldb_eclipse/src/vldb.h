#pragma once

#include <vector>
#include <float.h>
#include "matrix.h"

#define MINVAL DBL_MAX;
#define MAXVAL -DBL_MAX;

class Vldb
{
public:
	Vldb(const std::vector<bool>& labelList, int p, double delta, int maxIterations);

	matrix initOLSParams(const_matrix xMatrix, const_matrix yMatrix);

	matrix initICParams(const_matrix A, const_matrix B);

	matrix combine(const_matrix phi, const_matrix xMatrix);

	int repair(const_matrix yhatMatrix, const_matrix yMatrix);

	std::vector<bool> labelList; // whether the point is labeled
	int p; // AR(p) model
	double delta; // converge
	int maxIterations; // max iteration number
};