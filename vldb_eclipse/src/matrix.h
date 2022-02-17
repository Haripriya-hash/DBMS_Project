#pragma once

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

typedef ublas::matrix<double> matrix;
typedef const matrix& const_matrix;
typedef const matrix* const_matrix_ptr;

matrix inverse(const_matrix m);
