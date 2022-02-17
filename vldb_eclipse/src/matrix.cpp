#include "matrix.h"

matrix inverse(const_matrix m)
{
	typedef ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix A(m);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0) {
		return m;
	}

	matrix result(ublas::identity_matrix<double>(A.size1()));
	//std::cout << result << std::endl << std::endl;

	// backsubstitute to get the inverse
	lu_substitute(A, pm, result);
	
	return result;
}
