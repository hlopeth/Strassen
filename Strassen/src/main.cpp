#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include "Matrix.h"
#include <chrono>
#include <iostream>
#include <cassert>

Matrix strassenMult(Matrix m1, Matrix m2);
Matrix _strassenMult(const Matrix& m1, const Matrix& m2);

Matrix meatherNaiveMatrixMult(const Matrix& m1, const Matrix& m2)
{
	auto tStart = std::chrono::high_resolution_clock::now();
	Matrix m3(m1 * m2);
	auto tEnd = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count();
	std::cout << "duration of naive matrix multiplication = " << duration << " microseconds" << std::endl;
	return m3;
}

Matrix meatherStrassenMatrixMult(const Matrix& m1, const Matrix& m2)
{
	auto tStart = std::chrono::high_resolution_clock::now();
	Matrix m3 = strassenMult(m1, m2);
	auto tEnd = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count();
	std::cout << "duration of strassen matrix multiplication = " << duration << " microseconds" << std::endl;
	return m3;
}

bool IsPowerOfTwo(uint32_t x)
{
	return (x != 0) && ((x & (x - 1)) == 0);
}

uint32_t nextPowerOfTwo(uint32_t v)
{
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v++;
	return v;
}

int main(int argc, char* argv[]) {
	MPI_Status status;
	int procNum;
	int procRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	printf("rank %i, num %i \n", procRank, procNum);
	//meatherNaiveMatrixMult();

	uint32_t n = 256;
	const Matrix m1 = randomMatrix(n,n,-100,100);
	const Matrix m2 = randomMatrix(n,n,-100,100);
	Matrix r2 = meatherNaiveMatrixMult(m1, m2);
	Matrix r1 = meatherStrassenMatrixMult(m1, m2);
	/*printf(to_string(m1).c_str());
	printf(to_string(m2).c_str());
	printf(to_string(r1).c_str());
	printf(to_string(r2).c_str());*/

	printf("same %i", r1 == r2);

	MPI_Finalize();
}


Matrix _strassenMult(const Matrix& m1, const Matrix& m2)
{
	if (m1.colls() <= 32)
	{
		return m1 * m2;
	}

	Matrix A11, A12, A21, A22;
	Matrix::split(m1, A11, A12, A21, A22);

	Matrix B11, B12, B21, B22;
	Matrix::split(m2, B11, B12, B21, B22);

	Matrix P1 = _strassenMult((A11 + A22), (B11 + B22));
	Matrix P2 = _strassenMult((A21 + A22), B11);
	Matrix P3 = _strassenMult(A11, (B12 - B22));
	Matrix P4 = _strassenMult(A22, (B21 - B11));
	Matrix P5 = _strassenMult((A11 + A12), B22);
	Matrix P6 = _strassenMult((A21 - A11), (B11 + B12));
	Matrix P7 = _strassenMult((A12 - A22), (B21 + B22));

	Matrix C11 = P1 + P4 - P5 + P7;
	Matrix C12 = P3 + P5;
	Matrix C21 = P2 + P4;
	Matrix C22 = P1 - P2 + P3 + P6;

	Matrix C = Matrix::merge(C11, C12, C21, C22);
	return C;
}

Matrix strassenMult(Matrix m1, Matrix m2)
{
	assert(m1.colls() == m2.rows());
	assert(m1.rows() == m1.colls());
	assert(m2.rows() == m2.colls());

	if (!IsPowerOfTwo(m1.colls()))
	{
		uint32_t newSize = nextPowerOfTwo(m1.colls());
		m1.recize(newSize, newSize);
	}
	if (!IsPowerOfTwo(m2.colls()))
	{
		uint32_t newSize = nextPowerOfTwo(m2.colls());
		m2.recize(newSize, newSize);
	}
	return _strassenMult(m1,m2);
}