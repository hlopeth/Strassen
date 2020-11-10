#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include "Matrix.h"
#include <chrono>
#include <iostream>
#include <cassert>

#define STRASSEN_PARALEL 1

Matrix strassenMult(Matrix m1, Matrix m2);
Matrix _strassenMult(const Matrix& m1, const Matrix& m2);
Matrix mpi_strassenMult(const Matrix& m1, const Matrix& m2);
void master(int procRank, int procNum);
void slave(int procRank, int procNum);

static uint32_t _N = 1024;

static const int INIT_FLAG_READY = 1;
static const int INIT_FLAG_TERMINATE = 2;

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

void mpiSend(int dest, const Matrix& m)
{
	MPI_Send(m.data(), m.colls() * m.rows(), MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
}

Matrix mpiRecv(int source, uint32_t size)
{
	MPI_Status mpiStatus;
	Matrix m(size, size);
	MPI_Recv(m.data(), size * size, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &mpiStatus);
	return m;
}

int main(int argc, char* argv[]) {
	MPI_Status status;
	int procNum;
	int procRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	if (procNum != 7)
	{
		printf("only 7 pocesses supported");
		MPI_Finalize();
		return 0;		
	}
		
	switch (procRank)
	{
	case 0:
		master(procRank, procNum);
		break;
	default:
		slave(procRank, procNum);
		break;
	}

	MPI_Finalize();
	return 0;
}

void master(int procRank, int procNum)
{
	for (int i = 1; i < 7; i++)
	{
		int flag = (STRASSEN_PARALEL) ? INIT_FLAG_READY : INIT_FLAG_TERMINATE;
		MPI_Send( &flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	const Matrix m1 = randomMatrix(_N, _N, -100, 100);
	const Matrix m2 = randomMatrix(_N, _N, -100, 100);
	Matrix r2 = meatherNaiveMatrixMult(m1, m2);
	Matrix r1 = meatherStrassenMatrixMult(m1, m2);
	printf("same %i", r1 == r2);
}

void slave(int procRank, int procNum)
{
	MPI_Status mpiStatus;
	int initFlag;
	MPI_Recv(&initFlag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpiStatus);
	if (initFlag == INIT_FLAG_READY)
	{
		printf("proc %i ready", procRank);
		Matrix m1 = mpiRecv(0, _N / 2);
		Matrix m2 = mpiRecv(0, _N / 2);
		Matrix res = _strassenMult(m1, m2);
		mpiSend(0, res);
	}
	else 
	{
		printf("proc %i terminate", procRank);
	}
}


Matrix mpi_strassenMult(const Matrix& m1, const Matrix& m2)
{
	Matrix A11, A12, A21, A22;
	Matrix::split(m1, A11, A12, A21, A22);

	Matrix B11, B12, B21, B22;
	Matrix::split(m2, B11, B12, B21, B22);

#if STRASSEN_PARALEL
	//send P1
	//Matrix P1 = _strassenMult((A11 + A22), (B11 + B22));
	mpiSend(1, A11 + A22);
	mpiSend(1, B11 + B22);

	//P2
	//Matrix P2 = _strassenMult((A21 + A22), B11);
	mpiSend(2, A21 + A22);
	mpiSend(2, B11);

	//P3
	//Matrix P3 = _strassenMult(A11, (B12 - B22));
	mpiSend(3, A11);
	mpiSend(3, B12 - B22);

	//P4
	//Matrix P4 = _strassenMult(A22, (B21 - B11));
	mpiSend(4, A22);
	mpiSend(4, B21 - B11);

	//P5
	//Matrix P5 = _strassenMult((A11 + A12), B22);
	mpiSend(5, (A11 + A12));
	mpiSend(5, B22);

	//P6
	//Matrix P6 = _strassenMult((A21 - A11), (B11 + B12));
	mpiSend(6, A21 - A11);
	mpiSend(6, B11 + B12);

	Matrix P7 = _strassenMult((A12 - A22), (B21 + B22));

	uint32_t size = A11.colls();
	Matrix P1 = mpiRecv(1, size);
	Matrix P2 = mpiRecv(2, size);
	Matrix P3 = mpiRecv(3, size);
	Matrix P4 = mpiRecv(4, size);
	Matrix P5 = mpiRecv(5, size);
	Matrix P6 = mpiRecv(6, size);
#else
	Matrix P1 = _strassenMult((A11 + A22), (B11 + B22));
	Matrix P2 = _strassenMult((A21 + A22), B11);
	Matrix P3 = _strassenMult(A11, (B12 - B22));
	Matrix P4 = _strassenMult(A22, (B21 - B11));
	Matrix P5 = _strassenMult((A11 + A12), B22);
	Matrix P6 = _strassenMult((A21 - A11), (B11 + B12));
	Matrix P7 = _strassenMult((A12 - A22), (B21 + B22));
#endif

	Matrix C11 = P1 + P4 - P5 + P7;
	Matrix C12 = P3 + P5;
	Matrix C21 = P2 + P4;
	Matrix C22 = P1 - P2 + P3 + P6;

	Matrix C = Matrix::merge(C11, C12, C21, C22);
	return C;
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
	return mpi_strassenMult(m1, m2);
}

