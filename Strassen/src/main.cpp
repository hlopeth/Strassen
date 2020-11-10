#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include "Matrix.h"
#include <chrono>
#include <iostream>

void meatherNaiveMatrixMult()
{
	Matrix m1 = randomMatrix(100, 100, -1000, 1000);
	Matrix m2 = randomMatrix(100, 100, -1000, 1000);
	auto tStart = std::chrono::high_resolution_clock::now();
	Matrix m3(m1 * m2);
	auto tEnd = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count();
	std::cout << "duration of naive matrix multiplication = " << duration << " microseconds" << std::endl;
	return;
}

int main(int argc, char* argv[]) {
	MPI_Status status;
	int procNum;
	int procRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	printf("rank %i, num %i \n", procRank, procNum);
	meatherNaiveMatrixMult();
	/*
	Matrix m1(3, 5);
	Matrix m2(2, 3);
	for (uint32_t i = 0; i < 3 * 5; i++)
	{
		m1(i) = i + 1;
	}
	for (uint32_t i = 0; i < 2 * 3; i++)
	{
		m2(i) = i + 1;
	}
	printf(to_string(m1).c_str());
	printf(to_string(m2).c_str());
	printf(to_string(m1 * m2).c_str());*/

	MPI_Finalize();
}