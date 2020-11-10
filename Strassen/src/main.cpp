#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include "math.h"
#include <random>
#include "Matrix.h"

int main(int argc, char* argv[]) {
	MPI_Status status;
	int procNum;
	int procRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	printf("rank %i, num %i \n", procRank, procNum);

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
	printf(to_string(m1 * m2).c_str());

	MPI_Finalize();
}