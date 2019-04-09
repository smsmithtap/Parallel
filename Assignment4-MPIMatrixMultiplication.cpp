//================Assignment4-MPIMatrixMulitplication.cpp================
//Steven Smith
//CSC 6740
//
//Complilation:
//mpic++ -std=c++11 -o P Assignment4-MPIMatrixMulitplication.cpp
//
//Execute:
//mpiexec -np <NumberOfProcessors> P <MatrixSize>
//=======================================================================

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
using namespace std;

int main (int args, char* argv[])
{
	double runTimeStart = MPI_Wtime();

	//Check for correct number of arguments.  Should call with at least 1 argument for matrix size.
	if(args < 2)
	{
		fprintf(stderr, "Usage: mpirun -np <NumberofThreads> ./program <MatrixSize>\n");
		return -1;
	}
	
	//Read in matrix size.
	int matrixSize = atoi(argv[1]);
	
	//Initialize MPI.
	MPI_Init(NULL, NULL);

	int totalThreadCount, threadRank, i, j, k;
	MPI_Comm_size(MPI_COMM_WORLD, &totalThreadCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &threadRank);
	
	//Instantiate the random number engine.
	std::default_random_engine re;
	std::uniform_real_distribution<double> randNum(1, 10);
	
	//Allocate contiguous blocks of memory for Matrix 2 for sending and receiving with MPI.
	double **matrix2 = (double **)malloc(matrixSize*sizeof(double));
	double *matrix2Buffer = (double *)malloc(matrixSize*matrixSize*sizeof(double));

	for(i=0; i<matrixSize; i++)
	{
		matrix2[i] = &matrix2Buffer[i*matrixSize];
	}
	
	
	//Total size of matrix 2 for sending and receiving.
	int matrix2TotalSize = (matrixSize*matrixSize);
	
	//Calculate chunk size.
	int chunk = matrixSize/totalThreadCount;
	
	//Calculate last value for current split matrix.
	int last = chunk-1;
	int row = 0;
	int col = 0;
	
	int cellValue = 0;
	
	int matrixSplitSize = chunk;
	
	//Assign any leftover rows to the last thread.
	if(threadRank == (totalThreadCount-1))
	{
		int chunkRemainder = matrixSize % totalThreadCount;
		last = last+chunkRemainder;
		matrixSplitSize = chunk+chunkRemainder;
	}
	
	if(threadRank == 0)
	{
		//Allocate memory for Matrix 1 and the final Matrix.
		double **matrix1 = (double **)malloc(matrixSize*sizeof(double));
		double *matrix1Buffer = (double *)malloc(matrixSize*matrixSize*sizeof(double));
		
		double **finalMatrix = (double **)malloc(matrixSize*sizeof(double));
		double *finalMatrixBuffer = (double *)malloc(matrixSize*matrixSize*sizeof(double));
		
		//Allocate memory for Matrix 1.
		for(i=0; i<matrixSize; i++)
		{
			matrix1[i] = &matrix1Buffer[i*matrixSize];
			finalMatrix[i] = &finalMatrixBuffer[i*matrixSize];
		}
		
		//Assign random values to matrix 1 and 2.
		for (i =0; i<matrixSize; i++)
		{
			for (j=0; j<matrixSize; j++)
			{
				matrix1[i][j] = randNum(re);
				matrix2[i][j] = randNum(re);
			}
		}
		
		//First and last rows to send.
		int first;
		int sendFirst;
		int sendLast;
		int rowSize = chunk;
		int totalSplitMatrixSize;
		int curRow = 0;
		int chunkRemainder = matrixSize % totalThreadCount;
		
		for(i=1; i<totalThreadCount; i++)
		{
			sendFirst = i * chunk;
			sendLast = sendFirst+chunk-1;
			rowSize = chunk;
			curRow = 0;
			
			printf("Sending to thread %d with first row %d and last row %d.\n", i, sendFirst, sendLast);
			
			if(i == (totalThreadCount-1))
			{
				sendLast = sendLast+chunkRemainder;
				rowSize = rowSize+chunkRemainder;
			}
			//Total number of elements in the matrix to send.
			totalSplitMatrixSize = (rowSize*matrixSize);
			
			//Instantiate array for sending to thread rank i.
			//Allocate memory for send matrix using buffer to keep it contiguous in memory for sending with MPI.
			double **sendMatrix = (double **)malloc(rowSize*sizeof(double));
			double *sendMatrixBuffer = (double *)malloc(rowSize*matrixSize*sizeof(double));
			
			for(k=0; k<rowSize; k++)
			{
				sendMatrix[k] = &sendMatrixBuffer[k*matrixSize];
			}
			
			//Split matrix1 by row into smaller matrices to send to the other threads.
			for(j = sendFirst; j<=sendLast; j++)
			{
				for(k = 0; k<matrixSize; k++)
				{
					sendMatrix[curRow][k] = matrix1[j][k];
				}
				curRow++;
			}
			
			MPI_Send(&sendMatrix[0][0], totalSplitMatrixSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&matrix2[0][0], matrix2TotalSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		
		printf("Thread %d first row 0 last row %d.\n", threadRank, last);
		
		//Perform matrix multiplication for thread 0's assigned rows.
		for (int row=0; row<=last; row++)
		{
			for(int col=0; col<matrixSize; col++)
			{
				cellValue = 0;
				for (int k=0; k<matrixSize;k++)
				{
					cellValue += matrix1[row][k]*matrix2[k][col];
				}
				finalMatrix[row][col] = cellValue;
			}
		}
		
		for(i=1; i<totalThreadCount; i++)
		{
			//First and last rows in the final matrix the current thread rank i calculated.
			first = i*chunk;
			last = first+chunk;
			
			matrixSplitSize = chunk;
			
			if(i == (totalThreadCount-1))
			{
				last = last+chunkRemainder;
				matrixSplitSize = matrixSplitSize+chunkRemainder;
			}
			
			printf("Thread %d first row %d last row %d.\n", i, first, last);
			
			//Allocate memory for Matrix 2 using buffer to keep it contiguous in memory for sending with MPI.
			double **finalMatrixSplit = (double **)malloc(matrixSplitSize*sizeof(double));
			double *finalMatrixSplitBuffer = (double *)malloc(matrixSplitSize*matrixSize*sizeof(double));
			
			for(k=0; k<matrixSplitSize; k++)
			{
				finalMatrixSplit[k] = &finalMatrixSplitBuffer[k*matrixSize];
			}
			
			curRow = 0;
			//Total number of elements in the matrix to receive.
			totalSplitMatrixSize = (matrixSplitSize*matrixSize);
			
			MPI_Recv(&finalMatrixSplit[0][0], totalSplitMatrixSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//Add received portions of final matrix into the full matrix.
			for(row = first; row<last; row++)
			{
				for(col = 0; col<matrixSize; col++)
				{
					finalMatrix[row][col] = finalMatrixSplit[curRow][col];
				}
				curRow++;
			}
		}
		
		printf("Thread %d has successfully received and combined the portions of the final matrix.\n", threadRank);
	}
	//Perform matrix multiplication for the rows assigned to threads other than 0, then send the results back to thread 0.
	else
	{
		//Allocate memory for split matrices.
		double **matrix1Split = (double **)malloc(matrixSplitSize*sizeof(double));
		double *matrix1SplitBuffer = (double *)malloc(matrixSplitSize*matrixSize*sizeof(double));
		
		double **finalMatrixSplit = (double **)malloc(matrixSplitSize*sizeof(double));
		double *finalMatrixSplitBuffer = (double *)malloc(matrixSplitSize*matrixSize*sizeof(double));
		
		for(i=0; i<matrixSplitSize; i++)
		{
			matrix1Split[i] = &matrix1SplitBuffer[i*matrixSize];
			finalMatrixSplit[i] = &finalMatrixSplitBuffer[i*matrixSize];
		}
		
		int totalSplitArraySize = (matrixSplitSize*matrixSize);
		
		MPI_Recv(&matrix1Split[0][0], totalSplitArraySize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&matrix2[0][0], matrix2TotalSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		//Perform matrix multiplication for this thread's assigned rows.
		for (int row=0; row<matrixSplitSize; row++)
		{
			for(int col=0; col<matrixSize; col++)
			{
				cellValue = 0;
				for (int k=0; k<matrixSize;k++)
				{
					cellValue += matrix1Split[row][k]*matrix2[k][col];
				}
				finalMatrixSplit[row][col] = cellValue;
			}
		}

		MPI_Send(&finalMatrixSplit[0][0], totalSplitArraySize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
	double runTimeEnd = MPI_Wtime();
	
	//Output Run Time only in thread 0.
	if(threadRank == 0)
	{
		printf("\nRun Time Start:  %lf\n", runTimeStart);
		
		printf("Run Time End:  %lf\n", runTimeEnd);
		
		double totalRunTime = runTimeEnd-runTimeStart;
		
		printf("\nTotal Run Time:  %lf\n", totalRunTime);
	}
}