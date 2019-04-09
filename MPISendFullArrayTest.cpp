//================Assignment4-MPIMatrixMulitplication.c================
//Steven Smith
//CSC 6740
//
//=====================================================================

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int args, char* argv[])
{
	//Check for correct number of arguments.  Should call with at least 3 arguments for number of threads, string, 
	//and integer.  Exit if too few arguments.
	if(args < 2)
	{
		fprintf(stderr, "Usage: mpirun -np <NumberofThreads> ./program <MatrixSize>\n");
		return -1;
	}
	
	//Read in string and integer arguments.
	int matrixSize = atoi(argv[1]);
	
	//Initialize MPI.
	MPI_Init(NULL, NULL);

	int totalThreadCount, threadRank, i, j, k;
	MPI_Comm_size(MPI_COMM_WORLD, &totalThreadCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &threadRank);
	
	
	//Validate thread count to be greater than 1.  Exit if not.
	if(totalThreadCount < 2)
	{
		fprintf(stderr, "Error:  Number of threads must be greater than 2.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	double **matrix2 = (double **)malloc(matrixSize*sizeof(double));
	double *matrix2Buffer = (double *)malloc(matrixSize*matrixSize*sizeof(double));
	
	//Allocate memory for Matrix 1.
	for(i=0; i<matrixSize; i++)
	{
		matrix2[i] = &matrix2Buffer[i*matrixSize];
	}
	
	
	//Total size of matrix 1 for sending and receiving.
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
		//Allocate memory for Receive Matrix.
		printf("Thread %d allocating memory for matrix 1.\n", threadRank);
		double **matrix1 = (double **)malloc(matrixSize*sizeof(double));
		double *matrix1Buffer = (double *)malloc(matrixSize*matrixSize*sizeof(double));
		
		//Allocate memory for Matrix 1.
		for(i=0; i<matrixSize; i++)
		{
			matrix1[i] = &matrix1Buffer[i*matrixSize];
		}
		
		double **finalMatrix = (double **)malloc(matrixSize*sizeof(double));
		double *finalMatrixBuffer = (double *)malloc(matrixSize*matrixSize*sizeof(double));
		
		
		for(i=0; i<matrixSize; i++)
		{
			finalMatrix[i] = &finalMatrixBuffer[i*matrixSize];
		}
		
		//First and last rows to send.
		int first;
		int sendFirst;
		int sendLast;
		int rowSize = chunk;
		int totalSplitMatrixSize;
		int curRow = 0;
		int chunkRemainder = matrixSize % totalThreadCount;
		
		printf("Setting values for matrix 1 in thread %d.\n", threadRank);
		
		//Set value for each row to the row num+1.
		for(i=0; i < matrixSize; i++)
		{
			for(j=0; j<matrixSize; j++)
			{
				matrix1[i][j] = i+1;
			}
		}
		
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
					//printf("Value adding to send matrix at row %d col %d, from row %d col %d: %lf.\n", curRow, k, j, k, matrix1[j][k]);
					sendMatrix[curRow][k] = matrix1[j][k];
				}
				curRow++;
			}
			
			printf("Thread %d now sending portion of matrix 1 with total size of %d to thread %d.\n", threadRank, totalSplitMatrixSize, i);
			MPI_Send(&sendMatrix[0][0], totalSplitMatrixSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		
		//Double values in thread 0's assigned rows.
		for(row=0; row<=last; row++)
		{
			for(col=0; col<matrixSize; col++)
			{
				finalMatrix[row][col] = 2*matrix1[row][col];
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
			
			printf("Thread %d final matrix - first row %d last row %d.\n", i, first, last);
			
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
			
			printf("Thread %d receiving final matrix from thread %d.\n", threadRank, i);
			MPI_Recv(&finalMatrixSplit[0][0], totalSplitMatrixSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			printf("Thread %d received final matrix from thread %d.\n", threadRank, i);
			
			//Add received portions of final matrix into the full matrix.
			for(row = first; row<last; row++)
			{
				for(col = 0; col<matrixSize; col++)
				{
					printf("Adding row %d col %d from split matrix row %d col %d.\n", row, col, curRow, col);
					fflush(stdout);
					finalMatrix[row][col] = finalMatrixSplit[curRow][col];
				}
				curRow++;
			}
		}
		
		printf("Thread %d has finished receiving and combining the portions of the final matrix with size %d.\n", threadRank, matrixSize);
		
		for(row = 0; row<matrixSize; row++)
		{
			for(col = 0; col<matrixSize; col++)
			{
				printf("Final matrix row %d col %d value:  %lf.\n", row, col, finalMatrix[row][col]);
			}
		}
	}
	else
	{
		//Allocate memory for full size matrices.
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
		
		printf("Waiting to receive portion of matrix 1 from thread 0 in thread %d.\n", threadRank);
		
		MPI_Recv(&matrix1Split[0][0], totalSplitArraySize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		printf("Thread %d received portion of Matrix 1 from thread 0.\n", threadRank);

		for(i = 0; i<matrixSplitSize; i++)
		{
			//printf("Matrix 1 Split Row %d in thread %d.\n", i, threadRank);
			for(j = 0; j<matrixSize; j++)
			{
				//printf("Split Matrix 1 thread %d row %d column %d Value:  %lf.\n", threadRank, i, j, matrix1Split[i][j]);
				matrix1Split[i][j] = 2*matrix1Split[i][j];
			}
		}
		
		printf("Thread %d sending matrix 1 back to thread 0.\n", threadRank);
		MPI_Send(&matrix1Split[0][0], totalSplitArraySize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
}