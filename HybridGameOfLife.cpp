//================Assignment6-GameOfLife-Serial.cpp================
//Steven Smith
//CSC 6740
//
//Complilation:
//mpic++ -fopenmp -o HybridProg HybridGameOfLife.cpp
//
//Execute:
//mpirun -pernode ./HybridProg <input filename> <Number of MPI Processes> <Number of Threads Per Process> <number of generations> <output filename>
//
//Example:
//mpirun -pernode ./HybridProg 10000by10000_0.txt 4 2 600 outputFile.txt
//=======================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <omp.h>
using namespace std;

//Create output file stream handler for use globally.
ofstream outFile;

//Create halo for the passed grid with rowSize and columnSize including halo rows and columns.
void createHalo(int** grid, int rowSize, int columnSize)
{
	int i = 0;
	//First copy last row of main grid to the top row of the halo.
	for(i = 1; i< columnSize-1; i++)
	{
		grid[0][i] = grid[rowSize-2][i];
	}
	
	//Second copy first row of main grid to bottom row in the halo.
	for(i=1; i< columnSize-1; i++)
	{
		grid[rowSize-1][i] = grid[1][i];
	}
	
	//Third Copy last column of main grid to left column in the halo.
	for(i = 1; i < rowSize - 1; i++)
	{
		grid[i][0] = grid[i][columnSize-2];
	}
	
	//Fourth copy the first column of the main grid to the right column in the halo.
	for(i = 1; i< rowSize-1; i++)
	{
		grid[i][columnSize-1] = grid[i][1];
	}
	
	
	//Finally copy the corners of the halo.
	
	//Top left cell.
	grid[0][0] = grid[rowSize-2][columnSize-2];
	
	//Top right cell.
	grid[0][columnSize-1] = grid[rowSize-2][1];
	
	//Bottom left cell.
	grid[rowSize-1][0] = grid[1][columnSize-2];
	
	//Bottom right cell.
	grid[rowSize-1][columnSize-1] = grid[1][1];
}

//Output grid to outFile starting from the specified start row and column to the specified last row and column.
void outputGrid(int** grid, int startRow, int startColumn, int lastRow, int lastColumn)
{
	int i =0;
	int j = 0;
	
	for(i = startRow; i <= lastRow; i++)
	{
		for(j = startColumn; j <= lastColumn; j++)
		{
			outFile << grid[i][j];
			
			if(j < lastColumn)
			{
				outFile << " ";
			}
		}
		outFile << "\n";
	}
}

/*Calculate the value for the current cell in the next generation, depending on the state of its neighbors (adjacent cells).
Rules:
1. Any live cell with fewer than two live neighbors dies, as if caused by under-population. 
2. Any live cell with two or three live neighbors’ lives on to the next generation. 
3. Any live cell with more than three live neighbors dies, as if by overcrowding. 
4. Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction
*/
int calcNextGenValue(int** grid, int row, int col)
{
	//Living neighbor count.
	int liveCount = 0;
	
	//First, calculate the number of live neighbors.
	
	//Upper left neighbor.
	if(grid[row-1][col-1] == 1)
	{
		liveCount++;
	}
	
	//Left neighbor.
	if(grid[row][col-1] == 1)
	{
		liveCount++;
	}
	
	//Lower left neighbor.
	if(grid[row+1][col-1] == 1)
	{
		liveCount++;
	}
	
	//Bottom neighbor.
	if(grid[row+1][col] == 1)
	{
		liveCount++;
	}
	
	//Top neighbor.
	if(grid[row-1][col] == 1)
	{
		liveCount++;
	}
	
	//Upper Right Neighbor.
	if(grid[row-1][col+1] == 1)
	{
		liveCount++;
	}
	
	//Right neighbor.
	if(grid[row][col+1] == 1)
	{
		liveCount++;
	}
	
	//Lower right neighbor.
	if(grid[row+1][col+1] == 1)
	{
		liveCount++;
	}
	
	//Next, determine if the current cell is live or dead in the next generation.
	//Currently Live cell
	if(grid[row][col] == 1)
	{
		//Dead if less than 2 live neighbors nearby.
		if(liveCount < 2)
		{
			return 0;
		}
		//Lives on if 2 or three neighbors are nearby.
		else if(liveCount < 4)
		{
			return 1;
		}
		//Dead if the number of live neighbors is 4 or more.
		else
		{
			return 0;
		}
	}
	else
	{
		//A dead cell lives in the next generation if the number of live neighbors is exactly 3.
		if(liveCount == 3)
		{
			return 1;
		}
		
		//Otherwise, the cell will stay dead.
		return 0;
	}
	
}

int main (int args, char* argv[])
{
	//Check for correct number of arguments.  Should call with at least 3 arguments for input filename, number of generations, and output filename.
	if(args < 3)
	{
		fprintf(stderr, "Usage: mpirun -pernode ./HybridProg <input filename> <Number of MPI Processes> <Number of Threads Per Process> <number of generations> <output filename>");
		return -1;
	}
	
	//Read in input file.
	char* inputFile = argv[1];
	
	//Read in number of MPI Processes.
	int numMPIProcesses = atoi(argv[2]);
	
	//Read in the number of threads per process.
	int numThreadsPerProcess = atoi(argv[3]);
	
	//Read in number of generations.
	int numGenerations = atoi(argv[4]);
	
	//Read in output file.
	char* outputFile = argv[5];
	
	//Error flag to indicate if any parameters are invalid.
	int error = 0;
	
	//Row and Column sizes for the main 2D grid.
	int startRowSize, startColumnSize;
	
	if(strlen(inputFile) < 1)
	{
		fprintf(stderr, "Input file name is required.\n");
		error = 1;
	}
	
	if(numMPIProcesses < 1)
	{
		fprintf(stderr, "Number of MPI Processes must be 1 or greater.\n");
		error = 1;
	}
	
	if(numThreadsPerProcess < 1)
	{
		fprintf(stderr, "Number of threads per process must be 1 or greater.\n");
		error = 1;
	}
	
	if(numGenerations < 1)
	{
		fprintf(stderr, "Number of generations must be 1 or greater.\n");
		error = 1;
	}
	
	if(strlen(outputFile) < 1)
	{
		fprintf(stderr, "Output file name is required.\n");
		error = 1;
	}
	
	if(error > 0)
	{
		return -1;
	}
	
	//Initialize MPI.
	MPI_Init(NULL, NULL);
	
	//Get start time of program.
	double runTimeStart = MPI_Wtime();

	//Set total MPI Process Count in MPI_COMM_WORLD and current process rank in MPI_COMM_WORLD
	int totalProcessCount, processRank;
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcessCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	
	//printf("Process %d Input File Name:  %s\nNum Generations:  %i\nOutput File Name %s\n", processRank, inputFile, numGenerations, outputFile);
	
	//Open input file for reading.
	ifstream inFile(inputFile);
	
	//Read in main grid row size and column size.
	inFile >> startRowSize >> startColumnSize;
	
	//Close input file here for all processes besides process 0.
	if(processRank != 0)
	{
		//Close input file.
		inFile.close();
	}
	
	//Set full row size and full column size by adding 2 to each for including the halo rows and columns.
	int fullRowSize = startRowSize+2;
	int fullColumnSize = startColumnSize+2;
	
	int currentRow = 1;
	int currentColumn = 1;
	int gridValue;
	int i = 0;
	int j = 0;
	int k = 0;
	
	//Set row chunk size for each MPI process.
	int chunkSize = startRowSize/totalProcessCount;
	
	//Set chunk remainder to distribute to the last process if it exists.
	int chunkRemainder = startRowSize % totalProcessCount;
	
	//First row to send and buffer size to send for each process.
	int sendFirstRow;
	int sendRowSize;
	int sendSize;
	
	//Last row for each process to handle processing the grid for only their portion.  Especially needed for process 0 as it has the full grid but only needs to process its 
	//portion.
	int lastRow;
	
	//printf("Full Row Size:  %d, Full Column Size:  %d\n", fullRowSize, fullColumnSize);
	
	if(totalProcessCount > startRowSize)
	{
		fprintf(stderr, "Error:  The number of rows in the grid must be greater than the number of MPI processes.\n");
		return -1;
	}
	
	//Allocate contiguous blocks of memory for the grid.
	int **fullGrid;
	int **fullGridNextGeneration;
	
	//First, process 0 needs to read the input file to create the initial grid, create the initial halo, and send chunks of the full grid (including their halos) out to each 
	//other process.
	if(processRank == 0)
	{
		//Allocate memory for the full grid and next generation grid.
		fullGrid = new int*[fullRowSize];
		fullGridNextGeneration = new int*[fullRowSize];
		
		for(i = 0; i<fullRowSize; i++)
		{
			fullGrid[i] = new int[fullColumnSize];
			fullGridNextGeneration[i] = new int[fullColumnSize];
		}
		
		//printf("Reading input grid.\n");
		while(inFile >> gridValue)
		{
			//Increment row when reaching the end of a column in the file.
			if(currentColumn > startColumnSize)
			{
				currentRow++;
				currentColumn = 1;
			}
			
			//Exit loop and stop reading values from the file if the file values go beyond the designated size of the grid at the beginning of the file.  This is to prevent 
			//errors caused by assigning values past the size of the grid in memory if the input file happens to have too many values.
			if((currentRow > startRowSize) || (currentColumn > startColumnSize))
			{
				continue;
			}
			
			fullGrid[currentRow][currentColumn] = (int)gridValue;
			
			currentColumn++;
		}
		
		//Close input file.
		inFile.close();
		
		//printf("Finished reading file.\n");
		fflush(stdout);
		
		//Create initial halo.
		createHalo(fullGrid, fullRowSize, fullColumnSize);
		
		//Distribute initial row chunks of the grid out to the other processes (1- (totalMPIProcesses-1)).  0 does not need to distribute its own as it already has the full 
		//grid.
		for(i = 1; i< totalProcessCount; i++)
		{
			
			//Set sendFirstRow to i*chunkSize-1.  The -1 is to also get the top halo row for the process.
			sendFirstRow = i*chunkSize-1;
			
			//Add +2 to the chunk size to get the halo rows.
			sendRowSize = chunkSize+2;
			
			//Give the last process the remaining rows if the total rows are not fully divisible by the number of processes.
			if((i == totalProcessCount-1) && (chunkRemainder > 0))
			{
				sendRowSize = sendRowSize+chunkRemainder;
			}
			
			
			//Calculate the size of the send buffer as the amount of rows to send multiplied by the number of columns to send.
			sendSize = sendRowSize*fullColumnSize;
			
			printf("Process %d send first row %d, send Row Size %d, send buffer Size %d.\n", i, sendFirstRow, sendRowSize, sendSize);
			
			fflush(stdout);
			
			//printf("Sending portion of grid to process %d, first value %d, second row first value %d, last row second to last value %d, last row last at location %d,%d value %d.\n", 
			//       i, fullGrid[sendFirstRow][0], fullGrid[sendFirstRow+1][0], fullGrid[sendFirstRow+sendRowSize][fullColumnSize-2], 
			//	   (sendFirstRow+sendRowSize), (fullColumnSize-1), fullGrid[sendFirstRow+sendRowSize][fullColumnSize-1]);
			
			//fflush(stdout);
			
			//Send the grid row chunks to each other process including their halos.
			MPI_Send(&fullGrid[sendFirstRow][0], sendSize, MPI_INT, i, 0, MPI_COMM_WORLD);
			printf("Sent portion of grid to process %d from process 0.\n", i);
			fflush(stdout);
		}
		
		//Set sendRowSize to chunkSize + 2 for halo rows as this will be the number of rows that process 0 will need to process.
		sendRowSize = chunkSize+2;
	}
	//Initialize fullGrid and fullGridNextGeneration to the size for each processes' partial grid, then receive the process's portion of the grid rows from process 0.
	else
	{
		//Add +2 to the chunk size to get the halo rows.
		sendRowSize = chunkSize+2;
		
		//Give the last process the remaining rows if the total rows are not fully divisible by the number of processes.
		if((processRank == totalProcessCount-1) && (chunkRemainder > 0))
		{
			sendRowSize = sendRowSize+chunkRemainder;
		}
		
		//Calculate the size of the send buffer as the amount of rows to send multiplied by the number of columns to send.
		sendSize = sendRowSize*fullColumnSize;
		
		int *localGrid = new int[sendSize];
		int localIndex = 0;
		
		printf("Process %d receive Row Size %d, receive buffer Size %d.\n", processRank, sendRowSize, sendSize);
		
		//Allocate memory for the process's grid and next generation grid.
		fullGrid = new int*[sendRowSize];
		fullGridNextGeneration = new int*[sendRowSize];
		
		for(i = 0; i<sendRowSize; i++)
		{
			fullGrid[i] = new int[fullColumnSize];
			fullGridNextGeneration[i] = new int[fullColumnSize];
		}
		
		//printf("Waiting to receive portion of grid from process 0 for process %d.\n", processRank);
		
		fflush(stdout);
		
		//Receive in 1D array to prevent memory issues.
		MPI_Recv(localGrid, sendSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		//Assign 1D array to 2D grid.
		for(i = 0; i < sendRowSize; i++)
		{
			for(j= 0; j < fullColumnSize; j++)
			{
				fullGrid[i][j] = localGrid[localIndex];
				localIndex++;
			}
		}
		
		//printf("Received portion of grid for process %d, first value %d, second row first value %d, last row second to last value %d, last row last value %d.\n", 
		//processRank, fullGrid[0][0], fullGrid[1][0], fullGrid[sendRowSize-1][fullColumnSize-2], fullGrid[sendRowSize-1][fullColumnSize-1]);
		
		//fflush(stdout);
	}
	
	//Last row (including halo) for processing should be same as the row size as the grid portion will start from index 0 for all processes.
	lastRow = sendRowSize;
	
	printf("Finished sending and receiving initial grid portions for process %d.\n", processRank);
	fflush(stdout);
	
	//Wait for all processes to finish receiving their portions.
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(processRank == 0)
	{
		printf("All processes have received their portions of the initial grid.\n");
		fflush(stdout);
	}
	
	printf("Process %d closed input file.\n", processRank);
	fflush(stdout);
	
	//MPI_Finalize();
	//return 0;
	
	//Open output file for writing.
	//outFile.open(outputFile);
	
	//Create the halo for the new generation's grid to prepare for the next generation.
	//COMMENT OUT
	//createHalo(fullGrid, fullRowSize, fullColumnSize);
	
	//COMMENT OUT
	//outFile << endl << "Generation 1" << endl;
	//outputGrid(fullGrid, 0, 0, fullRowSize-1, fullColumnSize-1);
	
	//Set values for chunk size, chunk remainder, start row, and last row for openmp thread per process portion.
	int threadChunkSize = sendRowSize/numThreadsPerProcess;
	int threadChunkRemainder = sendRowSize % numThreadsPerProcess;
	int threadStartRow = 0;
	int threadLastRow = 0;

	//Create the halo for the first generation's grid to prepare for the next generation.
	//createHalo(fullGrid, sendRowSize, fullColumnSize);
	//Calculate the generations.
	for(i =1; i< numGenerations; i++)
	{
		#pragma omp parallel firstprivate(threadStartRow, threadLastRow) num_threads(numThreadsPerProcess)
		{
			//Set current thread rank in openmp.
			int ompThreadRank = omp_get_thread_num();
			
			threadStartRow = ompThreadRank*threadChunkSize;
			threadLastRow = threadStartRow+threadChunkSize;
			
			//Assign the extra rows to the last openmp thread.
			if(ompThreadRank == numThreadsPerProcess-1)
			{
				threadLastRow = threadLastRow+threadChunkRemainder;
			}
			
			for(j= threadStartRow+1; j <= (threadLastRow-2); j++)
			{
				for(k = 1; k <= (fullColumnSize-2); k++)
				{
					fullGridNextGeneration[j][k] = calcNextGenValue(fullGrid, j, k);
				}
			}
		}
		
		//Copy new grid to current grid.
		for(j = 1; j <= sendRowSize-2; j++)
		{
			for(k = 1; k<= fullColumnSize-2; k++)
			{
				fullGrid[j][k] = fullGridNextGeneration[j][k];
			}
		}
		
		//Generate halo for the new generation.
		createHalo(fullGrid, sendRowSize, fullColumnSize);
		
		//If odd process, Send top main row to lower ranks for odd processes.
		if(processRank%2 == 1)
		{
			//printf("Odd process %d Sending top main row to previous process %d.\n", processRank, (processRank-1));
			//fflush(stdout);
			
			MPI_Send(&fullGrid[1][0], fullColumnSize, MPI_INT, processRank-1, 0, MPI_COMM_WORLD);
			
			//printf("Odd process %d Sent top main row to previous process %d.\n", processRank, (processRank-1));
			//fflush(stdout);
		}
		
		//If even process, receive top main row from odd process and assign it as the bottom halo row.
		if((processRank % 2 == 0) || (processRank == 0))
		{
			
			//If process +1 exists, receive its top main row and assign it as the bottom halo row for the current process.
			if((processRank+1) < totalProcessCount)
			{
				//printf("Even Process %d receiving top main row from process %d and assigning it as the bottom halo row.\n", processRank, (processRank+1));
				//fflush(stdout);
			
				MPI_Recv(&fullGrid[sendRowSize-1][0], fullColumnSize, MPI_INT, processRank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				//printf("Even Process %d received top main row from process %d and assigned it as the bottom halo row.\n", processRank, (processRank+1));
				//fflush(stdout);
			}
		}
		
		//If even process and not the last process, send bottom main row to the next process.
		if((processRank % 2) == 0 || (processRank == 0) && (processRank != totalProcessCount-1))
		{
			//If the next process exists, send the bottom main row to it.
			if((processRank+1) < totalProcessCount)
			{
				//printf("Even process %d sending bottom main row to next process %d.\n", processRank, (processRank+1));
				//fflush(stdout);
				
				MPI_Send(&fullGrid[sendRowSize-2][0], fullColumnSize, MPI_INT, processRank+1, 0, MPI_COMM_WORLD);
				
				//printf("Even process %d sent bottom main row to next process %d.\n", processRank, (processRank+1));
				//fflush(stdout);
			}
		}
		
		//Receive bottom row from the previous process and assign as the top halo row for odd processes.
		if(processRank %2 == 1)
		{
			//printf("Odd process %d receiving bottom main row from process %d and assigning it as the top halo row.\n", processRank, (processRank-1));
			//fflush(stdout);
			
			MPI_Recv(&fullGrid[0][0], fullColumnSize, MPI_INT, processRank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//printf("Odd process %d received bottom main row from process %d and assigned it as the top halo row.\n", processRank, (processRank-1));
			//fflush(stdout);
		}
		
		//If the last process, send bottom main row to process 0.
		if(processRank == totalProcessCount-1)
		{
			//printf("Last process %d sending bottom main row to process 0.\n", processRank);
			//fflush(stdout);
			
			MPI_Send(&fullGrid[sendRowSize-2][0], fullColumnSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
			
			//printf("Last process %d sent bottom main row to process 0.\n", processRank);
			//fflush(stdout);
		}
		
		//If process 0, receive the bottom main row from the last process and assign it as the top halo row.
		if(processRank == 0)
		{
			//printf("Process 0 receiving bottom main row from last process %d and assigning it as the top halo row.\n", (totalProcessCount-1));
			//fflush(stdout);
			
			MPI_Recv(&fullGrid[0][0], fullColumnSize, MPI_INT, totalProcessCount-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//printf("Process 0 received bottom main row from last process %d and assigned it as the top halo row.\n", (totalProcessCount-1));
			//fflush(stdout);
		}
		
		//If odd process and not the last process, Send bottom row from odd process to next process.
		if((processRank %2 == 1) && (processRank != totalProcessCount-1))
		{
			if((processRank+1) < totalProcessCount)
			{
				//printf("Odd process %d sending bottom main row to next process %d.\n", processRank, (processRank+1));
				//fflush(stdout);
				
				MPI_Send(&fullGrid[sendRowSize-2][0], fullColumnSize, MPI_INT, processRank+1, 0, MPI_COMM_WORLD);
				
				//printf("Odd process %d sent bottom main row to next process %d.\n", processRank, (processRank+1));
				//fflush(stdout);
			}
		}
		
		//If even process and not process 0, Receive bottom main row from previous process and assign as the top halo row.
		if((processRank %2 == 0) && (processRank != 0))
		{
			//printf("Even process %d receiving bottom main row from process %d and assigning it to the top halo row.\n", processRank, (processRank-1));
			//fflush(stdout);
			
			MPI_Recv(&fullGrid[0][0], fullColumnSize, MPI_INT, processRank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//printf("Even process %d received bottom main row from process %d and assigned it to the top halo row.\n", processRank, (processRank-1));
			//fflush(stdout);
		}
		
		//If even process and not process 0, Send top main row to previous process.
		if((processRank %2 == 0) && (processRank != 0))
		{
			//printf("Even process %d sending top main row to previous process %d.\n", processRank, (processRank-1));
			//fflush(stdout);
			
			MPI_Send(&fullGrid[1][0], fullColumnSize, MPI_INT, processRank-1, 0, MPI_COMM_WORLD);
			
			//printf("Even process %d sent top main row to previous process %d.\n", processRank, (processRank-1));
			//fflush(stdout);
		}
		
		//If odd process and not the last process, receive top main row from next process and set as the bottom halo row.
		if(processRank % 2 == 1)
		{
			if((processRank+1) < totalProcessCount)
			{
				//printf("Odd process %d receiving top main row from next process %d and assigning it to the bottom halo row.\n", processRank, (processRank+1));
				//fflush(stdout);
				
				MPI_Recv(&fullGrid[sendRowSize-1][0], fullColumnSize, MPI_INT, processRank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				//printf("Odd process %d received top main row from next process %d and assigned it to the bottom halo row.\n", processRank, (processRank+1));
				//fflush(stdout);
			}
		}
		
		//If process 0, send top main row to last process.
		if(processRank == 0)
		{
			//printf("Process 0 sending top main row to last process %d.\n", (totalProcessCount-1));
			//fflush(stdout);
			
			MPI_Send(&fullGrid[1][0], fullColumnSize, MPI_INT, totalProcessCount-1, 0, MPI_COMM_WORLD);
			
			//printf("Process 0 sent top main row to last process %d.\n", (totalProcessCount-1));
			//fflush(stdout);
		}
		
		//If the last process, receive top main row from process 0 and assign as the top halo row.
		if(processRank == (totalProcessCount-1))
		{
			//printf("Last process %d receiving top main row from process 0 and assigning it to the top halo row.\n", (totalProcessCount-1));
			//fflush(stdout);
			
			MPI_Recv(&fullGrid[sendRowSize-1][0], fullColumnSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//printf("Last process %d receiving top main row from process 0 and assigning it to the top halo row.\n", (totalProcessCount-1));
			//fflush(stdout);
		}
		
		//printf("Generation %d finished.\n", i+1);
	}
	
	//Reassemble the grid in process 0.
	if(processRank == 0)
	{
		for(int i =1; i< totalProcessCount; i++)
		{
			sendFirstRow = i*chunkSize-1;
			sendRowSize = chunkSize+2;

			
			if(i == (totalProcessCount-1) && (chunkRemainder > 0))
			{
				sendRowSize = sendRowSize+chunkRemainder;
			}
			
			lastRow = sendFirstRow+sendRowSize;
			
			sendSize = sendRowSize*fullColumnSize;
			
			printf("Process 0 receiving portion of final grid from process %d.  Send first Row %d.  Send Row Size %d.  \n", i, sendFirstRow, sendRowSize);
			fflush(stdout);
			
			MPI_Recv(&fullGrid[sendFirstRow-1][0], sendSize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			printf("Process 0 received portion of final grid from process %d.\n", i);
			fflush(stdout);
		}
		
		printf("Outputting final grid.\n");
		
		//Open output file for writing.
		outFile.open(outputFile);
		
		//Output main portion of grid to file, excluding the halo rows and columns.
		outputGrid(fullGrid, 1, 1, fullRowSize-2, fullColumnSize-2);

		outFile.close();
	}
	else
	{
		printf("Process %d sending portion of final grid to process 0.\n", processRank);
		fflush(stdout);
		
		//Each process sends their portion of the grid to process 0.
		MPI_Send(&fullGrid[0][0], sendSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
		
		printf("Process %d sent portion of final grid to process 0.\n", processRank);
		fflush(stdout);
	}
	
	MPI_Finalize();
	
	//Get end time for program.
	double runTimeEnd = MPI_Wtime();
	
	//Output Run Time only in process 0.
	if(processRank == 0)
	{
		printf("\nRun Time Start:  %lf\n", runTimeStart);
		
		printf("Run Time End:  %lf\n", runTimeEnd);
		
		double totalRunTime = runTimeEnd-runTimeStart;
		
		printf("\nTotal Run Time:  %lf seconds.\n", totalRunTime);
	}
}