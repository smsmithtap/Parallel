//================Assignment6-GameOfLife-Serial.cpp================
//Steven Smith
//CSC 6740
//
//Complilation:
//g++ -std=c++11 -o SerialProg GameOfLife-Serial.cpp
//
//Execute:
//./SerialProg <inputFilename> <Number of Generations> <outputFilename>
//
//Example:
//./SerialProg inputFile.txt 4 outputFile.txt 
//=======================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;
using namespace std::chrono;

ofstream outFile;

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
		fprintf(stderr, "Usage: ./SerialProg <input filename> <number of generations> <output filename>");
		return -1;
	}
	
	//Read in input file.
	char* inputFile = argv[1];
	
	//Read in number of generations.
	int numGenerations = atoi(argv[2]);
	
	//Read in output file.
	char* outputFile = argv[3];
	
	//Error flag to indicate if any parameters are invalid.
	int error = 0;
	
	//Row and Column sizes for the main 2D grid.
	int startRowSize, startColumnSize;
	
	if(strlen(inputFile) < 1)
	{
		fprintf(stderr, "Input file name is required.\n");
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
	
	//Get start time.
	auto startTime = high_resolution_clock::now();
	
	printf("Input File Name:  %s\nNum Generations:  %i\nOutput File Name %s\n", inputFile, numGenerations, outputFile);
	
	//Open input file for reading.
	ifstream inFile(inputFile);
	
	//Read in main grid row size and column size.
	inFile >> startRowSize >> startColumnSize;
	
	//Set full row size and full column size by adding 2 to each for including the halo rows and columns.
	int fullRowSize = startRowSize+2;
	int fullColumnSize = startColumnSize+2;
	
	int currentRow = 1;
	int currentColumn = 1;
	int gridValue;
	int i = 0;
	int j = 0;
	int k = 0;
	
	printf("Full Row Size:  %d, Full Column Size:  %d\n", fullRowSize, fullColumnSize);
	
	//Allocate contiguous blocks of memory for the grid.
	int **fullGrid;
	int **fullGridNextGeneration;
	
	fullGrid = new int*[fullRowSize];
	fullGridNextGeneration = new int*[fullRowSize];
	
	for(int i = 0; i<fullRowSize; i++)
	{
		fullGrid[i] = new int[fullColumnSize];
		fullGridNextGeneration[i] = new int[fullColumnSize];
	}
	
	printf("Reading input grid.\n");
	while(inFile >> gridValue)
	{
		//printf("Read in character %d.\n", gridValue);
		//fflush(stdout);
		//Increment row when reaching the end of a column in the file.
		if(currentColumn > startColumnSize)
		{
			currentRow++;
			currentColumn = 1;
		}
		
		//Exit loop and stop reading input if the file values go beyond the designated size of the grid at the beginning of the file to prevent errors in memory.
		if((currentRow > startRowSize) || (currentColumn > startColumnSize))
		{
			continue;
		}

		//printf("Assigning cell value %d for grid at row %d column %d.\n", gridValue, currentRow, currentColumn);
		
		//fflush(stdout);
		
		fullGrid[currentRow][currentColumn] = (int)gridValue;
		
		//printf("File char for [%d][%d]: %d\n", currentRow, currentColumn, gridValue);
		//fflush(stdout);
		
		currentColumn++;
	}
	
	printf("Finished reading file.\n");
	fflush(stdout);
	
	inFile.close();
	
	//Calculate the generations.
	for(i =1; i< numGenerations; i++)
	{
		//Create the halo for the current generation's grid.
		createHalo(fullGrid, fullRowSize, fullColumnSize);
		for(j= 1; j <= (fullRowSize-2); j++)
		{
			for(k = 1; k <= (fullColumnSize-2); k++)
			{
				fullGridNextGeneration[j][k] = calcNextGenValue(fullGrid, j, k);
				
			}
		}
		
		//Copy new grid to current grid.
		for(j = 1; j <= fullRowSize-2; j++)
		{
			for(k = 1; k<= fullColumnSize-2; k++)
			{
				fullGrid[j][k] = fullGridNextGeneration[j][k];
			}
		}
	}
	
	//Open output file for writing.
	outFile.open(outputFile);
	
	//Output grid to file.
	outputGrid(fullGrid, 1, 1, fullRowSize-2, fullColumnSize-2);

	outFile.close();
	
	//Get end time.
	auto endTime = high_resolution_clock::now();
	
	//Output the execution time in microseconds.
	auto duration = duration_cast<microseconds>(endTime-startTime).count();
	
	cout << "Total Runtime:  " << duration << " microseconds." << endl;
	
	return 0;
}