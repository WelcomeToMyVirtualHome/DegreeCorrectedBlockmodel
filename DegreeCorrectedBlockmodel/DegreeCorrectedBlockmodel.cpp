// DegreeCorrectedBlockmodel.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "Graph.h"
#include <ctime>

int main()
{
	srand(time(NULL));
	clock_t start;
	double duration;
	start = clock();
	Graph graph = Graph();
	graph.numberOfBlocks = 8;
	graph.blockSize = 256;
	graph.Ers = 3; 
	graph.size = graph.numberOfBlocks*graph.blockSize;
	graph.NT = 100000000; //180000000, 8x512;	
	graph.snap = 100;
	graph.blockXMins.resize(graph.numberOfBlocks, 1);
	graph.blockExponents.resize(graph.numberOfBlocks, -3);
	graph.setDegreeSequence();
	graph.setLagrangeMultipliers();
	graph.monteCarlo();
	graph.writeParams("outputParams.dat");
	graph.writeVector(graph.nodeDegrees, "outputInternalND.dat");
	graph.writeVector(graph.internalNodeDegreesSym, "outputInternalNDSym.dat");
	graph.writeVector(graph.H, "outputH.dat");
	graph.writeVector(graph.internalTheta, "outputInternalTheta.dat");
	graph.writeVector(graph.externalTheta, "outputExternalTheta.dat");
	graph.writeEdgeSumMatrixes("outputEdgeSumMatrix.dat", "outputEdgeSumMatrixSym.dat");
	graph.writeCytoscapeMatrix("outputMatrix.sif");
	//graph.writeMatrix("outputMatrix.dat");
	//graph.writeMatrixTheta("outputMatrixTheta.dat");
	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "Elapsed time: " << duration << endl;
	ofstream myfile;
	myfile.open(graph.Log, ofstream::app);
	myfile << "Elapsed time: " << duration << endl;
	myfile.close();
	system("PAUSE");
	return 0;
}
	
