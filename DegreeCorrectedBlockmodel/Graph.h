#pragma once
#ifndef _GRAPH_H_ 
#define _GRAPH_H_
#include <vector>
#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <random>
#include <chrono>
using namespace std;

class Graph
{
public:
	size_t NT;
	vector<double> blockXMins;
	vector<double> blockExponents;
	vector<double> nodeDegrees;
	vector<double> internalTheta;
	vector<double> externalTheta;
	vector<double> internalNodeDegreesSym;
	vector<vector<double>> matrixESym;
	vector<vector<double>> matrixE;
	vector<vector<double> > matrixTheta;
	vector<vector<int> > graph;
	vector<double> Err;
	vector<double> H;
	string Log;
	double Ers;
	int interblockConnections;
	size_t numberOfBlocks;
	size_t blockSize;
	size_t size;
	size_t snap;
	Graph();
	void setDegreeSequence();
	void setLagrangeMultipliers();
	void monteCarlo(); 
	void monteCarlo(float NTfraction);
	double drand48();
	double vectorSum(vector<double> &vect);
	double vectorSum(vector<double> &vect, int N0, int N);
	void writeVector(vector<double> &vect, string file);
	void writeParams(string file);
	void writeMatrix(string file);
	void writeCytoscapeMatrix(string file);
	void writeEdgeSumMatrixes(string file, string fileSym);
	void writeMatrixTheta(string file);
	~Graph();
};

#endif