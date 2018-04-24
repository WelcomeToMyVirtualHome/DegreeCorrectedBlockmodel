#include "stdafx.h"
#include "Graph.h"

#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC ((znew<<16)+wnew )
#define SHR3 (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define CONG (jcong=69069*jcong+1234567)
#define KISS ((MWC^CONG)+SHR3)
#define UNI (KISS*2.328306e-10)
typedef unsigned long UL;
static UL z = 362436069, w = 521288629, jsr = 123456789, jcong = 380116160;
static UL a = 224466889, b = 7584631, t[256];

using namespace std;

Graph::Graph()
{
	Log = "log.dat";
}
double Graph::drand48() 
{
	return rand() / ((double)RAND_MAX + 1.0);
}
double Graph::vectorSum(vector<double> &vect)
{
	double sum = 0;
	for (int i = 0; i < vect.size(); i++)
	{
		sum += vect[i];
	}
	return sum;
}
double Graph::vectorSum(vector<double> &vect, int N0, int N)
{
	double sum = 0;
	for (int i = N0; i < N; i++) {
		sum += vect[i];
	}
	return sum;
}
void Graph::setDegreeSequence()
{
	cout << "Getting node degrees..." << endl;
	nodeDegrees.resize(size);
	for (int i = 0; i < numberOfBlocks; i++) 
	{
		float x1 = pow((float)blockSize, (float)(blockExponents[i] + 1));
		float x0 = pow((float)blockXMins[i], (float)(blockExponents[i] + 1));
		for (int j = i * blockSize; j < (i + 1)*blockSize; j++)
		{
			nodeDegrees[j] = pow((x1 - x0) * UNI  + x0, (float)(1 / (blockExponents[i] + 1)));
		}
	}
	cout << "..done" << endl;
}
void Graph::setLagrangeMultipliers()
{
	cout << "Calculating multipliers..." << endl;
	internalTheta.resize(size);
	Err.resize(numberOfBlocks,0);
	externalTheta.resize(size);
	for (int i = 0; i < numberOfBlocks; i++)
	{
		for (int j = i * blockSize; j < (i + 1)*blockSize; j++)
		{
			Err[i] += nodeDegrees[j]/2;
		}
		for (int j = i * blockSize; j < (i + 1)*blockSize; j++)
		{
			internalTheta[j] = log(nodeDegrees[j] / sqrt(2*Err[i]));
			externalTheta[j] = log(nodeDegrees[j] * sqrt(Ers) / (2 * Err[i]));
		}
	}
	matrixTheta.resize(size, vector<double>(size, 0));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i / blockSize == j / blockSize)
			{
				matrixTheta[i][j] = internalTheta[i] + internalTheta[j];
				matrixTheta[j][i] = internalTheta[i] + internalTheta[j];
			}
			else
			{
				matrixTheta[i][j] = externalTheta[i] + externalTheta[j];
				matrixTheta[j][i] = externalTheta[i] + externalTheta[j];
			}
		}
	}
	cout << "...done" << endl;
}
void Graph::monteCarlo()
{
	ofstream myfile;
	myfile.open(Log, ios::trunc);
	graph.resize(size, vector<int>(size, 1));
	matrixESym.resize(numberOfBlocks, vector<double>(numberOfBlocks, 0));
	internalNodeDegreesSym.resize(size);
	H.resize(snap);
	float theta_ij = 0, nodeDegree = 0, sum = 0, ETA = 0;
	int i, j, dH, k, it;
	float snapTemp = (float)NT / snap;
	dH = size * size / 2;
	k = 0;
	clock_t start = 0;
	for (it = 0; it < NT; it++)
	{
		i = (int)(UNI * size);
		j = (int)(UNI * size);
		theta_ij = matrixTheta[i][j];
		if (graph[i][j] == 1)
		{
			graph[i][j] = 0;
			graph[j][i] = 0;
			dH -= 1;
		}
		else
		{
			if (UNI < exp(theta_ij))
			{
				graph[i][j] = 1;
				graph[j][i] = 1;
				dH += 1;
			}
		}
		if (fmod(it, snapTemp) == 0)
		{
			H[k] = dH;
			printf("E = %1.0f, %1.0f%c done", H[k], (float)it / NT * 100, '%');
			myfile << "E= " << H[k] << " " << (float)it / NT * 100 << "% done" << " ETA: " << ETA << "s" << endl;
			k++;
			cout << " ETA: " << ETA << "s" << endl;
			ETA = (clock() - start) / (double)CLOCKS_PER_SEC * (NT - it) / snapTemp;
			start = clock();

		}
	}
	for (int i = 0; i < numberOfBlocks; i++)
	{
		for (int j = 0; j < numberOfBlocks; j++)
		{
			sum = 0;
			for (int p = i * blockSize; p < (i + 1) * blockSize; p++)
			{
				for (int k = j * blockSize; k < (j + 1) * blockSize; k++)
				{
					sum += graph[k][p];
					if (i == j)
					{
						nodeDegree += graph[k][p];
					}
				}
				if (i == j)
				{
					internalNodeDegreesSym[p] = nodeDegree;
					nodeDegree = 0;
				}
			}
			if (i == j)
			{
				matrixESym[i][j] = sum/2;
			}
			else
			{
				matrixESym[i][j] = sum;
			}
		}
	}
	myfile.close();
	myfile << "100% done" << endl;
	cout << "100% done" << endl;
}
void Graph::monteCarlo(float NTfraction)
{
	ofstream myfile;
	myfile.open(Log, ios::trunc);
	graph.resize(size, vector<int>(size, 1));
	matrixESym.resize(numberOfBlocks, vector<double>(numberOfBlocks, 0));
	internalNodeDegreesSym.resize(size, 0);
	H.resize(snap);
	float theta_ij = 0, nodeDegree = 0, sum = 0, ETA = 0;
	int i, j, dH, k, it;
	float snapTemp = (float)NT / snap;
	dH = size * size / 2;
	k = 0;
	clock_t start = 0;
	for (it = 0; it < NT; it++)
	{
		i = (int)(UNI * size);
		j = (int)(UNI * size);
		theta_ij = matrixTheta[i][j];
		if (graph[i][j] == 1)
		{
			graph[i][j] = 0;
			graph[j][i] = 0;
			dH -= 1;
		}
		else
		{
			if (UNI < exp(theta_ij))
			{
				graph[i][j] = 1;
				graph[j][i] = 1;
				dH += 1;
			}
		}
		if (fmod(it, snapTemp) == 0)
		{
			H[k] = dH;
			printf("E = %1.0f, %1.0f%c done", H[k], (float)it / NT * 100, '%');
			myfile << "E= " << H[k] << " " << (float)it / NT * 100 << "% done" << " ETA: " << ETA << "s" << endl;
			k++;
			cout << " ETA: " << ETA << "s" << endl;
			ETA = (clock() - start) / (double)CLOCKS_PER_SEC * (NT - it) / snapTemp;
			start = clock();
			if (it > (float)NT*NTfraction )
			{
				for (int i = 0; i < numberOfBlocks; i++)
				{
					for (int j = 0; j < numberOfBlocks; j++)
					{
						sum = 0;
						for (int p = i * blockSize; p < (i + 1) * blockSize; p++)
						{
							for (int k = j * blockSize; k < (j + 1) * blockSize; k++)
							{
								sum += graph[k][p];
								if (i == j)
								{
									nodeDegree += graph[k][p];
								}
							}
							if (i == j)	
							{
								internalNodeDegreesSym[p] += nodeDegree / ((1 - NTfraction)*snap);
								nodeDegree = 0;
							}
						}
						if (i == j)
						{
							matrixESym[i][j] += sum / 2 / ((1 - NTfraction)*snap);
						}
						else
						{
							matrixESym[i][j] = sum  / ((1 - NTfraction)*snap);
						}
					}
				}
				
			}		
		}
	}
	myfile.close();
	myfile << "100% done" << endl;
	cout << "100% done" << endl;
}

Graph::~Graph()
{
}
void Graph::writeMatrix(string file)
{
	ofstream myfile;
	myfile.open(file, ios::trunc);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			myfile << graph[i][j] << " ";
		}
		myfile << endl;
	}
	myfile.close();
}
void Graph::writeMatrixTheta(string file)
{
	ofstream myfile;
	myfile.open(file, ios::trunc);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			myfile << matrixTheta[i][j] << " ";
		}
		myfile << endl;
	}
	myfile.close();
}
void Graph::writeCytoscapeMatrix(string file)
{
	ofstream myfile;
	myfile.open(file, ios::trunc);
	myfile << "source \t target" << endl;
	for (int j = 0; j < size; j++)
	{
		for (int k = j + 1; k < size; k++)
		{
			if (graph[j][k] == 1)
			{
				myfile << j << "\t" << k << endl;
			}
		}
	}
	myfile.close();
}
void Graph::writeEdgeSumMatrixes(string file, string fileSym)
{
	ofstream myfile;
	myfile.open(fileSym, ios::trunc);
	for (int i = 0; i < numberOfBlocks; i++)
	{
		for (int j = 0; j < numberOfBlocks; j++)
		{
			myfile << matrixESym[i][j] << '\t';
		}
		myfile << endl;
	}
	myfile.close();
	myfile.open(file, ios::trunc);
	for (int i = 0; i < numberOfBlocks; i++)
	{
		for (int j = 0; j < numberOfBlocks; j++)
		{
			if (i == j)
			{
				myfile << Err[i] << '\t';
			}
			else
			{
				myfile << Ers << '\t';
			}
		}
		myfile << endl;
	}
	myfile.close();
}
void Graph::writeVector(vector<double> &vect, string file)
{
	ofstream myfile;
	myfile.open(file, ios::trunc);
	for (int i = 0; i < vect.size(); i++)
		myfile << vect[i] << endl;;
	myfile.close();
}
void Graph::writeParams(string file)
{
	ofstream myfile;
	myfile.open(file, ios::trunc);
	myfile << numberOfBlocks << '\n' << blockSize << endl;;
	for (int i = 0; i < numberOfBlocks; i++) {
		myfile << blockXMins[i] << endl;;
	}
	for (int i = 0; i < numberOfBlocks; i++) {
		myfile << blockExponents[i] << endl;
	}
	myfile.close();
}
