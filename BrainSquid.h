#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <complex>
using namespace std;

class Squid {

private:
int count;
double** x;
complex<double>*** y;
double*** periodogram;
double ** compressedpsd;
double* NumericallyValuedData;
char* Buffer;

public:
	Squid(void);
	double stringToDouble(string temp);
	void load(void);
	void printData(void);
	void printFFT(void); // NEW ADDED 3/13
	//void readIn(void);
	//void write(void);
	void write0(void);
	double dot(double A[], double B[]);
	double magnitude(double array[]);
	double similarity(double A[], double B[]);
	//void init(double*** x, complex<double>*** y);
	//void init(double*** x, complex<double>**** y);
	void init0(void); // NEW ADDED 3/13
	void readIn0(void); // NEW ADDED 3/13
	//void psd(complex<double>*** y, double*** periodogram);
	//void fft(double** x, complex<double>** y);
	//void fft(double** x, complex<double>*** y,int inc, int offset);
	void fft0(int inc, int offset); // NEW ADDED 3/13
	void psd0(void); // NEW ADDED 3/13
	void DoFFT(void); // NEW ADDED 3/13
	//void insertionSort(void);
	void doMedian(void); // NEW ADDED 3/13
	void insertionSort0(void); // NEW ADDED 3/13
	~Squid(void);
};

