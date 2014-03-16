#pragma once
#include "Squid2.h"
using namespace std;
//C:\Users\Frazer\Documents\Visual Studio 2010\Projects\fft\fft
#define INPUT_FILE "C:\\Users\\Frazer\\Documents\\Visual Studio 2010\\Projects\\fft\\fft\\sinewave0.csv"
#define OUTPUT_FILE "C:\\Users\\Frazer\\Documents\\Visual Studio 2010\\Projects\\fft\\fft\\SquidOut.csv"
#define OUTPUT_FFT "C:\\Users\\Frazer\\Documents\\Visual Studio 2010\\Projects\\fft\\fft\\FFTOut.csv"
#define PI 3.14159265359
#define SIZE 640
#define CHANNELS 14
#define SAMPLING_FREQUENCY_Hz 128 
#define FFT_SAMPLE_SIZE_N 256 
#define FFT_STEP 16 
#define MAX_DATA 1920 // 1920 = 128 * 15 (15 seconds of sampling)
#define complexi complex<double>(0,1)


Squid::Squid() {
	Squid::count = 0;
	Squid::NumericallyValuedData = new double[SIZE];
	Squid::Buffer = new char[SIZE];
	Squid::x = new double* [CHANNELS];
	Squid::y = new complex<double>** [CHANNELS]; 
	Squid::periodogram = new double** [CHANNELS];
	Squid::compressedpsd = new double* [CHANNELS];
}
/*
void Squid::psd(complex<double>*** y, double*** periodogram){
	double T = FFT_SAMPLE_SIZE_N/SAMPLING_FREQUENCY_Hz;
	double deltaT = 1/SAMPLING_FREQUENCY_Hz;
	for(int i = 0; i < CHANNELS; i++){
		for(int j = 0; j < (SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP; j++){
			for(int k = 0; k < FFT_SAMPLE_SIZE_N; k++){
				periodogram[i][j][k] = ((deltaT*deltaT)/T) * abs(y[i][j][k])*abs(y[i][j][k]);
			}
		}
	}
}
*/

//Where (2*pi*k)/N is omega
void Squid::psd0(){
	double T = (double)FFT_SAMPLE_SIZE_N/(double)SAMPLING_FREQUENCY_Hz;
	double deltaT = 1.0/(double)SAMPLING_FREQUENCY_Hz;
	for(int i = 0; i < CHANNELS; i++){
		for(int j = 0; j < (SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP; j++){
			for(int k = 0; k < FFT_SAMPLE_SIZE_N; k++){
				periodogram[i][j][k] = ((deltaT*deltaT)/T) * abs(y[i][j][k])*abs(y[i][j][k]);
			}
		}
	}
}
/*
void Squid::init(double*** x, complex<double>*** y){
	*y = new complex<double>* [CHANNELS];
	*x = new double* [CHANNELS];
	for(int j = 0; j < CHANNELS; j++){
		(*y)[j] = new complex<double>[SIZE];
		(*x)[j] = new double[SIZE];
	}
}

void Squid::init(double*** x, complex<double>**** y){
	*y = new complex<double>** [CHANNELS];
	*x = new double* [CHANNELS];
	for(int j = 0; j < CHANNELS; j++){
		(*y)[j] = new complex<double>* [FFT_SAMPLE_SIZE_N/FFT_STEP];
		(*x)[j] = new double[SIZE];
		for(int i = 0; i < (FFT_SAMPLE_SIZE_N/FFT_STEP); i++){
			(*y)[j][i] = new complex<double>[FFT_SAMPLE_SIZE_N];
		}
	}
}
*/
void Squid::init0(){
	for(int j = 0; j < CHANNELS; j++){
		y[j] = new complex<double>* [(SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP];
		periodogram[j] = new double* [(SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP];
		x[j] = new double[MAX_DATA];
		compressedpsd[j] = new double[FFT_SAMPLE_SIZE_N];
		for(int i = 0; i < ((SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP); i++){
			y[j][i] = new complex<double>[FFT_SAMPLE_SIZE_N];
			periodogram[j][i] = new double[FFT_SAMPLE_SIZE_N];
		}
	}
}
/*
void Squid::fft(double** x, complex<double>** y){
	complex<double> sum(0,0);
	double z = 0;
	for(int i = 0 ; i < CHANNELS; i++){
		for(int k = 0; k < SIZE; k++){
			sum = complex<double>(0,0);
			for(int n = 0; n < SIZE; n++){
				z = ((-1) * 2 * PI * n * k) / SIZE;
				sum += x[i][n] * exp(z);
			}
			y[i][k] = sum;
		}
	}
}
*/
/*
void Squid::fft(double** x, complex<double>*** y, int inc, int offset){
	complex<double> sum(0,0);
	double z = 0;
	for(int i = 0 ; i < CHANNELS; i++){
		for(int k = 0; k < FFT_SAMPLE_SIZE_N; k++){
			sum = complex<double>(0,0);
			for(int n = inc; n < (inc+FFT_SAMPLE_SIZE_N); n++){
				z = ((-1) * 2 * PI * n * k) / SIZE;
				sum += x[i][n] * exp(z);
			}
			y[i][(inc-offset)/FFT_STEP][k] = sum;
		}
	}
}
*/
void Squid::DoFFT(){
	int offset;
	int inc = 0;
	if(count < SIZE){
		cout << "Error: Data Sample too small" << endl;
		return;
	}
	else if(count == SIZE){
		offset = 0;
	}
	else{
		offset = (count/2) - (SIZE/2);
	}
	for(inc = offset; inc < offset+SIZE-FFT_SAMPLE_SIZE_N; inc+=FFT_STEP){
		fft0(inc, offset);
	}
}
void Squid::fft0(int inc, int offset){
	complex<double> sum(0,0);
	double z = 0;
	for(int i = 0 ; i < CHANNELS; i++){
		for(int k = 0; k < FFT_SAMPLE_SIZE_N; k++){
			sum = complex<double>(0,0);
			for(int n = inc; n < (inc+FFT_SAMPLE_SIZE_N); n++){
				z = ((-1) * 2 * PI * n * k) / FFT_SAMPLE_SIZE_N;
				sum += x[i][n] * exp(z*complexi);
			}
			y[i][(inc-offset)/FFT_STEP][k] = sum;
		}
	}
}

/*
void Squid::load(){
	for(int i = 0 ; i < CHANNELS; i++){
		for(int j = 0; j < SIZE; j++){
			x[i][j] = 1;
			y[i][j] = complex<double>(0,0);
		}
	}
}
*/
double Squid::similarity(double A[], double B[]){
	double dotAB = dot(A,B);
	double magA = magnitude(A);
	double magB = magnitude(B);
	return dotAB/(magA*magB);
}

double Squid::dot(double A[], double B[]){
	double sum = 0;
	for(int i=0; i<sizeof(A)/sizeof(double);i++){
		sum += A[i]*B[i];
	}
	return sum;
}
double Squid::magnitude(double array[]){
	double sum = 0;
	for(int i=0; i<sizeof(array)/sizeof(double);i++){
		sum += array[i]*array[i];
	}
	return sqrt(sum);
}

double Squid::stringToDouble(string temp) {
	return atof(temp.c_str());
}
/*
void Squid::insertionSort(void){
		double temp = 0;
		double* tempArray = new double[count];
		for(int i=0; i<count;i++){
			tempArray[i]=NumericallyValuedData[i];
		}
		NumericallyValuedData = NULL;
		NumericallyValuedData = tempArray;
		for (int i = 0; i < count; i++) {
			for (int j = 0; j < count; j++) {
				if (NumericallyValuedData[j] < NumericallyValuedData[i]) {
					temp = NumericallyValuedData[j];
					NumericallyValuedData[j] = NumericallyValuedData[i];
					NumericallyValuedData[i] = temp;
				}
			}
		}
}
*/
void Squid::doMedian(void){
	insertionSort0();
	for(int i = 0; i < CHANNELS; i++){
		for(int k = 0; k < FFT_SAMPLE_SIZE_N; k++){
			compressedpsd[i][k] = periodogram[i][((SIZE-FFT_SAMPLE_SIZE_N)/(2*FFT_STEP))][k];
		}
	}
}
void Squid::insertionSort0(void){
	double temp = 0;
	/*double* tempArray = new double[count];
	for(int i=0; i<count;i++){
		tempArray[i]=NumericallyValuedData[i];
	}
	NumericallyValuedData = NULL;
	NumericallyValuedData = tempArray;
	*/
	for (int i = 0; i < CHANNELS; i++) {
		for (int k = 0; k < FFT_SAMPLE_SIZE_N; k++) {
			for(int j = 0; j < (SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP; j++){
				for(int n = 0; n < (SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP; n++){
					if (periodogram[i][n][k] < periodogram[i][j][k]) {
						temp = periodogram[i][j][k];
						periodogram[i][j][k] = periodogram[i][n][k];
						periodogram[i][n][k] = temp;
					}
				}
			}
		}
	}
}
/*
void Squid::readIn(void) {
	ifstream brainWaveInputStream(INPUT_FILE,std::ifstream::in);
	string inputStringBuffer;
	
	while (getline(brainWaveInputStream, inputStringBuffer)) {
		stringstream inputStringStream(inputStringBuffer);
		string tempString;
		while (getline(inputStringStream, tempString, ',')) {
			NumericallyValuedData[count++] = stringToDouble(tempString);
		}
	}
	brainWaveInputStream.close();
}
*/
void Squid::readIn0(void) {
	int channelNumber = 0;
	ifstream brainWaveInputStream(INPUT_FILE,std::ifstream::in);
	string inputStringBuffer;
	getline(brainWaveInputStream, inputStringBuffer); //Gets first information line
	//can add stuff here to see when the data was collected, etc.
	while (getline(brainWaveInputStream, inputStringBuffer)) {
		stringstream inputStringStream(inputStringBuffer);
		string tempString;
		getline(inputStringStream, tempString, ','); // Emotiv Data, for sample number
		getline(inputStringStream, tempString, ','); // Emotiv Data for DC Offset
		while (getline(inputStringStream, tempString, ',') && channelNumber < CHANNELS) {
			x[channelNumber++][count] = stringToDouble(tempString);
		}
		channelNumber = 0;
		count++;
		if(count == MAX_DATA){
			break;
		}
	}
	brainWaveInputStream.close();
}
/*
void Squid::write(){
	ofstream brainWaveOutputStream(OUTPUT_FILE,std::ofstream::out);
	for(int i =0; i<sizeof(NumericallyValuedData)/sizeof(double);i++){
		brainWaveOutputStream<<NumericallyValuedData[i]<<',';
	}
	brainWaveOutputStream.close();
}
*/

void Squid::write0(){
	ofstream brainWaveOutputStream(OUTPUT_FILE,std::ofstream::out);
	for(int j = 0; j < FFT_SAMPLE_SIZE_N; j++){
		for(int i = 0; i < CHANNELS; i++){
			brainWaveOutputStream << compressedpsd[i][j] << ',';
		}
		brainWaveOutputStream << endl;
	}
	brainWaveOutputStream.close();
}

/*
void Squid::printData(void) {
	for (int i = 0; i < count; i++) {
		printf("%.10f \n", NumericallyValuedData[i]);
	}
}
*/
/*
	FUNCTION: printFFT
	PURPOSE: save fft information as csv file. In this case, it takes the 14 channel FFT,
		and saves both real and imaginary values delimited by commas.
		Each line represents a new sample number, each other pair of values in a line
		represent a separate channel.
void Squid::printFFT(void) {
	ofstream brainWaveOutputStream(OUTPUT_FFT,std::ofstream::out);
	for (int i = 0; i < CHANNELS; i++) {
		for(int k = 0; k < FFT_SAMPLE_SIZE_N; k++){
			for(int j = 0; j < (SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP; j++){
				brainWaveOutputStream << y[i][j][k].real() << "," << y[i][j][k].imag() << ",";
			}
			brainWaveOutputStream << endl;
		}
		brainWaveOutputStream << endl;
	}
	brainWaveOutputStream.close();
}
*/
Squid::~Squid() {
	for(int j = 0; j < CHANNELS; j++){
		for(int i = 0; i < ((SIZE-FFT_SAMPLE_SIZE_N)/FFT_STEP); i++){
			delete y[j][i];
			delete periodogram[j][i];
		}
		delete compressedpsd[j];
		delete x[j];
		delete periodogram[j];
		delete y[j];
	}
	delete [] x;
	delete [] y;
	delete [] periodogram;
	delete [] compressedpsd;
	delete [] NumericallyValuedData;
	delete [] Buffer;
}