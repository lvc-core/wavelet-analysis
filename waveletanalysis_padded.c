#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <complex.h>
#include <time.h>

/*
 * Parameters which should be used for the wavelet
 */
typedef struct Parameter_t{
    double freq_min;
    double freq_max;
    double delta_freq;

    double time_min;
    double time_max;
    double delta_time;

    double b, k, omega_0;
} Parameter_t;

/*
 * Regarding the data itself
 */
typedef struct Data_t{
    double dt_data;
    int numLines;
    int padding;
    int numLines_padded;

    double* dataArray;
    double complex* waveletArray;
    double complex* productArray;
} Data_t;

void allocate_memory_for_data(Data_t* data);
int getNumberOfLines(const char* filename, Data_t* data);
int readDataFromFile(const char* filename, Data_t* data, double *targetArray);
void createMorletWavelet(double complex* array, int numElements, double dt_morlet, double k, double omega0, double a, double b);
double complex calculateConvolutionValue(double* dataArray, double complex* waveletArray, double complex* productArray, int numLines, double h);
void analyzeWaveletTD(Parameter_t* parameters, Data_t* data);
void analyzeWaveletFD(Parameter_t* parameters, Data_t* data);
void createDummyData(double dt_data);
void freeData(Data_t* data);



int main(int argc, char* argv[]){
	Parameter_t parameters = {
	    0.1,	// freq_min
	    20.0,	// freq_max
	    0.01,	// delta_freq
	    0.0,	// time_min
	    20.0,	// time_max
	    0.01,	// delta_time
	    0.0,	// b
	    0.0,	// k
	    atof(argv[1])	// omega_0
	};

	Data_t data = {
	    0.01,	// dt_data
	    -1,		// numLines
	    -1,		// padding
	    -1,		// numLines_padded
	    NULL,	// *dataArray
	    NULL,	// *waveletArray
	    NULL	// *productArray
	};

	// Create artificial data (here a linear chirp function / a sin with constant frequency)
	createDummyData(data.dt_data);


	// Get length of file/data
	getNumberOfLines("data_input/hellothere.dat", &data);
	

	// Allocate memory for arrays in data-structure
	allocate_memory_for_data(&data);


	// Read data from file into data-array
	readDataFromFile("data_input/hellothere.dat", &data, data.dataArray);


	// Perform wavelet calculation depending on the mode set (frequency and/or time domain)
#ifdef TD
	analyzeWaveletTD(&parameters, &data);
#endif
#ifdef FD
	analyzeWaveletFD(&parameters, &data);
#endif
	freeData(&data);
}

void createDummyData(double dt_data){
	double t;
	FILE *fp = fopen("data_input/hellothere.dat", "w");
	for(t=0.0; t<20.0; t+=dt_data){
	    if(t<10.0)
		fprintf(fp, "%.2f\n", sin(2*M_PI*(0.5*t*t + 1.0*t)));
	    else
		//fprintf(fp, "%.2f\n", sin(2*M_PI*(0.5*t*t + 1.0*t)));
		fprintf(fp, "%.2f\n", sin(2*M_PI*(5.0*t)));

	}
	fclose(fp);
}

int getNumberOfLines(const char* filename, Data_t* data){
	FILE *fp = fopen(filename, "r");	
	if(fp==NULL){
		printf("Could not open file!\n");
		return -1;
	}

	int lines = 0;
	char buffer[1024];

	while(fgets(buffer, sizeof(buffer), fp) != NULL){
		lines++;
	}

	fclose(fp);

	data->numLines = lines;

	// add zero-padding
	int padding = ceil(lines / 4.0);
	printf("numLines: %d\n", lines);
	printf("padding: %d\n", padding);
	data->padding = padding;
	data->numLines_padded = lines + 2*padding;
	printf("numLines_padded: %d\n", data->numLines_padded);

	return 0;
}

int readDataFromFile(const char* filename, Data_t* data, double *targetArray){
	FILE *fp = fopen(filename, "r");

	char line[10];
	int i;
	char *endptr;

	for(i=0; i<data->padding; i++){
	    targetArray[i] = 0.0;
	}

	for(i=0; i<data->numLines; i++){
		if(!fgets(line, sizeof(line), fp)){
			fprintf(stderr, "Could not read line %d\n", i+1);
			fclose(fp);
			return EXIT_FAILURE;
		}	

		targetArray[i + data->padding] = strtod(line, &endptr);
		if(endptr == line){
			fprintf(stderr, "Unknown value in line %d\n", i+1);
			fclose(fp);
			return EXIT_FAILURE;
		}
	}
	fclose(fp);

	for(i=data->numLines + data->padding + 1; i<data->numLines_padded; i++){
	    targetArray[i] = 0.0;
	}

	return EXIT_SUCCESS;
}

void allocate_memory_for_data(Data_t* data){
	double *dataArray = malloc(data->numLines_padded * sizeof(double));
	if(dataArray == NULL){
		printf("Could not allocate memory for data array!\n");
	}


	double complex *waveletArray = malloc(data->numLines_padded * sizeof(double complex));
	if(waveletArray == NULL){
		printf("Could not allocate memory for wavelet array!\n");
	}


	double complex *productArray = malloc(data->numLines_padded * sizeof(double complex));
	if(productArray == NULL){
		printf("Could not allocate memory for product array!\n");
	}


	data->dataArray = dataArray;
	data->waveletArray = waveletArray;
	data->productArray = productArray;
}


void createMorletWavelet(double complex* array, int numElements, double dt_morlet, double k, double omega0, double a, double b){
	// psi = k_0 * cos(wt) * exp(-t²/2)
	// psi = k * exp(i*w_0*t) * exp(-t²/2)	
	
	int i;
	double tau;
	for(i=0; i<numElements; i++){
		tau = (i*dt_morlet - a)/b;
		if (fabs(tau) > 6.0){
		    array[i] = 0.0;
		}else{
		    array[i] = k * cexp(I * omega0 * tau) * exp(-pow(tau, 2) / 2.0);
		}
	}
}

/*
 *
 * ===========================================================================================================
 * 					TIME-DOMAIN
 * ===========================================================================================================
 */

double complex calculateConvolutionValue(double* dataArray, double complex* waveletArray, double complex* productArray, int numLines, double h){
	int i;
	double complex integralValue = 0.0;

	for(i=0; i<numLines; i++){
		productArray[i] = dataArray[i] * waveletArray[i]; 
	}


	// Simpson's integration formula + 3/8-rule
	if(numLines % 2 == 0){	
		integralValue += productArray[0];
		for(i=1; i < numLines-3; i+=2){
			integralValue += 4.0 * productArray[i] + 2.0 * productArray[i+1];
		}
		integralValue *= h/6.0;

		integralValue += h/8.0 * (productArray[numLines-4] + 3.0*productArray[numLines-3] + 3.0*productArray[numLines-2] + productArray[numLines-1]);
	
	}else{ // Simpson's integration formula
		integralValue += productArray[0] + 4.0*productArray[numLines - 2] + productArray[numLines - 1];
		for(i=1; i < numLines-3; i+=2){
			integralValue += 4.0 * productArray[i] + 2.0 * productArray[i+1];
		}
		integralValue *= h/6.0;
	}
	
	return integralValue;	
}

void analyzeWaveletTD(Parameter_t* parameters, Data_t* data){
	// *** TIME DOMAIN ***
	double freq, time;
	double complex value;

	parameters->k = 1.0;

	FILE *results_TD = fopen("data_results/results_TD.dat", "w");
	clock_t start_TD = clock();
	for(freq=parameters->freq_min; freq <= parameters->freq_max; freq += parameters->delta_freq){
		for(time=parameters->time_min; time <= parameters->time_max; time += parameters->delta_time){
			parameters->b = parameters->omega_0/(2.0*M_PI*freq);
			parameters->k = 1/sqrt(parameters->b); 		// for energy conservation
			createMorletWavelet(data->waveletArray, data->numLines_padded, data->dt_data, parameters->k, parameters->omega_0, time, parameters->b);
			value = calculateConvolutionValue(data->dataArray, data->waveletArray, data->productArray, data->numLines_padded, data->dt_data);
			fprintf(results_TD, "%.4e %.4e %.4e\n", time, freq, cabs(value));
		}	
	}
	clock_t end_TD = clock();
	printf("Time spent for time-domain calculation: %lf s\n", (double)(end_TD - start_TD) / CLOCKS_PER_SEC);
	fclose(results_TD);
	
}

/*
 *
 * ===========================================================================================================
 * 					FREQUENCY-DOMAIN
 * ===========================================================================================================
 */

void analyzeWaveletFD(Parameter_t* parameters, Data_t* data){
	// *** FREQUENCY DOMAIN ***
	clock_t start_FD = clock();
	int numLines = data->numLines_padded;

	// Prepare Fourier-transformation
	fftw_complex *data_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numLines);
	fftw_complex *data_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numLines);
	fftw_plan plan_data = fftw_plan_dft_1d(numLines, data_in, data_out, FFTW_FORWARD, FFTW_MEASURE);

	fftw_complex *psi_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numLines);
	fftw_complex *psi_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numLines);
	fftw_plan plan_psi = fftw_plan_dft_1d(numLines, psi_in, psi_out, FFTW_FORWARD, FFTW_MEASURE);

	fftw_complex *result_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numLines);
        fftw_complex *result_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numLines);
        fftw_plan plan_result = fftw_plan_dft_1d(numLines, result_in, result_out, FFTW_BACKWARD, FFTW_MEASURE);

	// Copy data into real-part of data_in array and set complex part equal to 0.0
	int i;
	for(i=0; i<numLines; i++){
		data_in[i][0] = data->dataArray[i];
		data_in[i][1] = 0.0;
	}
	// Fourier-transformation
	fftw_execute(plan_data);

	int shifted_index;
	double A, B, C, D;
	double freq, time, value;
	FILE *results_FD = fopen("data_results/results_FD.dat", "w");

	parameters->k = 1.0;
	for(freq=parameters->freq_min; freq <= parameters->freq_max; freq += parameters->delta_freq){
		parameters->b = parameters->omega_0/(2.0*M_PI*freq);

		// set k = 1/sqrt(b) for energy conservation
		parameters->k = 1.0 / sqrt(parameters->b);

		createMorletWavelet(data->waveletArray, numLines, data->dt_data, parameters->k, parameters->omega_0, data->numLines/2.0 * data->dt_data, parameters->b);
		
		/*
		 * Here the wavelet needs to be sliced in half/shifted so that its origin is periodically 
		 * centered around 0 to avoid anomalies at the center of the wavelet.
		 */
		for(i=0; i<numLines; i++){
			shifted_index = (i + numLines/2) % numLines;
			psi_in[i][0] = creal(data->waveletArray[shifted_index]);
			psi_in[i][1] = cimag(data->waveletArray[shifted_index]);
			//psi_in[i][0] = creal(waveletArray[i]);
			//psi_in[i][1] = cimag(waveletArray[i]);
		}	
		// Fourier-transformation of the wavelet
		fftw_execute(plan_psi);
		
		// Calculate convolution via the product in frequency-space
		for(i=0; i<numLines; i++){
			A = data_out[i][0];
			B = data_out[i][1];
			C = psi_out[i][0];
			D = -psi_out[i][1];
			
			result_in[i][0] = A*C - B*D;
			result_in[i][1] = A*D + B*C;
		}
		// Fourier-transformation of the results
		fftw_execute(plan_result);

		// Write results to file but taking zero-padding into account
		for(i=data->padding; i<data->numLines+data->padding; i++){
			//time = (i - numLines / 2) * data->dt_data;
			time = (i - data->padding) * (data->dt_data);
			value = sqrt(result_out[i+data->padding][0]*result_out[i+data->padding][0] + result_out[i+data->padding][1]*result_out[i+data->padding][1]) / numLines;
			
			//double threshold = 20.0;
			//if (value < threshold) value = 0.0;
			
			fprintf(results_FD, "%.4e %.4e %.4e\n", time, freq, value);
		}
		fprintf(results_FD, "\n");
	}
	clock_t end_FD = clock();
	printf("Time spent for frequency-domain calculation: %lf s\n", (double)(end_FD - start_FD) / CLOCKS_PER_SEC);
	fclose(results_FD);
	
	fftw_free(data_in);
	fftw_free(data_out);
	fftw_free(psi_in);
	fftw_free(psi_out);
	fftw_free(result_in);
	fftw_free(result_out);
}



void freeData(Data_t* data){
    free(data->dataArray);
    free(data->waveletArray);
    free(data->productArray);
}







