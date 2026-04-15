#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <complex.h>
#include <time.h>
#include <string.h>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
 * TODO:
 * - create separate file for each atom, then average them
 * - overlap save stuff such as zero padding and deleting parts
 * - logarithmic scaling for the tested frequencies
 */
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*
 * here we assume that the wavelet is relevant at +-4sigma
 * further more we add 1 so that the wavelet is constructed with an uneven number of sampling points
 */
#define SIGMA 4.0

typedef struct WaveletParameter_s{
    int M;

    double freq_min;
    double freq_max;
    double delta_freq;

    double scale, d, omega_0;
} WaveletParameter_s;

typedef struct DataParameters_s{
    double dt_data;
    int numAtoms;
    int numTimesteps;

    int chunkSize;	// N
    int L;
    int numberOfBlocks;
} DataParameters_s;

typedef struct ArrayContainer{
    double* dataArray;
    double complex* waveletArray;
} ArrayContainer;

typedef struct FFT_wrapper{
    fftw_complex* in;
    fftw_complex* out;

    fftw_plan plan;
} FFT_wrapper;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void setFilenames(char* read_filename);
void initWaveletParameters(double omega_0);
void initDataParameters();
void initArrayContainer();
void createArtificialDataset(double dt_data);
void calcAbsVelocity(char* filename);
int getNumberOfElements(const char* filename);
void getNumberOfTimestepsBinary(const char* write_filename_binary);
int getNumberOfAtoms(const char* filename);
int getDeltaTime(const char* filename);
void setChunkSize();
void changeFileFormat(char* filename);
void copyChunkFromFileToArray(int atom, int chunk, const char* filename, double *targetArray);
void createMorletWavelet(int M);
void WaveletAnalysisFD(int currentBlock, float currentFrequency, FILE* fp);
void prepareData(int chunkSize);
void prepareWavelet(int chunkSize);
void prepareFFT(FFT_wrapper* data, FFT_wrapper* wavelet, FFT_wrapper* result, int chunkSize);
void freeData();

static WaveletParameter_s wp;
static DataParameters_s dp;
static ArrayContainer ac;

FFT_wrapper data, wavelet, result;

/*
 * Here not all files will be needed. A few are just for debugging
 */
char* read_filename;
char* write_filename;
char* write_filename_binary = "abs_trajectory_bin";
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char* argv[]){
    createArtificialDataset(1.0);

    setFilenames(argv[1]);
    initWaveletParameters(atof(argv[2]));
    initDataParameters();

    calcAbsVelocity(read_filename);
    getNumberOfTimestepsBinary(write_filename_binary);
    changeFileFormat(write_filename_binary);

    setChunkSize();
    initArrayContainer();
    prepareFFT(&data, &wavelet, &result, dp.chunkSize);


    int atom, block;
    float freq;
    for(atom=0; atom<dp.numAtoms; atom++){
	clock_t start_FD = clock();

	FILE *results_FD = fopen("data_results/results_FD.dat", "w");

	for(freq=wp.freq_min; freq <= wp.freq_max; freq += wp.delta_freq){
	    wp.scale = wp.omega_0 / (2.0*M_PI*freq*dp.dt_data);

	    wp.M = 2.0 * SIGMA * (int)(wp.scale+1) + 1;	
	    //printf("M: %d\n", wp.M);

	    /*
	     * here the actual number of useable data is calculated
	     */
	    dp.L = dp.chunkSize - (wp.M - 1);

	    int numberOfBlocks = (dp.numTimesteps + dp.L - 1) / dp.L;

	    for(block=0; block<numberOfBlocks; block++){
		copyChunkFromFileToArray(atom, block, write_filename_binary, ac.dataArray);
		prepareData(dp.chunkSize);

		createMorletWavelet(wp.M);
		prepareWavelet(dp.chunkSize);
		WaveletAnalysisFD(block, freq, results_FD);
	    }

	}
	clock_t end_FD = clock();
	printf("Time spent for atom %d: %lf s\n", atom, (double)(end_FD - start_FD) / CLOCKS_PER_SEC);
	fclose(results_FD);
    }

    freeData();
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void setFilenames(char* filename){ 
    read_filename = filename;

    int length = strlen(filename);
    write_filename = malloc((length + 4 + 1) * sizeof(char));
    strcpy(write_filename, "abs_");
    strcat(write_filename, filename);
    printf("read_file: %s\n", read_filename);
    printf("\n");
}

void initWaveletParameters(double omega_0){
    wp.freq_min = 20.0e-3;
    wp.freq_max = 100.0e-3;	
    wp.delta_freq = 0.1e-3; // here include logarithmic scaling not linear
    wp.scale = -1.0;
    wp.d = 1.0;
    wp.omega_0 = omega_0;
}

void initDataParameters(){
    dp.numAtoms = getNumberOfAtoms(read_filename);
    dp.dt_data = getDeltaTime(read_filename);
    dp.numTimesteps = -1;

    dp.L = -1;
    dp.chunkSize = -1;
    dp.numberOfBlocks = -1;
}


void getNumberOfTimestepsBinary(const char* filename){
    // TODO: only valid for one atom?
    FILE* fp = fopen(filename, "rb"); 

    fseek(fp, 0, SEEK_END);
    long filesize = ftell(fp);
    rewind(fp);

    int numElements = filesize / sizeof(double);
    int numTimesteps = numElements / dp.numAtoms;
    dp.numTimesteps = numTimesteps;
    printf("timesteps: %d\n\n", dp.numTimesteps);

    fclose(fp);
}

int getNumberOfAtoms(const char* filename){
    FILE* fp = fopen(filename, "r");
    if (fp == NULL){ 
	fprintf(stderr, "Could not open file! (blocksize func)\n");
	return -1;
    }

    int numberOfAtoms = -1;
    char buffer[256];

    fgets(buffer, sizeof(buffer), fp);
    numberOfAtoms = strtod(buffer, NULL);
    printf("reading number of atoms as: %d\n", numberOfAtoms);
    rewind(fp);

    fclose(fp);

    return numberOfAtoms;
};

/**
 * This function reads dt from the file
 */
int getDeltaTime(const char* filename){
    FILE* fp = fopen(filename, "r");
    if (fp == NULL){ 
	fprintf(stderr, "Could not open file! (deltatime func)\n");
	return -1;
    }

    char buffer[256];
    int line;
    double t1, t2;
    char* unit;
    char* ptr;

    for(line=1; line<=dp.numAtoms + 4; line++){
	fgets(buffer, sizeof(buffer), fp);

	if ( (line == 2) | (line == dp.numAtoms + 4)){
	    ptr = strtok(buffer, ",");
	    ptr = strtok(NULL, ",");
	    //printf("%s\n", ptr);

	    if (line == 2){
		t1 = strtod(strstr(ptr, "=") + 1, &unit);
	    }else{
		t2 = strtod(strstr(ptr, "=") + 1, &unit);
		printf("reading timestep as: %f%s\n", t2-t1, unit);
		printf("\n");
	    }
	}	
    }

    fclose(fp);
    return t2 - t1;
}

/**
 * The initial velocities are given in vel.x vel.y vel.z.
 * This function creates a new file which contains the absolute velocity of each atom at each timestep. 
 * So three columns are reduced to one column.
 */
void calcAbsVelocity(char* filename){

    FILE* fp1 = fopen(filename, "r");
    FILE* fp2 = fopen(write_filename, "w");
    FILE* fp3 = fopen(write_filename_binary, "wb");

    if (fp1 == NULL){ 
	fprintf(stderr, "Could not open copy-file!\n");
    }
    if (fp2 == NULL){ 
	fprintf(stderr, "Could not open paste-file!\n");
    }
    if (fp3 == NULL){ 
	fprintf(stderr, "Could not open binary-paste-file!\n");
    }

    char buffer[256];
    char* ptr;
    int line = 1;

    double velocities[3];
    double absVel;

    while(fgets(buffer, sizeof(buffer), fp1) != NULL){
	if((line % (dp.numAtoms + 2) != 1) && (line % (dp.numAtoms + 2) != 2)){
	    ptr = strtok(buffer, " ");
	    ptr = strtok(NULL, " ");
	    //printf("%s\n", ptr);

	    velocities[0] = strtod(ptr, NULL);
	    ptr = strtok(NULL, " ");
	    velocities[1] = strtod(ptr, NULL);
	    ptr = strtok(NULL, " ");
	    velocities[2] = strtod(ptr, NULL);


	    /*
	     *
	     *
	     *
	     *
	     * 				CHANGE THIS LATER ON IF TESTS ARE PASSED !!!!
	     *
	     *
	     *
	     *
	     *
	     */
	    //absVel = sqrt(velocities[0]*velocities[0] + velocities[1]*velocities[1] + velocities[2]*velocities[2]);
	    absVel = velocities[0];
	    fprintf(fp2, "%f\n", absVel);
	    fwrite(&absVel, sizeof(double), 1, fp3);
	}
	line++;
    }

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
}
/*
 * Until now the only the .xyz file-type-specific lines have been read and removed. At this point each block
 * represents all atoms at time t. 
 * Now this should be changed so that one block contains the trajectory of a certain atom for all times.
 *
 * 		------------    	------------
 * 		|atom1	t1 |		|atom1	t1 |
 * 		|atom2	t1 |	--->	|atom1	t2 |
 * 		|atom3	t1 |		|atom1	t3 |
 * 		|atom4	t1 |		|atom1	t4 |
 * 		------------		-----------
 *
 * This makes the readin later on a bit easier.
 */
void changeFileFormat(char* filename){
    FILE* in = fopen(filename, "rb");
    FILE* out = fopen("abs_trajectory_bin_formatted", "wb");

    double value;
    for(int timestep=0; timestep<dp.numTimesteps; timestep++){
	for(int atom=0; atom<dp.numAtoms; atom++){

	    fread(&value, sizeof(double), 1, in);
	    long offset = (atom * dp.numTimesteps + timestep) * sizeof(double); 
	    fseek(out, offset, SEEK_SET);
	    fwrite(&value, sizeof(double), 1, out);
	}
    }

    fclose(out);


    /*
     * For testing
     *
     
    out = fopen("abs_trajectory_bin_formatted", "rb");
    FILE* new = fopen("test_file.dat", "w");
    rewind(out);
    for(int i=0; i<401; i++){
	fread(&value, sizeof(double), 1, out);
	fprintf(new, "%f\n", value);
    }
    fclose(new);
    fclose(out);
    */

    fclose(in);
}

/*
 * Allocate memory for the signal and the wavelet. 		
 */
void initArrayContainer(){
    ac.dataArray = malloc(dp.chunkSize * sizeof(double));
    if(ac.dataArray == NULL){
	printf("Could not allocate memory for data array!\n");
    }
    memset(ac.dataArray, 0, dp.chunkSize * sizeof(double));


    ac.waveletArray = malloc(dp.chunkSize * sizeof(double complex));
    if(ac.waveletArray == NULL){
	printf("Could not allocate memory for wavelet array!\n");
    }
    memset(ac.waveletArray, 0, dp.chunkSize * sizeof(double complex));
}

void setChunkSize(){
    /*
     * N (chunksize) = L + (M - 1) such that N is power of 2 if possible
     *
     * L is the desired block size and N is the whole chunk, so additionally including M-1 elements from the previous block.
     * Here L is chosen in such a way that N is a power of 2 so that the FFT can use a more efficient algorithm.
     */

    double scale1 = wp.omega_0/(2.0*M_PI*wp.freq_min*dp.dt_data); 
    double scale2 = wp.omega_0/(2.0*M_PI*wp.freq_max*dp.dt_data); 
    printf("scale1: %f (@ %f THz) --> M: %d\n", scale1, wp.freq_min/1e-3, 8*(int)(scale1+1)+1);
    printf("scale2: %f (@ %f THz) --> M: %d\n", scale2, wp.freq_max/1e-3, 8*(int)(scale2+1)+1);

    double scale = wp.omega_0/(2.0*M_PI*wp.freq_min*dp.dt_data);
    long M_max = 8 * (int)(scale+1) + 1;

    if(M_max > dp.numTimesteps){
	printf("Warning: This dataset cannot be used to analyse %f THz! Please choose a higher frequency!\n", wp.freq_min/1e-3);
	double min_freq = (4 * wp.omega_0) / (M_PI * dp.numTimesteps * dp.dt_data);
	printf("The lowest possible frequency would be: %f THz.\n", min_freq / 1e-3);
    }

    int i;
    long N=0, L=0;
    for(i=0; ; i++){
	/*
	 * Trying to find a suitable chunksize as a power of 2
	 */
	if(M_max > pow(2, i)){;	
	    continue;
	}else{
	    N = pow(2, i);
	    if(N <= dp.numTimesteps){
		L = N - M_max + 1;
		if((L<=M_max) && (2*N<dp.numTimesteps)){
		    N = pow(2, i+1);
		    L = N - M_max + 1;
		}
		printf("Setting chunk size to: %ld with a block length of L = %ld\n", N, L);
	    }

	    /*
	     * Choosing L as the first power of 10^j which is smaller than the number of timesteps.
	     */
	    if(N > dp.numTimesteps){
		printf("Not possible to choose a chunksize as a power of 2!\n");
		printf("Falling back to calculting L manually ...\n");

		int j;
		for(j=5;j>=0;j--){
		    L = pow(10, j);
		    N = L + M_max - 1;
		    if(N > dp.numTimesteps){
			continue;
		    }else{
			printf("Setting chunk size to: %ld with a block length of L = %ld\n", N, L);
			break;
		    }
		}
	    }

	    break;
	}
    }
    dp.chunkSize = N;

    printf("N: %ld\n", N);
    printf("L: %ld\n", L);
    printf("M_max: %ld\n", M_max);
}


void copyChunkFromFileToArray(int atom, int block, const char* filename, double *targetArray){
    /*
     * here the block itself and the M-1 preceding elements (which together make up a chunk) will be copied into the array
     *
     * There are a few cases which have to be taken care of:
     * 1) the first block should be copied (so there are no M-1 previous elements)
     * 		--> zeros are added
     * 		same for any block so that blockPos - (M-1) < 0
     * 2) the number of timesteps is not divisible by L (so the last block does not have length L)
     * 		 --> zeros are added such that the last block is of length L
     */

    // TODO: File descriptor is opened each time the function is called. Only open it once per atom.
    FILE *fp = fopen(filename, "rb");
    long atomStart = atom * dp.numTimesteps;
    long atomEnd = (atom+1) * dp.numTimesteps;
    
    long blockStart = atomStart + block * dp.L; 
    long chunkStart = blockStart - (wp.M - 1);
    long chunkEnd = chunkStart + dp.chunkSize;
    
    memset(targetArray, 0, dp.chunkSize * sizeof(double));

    long readStart = chunkStart;
    long readEnd = chunkEnd;

    if(readStart < atomStart){
	readStart = atomStart;
    }

    if(readEnd > atomEnd){
	readEnd = atomEnd;
    }

    long elementsToRead = readEnd - readStart;

    if(elementsToRead > 0){
	long shift = readStart - chunkStart;

	fseek(fp, readStart * sizeof(double), SEEK_SET);

	fread(targetArray + shift, sizeof(double), elementsToRead, fp);
    }

    fclose(fp);
}

void createMorletWavelet(int M){
    /*
     * M should be an uneven number so the wavelet can be centered around 0
     * This function creates a wavelet with M datapoints which are centered around 0 (therefore wrapped) 
     */

    // psi_mother = 1/(pi^0.25 * sqrt(d)) * exp(-t^2/(2*d^2)) * exp(i*w_0*t)
    // psi = 1/sqrt(s) * psi_mother( (t-t0) / s)
    // psi = 1/(pi^0.25 * sqrt(s) * sqrt(d)) * exp(-(t-t_0)^2/(2*s^2*d^2)) * exp(i*w_0*t)

    if(M % 2 == 0){
	printf("wavelet has been created with an even number of sampling points!\n");
    }

    memset(ac.waveletArray, 0, dp.chunkSize * sizeof(double complex));

    /* 
     * TODO: Why do those two block yield different results???
     */

    /*
     * Comment this line out if you want to execute the other block
     */
    #define BLOCK1

    #ifdef BLOCK1
	FILE* fp = fopen("wavelet_test.dat", "w");
	int half = wp.M / 2 + 1;
	for(int n = -half; n <= half; n++) {

	    double tau = n * dp.dt_data / wp.scale;

	    double complex psi =
		1.0 / (pow(M_PI, 0.25) * sqrt(wp.scale * wp.d))
		* exp(-tau*tau/(2.0*wp.d*wp.d))
		* cexp(I * wp.omega_0 * tau);

	    int idx;
	    if(n >= 0)
		idx = n;
	    else
		idx = dp.chunkSize + n;   // negative times at the end of the array

	    ac.waveletArray[idx] = psi;
	}
    #else
	int i;
	double tau;
	double a = M/2.0;

	FILE* fp = fopen("wavelet_test.dat", "w");

	for(i=0; i<M; i++){
	    // wraps the data around 0. (M+1)/2 on the left and (M-1)/2 on the right.
	    // this is done to eliminate artifacts in the fft plot later on.

	    int shifted = (i + M/2 + 1) % M;
	    //int shifted = i;
	    //printf("i: %d --> shifted index: %d\n", i, shifted);
	    a = 0;

	    tau = (shifted - a) * dp.dt_data / wp.scale;

	    ac.waveletArray[shifted] = 1.0 / (pow(M_PI, 1.0/4.0) 
		    * sqrt(wp.scale * wp.d)) 
		* exp(-pow(tau, 2) / (2.0*pow(wp.d, 2.0))) 
		* cexp(I * wp.omega_0 * tau);

	    fprintf(fp, "%f\n", creal(ac.waveletArray[shifted]));
	}
    #endif

    fclose(fp);
}


void prepareFFT(FFT_wrapper* data, FFT_wrapper* wavelet, FFT_wrapper* result, int chunkSize){
    /*
     * Prepare Fourier-transformation
     */
    data->in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunkSize);
    data->out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunkSize);
    data->plan = fftw_plan_dft_1d(chunkSize, data->in, data->out, FFTW_FORWARD, FFTW_MEASURE);

    wavelet->in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunkSize);
    wavelet->out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunkSize);
    wavelet->plan = fftw_plan_dft_1d(chunkSize, wavelet->in, wavelet->out, FFTW_FORWARD, FFTW_MEASURE);

    result->in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunkSize);
    result->out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunkSize);
    result->plan = fftw_plan_dft_1d(chunkSize, result->in, result->out, FFTW_BACKWARD, FFTW_MEASURE);
}

void FFT(FFT_wrapper* wrapper){
    fftw_execute(wrapper->plan);
}

void prepareData(int chunkSize){
    /*
     * Copy data into real-part of data_in array and set complex part equal to 0.0
     */
    int i;
    memset(data.in, 0, chunkSize * sizeof(double complex)); // TODO: needed?
    for(i=0; i<chunkSize; i++){
	data.in[i][0] = ac.dataArray[i];
	data.in[i][1] = 0.0;
    }
}

void prepareWavelet(int chunkSize){
    /* 
     * The created wavelet has length M however the array used for the FFT
     * has to has to have length N so the waveletarray is zero padded
     * Basically a 0-array is created and only the first M values are overwritten
     */
    int i;
    memset(wavelet.in, 0, chunkSize * sizeof(double complex));
    for(i=0; i<wp.M; i++){
	wavelet.in[i][0] = creal(ac.waveletArray[i]);
	wavelet.in[i][1] = cimag(ac.waveletArray[i]);
    }	
}

void WaveletAnalysisFD(int currentBlock, float currentFrequency, FILE* fp){

    int numAtoms = dp.numAtoms;
    int chunkSize = dp.chunkSize;

    FFT(&data);
    FFT(&wavelet);

    /*
     * Calculate convolution via the product in frequency-space
     */
    int i;
    double A, B, C, D;
    for(i=0; i<chunkSize; i++){
	A = data.out[i][0];		// Re{f}
	B = data.out[i][1];		// Im{f}
	C = wavelet.out[i][0];		// Re{Psi}
	D = wavelet.out[i][1];		// Im{Psi}

	/*
	 * Convolution = Psi_conj * f
	 */
	result.in[i][0] = A*C + B*D;
	result.in[i][1] = B*C - A*D;
    }
    // inverse Fourier-transformation of the results
    FFT(&result);



    /*
     * Write results to file
     */
    double time, value;
    int offset = 0;
    //for(i=chunkSize-dp.L; i<chunkSize; i++){
    for(i=wp.M-1; i<wp.M-1+dp.L; i++){
	time = (currentBlock * dp.L + offset) * dp.dt_data; 

	/*
	 * TODO: should not be neccessary
	 * discard values which are greater than the chunkSize (last block is padded with zeros) so gnuplot gets an even grid to plot
	 */
	if(currentBlock *dp.L + offset >= dp.numTimesteps){
	    break;
	}

	value = sqrt(result.out[i][0]*result.out[i][0] + result.out[i][1]*result.out[i][1]) / numAtoms; // TODO: why dividing by numatoms?

	fprintf(fp, "%.4e %.4e %.4e\n", time, currentFrequency, value);
	offset++;
    }
    fprintf(fp, "\n");
}


void freeData(){
    free(ac.dataArray);
    free(ac.waveletArray);

    free(write_filename);

    fftw_free(data.in);
    fftw_free(data.out);
    fftw_free(wavelet.in);
    fftw_free(wavelet.out);
    fftw_free(result.in);
    fftw_free(result.out);
}


void createArtificialDataset(double dt_data){
    /*
     * This function creates an artificial dataset.
     * Here a linear chirp is used for the first 200ps and a constant frequency for the last 200ps
     */
    double t;
    FILE *fp = fopen("hellothere.xyz", "w");
    for(t=0.0; t<400.0; t+=dt_data){
	//printf("t: %f\n", t);
	if(t<=200.0){
	    fprintf(fp, "1\n");
	    fprintf(fp, "# ORCA AIMD Position Step 0, t=%.2f fs, E_Pot=-4816.70201257 Hartree, Unit is Angstrom\n", t);
	    double value = sin(2*M_PI*(0.000117*t*t + 0.02*t));
	    //double value = sin(2*M_PI*(0.05*t));
	    double value2 = pow(3, -(1.0/2.0)) * value;
	    fprintf(fp, "S %.2f %.2f %.2f\n", value, value, value);
	}else{
	    fprintf(fp, "1\n");
	    fprintf(fp, "# ORCA AIMD Position Step 0, t=%.2f fs, E_Pot=-4816.70201257 Hartree, Unit is Angstrom\n", t);
	    //fprintf(fp, "%.2f\n", sin(2*M_PI*(0.5*t*t + 1.0*t)));
	    double value3 = sin(2*M_PI*(0.05*t));
	    double value4 = pow(3, -(1.0/2.0)) * value3;
	    fprintf(fp, "S %.2f %.2f %.2f\n", value4, value4, value4);
	}
    }
    fclose(fp);
}


