/********************************************************************************/
/* School of Chemical Engineering, University of Campinas, Brazil               */
/* Code to calculate viscosity using Green-Kubo equation   			*/
/* Developer: Dr. Juliane Fiates						*/
/* Paralelization: Leandro Negrini Zanotto					*/
/********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <new>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cuda.h>
using namespace std;

__global__ void integrate_kernel(const float *savg, const int *time,
	float *vis_vec, const int n_max, const float cons) {

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int step;
	float vis = savg[0];
	if (id < n_max) {		
		for (step = 1; step < id - 1; step++)
			vis += 2.0 * savg[step];			
		
		vis_vec[id] = cons * (float) (time[step] - time[0])
				* (vis + savg[step]) / (float) (id - 1);
	}
}

__global__ void stress_acf_kernel(const float *presxy, const float *presxz,
		const float *presyz, const float *presx1, const float *presy1,
		const float *presz1, const int half_n_steps, float *savg) {

	// Get our global thread ID
	int id = blockIdx.x * blockDim.x + threadIdx.x;

    float sxy_n, sxz_n, syz_n, sx1_n, sy1_n, sz1_n;
    float sxy = 0.0, sxz = 0.0, syz = 0.0, sx1 = 0.0, sy1 = 0.0, sz1 = 0.0;
    float f_half_n_steps = float(half_n_steps);
	// Make sure we do not go out of bounds
	if (id < half_n_steps) {
		for (int step = 0; step < half_n_steps; step++) {
			sxy += presxy[step] * presxy[step + id];
			sxz += presxz[step] * presxz[step + id];
			syz += presyz[step] * presyz[step + id];
			sx1 += presx1[step] * presx1[step + id];
			sy1 += presy1[step] * presy1[step + id];
			sz1 += presz1[step] * presz1[step + id];

		}
		sxy_n = sxy / f_half_n_steps;
		sxz_n = sxz / f_half_n_steps;
		syz_n = syz / f_half_n_steps;
		sx1_n = sx1 / f_half_n_steps;
		sy1_n = sy1 / f_half_n_steps;
		sz1_n = sz1 / f_half_n_steps;
		savg[id] =  (sxy_n + sxz_n + syz_n + sx1_n + sy1_n + sz1_n) / 6.0;
	}
	//__syncthreads();

   //integrate_kernel(savg, time, vis_vec, n_max, cons);

}


int main(int argc, char **argv) {
	/************************************************************************/
	/* Variables declaration                                                */
	/************************************************************************/
	/* Counter 								*/
	//int k;
	/* Counter for data points 						*/
	int step;
	/* Total number of steps                                                */
	int n_steps;
	/* Maximum number of steps for integration				*/
	int n_max;
	/* Counter for pressure tensor components 				*/
	int p;
	/* Boltzmann constant in m2.kg/(s2.K) 					*/
	float kb;
	/* Box volume in m3 							*/
	float volume, density;
	/* Absolute temperature in K 						*/
	float t;
	/* CONST = V/kB/T 							*/
	/* Viscosity in cP = mPa.s 						*/
	/* Sum of each pressure tensor component 				*/
	float sum[6] = { };
	/* Average of each pressure tensor component 				*/
	float avg[6] = { };
	/* Auxiliary string for comments in the input file                      */
	char com[1000];
	/* Input file with the tensor components                                */
	char p_file[400];
	/* Output file with the stress auto-correlation function                */
	char out_sacf_file[500];
	/* Output file with the viscosity coefficient as function of time       */
	char vis_file[500];

	double start = omp_get_wtime();

	FILE *in, *out;

	/************************************************************************/
	/* Physical constants                                                   */
	/************************************************************************/
	kb = 1.38064852e-23;

	/************************************************************************/
	/* Reading input file                                                   */
	/************************************************************************/
	in = fopen(argv[1], "r");
	if ((in = fopen(argv[1], "r")) == NULL) {
		cout << "No such file" << "\n" << argv[1];
		exit(1);
	}

	fscanf(in, "%s", com);
	fscanf(in, "%f", &t);
	fscanf(in, "%s", com);
	fscanf(in, "%d", &n_steps);
	fscanf(in, "%s", com);
	fscanf(in, "%d", &n_max);
	fscanf(in, "%s", com);
	fscanf(in, "%s", p_file);
	fscanf(in, "%s", com);
	fscanf(in, "%s", out_sacf_file);
	fscanf(in, "%s", com);
	fscanf(in, "%s", vis_file);
	fclose(in);
	int half_n_steps = n_steps / 2;

	if (n_max > half_n_steps) {
		printf(" Error! n_max must be lower than n_steps/2!");
		return 0;
	}

	/************************************************************************/
	/* Memory allocation                                                    */
	/************************************************************************/
	/* Time in fs 								*/
	int *time = new int[n_steps];

	/* Pressure tensor components on CPU 	    		                */
	float presxx, presyy, preszz;
	float *presxy = new float[n_steps]();
	float *presxz = new float[n_steps]();
	float *presyz = new float[n_steps]();
	float *presx1 = new float[n_steps]();
	float *presy1 = new float[n_steps]();
	float *presz1 = new float[n_steps]();
	/* Average of the stress correlation function 				*/
	float *savg = new float[half_n_steps];
	float *vis_vec = new float[n_max];

	// Size, in bytes, of each vector
	size_t bytes = n_steps * sizeof(float);
	/* Pressure tensor components on GPU 	    		                */
	float *presxy_d = NULL;
	float *presxz_d = NULL;
	float *presyz_d = NULL;
	float *presx1_d = NULL;
	float *presy1_d = NULL;
	float *presz1_d = NULL;
	float *vis_vec_d = NULL;
	int *time_d = NULL;
	/* Average of the stress correlation function 				*/
	float *savg_d;

	/************************************************************************/
	/* Reading pressure tensor components file 				                */
	/************************************************************************/
	in = fopen(p_file, "r");
	cout << "\nReading input file...\n";
	for (int step = 0; step < n_steps; step++) {
		fscanf(in, "%d %f %f %f %f %f %f %f %f", &time[step], &presxx,
				&presxy[step], &presxz[step], &presyy, &presyz[step], &preszz,
				&density, &volume);
		presx1[step] = 0.5 * (presxx - presyy);
		presy1[step] = 0.5 * (presyy - preszz);
		presz1[step] = 0.5 * (presxx - preszz);
	}
	fclose(in);
	volume = volume * 1e-30;
	float cons = (0.010266755 * volume / kb / t) * 0.5;
	
	/************************************************************************/
	/* Normalization of the pressure tensor components 	 		*/
	/************************************************************************/
	for (int step = 0; step < n_steps; step++) {
		sum[0] += presx1[step];
		sum[1] += presxy[step];
		sum[2] += presxz[step];
		sum[3] += presy1[step];
		sum[4] += presyz[step];
		sum[5] += presz1[step];
	}

	for (p = 0; p < 6; p++)
		avg[p] = sum[p] / (float) (n_steps);

	for (int k = 0; k < n_steps; k++) {
		presx1[k] -= avg[0];
		presxy[k] -= avg[1];
		presxz[k] -= avg[2];
		presy1[k] -= avg[3];
		presyz[k] -= avg[4];
		presz1[k] -= avg[5];
	}

	int deviceCount = 0;
        cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess) {
	    printf("cudaGetDeviceCount returned %d\n-> %s\n",
            static_cast<int>(error_id), cudaGetErrorString(error_id));
	    printf("Result = FAIL\n");
	    exit(EXIT_FAILURE);
	}

 	// This function call returns 0 if there are no CUDA capable devices.
	if (deviceCount == 0) {
	    printf("There are no available device(s) that support CUDA the Code will run on CPU\n");
	} else {
	    printf("Detected %d CUDA Capable device(s)\n", deviceCount);
        }
	
	if (deviceCount != 0) {
		

		// Allocate memory for each vector on GPU
		cudaMalloc(&presxy_d, bytes);
		cudaMalloc(&presxz_d, bytes);
		cudaMalloc(&presyz_d, bytes);
		cudaMalloc(&presx1_d, bytes);
		cudaMalloc(&presy1_d, bytes);
		cudaMalloc(&presz1_d, bytes);
		cudaMalloc(&vis_vec_d, n_max * sizeof(float));
		cudaMalloc(&time_d, n_steps * sizeof(int));
		cudaMalloc(&savg_d, half_n_steps * sizeof(float));


		/************************************************************************/
		/* Copying Arrays from host to device			 		*/
		/************************************************************************/
		cudaMemcpy(presxy_d, presxy, bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(presxz_d, presxz, bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(presyz_d, presyz, bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(presx1_d, presx1, bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(presy1_d, presy1, bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(presz1_d, presz1, bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(time_d, time, n_steps * sizeof(int), cudaMemcpyHostToDevice);

		/************************************************************************/
		/* Stress Auto-Correlation Function calculation                         */
		/************************************************************************/

		cout << "Calculating Stress ACF on GPU and Integrating on GPU...\n";

		int threadsPerBlock = 1024;
		int blocksPerGrid = (half_n_steps + threadsPerBlock - 1) / threadsPerBlock;

		// Execute the kernel
		stress_acf_kernel<<<blocksPerGrid, threadsPerBlock>>>(presxy_d, presxz_d, presyz_d,
				presx1_d, presy1_d, presz1_d, half_n_steps, savg_d);

		blocksPerGrid = (n_max + threadsPerBlock - 1) / threadsPerBlock;
		integrate_kernel<<<blocksPerGrid, threadsPerBlock>>>(savg_d, time_d, vis_vec_d, n_max, cons);


		cudaDeviceSynchronize();
		/** Copying the array back to write the file **/
		cudaMemcpy(savg, savg_d, half_n_steps * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(vis_vec, vis_vec_d, n_max * sizeof(float), cudaMemcpyDeviceToHost);
	} else {


	//Run on CPU in Parallel
	//Stress correlation functions
	int numthreads = omp_get_num_threads();
	cout << "Calculating Stress ACF on CPU using " << "\t" << numthreads << "\n";
	float f_half_n_steps = (float) half_n_steps;
	float sxy_n = 0.0, sxz_n = 0.0, syz_n = 0.0, sx1_n = 0.0, sy1_n = 0.0, sz1_n = 0.0;
		 #pragma omp parallel for private(sxy_n,sxz_n,syz_n,sx1_n,sy1_n,sz1_n)
		 for (int k = 0; k < half_n_steps; k++){
			 float sxy = 0.0, sxz = 0.0, syz = 0.0, sx1 = 0.0, sy1 = 0.0, sz1 = 0.0;
			 for (int step = 0; step < half_n_steps; step++){
				 sxy += presxy[step] * presxy[step + k];
				 sxz += presxz[step] * presxz[step + k];
				 syz += presyz[step] * presyz[step + k];
				 sx1 += presx1[step] * presx1[step + k];
				 sy1 += presy1[step] * presy1[step + k];
				 sz1 += presz1[step] * presz1[step + k];
			 }
		 sxy_n = sxy / f_half_n_steps;
		 sxz_n = sxz / f_half_n_steps;
		 syz_n = syz / f_half_n_steps;
		 sx1_n = sx1 / f_half_n_steps;
		 sy1_n = sy1 / f_half_n_steps;
		 sz1_n = sz1 / f_half_n_steps;
		 savg[k] = (sxy_n + sxz_n + syz_n + sx1_n + sy1_n + sz1_n)/6.0;
	 }


	/************************************************************************/
	/* Integration of Stress Auto-Correlation Function 		 	*/
	/* Green-Kubo equation for viscosity                                    */
	/************************************************************************/
	cout << "Integrating Stress ACF to calculate viscosity on CPU...\n";

	 // Run on CPU in Parallel	  
	 #pragma omp parallel for private(step)
	 for (int k = 0; k < n_max; k++){
		 float vis = savg[0];
		 for (step = 1; step < k-1; step++){
			 vis += 2.0 * savg[step];			 
		 }
	 	 vis_vec[k] = cons * (float)(time[step] - time[0]) * (vis + savg[step]) / (float)(k-1);
	 }
	}	

 	cout << "Writing the Files...\n";
	out = fopen(out_sacf_file, "w");
	for (int k = 0; k < half_n_steps; k++)
	    fprintf(out, "%d %f %f\n", time[k]-time[0], savg[k], savg[k]/savg[0]);
	fclose(out);

	 out = fopen(vis_file, "w");
	 for (int k = 3; k < n_max; k++){
		 for (step = 1; step < k-1; step++){}
		 fprintf(out, "%d %f\n", time[step]-time[0], vis_vec[k]);
	 }
	 fclose(out);

	//Free Memory from CPU and GPU
	delete[] presxy;
	delete[] presxz;
	delete[] presyz;
	delete[] presx1;
	delete[] presy1;
	delete[] presz1;
	delete[] savg;
	delete[] vis_vec;

	if (deviceCount != 0) {
		cudaFree(presxy_d);
		cudaFree(presxz_d);
		cudaFree(presyz_d);
		cudaFree(presx1_d);
		cudaFree(presy1_d);
		cudaFree(presz1_d);
		cudaFree(savg_d);
		cudaFree(vis_vec_d);
		cudaFree(time_d);
	}
	double end = omp_get_wtime();
	cout << "The calculation is ended...\n";
	cout << "Elapsed Time: " << (end - start) << "s";

	return 0;

}

