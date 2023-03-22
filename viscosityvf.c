/********************************************************************************/
/* School of Chemical Engineering, University of Campinas, Brazil               */
/* Code to calculate viscosity using Green-Kubo equation   			*/
/* Developer: Dr. Juliane Fiates						*/
/********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main(){
	/************************************************************************/
	/* Variables declaration                                                */
	/************************************************************************/	
	/* Counter 								*/
	int k;
	/* Counter for data points 						*/
	int step;
        /* Total number of steps                                                */
        int n_steps;
        /* Maximum number of steps for integration				*/
        int n_max;
	/* Counter for pressure tensor components 				*/
	int  p;
	/* Boltzmann constant in m2.kg/(s2.K) 					*/
	float kb;
	/* Box volume in m3 							*/
	float vol, volume, density;
	/* Absolute temperature in K 						*/
	float t; 
	/* CONST = V/kB/T 							*/
	float cons; 
	/* Stress correlation functions 					*/
	float sxy, sxz, syz, sx1, sy1, sz1;
	/* Viscosity in cP = mPa.s 						*/
	float vis;
	/* Sum of each pressure tensor component 				*/
	float sum[6];
	/* Average of each pressure tensor component 				*/
	float avg[6];
        /* Input file                                                           */ 
        char inpfile[400];
        /* Auxiliary string for comments in the input file                      */ 
        char com[1000];
        /* Input file with the tensor components                                */ 
        char p_file[400];
        /* Output file with the stress auto-correlation function                */ 
        char out_sacf_file[500];
        /* Output file with the viscosity coefficient as function of time       */ 
        char vis_file[500]; 

	clock_t start, end;
	double cpu_time_used;

	FILE *in, *out;
	
	/************************************************************************/	
	/* Physical constants                                                   */
	/************************************************************************/	
	kb = 1.38064852e-23;
	

	/************************************************************************/
	/* Reading input file                                                   */
	/************************************************************************/
	start = clock();
	scanf("%s", inpfile);
        in = fopen(inpfile, "r");
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
        int half_n_steps = n_steps/2;

        if (n_max > half_n_steps){
            printf(" Error! n_max must be lower than n_steps/2!");
        }
        else{

        /************************************************************************/
        /* Memory allocation                                                    */
        /************************************************************************/
        int r = n_steps;
	/* Time in fs 								*/
	int   *time   = calloc (r, sizeof (int));
	/* Pressure tensor components in atm 					*/
	float presxx, presyy, preszz;  
        float *presxy   = calloc (r, sizeof (float));
        float *presxz   = calloc (r, sizeof (float));
        float *presyz   = calloc (r, sizeof (float));
        float *presx1   = calloc (r, sizeof (float));
        float *presy1   = calloc (r, sizeof (float));
        float *presz1   = calloc (r, sizeof (float));
	/* Normalized pressure tensor components 				*/
        float *presxy_n = calloc (r, sizeof (float));
        float *presxz_n = calloc (r, sizeof (float));
        float *presyz_n = calloc (r, sizeof (float));
        float *presx1_n = calloc (r, sizeof (float));
        float *presy1_n = calloc (r, sizeof (float));
        float *presz1_n = calloc (r, sizeof (float));
	/* Average of the stress correlation function 				*/
        float *savg     = calloc (half_n_steps, sizeof (float));
  

        /************************************************************************/
        /* Reading pressure tensor components file 		                */
        /************************************************************************/
	in = fopen(p_file, "r");
	printf("\nReading input file...\n");
        for (step = 0; step < n_steps; step++){
	     fscanf(in, "%d %f %f %f %f %f %f %f %f", &time[step], &presxx, &presxy[step], &presxz[step], &presyy, &presyz[step], &preszz, &density, &volume);
	     presx1[step] = 0.5*(presxx - presyy); 	
	     presy1[step] = 0.5*(presyy - preszz); 	
	     presz1[step] = 0.5*(presxx - preszz);
	}
	fclose(in);
	vol   = volume*1e-30;
	cons  = 0.010266755*vol/kb/t;
        /************************************************************************/
        /* Normalization of the pressure tensor components 	 		*/
        /************************************************************************/
	for (p = 0; p < 6; p++){
	     sum[p] = 0.0;
	     avg[p] = 0.0;
	}	
        for (step = 0; step < n_steps; step++){
	     sum[0] = sum[0]+presx1[step];
	     sum[1] = sum[1]+presxy[step];
	     sum[2] = sum[2]+presxz[step];
	     sum[3] = sum[3]+presy1[step];
	     sum[4] = sum[4]+presyz[step];
	     sum[5] = sum[5]+presz1[step];
	}
	for (p = 0; p < 6; p++)
	     avg[p] = sum[p]/(float)(n_steps);	
	for (k = 0; k < n_steps; k++){
	     presx1_n[k] = presx1[k]-avg[0];
	     presxy_n[k] = presxy[k]-avg[1];
	     presxz_n[k] = presxz[k]-avg[2];
	     presy1_n[k] = presy1[k]-avg[3];
	     presyz_n[k] = presyz[k]-avg[4];
	     presz1_n[k] = presz1[k]-avg[5];
	}

        /************************************************************************/
        /* Memory deallocation	                                                */
        /************************************************************************/
        free (presx1); 
        free (presxy);
        free (presxz);
        free (presy1); 
        free (presyz);
        free (presz1);

	/************************************************************************/
	/* Stress Auto-Correlation Function calculation                         */
	/************************************************************************/
	printf("Calculating Stress ACF...\n");	
	out = fopen(out_sacf_file, "w");
	for (k = 0; k < half_n_steps; k++){
	     sxy = 0.0;
	     sxz = 0.0;
	     syz = 0.0;
	     sx1 = 0.0;
	     sy1 = 0.0;
	     sz1 = 0.0;
	     for (step = 0; step < half_n_steps; step++){
		  sxy = sxy+presxy_n[step]*presxy_n[step+k];
		  sxz = sxz+presxz_n[step]*presxz_n[step+k];
		  syz = syz+presyz_n[step]*presyz_n[step+k];
		  sx1 = sx1+presx1_n[step]*presx1_n[step+k];
		  sy1 = sy1+presy1_n[step]*presy1_n[step+k];
		  sz1 = sz1+presz1_n[step]*presz1_n[step+k];
	     }
	     sxy = sxy/(float)(half_n_steps); 
	     sxz = sxz/(float)(half_n_steps); 
	     syz = syz/(float)(half_n_steps); 
	     sx1 = sx1/(float)(half_n_steps); 
	     sy1 = sy1/(float)(half_n_steps); 
	     sz1 = sz1/(float)(half_n_steps); 
	     savg[k] = (sxy+sxz+syz+sx1+sy1+sz1)/6.0;
	     fprintf(out, "%d %f %f\n", time[k]-time[0], savg[k], savg[k]/savg[0]);
	}
	fclose(out);

	/************************************************************************/
	/* Integration of Stress Auto-Correlation Function 		 	*/
	/* Green-Kubo equation for viscosity                                    */
        /************************************************************************/
	printf("Integrating Stress ACF to calculate viscosity...\n");	
	out = fopen(vis_file, "w");
	for (k = 0; k < n_max; k++){
	     vis = savg[0];
	     for (step = 1; step < k-1; step++)
		  vis = vis+2.0*savg[step];
	     vis = cons*0.5*(float)(time[step]-time[0])*(vis+savg[step])/(float)(k-1);
	     if (step > 1) 
		 fprintf(out, "%d %f\n", time[step]-time[0], vis);
	}
	fclose(out);
	end           = clock();
	cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;
	printf("The calculation is ended...\n\n");	
	printf("Elapsed time = %lf seconds...\n\n", cpu_time_used);	
	}
	
	return 0;

}
	
