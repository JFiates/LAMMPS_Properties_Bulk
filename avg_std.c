/*********************************************************************************/
/* School of Chemical Engineering, University of Campinas, Brazil    		 */
/* Code to calculate the average and standard deviation of properties	 	 */
/* Developer: Dr. Juliane Fiates						 */
/*********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(){
        /*************************************************************************/
        /* Variables declaration                                                 */
        /*************************************************************************/
	/* Counters								 */
        int i,k,step;
	float ppt, ppt2, idx;
	/* Number of trajectories						 */
	int n;
        /* Total number of steps                                                 */
        int n_steps;
        /* Input file                                                            */ 
        char inpfile[100];
        /* Auxiliary string for comments in the input file                       */ 
        char com[100];
        /* Input files with the properties                                       */
        char filename[100];
        /* Output file with the average and standard deviation	                 */ 
        char out_file[100];


        FILE *in, *out;

        /*************************************************************************/
        /* Reading input file                                                    */
        /*************************************************************************/
        scanf("%s", inpfile);
        in = fopen(inpfile, "r");
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_steps);
        fscanf(in, "%s", com);
        fscanf(in, "%s", out_file); 
	fclose(in);

        /*************************************************************************/
        /* Memory allocation                                                     */
        /*************************************************************************/
        int r = n;
	int x = n_steps;
	float *sum1 = calloc (x, sizeof(float));
	float *sum2 = calloc (x, sizeof(float));
	float *avg  = calloc (x, sizeof(float));
	float *var  = calloc (x, sizeof(float));
	float *std  = calloc (x, sizeof(float));
        float *ppt_n[r];
        for (i = 1; i <= r; i++)
             ppt_n[i] = (float*) calloc (x, sizeof(float));

        /*************************************************************************/
        /* Reading property input file 		                                 */
        /*************************************************************************/	
        printf("Reading input file...\n");
	for (i = 1; i <= n; i++){
	     sprintf(filename, "prop%d.data", i);
             in = fopen(filename, "r");
	     if (in == NULL){
                printf("Change the information of the input file in line 69! \n");
                exit(1);
	     } 
             for (step = 0; step < n_steps; step++){
		  fscanf(in, "%f %f %f", &idx, &ppt, &ppt2);
		  //printf("%d %f %f \n", i, idx, vis);
                  ppt_n[i][step] = ppt;
		  //printf("%d %d %f \n", i, idxx, vis_n[i][idxx]);
	     }


	}
	fclose(in);

        /*************************************************************************/
        /* Calculation of the average, variance and standard deviation           */
        /*************************************************************************/
        printf("Calculating average...\n");
        out = fopen(out_file, "w");
        for (step = 0; step < n_steps; step++){	
	     for (i = 1; i <= n; i++){
		  sum1[step] = sum1[step] + ppt_n[i][step];
	     }
	avg[step] = sum1[step]/(float)(n);
	}
        printf("Calculating variance and standard deviation...\n");
        for (step = 0; step < n_steps; step++){	
	     for (i = 1; i <= n; i++){
		  sum2[step] = sum2[step] + pow((ppt_n[i][step]-avg[step]), 2);
	     }
	var[step] = sum2[step]/(float)(n-1);
	std[step] = sqrt(var[step]);
	//printf("%d %f %f \n", step, avg[step], std[step]);
	fprintf(out, "%f %f %f \n", (float)(step)/10.0, avg[step], std[step]);	
        }
        fclose(out);
		


        return 0;
}
