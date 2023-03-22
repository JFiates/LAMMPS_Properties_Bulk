/*********************************************************************************/
/* School of Chemical Engineering, University of Campinas, Brazil                */
/* Code to calculate the self-diffusion coefficient using Green-Kubo equation    */
/* Developer: Dr. Juliane Fiates						 */
/*********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main(){
        /*************************************************************************/
        /* Variables declaration                                                 */
        /*************************************************************************/
	/* Counters								 */
        int i,j,k,step;
	/* Molecules per frame							 */ 
	int af,idx,type;
        /* Box size		 						 */
        float  bxl,bxu;
        float  byl,byu;
        float  bzl,bzu;
	/* Unit conversion 							 */
	double unit_conv;
        /* Total number of steps                                                 */
        int n_steps;
	/* ID number of the first atom in Dump file				 */
	int st_atom; 
        /* Maximum number of steps for integration                               */
        int n_max;
        /* Total number of molecules                                             */ 
        int n_mol;
        /* Number of atoms per molecule                                          */ 
        int n_atoms; 
        /* Masses                                                                */ 
	float mass;
	/* Total mass								 */
	float totmass;
        /* Velocity auto-correlation funcion                                     */ 
        double vcfx,vcfy,vcfz,vx,vy,vz;
        /* Self-diffusion coefficient                                            */ 
        double diff;
        /* Input file                                                            */ 
        char inpfile[400];
        /* Auxiliary string for comments in the input file                       */ 
        char com[1000];
        /* Input file with velocities                                            */ 
        char vel_file[400];
        /* Output file with the velocity auto-correlation function               */ 
        char out_vacf_file[500];
        /* Output file with the self-diffusion coefficient as function of time   */ 
        char diff_file[500]; 

	clock_t start, end;
	double cpu_time_used;

        FILE *in, *out;

        /*************************************************************************/
        /* Reading input file                                                    */
        /*************************************************************************/
	start = clock();
        scanf("%s", inpfile);
        in = fopen(inpfile, "r");
        fscanf(in, "%s", com);
        fscanf(in, "%lf", &unit_conv);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &st_atom);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_steps); 
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_max); 
        fscanf(in, "%s", com);	
        fscanf(in, "%d", &n_mol);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_atoms);
        fscanf(in, "%s", com);
	fscanf(in, "%f", &totmass);
        fscanf(in, "%s", com);
        fscanf(in, "%s", vel_file);
        fscanf(in, "%s", com);
        fscanf(in, "%s", out_vacf_file); 
        fscanf(in, "%s", com);
        fscanf(in, "%s", diff_file); 
        fclose(in);
        int half_n_steps = n_steps/2;
        /*************************************************************************/

        if (n_max > half_n_steps){
            printf(" Error! n_max must be lower than n_steps/2!");
        }
        else{

        /*************************************************************************/
        /* Memory allocation                                                     */
        /*************************************************************************/
        int r = n_steps;
        int c = n_mol*n_atoms;
	int x = n_mol;
        double *vx_n     = calloc (c, sizeof (double));
        double *vy_n     = calloc (c, sizeof (double));
        double *vz_n     = calloc (c, sizeof (double));
	int    *time     = calloc (r, sizeof (int));
	double *vacf_avg = calloc (half_n_steps, sizeof (double));                    
        double *vcx[r], *vcy[r], *vcz[r];
        for (i = 0; i < r; i++){
             vcx[i] = (double*) calloc (x, sizeof(double));
             vcy[i] = (double*) calloc (x, sizeof(double));
             vcz[i] = (double*) calloc (x, sizeof(double));
        }

        /*************************************************************************/
        /* Reading velocity input file 		                                 */
        /*************************************************************************/	
        printf("Reading input file...\n");
        in = fopen(vel_file, "r");
        for (step = 0; step < n_steps; step++){
             fscanf(in, "%d", &time[step]);
	     fscanf(in, "%d", &af);
             fscanf(in, "%e %e", &bxl, &bxu);
	     fscanf(in, "%e %e", &byl, &byu);
	     fscanf(in, "%e %e", &bzl, &bzu);
             for (j = 0; j < (n_mol*n_atoms); j++){
                  fscanf(in, "%d %d %f %lf %lf %lf", &idx, &type, &mass,  &vx, &vy, &vz);
	     	  vx_n[idx]=(vx*mass)/totmass;
	     	  vy_n[idx]=(vy*mass)/totmass;
	     	  vz_n[idx]=(vz*mass)/totmass;
	     }

             /********************************************************************/
             /* Calculation of the velocity of the center of mass of the molecule*/
             /********************************************************************/
	     i=0;
             for(idx = st_atom; idx < (st_atom+(n_mol*n_atoms)); idx=idx+n_atoms){
                 for (k = 0; k < n_atoms; k++){
                      vcx[step][i] = vcx[step][i]+vx_n[idx+k];
                      vcy[step][i] = vcy[step][i]+vy_n[idx+k];
                      vcz[step][i] = vcz[step][i]+vz_n[idx+k];
	           }
	     i++; 
	     }
	   
        }
        fclose(in);

        /*************************************************************************/
        /* Velocity Auto-Correlation Function (VACF) calculation                 */
        /*************************************************************************/	
        printf("Calculating VACF...\n");
        out = fopen(out_vacf_file, "w");
        for (k = 0; k < half_n_steps; k++){
             vcfx = 0.0;
             vcfy = 0.0;
             vcfz = 0.0;
             for (step = 0; step < half_n_steps; step++){
                  for (i = 0; i < n_mol; i++){
                       vcfx = vcfx+vcx[step][i]*vcx[step+k][i];
                       vcfy = vcfy+vcy[step][i]*vcy[step+k][i];
                       vcfz = vcfz+vcz[step][i]*vcz[step+k][i];
                  }
             }
             vacf_avg[k] = (vcfx+vcfy+vcfz)/(double)(n_mol)/(double)(half_n_steps)/3.0;
             fprintf(out,"%d %e %e\n", time[k]-time[0], vacf_avg[k], vacf_avg[k]/vacf_avg[0]);
        }
        fclose(out);


        /*************************************************************************/
        /* Numerical Integration of Velocity Auto-Correlation Function           */
        /* Green-Kubo equation for self-diffusion coefficient                    */
        /*************************************************************************/
        printf("Integrating VACF...\n");
        out = fopen(diff_file, "w");
        for (k = 0; k < n_max; k++){
             diff = vacf_avg[0];
             for (step = 1; step < k-1; step++)
                  diff = diff+2.0*vacf_avg[step];
             diff = unit_conv*0.5*(double)(time[step]-time[0])*(diff+vacf_avg[step])/(double)(k-1);
             if (step > 1)
                 fprintf(out, "%d %e \n", time[step]-time[0], diff);	
        }
        fclose(out);
	end           = clock();
	cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;
	printf("The calculation is ended...\n\n");	
	printf("Elapsed time = %lf seconds...\n\n", cpu_time_used);			
        }

        return 0;
	
}
