/******************************************************************************************/
/* School of Chemical Engineering, University of Campinas, Brazil               	  */
/* Code to calculate the ionic conductivity coefficient using Green-Kubo equation         */
/* Developer: Dr. Juliane Fiates							  */
/******************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main(){
        /*********************************************************************************/
        /* Variables declaration                                                	 */
        /*********************************************************************************/
	/* Counters								 	 */
        int i,j,k,l,step;
	/* Boltzmann constant in m2.kg/(s2.K)						 */
	double kb;
	/* Elemental electric charge in Coulomb						 */
	double ee;
	/* Temperature in K							 	 */
	float t;
	/* Box volume in m3								 */
	float vol;
	/* Constant									 */
	double cons;
	/* Molecules per frame								 */ 
	int af,idx,type;
        /* Box size		 							 */
        float  bxl, bxu;
        float  byl, byu;
        float  bzl, bzu;
	/* Unit conversion 								 */
	double unit_conv;
	/* Anion electrical charge							 */
	float qa;
	/* Cation electrical charge						 	 */
	float qc;
        /* Total number of steps							 */
        int n_steps;
        /* Maximum number of steps for integration					 */
        int n_max;
        /* Total number of molecules							 */ 
        int n_mol;
        /* Masses                                                                	 */ 
	float mass;
	/* Total mass								 	 */
	float totmass;
        /* Number of atoms per molecule of anion 					 */ 
        int n_anion;
	/* ID number of the first anion atom in Dump file				 */
	int st_at_a;  
        /* Number of atoms per molecule of cation					 */ 
        int n_cation; 
	/* ID number of the first cation atom in Dump file				 */
	int st_at_c; 
        /* Velocities				                                         */ 
        double vx,vy,vz; 
        /* Electrical current auto-correlation funcion					 */ 
        double jcfx,jcfy,jcfz;
        /* Electrical conductivity coefficient in Si/m					 */ 
        double cond;
        /* Input file									 */ 
        char inpfile[400];
        /* Auxiliary string for comments in the input file				 */ 
        char com[1000];
        /* Input anion file with velocities						 */ 
        char anion_vel_file[400];
        /* Input cation file with velocities 						 */ 
        char cation_vel_file[400];
        /* Output file with the electrical auto-correlation function			 */ 
        char out_jacf_file[500];
        /* Output file with the electrical conductivity coefficient as function of time  */ 
        char cond_file[500]; 

	clock_t start, end;
	double cpu_time_used;

        FILE *in, *out;

        /*********************************************************************************/
        /* Reading input file                                                    	 */
        /*********************************************************************************/
        scanf("%s", inpfile);
        in = fopen(inpfile, "r");
        fscanf(in, "%s", com);
        fscanf(in, "%lf", &unit_conv); 
        fscanf(in, "%s", com);
        fscanf(in, "%f", &t); 
        fscanf(in, "%s", com);
        fscanf(in, "%f", &vol); 
        fscanf(in, "%s", com);
        fscanf(in, "%f", &qa); 
        fscanf(in, "%s", com);
        fscanf(in, "%f", &qc); 
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_steps); 
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_max); 
        fscanf(in, "%s", com);	
        fscanf(in, "%d", &n_mol);
        fscanf(in, "%s", com);
	fscanf(in, "%f", &totmass);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_anion);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &st_at_a);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_cation);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &st_at_c);
        fscanf(in, "%s", com);
        fscanf(in, "%s", anion_vel_file);
        fscanf(in, "%s", com);
        fscanf(in, "%s", cation_vel_file);
        fscanf(in, "%s", com);
        fscanf(in, "%s", out_jacf_file); 
        fscanf(in, "%s", com);
        fscanf(in, "%s", cond_file); 
        fclose(in);
        int half_n_steps = n_steps/2;
        /*********************************************************************************/

	kb   = 1.38064852e-23;
	ee   = 1.60217657e-19;
	vol  = vol * 1e-30;
	cons = unit_conv*ee*ee/kb/t/vol;

        if (n_max > half_n_steps){
            printf(" Error! n_max must be lower than n_steps/2!");
        }
        else{

        /********************************************************************************/
        /* General memory allocation                                                  	*/
        /********************************************************************************/
        int r = n_steps;
	int x = n_mol;

	int    *time  = calloc (r, sizeof (int));
	double *j_avg = calloc (half_n_steps, sizeof (double));
        double *jx[r], *jy[r], *jz[r];
             for (i = 0; i < r; i++){
             jx[i] = (double*) calloc (x, sizeof(double));
             jy[i] = (double*) calloc (x, sizeof(double));
             jz[i] = (double*) calloc (x, sizeof(double));
        }  


        /********************************************************************************/
        /* Memory allocation to anion                                                  	*/
        /********************************************************************************/
        int a = n_mol*n_anion;

        double *vx_na = calloc (a, sizeof (double));
        double *vy_na = calloc (a, sizeof (double));
        double *vz_na = calloc (a, sizeof (double));
        double *vcxa  = calloc (a, sizeof (double));
        double *vcya  = calloc (a, sizeof (double));
        double *vcza  = calloc (a, sizeof (double));
        double *jxa[r], *jya[r], *jza[r];
             for (i = 0; i < r; i++){
             jxa[i] = (double*) calloc (x, sizeof(double));
             jya[i] = (double*) calloc (x, sizeof(double));
             jza[i] = (double*) calloc (x, sizeof(double));
        }

	
        /********************************************************************************/
        /* Reading velocity input file of the anion                              	*/
        /********************************************************************************/	
        printf("Reading anion input file...\n");
        in = fopen(anion_vel_file, "r");
        for (step = 0; step < n_steps; step++){
             fscanf(in, "%d", &time[step]);
	     //printf("%d\n", time[step]);
	     fscanf(in, "%d", &af);
             fscanf(in, "%e %e", &bxl, &bxu);
	     fscanf(in, "%e %e", &byl, &byu);
	     fscanf(in, "%e %e", &bzl, &bzu);
	     //printf("%e %e\n", bzl, bzu);
             for (j = 0; j < (n_mol*n_anion); j++){
                  fscanf(in, "%d %d %f %lf %lf %lf", &idx, &type, &mass,  &vx, &vy, &vz);
                 // printf("%d %d %f  %lf %lf %lf\n", idx, type, mass,vx, vy, vz);
	     	  vx_na[idx]=(vx*mass)/totmass;
	     	  vy_na[idx]=(vy*mass)/totmass;
	     	  vz_na[idx]=(vz*mass)/totmass;
	     }

             /********************************************************************/
             /* Calculation of the velocity of the center of mass of the molecule*/
             /********************************************************************/
	     i=0;
             for(idx = st_at_a; idx < (st_at_a+(n_mol*n_anion)); idx=idx+n_anion){
                 for (k = 0; k < n_anion; k++){
                      vcxa[i] = vcxa[i]+vx_na[idx+k];
                      vcya[i] = vcya[i]+vy_na[idx+k];
                      vcza[i] = vcza[i]+vz_na[idx+k];
	              //printf("%d %d %lf %lf %lf\n",step, i,vcx[step][i], vcy[step][i], vcz[step][i]);
             	/************************************************************************/
             	/* Calculation of the electrical current of the anion		        */
             	/************************************************************************/
	     	      jxa[step][i] = qa * vcxa[i];
		      jya[step][i] = qa * vcya[i];
		      jza[step][i] = qa * vcza[i]; 
	         }
	     i++;
	     }  
        }
        fclose(in);

        /********************************************************************************/
        /* Memory deallocation to anion                                                  */
        /********************************************************************************/
        free (vx_na); 
        free (vy_na);
        free (vz_na);
        free (vcxa); 
        free (vcya);
        free (vcza);



        /********************************************************************************/
        /* Memory allocation to cation                                                 	*/
        /********************************************************************************/
        int c = n_mol*n_cation;

        double *vx_nc = calloc (c, sizeof (double));
        double *vy_nc = calloc (c, sizeof (double));
        double *vz_nc = calloc (c, sizeof (double));
        double *vcxc  = calloc (c, sizeof (double));
        double *vcyc  = calloc (c, sizeof (double));
        double *vczc  = calloc (c, sizeof (double));
        double *jxc[r], *jyc[r], *jzc[r];
             for (i = 0; i < r; i++){
             jxc[i] = (double*) calloc (x, sizeof(double));
             jyc[i] = (double*) calloc (x, sizeof(double));
             jzc[i] = (double*) calloc (x, sizeof(double));
        }

        /********************************************************************************/
        /* Reading velocity input file of the cation                              	*/
        /********************************************************************************/	
        printf("Reading cation input file...\n");
        in = fopen(cation_vel_file, "r");
        for (step = 0; step < n_steps; step++){
             fscanf(in, "%d", &time[step]);
	     //printf("%d\n", time[step]);
	     fscanf(in, "%d", &af);
             fscanf(in, "%e %e", &bxl, &bxu);
	     fscanf(in, "%e %e", &byl, &byu);
	     fscanf(in, "%e %e", &bzl, &bzu);
	     //printf("%e %e\n", bzl, bzu);
             for (j = 0; j < (n_mol*n_cation); j++){
                  fscanf(in, "%d %d %f %lf %lf %lf", &idx, &type, &mass,  &vx, &vy, &vz);
                 // printf("%d %d %f  %lf %lf %lf\n", idx, type, mass,vx, vy, vz);
	     	  vx_nc[idx]=(vx*mass)/totmass;
	     	  vy_nc[idx]=(vy*mass)/totmass;
	     	  vz_nc[idx]=(vz*mass)/totmass;
	     }

             /********************************************************************/
             /* Calculation of the velocity of the center of mass of the molecule*/
             /********************************************************************/
	     i=0;
             for(idx = st_at_c; idx < (st_at_c+(n_mol*n_cation)); idx=idx+n_cation){
                 for (k = 0; k < n_cation; k++){
                      vcxc[i] = vcxc[i]+vx_nc[idx+k];
                      vcyc[i] = vcyc[i]+vy_nc[idx+k];
                      vczc[i] = vczc[i]+vz_nc[idx+k];
	              //printf("%d %d %lf %lf %lf\n",step, i,vcx[step][i], vcy[step][i], vcz[step][i]);
             	/************************************************************************/
             	/* Calculation of the electrical current of the cation		        */
             	/************************************************************************/
	     	      jxc[step][i] = qc * vcxc[i];
		      jyc[step][i] = qc * vcyc[i];
		      jzc[step][i] = qc * vczc[i]; 
	         }
	     i++;
	     }  
        }
        fclose(in);

        /********************************************************************************/
        /* Memory deallocation to cation                                                */
        /********************************************************************************/
        free (vx_nc); 
        free (vy_nc);
        free (vz_nc);
        free (vcxc); 
        free (vcyc);
        free (vczc);


        /********************************************************************************/
        /* Calculation of the electrical current	                             	*/
        /********************************************************************************/
        for (step = 0; step < n_steps; step++){
 		for (i = 0; i < n_mol; i++){
        	     jx[step][i] = jxa[step][i]+jxc[step][i];
                     jy[step][i] = jya[step][i]+jyc[step][i];
                     jz[step][i] = jza[step][i]+jzc[step][i];
		}
	}

        /********************************************************************************/

        /********************************************************************************/
        /* Electrical Current Auto-Correlation Function calculation       		*/
        /********************************************************************************/	
        printf("Calculating JACF...\n");
        out = fopen(out_jacf_file, "w");
        for (k = 0; k < half_n_steps; k++){
             jcfx = 0.0;
             jcfy = 0.0;
             jcfz = 0.0;
             for (step = 0; step < half_n_steps; step++){
                  for (i = 0; i < n_mol; i++){
                       jcfx = jcfx+jx[step][i]*jx[step+k][i];
                       jcfy = jcfy+jy[step][i]*jy[step+k][i];
                       jcfz = jcfz+jz[step][i]*jz[step+k][i];
                  }
             }
             j_avg[k] = (jcfx+jcfy+jcfz)/(double)(half_n_steps)/3.0;
             fprintf(out,"%d %e %e\n", time[k]-time[0], j_avg[k], j_avg[k]/j_avg[0]);
        }
        fclose(out);
        /********************************************************************************/	


        /********************************************************************************/
        /* Numerical Integration of Electrical Current Auto-Correlation Function        */
        /* Green-Kubo equation for electrical conductivity coefficient                  */
        /********************************************************************************/
        printf("Integrating JACF...\n");
        out = fopen(cond_file, "w");
        for (k = 0; k < n_max; k++){
             cond = j_avg[0];
             for (step = 1; step < k-1; step++)
                  cond = cond+2.0*j_avg[step];
             cond = cons*0.5*(double)(time[step]-time[0])*(cond+j_avg[step])/(double)(k-1);
             if (step > 1)
                 fprintf(out, "%d %e \n", time[step]-time[0], cond);	
        }
        fclose(out);	
	end           = clock();
	cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;
	printf(" The calculation is ended...\n\n");	
	printf(" Elapsed time = %lf seconds...\n\n", cpu_time_used);		
        /********************************************************************************/
        }

        return 0;
	
}
