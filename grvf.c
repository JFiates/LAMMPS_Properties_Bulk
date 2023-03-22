/**********************************************************************************/
/* School of Chemical Engineering, University of Campinas, Brazil                 */
/* Code to calculate the radial distribution function and the coordinate number	  */
/* at the bulk system based on atoms type					  */
/* Code to unwrapped coordinates						  */
/* Developer: Dr. Juliane Fiates						  */
/**********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>



int main(){
        /*************************************************************************/
        /* Variables declaration                                                 */
        /*************************************************************************/
	/* Counters								 */
        int i,j,k,bin,step,idx;
	/* Molecules per frame							 */ 
	int af,type;
        /* Box size fom dump file 						 */
        float bxl,bxu;
        float byl,byu;
        float bzl,bzu;
	/* Time 								 */
	int time;
        /* Total number of steps                                                 */
        int n_steps;
        /* Total number of molecules                                             */ 
        int n_molt;
	/* Atom types used to compute the distribution function			 */
	int tp_1,tp_2; 
	/* Number of atoms of each type 					 */
	int cont1,cont2;
        /* Coordinates 			                                         */ 
        float x,y,z;
	/* Width of the bin							 */
	float delr;
	/* The lower value of radius						 */
	float rlower; 
	/* The upper value of radius						 */
	float rupper;
	/* Number of ideal gas particles at same density			 */
	float nideal;  
	/* Number of bins							 */
	int maxbin;
	/* Box size 	 	 						 */
	float boxl; 
	/* Volume	 	 						 */
	float vol;
	/* Distance between the atoms						 */
 	float rxij,ryij,rzij,rijsq,rij;
	/* Coordination number 							 */
 	float n;
        /* Input file                                                            */ 
        char inpfile[50];
        /* Auxiliary string for comments in the input file                       */ 
        char com[500];
        /* Input file with the coordinates                                       */ 
        char cor_file[50];
        /* Output file with the RDF				                 */ 
        char g_file[50];
        /* Output file with the coordination number		                 */ 
        char n_file[50];


        FILE *in, *out;

        /*************************************************************************/
        /* Reading input file                                                    */
        /*************************************************************************/
        scanf("%s", inpfile);
        in = fopen(inpfile, "r");
        fscanf(in, "%s", com);
        fscanf(in, "%d", &n_steps); 
	//printf("%d\n", n_steps);
        fscanf(in, "%s", com);	
        fscanf(in, "%d", &n_molt);
        fscanf(in, "%s", com);
        fscanf(in, "%d", &tp_1);
        fscanf(in, "%s", com);	
        fscanf(in, "%d", &tp_2);
        fscanf(in, "%s", com);		
        fscanf(in, "%d", &maxbin);
        fscanf(in, "%s", com);	
        fscanf(in, "%f", &vol);
	//printf("%f\n", vol);
        fscanf(in, "%s", com);
        fscanf(in, "%s", cor_file);	
        fscanf(in, "%s", com);
        fscanf(in, "%s", g_file); 
        fscanf(in, "%s", com);
        fscanf(in, "%s", n_file); 
        fclose(in);

        /*************************************************************************/
        /* Memory allocation                                                     */
        /*************************************************************************/
	int c = n_molt;
        int a = maxbin;
        /* Coordinates of atom 1 and 2	                                         */ 
        float *rx_i = calloc (c, sizeof (float));
        float *ry_i = calloc (c, sizeof (float));
        float *rz_i = calloc (c, sizeof (float));
        float *rx_j = calloc (c, sizeof (float));
        float *ry_j = calloc (c, sizeof (float));
        float *rz_j = calloc (c, sizeof (float));
	/* Histogram and RDF							 */
	int   *hist = calloc (a, sizeof (int));
	float *gr   = calloc (a, sizeof (float));
	float *gr_n = calloc (a, sizeof (float));
	float *r    = calloc (a, sizeof (float));

        /*************************************************************************/
        /* Reading the coordinates input file	                                 */
        /*************************************************************************/
        printf("Reading input file...\n");
        in  = fopen(cor_file, "r");	
        for (step = 0; step < n_steps; step++){	
             fscanf(in, "%d", &time);
	     fscanf(in, "%d", &af);
             fscanf(in, "%e %e", &bxl, &bxu);
	     fscanf(in, "%e %e", &byl, &byu);
	     fscanf(in, "%e %e", &bzl, &bzu);
	     cont1 = 0;
	     cont2 = 0;
             for (i = 0; i < n_molt; i++){
                  fscanf(in, "%d %d %f %f %f", &idx, &type, &x, &y, &z);
		  //printf("%d %d %f %f %f\n", idx,type, x, y, z); 
		     if(type == tp_1){
		        rx_i[cont1] = x;
	     	        ry_i[cont1] = y;
	     	        rz_i[cont1] = z;
	     	     //printf(" %d %d  %f %f %f\n", cont1, type, rx_i[cont1], ry_i[cont1], rz_i[cont1]); 
		        cont1++;
		     }
		     if(type == tp_2){
		        rx_j[cont2] = x;
	     	        ry_j[cont2] = y;
	     	        rz_j[cont2] = z;
	     	     //printf(" %d %d %f %f %f\n", cont2, type, rx_j[cont2], ry_j[cont2], rz_j[cont2]);
		        cont2++; 
		     }

	     }
		//printf("%d %d\n", cont1,cont2);
        /*************************************************************************/
        /* Calculation of the histogram 	                                 */
        /*************************************************************************/
	boxl = pow(vol,(1.0/3.0));
	delr = (0.5*boxl)/(float)(maxbin);
             for (i = 0; i < cont1; i++){
		  for (j = 0; j < cont2; j++){
		       rxij  = rx_i[i]-rx_j[j];
		       ryij  = ry_i[i]-ry_j[j];
		       rzij  = rz_i[i]-rz_j[j];
		       //printf( "%d %f %f %f\n", j, rx_j[j], ry_j[j], rz_j[j]);
		       rxij  = rxij-boxl*round(rxij/boxl);
		       ryij  = ryij-boxl*round(ryij/boxl);
		       rzij  = rzij-boxl*round(rzij/boxl);
		       //printf( "%d %d %f %f %f\n", i, j, rxij, ryij, rzij);
		       rijsq = rxij*rxij+ryij*ryij+rzij*rzij;
		       rij   = pow(rijsq,(1.0/2.0));
		       //printf( "%f \n", rij);
		       bin   = (int) (rij/delr)+1;
		       if (bin <= maxbin)
		           hist[bin]=hist[bin]+1;
		  }
	     }
        }
        fclose (in);

	/*************************************************************************/
        /* Calculation of the Radial Distribution Function (RDF)                 */
        /*************************************************************************/
        printf("Calculating the RDF...\n");
	out = fopen(g_file, "w");
	float cons_g = (4.0*M_PI*((float)(cont1*cont2)/vol))/3.0;
	float cons_n = ((float)(cont2)/vol)*4.0*M_PI;
	//printf("%f", cons);
	for(bin = 0; bin < maxbin; bin++){
	    rlower    = (float)(bin)*delr;
	    rupper    = rlower + delr;
	    nideal    = cons_g * ((rupper*rupper*rupper)-(rlower*rlower*rlower));
	    gr[bin]   = (float)(hist[bin])/(float)(n_steps)/nideal;
	    gr_n[bin] = gr[bin]*(rlower*rlower);
	    //printf("%e %e \n", rlower, gr[bin]);		
	    fprintf(out,"%e %e \n", rlower, gr[bin]);
	}
	fclose (out);


	/*************************************************************************/
        /* Calculation of the coordination number (n)             		 */
        /*************************************************************************/
        printf("Integration of the RDF...\n");
	out = fopen(n_file, "w");
        for (bin = 0; bin < maxbin; bin++){
	     r[bin] = (float)(bin)*delr;
             n      = gr_n[0];
             for (k = 1; k < bin-1; k++)
                  n = n+2.0*gr_n[k];
             n = cons_n*0.5*(r[k]-r[0])*(n+gr_n[k])/(float)(bin-1);
             if (k > 1)
                 fprintf(out, "%e %e \n", r[k]-r[0], n);	
	}
	fclose (out);

	return 0;

}

