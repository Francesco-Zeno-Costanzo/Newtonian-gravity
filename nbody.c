#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/*Structure for a 2D vector*/
typedef struct{
    double x, y;
}vector;

/************************************************/
/***************Global variables*****************/
/************************************************/

int num_body, steps;     // Number of body, number of evolutions' steps
double *masses, G;       // Pointer for all masses, gravitational constant
double softening, dt;    // Softening parameter for potential, time step
vector *pos, *vel, *acc; // Pointer for positions, velocities and accelerations

/************************************************/
/***************Auxiliar functions***************/
/************************************************/

vector add(vector x1, vector x2){
    vector z = {x1.x + x2.x, x1.y + x2.y};
    return z;
}

vector sub(vector x1, vector x2){
    vector z = {x1.x - x2.x, x1.y - x2.y};
    return z;
}

vector prod(double s, vector x1){
    vector z = {x1.x * s, x1.y * s};
    return z;
}

double v_mod(vector x1){
    /*module of a vector plus softening to avoid divergence*/
    return sqrt(x1.x * x1.x + x1.y * x1.y + softening);
}

/************************************************/
/***************Initialize system****************/
/************************************************/

void init(char* fileName){
    /*Open file*/
	FILE* fp = fopen(fileName,"r");
	
    /* Read computational parameters*/
	if(fscanf(fp,"%lf%lf%lf%d%d",&G, &softening, &dt, &num_body, &steps));

    /*Define all vectors for evolution*/
	masses = (double*)malloc(num_body*sizeof(double));
	pos    = (vector*)malloc(num_body*sizeof(vector));
	vel    = (vector*)malloc(num_body*sizeof(vector));
	acc    = (vector*)malloc(num_body*sizeof(vector));
	
    /*Set initial condition*/
	for(int i=0; i<num_body; i++){
		if(fscanf(fp,"%lf%lf%lf%lf%lf",&masses[i], &pos[i].x, &pos[i].y, &vel[i].x, &vel[i].y));
	}
	/*Close file*/
	fclose(fp);
}

/******************************************************/
/***************Function for simulation****************/
/******************************************************/

void compute_acc(){

	/*Lopp over all bodies*/
	for(int i=0; i<num_body; i++){

		acc[i].x = 0;
		acc[i].y = 0;

		for(int j=0; j<num_body; j++){
            /*loop again and avoid autointeraction*/
			if(i != j){
				acc[i] = add(acc[i], prod(G*masses[j]/pow(v_mod(sub(pos[j], pos[i])), 3), sub(pos[j],pos[i])));
			}
		}
	}
}

void update(){
    /* Function to update the state of the system of a step od size dt
    */
    double l, w0, w1, c1, c2, c3, c4, d1, d2, d3;

    /*Some funny coefficents*/
    l  = pow(2.0, (1.0/3.0));
    w0 = -l/(2.0-l);
    w1 = 1.0/(2.0-l);
    /*Other funny coefficents*/
    c1 = c4 = w1/2.0;
    c2 = c3 = (w0 + w1)/2.0;
    d1 = d3 = w1;
    d2 = w0;

    /*firstly update positions*/
    for(int i=0; i<num_body; i++){
        pos[i] = add(pos[i], prod(c1*dt, vel[i]));
    }
    /*
    then compute the accelerations this
    calculation must be out of for because
    there is already a for inside compute_acc()
    */
    compute_acc();
    /*
    finally update velocities with the accelerations
    and position with the new velocities. And repeat...
    */
    for(int i=0; i<num_body; i++){
        vel[i] = add(vel[i], prod(d1*dt, acc[i]));
        pos[i] = add(pos[i], prod(c2*dt, vel[i]));
    }
    compute_acc();
    for(int i=0; i<num_body; i++){
        vel[i] = add(vel[i], prod(d2*dt, acc[i]));
        pos[i] = add(pos[i], prod(c3*dt, vel[i]));
    }
    compute_acc();
    for(int i=0; i<num_body; i++){
        vel[i] = add(vel[i], prod(d3*dt, acc[i]));
        pos[i] = add(pos[i], prod(c4*dt, vel[i]));
    }
}

/******************************************************/
/*****************Physical quantities******************/
/******************************************************/

double energy(){
    /* Compute the total energy */
    double K, V;
    K = V = 0;
 
    /*Kinetic term*/
    for(int i=0; i < num_body; i++){
        K += 0.5*masses[i]*(vel[i].x * vel[i].x + vel[i].y * vel[i].y);
    }

    /*Potential term*/
    for(int i=0; i<num_body; i++){
		for(int j=0; j<num_body; j++){
            if(i != j){
                V += -masses[i]*masses[j]*G/v_mod(sub(pos[j], pos[i]));
            }
        }
    }
    /* V/2 is to not overcount interactions*/
    return K + V/2.0;
}

double angular(){
    double L;
    L = 0;
    for(int i=0; i < num_body; i++){
        L += masses[i]*(pos[i].x * vel[i].y - pos[i].y * vel[i].x);
    }
    return L;
}

/******************************************************/
/******************Main of the code********************/
/******************************************************/

int main(int argC, char* argV[]){
    clock_t start = clock(); // Misura il tempo di esecuzione
    
    FILE* fpd = fopen(argV[2],"w");
    FILE* fpq = fopen(argV[3],"w");

	if(argC!=4){
		printf("Usage: %s <file containing initial configuration> <file to store output> <file for E & L> \n",argV[0]);
    }
    else{
        /*Initialize the system*/
		init(argV[1]);

        /*Save initail configuration*/
        for(int j=0; j<num_body; j++){
			fprintf(fpd, "%.16f\t%.16f\t%.16f\t%.16f\n",pos[j].x,pos[j].y,vel[j].x,vel[j].y);
        }
        fprintf(fpq, "%.16f\t%.16f\n", energy(), angular());

        /*Time evolution*/
    	for(int i=0; i<steps; i++){

			update(); //Update sistem

            /*Save configuration*/
			for(int j=0; j<num_body; j++){
				fprintf(fpd, "%.16f\t%.16f\t%.16f\t%.16f\n",pos[j].x,pos[j].y,vel[j].x,vel[j].y);
            }
            fprintf(fpq, "%.16f\t%.16f\n", energy(), angular());
		}
	}

    /*Close file*/
	fclose(fpd);
    fclose(fpq);

    clock_t end = clock();
    printf("Elapsed time =  %f secondi \n", ((double)(end - start)) / CLOCKS_PER_SEC);

	return 0;
}