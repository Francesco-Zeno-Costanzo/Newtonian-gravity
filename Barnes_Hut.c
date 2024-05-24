#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define THETA 0.5      // Parameter of the algorithm

int num_body, steps;   // Number of body, number of evolutions' steps
double  G;             // Gravitational constant
double softening, dt;  // Softening parameter for potential, time step

/*Structure for a body of the system*/
typedef struct{
    double x, y;    // Posistion
    double mass;    // mass of body
    double vx, vy;  // velocity
    double fx, fy;  // force on the body
} Body;

/*Structure for a Node of the quadtree*/
typedef struct Node{
    double x, y;                     // Position of ceter of mass 
    double mass;                     // Total mass of the node
    double xmin, ymin;               // Limit of the node
    double xmax, ymax;               // Limit of the node
    struct Node *NW, *NE, *SW, *SE;  // Pointer to sub node
    Body *body;                      // Body
} Node;

Node* create_node(double xmin, double ymin, double xmax, double ymax){
    /*Function for creation and initializzation of the node*/
    
    /*Create Node*/
    Node *node = (Node *)malloc(sizeof(Node));

    /*Initialize Node*/
    node->x = 0.0;
    node->y = 0.0;
    node->mass = 0.0;
    node->xmin = xmin;
    node->ymin = ymin;
    node->xmax = xmax;
    node->ymax = ymax;

    /*The node is initially empty and without subnode*/
    node->NW = node->NE = node->SW = node->SE = NULL;
    node->body = NULL;
    return node;
}

void insert_body(Node *node, Body *body){
    /*Function that recursivelly put all body in node*/

    if (node->body == NULL && node->NW == NULL){
        /*If the node is empty and not have subnodes, put here the body*/
        node->body = body;
        node->x    = body->x;
        node->y    = body->y;
        node->mass = body->mass;
    }
    else{
        if (node->NW == NULL){
            /*If the node has a body inside we subdivide the node*/
            /*Compute mid point*/
            double midx = (node->xmin + node->xmax) / 2.0;
            double midy = (node->ymin + node->ymax) / 2.0;
            /*create subnode*/
            node->NW = create_node(node->xmin, node->ymin, midx, midy);
            node->NE = create_node(midx, node->ymin, node->xmax, midy);
            node->SW = create_node(node->xmin, midy, midx, node->ymax);
            node->SE = create_node(midx, midy, node->xmax, node->ymax);
            
            /*Reinsert the existing body in the new subdivision of the space*/
            insert_body(node, node->body);
            /*Remove the body beacuse now is in the subnode*/
            node->body = NULL;
        }
        
        /* Update the center of mass and total mass of the node*/
        node->mass += body->mass;
        node->x = (node->x * (node->mass - body->mass) + body->x * body->mass) / node->mass;
        node->y = (node->y * (node->mass - body->mass) + body->y * body->mass) / node->mass;
        
        /* Insert the new body into the appropriate quadrant*/
        if (body->x < (node->xmin + node->xmax) / 2.0){
            if (body->y < (node->ymin + node->ymax) / 2.0){
                insert_body(node->NW, body);
            } 
            else{
                insert_body(node->SW, body);
            }
        }
        else{
            if (body->y < (node->ymin + node->ymax) / 2.0){
                insert_body(node->NE, body);
            }
            else{
                insert_body(node->SE, body);
            }
        }
    }
}

void compute_force(Body *body, Node *node){
    /*No body no force*/
    if (node == NULL || node->mass == 0){
        return;
    }
    
    /*Compute distance*/
    double force;
    double dx = node->x - body->x;
    double dy = node->y - body->y;
    double d  = sqrt(dx * dx + dy * dy);
    
    /*Avoid autointeraction*/
    if (node->body == body){
        return;
    }
    
    /*Compute interaction with the certer of mass of the node*/
    if ((node->xmax - node->xmin) / d < THETA || node->body != NULL){
        force = G * body->mass * node->mass / (d*d*d + softening);
        body->fx += force * dx;
        body->fy += force * dy;
    }
    else{
        /*Or compute force in the subnode*/
        compute_force(body, node->NW);
        compute_force(body, node->NE);
        compute_force(body, node->SW);
        compute_force(body, node->SE);
    }
}

void compute_acc(Node *root, Body bodies[]){
    for (int i = 0; i < num_body; i++) {
        bodies[i].fx = 0;
        bodies[i].fy = 0;
        compute_force(&bodies[i], root);
    }
}

void bounding(Body bodies[], double *min_x, double *min_y, double *max_x, double *max_y){
    /********************************************************************
    Function to compute the maximum extention of the tree.
    The particles are not bounded in a box so we must copute the maximum
    and the minimum so we can copute the quadtree.
    ********************************************************************/
    *min_x = bodies[0].x;
    *min_y = bodies[0].y;
    *max_x = bodies[0].x;
    *max_y = bodies[0].y;

    for (int i = 1; i < num_body; i++){
        if (bodies[i].x < *min_x) *min_x = bodies[i].x;
        if (bodies[i].y < *min_y) *min_y = bodies[i].y;
        if (bodies[i].x > *max_x) *max_x = bodies[i].x;
        if (bodies[i].y > *max_y) *max_y = bodies[i].y;
    }
}

Node* Quadtree(Body bodies[], double min_x, double min_y, double max_x, double max_y){
    /*Function to compute the quadtree*/
    Node *root;
    root = create_node(min_x, min_y, max_x - min_x, max_y - min_y);
    for (int i = 0; i < num_body; i++){
        insert_body(root, &bodies[i]);
    }
    return root;
}

void freeTree(Node *node) {
    /*********************************
    Function to delete the tree.
    The tree must be compute at each
    update of positions
    *********************************/
    if (node == NULL) {
        return;
    }
    freeTree(node->NW);
    freeTree(node->NE);
    freeTree(node->SW);
    freeTree(node->SE);
    free(node);
}

void update(Body bodies[]){
    /*Function to update the state of the system of a step od size dt*/
    double l, w0, w1, c1, c2, c3, c4, d1, d2, d3;
    double min_x, min_y, max_x, max_y;
    Node *root;

    /*Some funny coefficents*/
    l  = pow(2.0, (1.0/3.0));
    w0 = -l/(2.0-l);
    w1 = 1.0/(2.0-l);
    /*Other funny coefficents*/
    c1 = c4 = w1/2.0;
    c2 = c3 = (w0 + w1)/2.0;
    d1 = d3 = w1;
    d2 = w0;

    /*First step for poition*/
    for (int i = 0; i < num_body; i++){
        bodies[i].x += bodies[i].vx * c1 * dt;
        bodies[i].y += bodies[i].vy * c1 * dt;
    }

    /*Compute the bound of the system according the configuration of the particles*/
    bounding(bodies, &min_x, &min_y, &max_x, &max_y);
    /*Compute accelerations*/
    root = Quadtree(bodies, min_x, min_y, max_x, max_y);
    compute_acc(root, bodies);
    freeTree(root);

    // Primo aggiornamento delle velocità e secondo aggiornamento delle posizioni
    for (int i = 0; i < num_body; i++){
        bodies[i].vx += bodies[i].fx / bodies[i].mass * d1 * dt;
        bodies[i].vy += bodies[i].fy / bodies[i].mass * d1 * dt;
        bodies[i].x  += bodies[i].vx * c2 * dt;
        bodies[i].y  += bodies[i].vy * c2 * dt;
    }

    bounding(bodies, &min_x, &min_y, &max_x, &max_y);
    root = Quadtree(bodies, min_x, min_y, max_x, max_y);
    compute_acc(root, bodies);
    freeTree(root);

    for (int i = 0; i < num_body; i++){
        bodies[i].vx += bodies[i].fx / bodies[i].mass * d2 * dt;
        bodies[i].vy += bodies[i].fy / bodies[i].mass * d2 * dt;
        bodies[i].x  += bodies[i].vx * c3 * dt;
        bodies[i].y  += bodies[i].vy * c3 * dt;
    }

    bounding(bodies, &min_x, &min_y, &max_x, &max_y);
    root = Quadtree(bodies, min_x, min_y, max_x, max_y);
    compute_acc(root, bodies);
    freeTree(root);


    for (int i = 0; i < num_body; i++){
        bodies[i].vx += bodies[i].fx / bodies[i].mass * d3 * dt;
        bodies[i].vy += bodies[i].fy / bodies[i].mass * d3 * dt;
        bodies[i].x  += bodies[i].vx * c4 * dt;
        bodies[i].y  += bodies[i].vy * c4 * dt;
    }
}

void update_e(Body bodies[]){
    
    double min_x, min_y, max_x, max_y;
    Node *root;

    /*Compute the bound of the system according the configuration of the particles*/
    bounding(bodies, &min_x, &min_y, &max_x, &max_y);
    /*Compute accelerations*/
    root = Quadtree(bodies, min_x, min_y, max_x, max_y);
    compute_acc(root, bodies);

    /*First step for poition*/
    for (int i = 0; i < num_body; i++){
        bodies[i].vx += bodies[i].fx / bodies[i].mass * dt;
        bodies[i].vy += bodies[i].fy / bodies[i].mass * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
    freeTree(root);
}


void init(const char* filename, Body** bodies){
    FILE *file = fopen(filename, "r");
    if (!file){
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Leggi i parametri globali
    if(fscanf(file, "%lf %lf %lf %d %d", &G, &softening, &dt, &num_body, &steps));

    // Alloca memoria per i corpi
    *bodies = (Body *)malloc(num_body * sizeof(Body));
    if (!*bodies){
        perror("Failed to allocate memory for bodies");
        exit(EXIT_FAILURE);
    }

    // Leggi le proprietà dei corpi
    for (int i = 0; i < num_body; i++){
        if(fscanf(file, "%lf %lf %lf %lf %lf", &(*bodies)[i].mass, &(*bodies)[i].x, &(*bodies)[i].y, &(*bodies)[i].vx, &(*bodies)[i].vy));
        (*bodies)[i].fx = 0.0;
        (*bodies)[i].fy = 0.0;
    }

    fclose(file);
}

double energy(Body bodies[]){
    /* Compute the total energy */
    double K, V;
    double dx, dy, d;
    K = V = 0;
 
    /*Kinetic term*/
    for(int i=0; i < num_body; i++){
        K += 0.5*bodies[i].mass*(bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy);
    }

    /*Potential term*/
    for(int i=0; i<num_body; i++){
		for(int j=0; j<num_body; j++){
            if(i != j){
                dx = bodies[i].x - bodies[j].x;
                dy = bodies[i].y - bodies[j].y;
                d  = sqrt(dx * dx + dy * dy + softening/sqrt(dx * dx + dy * dy) );
                V += -bodies[i].mass*bodies[j].mass*G/d;
            }
        }
    }
    /* V/2 is to not overcount interactions*/
    return K + V/2.0;
}

double angular(Body bodies[]){
    double L;
    L = 0;
    for(int i=0; i < num_body; i++){
        L += bodies[i].mass*(bodies[i].x * bodies[i].vy - bodies[i].y * bodies[i].vx);
    }
    return L;
}

int main(int argC, char* argV[]){
    clock_t start = clock(); // Misura il tempo di esecuzione
    
    FILE* fpd = fopen(argV[2],"w");
    FILE* fpq = fopen(argV[3],"w");

    Body *bodies;

	if(argC!=4){
		printf("Usage: %s <file containing initial configuration> <file to store output> <file for E & L> \n",argV[0]);
    }
    else{
        /*Initialize the system*/
		init(argV[1], &bodies);

        /*Save initail configuration*/
        for(int i=0; i<num_body; i++){
			fprintf(fpd, "%.16f\t%.16f\t%.16f\t%.16f\n", bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
        }
        fprintf(fpq, "%.16f\t%.16f\n", energy(bodies), angular(bodies));

        /*Time evolution*/
    	for(int j=0; j<steps; j++){

			update(bodies); //Update sistem

            /*Save configuration*/
            for (int i = 0; i < num_body; i++) {
                fprintf(fpd, "%.16f\t%.16f\t%.16f\t%.16f\n", bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
            }
            fprintf(fpq, "%.16f\t%.16f\n", energy(bodies), angular(bodies));
		}
	}

    /*Close file*/
	fclose(fpd);
    fclose(fpq);

    clock_t end = clock();
    printf("Elapsed time =  %f secondi \n", ((double)(end - start)) / CLOCKS_PER_SEC);

	return 0;
}