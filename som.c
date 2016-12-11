
/*******************************************************************************
Submitter: Fida Mohammad Thoker
           Kunwar Abhinav adithya
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>

#define NUMPAT 1000     //change as per No of Pattersnss 
#define INP_DIM  2      //change as per Dimension of Input 
#define K_NEURONS 50         //change as per No of Rbf Neurons 
#define PARTNERS 4      //change as per No of Partner networks 
#define G_DIM  2      //change as per Dimension of Grid 
#define GAUSIAN_SIZE  3
#define INDIVUAL_PAT_SIZE  ((int)(K_NEURONS / PARTNERS))  // size of each som according to ( Partners and K_neurons)

#define rando() ((double)rand()/(RAND_MAX))




static double eucledian(double  *a,double  *b);
static void copy_input_to_centres(double  *a,double  *b);
void create_grid();
int find_winner_neuron(double *Input);
void AdjustWeights(double *centers, double *input,double LearningRate, double Influence);

void display_centres(double *c);
int find_Som_of_winner(int winner);

double Grid[K_NEURONS][G_DIM] = {0}; //  Create Array of G_Dimensional vectors
double MPATNERS[PARTNERS]= {0}; //  Create Array of G_Dimensional vectors

double centers[K_NEURONS][INP_DIM]= {0,};  //  Create Array of Input dimensional centeres

int main(int argc,char *argv[]) 
{
	int    i=0, j, k,z=0, p, np, op, ranpat[NUMPAT+1], epoch;
	int    NumPattern = NUMPAT;
	double Input[NUMPAT][INP_DIM] = {0};


	int seed = 1000;    /* to control the random number generator and randomize the input*/
	srand(seed);

	FILE* Error_file;
	Error_file = fopen("error_file", "w");
	FILE* output_file;
	output_file = fopen("output_file", "w");



	double sigma = GAUSIAN_SIZE;  // width of the gaussain

	double Error,  eta2 = 0.5;   // use eta2  output layer neurons as learning rates

	FILE* data;
	data = fopen(argv[1], "r");
	char *line = NULL;
	size_t size =0;
	char *string[20]; 
	char delimit[]=" \t";

	const int    Max_iterations  = 1000;

	//the value of the learning rate at the start of training
	const double Initial_LearningRate   = 0.1;
	const double Final_LearningRate   = 0.001;

	// Size  of the each som 
	const int const_size_som    = INDIVUAL_PAT_SIZE;

	double   MapRadius= 0;

	//used in the calculation of the neighbourhood width of influence
	double              TimeConstant= 0;

	//the number of training iterations
	int                 NumIterations= 0;

	//keeps track of what iteration each som has reached
	int                 IterationCount[PARTNERS];

	//the current size of the winning node's area of influence
	double              NeighbourhoodRadius [PARTNERS];

	//how much the learning rate is adjusted for nodes within
	//the area of influence
	double              Influence = 0;  // adjustment parameter

	double              LearningRate [PARTNERS];
	//this is the topological 'radius' each Som 
	MapRadius = const_size_som;

	//used in the calculation of the neighbourhood width of Influence
	TimeConstant = Max_iterations / log(MapRadius);



	for (i=0;i<NUMPAT;i++){  // Read the input and target output from the file 
		k=0;
		z=0;
		getline(&line,&size,data);

		string[k]=strtok(line,delimit);    
		while(string[k]!=NULL)                    
		{
			Input[i][k]=atof(string[k]);
			//		printf("Input=%f	",Input[i][k]);
			k++;
			string[k]=strtok(NULL,delimit);
		}
	}

	// create grid structure
	create_grid();


	// Initialize K centre neurons
	int c; 
	// select  centers from the input data (Input  subset selection)
	for (i=0;i<K_NEURONS;i++){

		c = random()%NUMPAT;
		copy_input_to_centres(centers[i],Input[c]);
		//printf(" centers %f	",centers[i][j]);
	}

	for (i=0;i<PARTNERS;i++){
		NeighbourhoodRadius[i]=0;
		LearningRate[i]=0.1;
		IterationCount[i] = 1;
	}

	int winner=0;
	double dist;
	int som =0;
	int start ;
	int end; 



	for (i=0;i<NumPattern;i++) {
		p=i;

			// find the winner neuron
			winner = find_winner_neuron(Input[p]);
			// find the respective som for the winner neuron
			som = find_Som_of_winner(winner);
			//start and end node of the winner som
			start = MPATNERS[som];
			end = MPATNERS[som+1];

		//calculate the width of the neighbourhood for this timestep for this SOM
		NeighbourhoodRadius[som] = MapRadius * exp(-(double)IterationCount[som]/TimeConstant);
		//printf("Neighbourhood Radius is %f \n",NeighbourhoodRadius);


		// Now adjust the weights  of the winner som
		for (k =start;k <end;k++){
			dist = eucledian(Grid[winner],Grid[k]);
				if( dist < NeighbourhoodRadius[som])
				{

					//calculate by how much its weights are adjusted
					sigma = NeighbourhoodRadius[som];
					Influence = exp(-(dist *dist) / (2*sigma*sigma));
	//	printf("Neighbourhood influence Rate is %f \n",Influence);

					AdjustWeights(centers[k],Input[p], LearningRate[som], Influence);
				}
		}
		//reduce the learning rate
		LearningRate[som] = Initial_LearningRate * exp(-(double)IterationCount[som]/Max_iterations );
		//printf("Learning Rate is %d \n",LearningRate);

		IterationCount[som] =IterationCount[som]+1;


	}
//}

for (i =0; i<K_NEURONS;i++){

	printf("Centre %d = {",i);
	display_centres(centers[i]);

}


return 0;
}

void display_centres(double *c)
{
int i =0;
for (i=0;i<INP_DIM;i++)
	printf("%f,  ",c[i]);
	printf("}\n");
}


static void copy_input_to_centres(double  *a,double  *b)
{

	int i;
	for (i =0;i<INP_DIM;i++) {
		 a[i] = b[i];
	}
	return ;
}

static void copy_indicies(double  *a,double  *b)
{
	int i;
	for (i =0;i<G_DIM;i++) {
		 a[i] = b[i];
	}
	return ;
}

/***************************************************
* Function to  implement the euclidean distance 
  between input and the center
****************************************************/
static double eucledian(double  *a,double  *b)
{
	double distance=0;
	int i;
	for (i =0;i<INP_DIM;i++) {
		distance +=  (a[i]-b[i]) *(a[i]-b[i]);

	}
	return sqrt(distance);
}

/***************************************************
* Function to  find the winner neuron based on 
  euclidean distance between input and the centers 
****************************************************/

 int find_winner_neuron(double *Input)
{
	double distance =0;
	int winner =0,i;
	double min =0;

	min = eucledian(Input,centers[0]);
	for (i=1;i<K_NEURONS;i++){
		distance = eucledian(Input,centers[i]);
		if(min > distance){
			min =distance;
			winner = i;
		}
	}
	return winner;
}

/***************************************************
* Function to apply adjustment parameter to the neighbour
  neuron of the winner som 
****************************************************/

void AdjustWeights(double *centers, double *input,double LearningRate, double Influence) 
{ 
	double delta =0;
	int i;
	for (i =0;i<INP_DIM;i++) {
		 delta = (LearningRate * (input[i] - centers[i])* Influence);
		 centers[i] = centers[i] +delta;
	}
	return;

}



/***************************************************
* Function to  Create the structure lattic of the 
  Grid represting multiple soms
****************************************************/


void create_grid()
{
	double temp_array[5];
	int count=0;
	int i,j,k,g,z,l;
	int no_of_indiv_k_neurons;
	switch(G_DIM)  // Create K grid indexs as per the grid dimension(Lattice structure) and copy them into GRID vectors
// Array of vector Grid[K_NEURONS] represents the whole grid which contains all the Partner soms
	{
		case 1:
			for (i =0,g=0;i<K_NEURONS;i++){
				temp_array[0]=j;
				copy_indicies(Grid[g],temp_array);
				count++;
				if(count == K_NEURONS-1)
					break;
			}
			break;
		case 2:
			for (i =0,g=0;i<K_NEURONS;i++){
				for(j=0;j<2;j++){
					temp_array[0]=i;
					temp_array[1]=j;
					copy_indicies(Grid[g],temp_array);
					g++;
					count++;
				}
				if(count == K_NEURONS-1)
					break;
			}
			break;
		case 3:
			for (i =0,g=0;i<K_NEURONS;i++){
				for(j=0;j<3;j++){
					for(k=0;k<3;k++){
						temp_array[0]=i;
						temp_array[1]=j;
						temp_array[2]=k;
						copy_indicies(Grid[g],temp_array);
						g++;
						count++;
						if(count == K_NEURONS-1)
							break;
					}
					if(count == K_NEURONS-1)
						break;
				}
				if(count == K_NEURONS-1)
					break;
			}
			break;
		case 4:
			for (i =0,g=0;i<K_NEURONS;i++){
				for(j=0;j<4;j++){
					for(k=0;k<4;k++){
						for(z=0;z<4;z++){
							temp_array[0]=i;
							temp_array[1]=j;
							temp_array[2]=k;
							temp_array[3]=z;
							copy_indicies(Grid[g],temp_array);
							g++;
							count++;
							if(count == K_NEURONS-1)
								break;
						}
						if(count == K_NEURONS-1)
							break;
					}
					if(count == K_NEURONS-1)
						break;
				}
				if(count == K_NEURONS-1)
					break;
			}
			break;
		case 5:
			for (i =0,g=0;i<K_NEURONS;i++){
				for(j=0;j<5;j++){
					for(k=0;k<5;k++){
						for(z=0;z<5;z++){
							for(l=0;l<5;l++){
								temp_array[0]=i;
								temp_array[1]=j;
								temp_array[2]=k;
								temp_array[3]=z;
								temp_array[4]=z;
								copy_indicies(Grid[g],temp_array);
								g++;
								count++;
								if(count == K_NEURONS-1)
									break;
							}
							if(count == K_NEURONS-1)
								break;
						}
						if(count == K_NEURONS-1)
							break;
					}
					if(count == K_NEURONS-1)
						if(count == K_NEURONS-1)
							break;
				}
				break;
			}
			break;

	}
	// calculate the number of k neurons in each partner form the totak K neurons 
	for (i =0,j=0;i<K_NEURONS;i =i+INDIVUAL_PAT_SIZE,j++)
		MPATNERS[j] = i;    // Each MPATNERS element contains the starting index of partner som in the FULL single grid
		

	return ;
}

/***************************************************
* Function to find the Som of the winner neuron
*
****************************************************/
int find_Som_of_winner(int winner)
{

	int i;
	for (i=0;i<PARTNERS-1;i++){
		if(winner>=MPATNERS[i] && winner <MPATNERS[i+1])
			return i;
	}

}
	
