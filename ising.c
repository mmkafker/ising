/* To compile, use "cc hsmc.c -o sim -lm" */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


void get1Dcoord(long int i,long int j,int N,long int * out) ;

void get2Dcoord(long int k, int N, long int *i, long int*j);

void getNeighbors(long int p, long int *neigh, int N);

int main(int argc, char *argv[])
{
    int N = 1000;
    long int Ns = N*N;
    int out[1], outx[1], outy[1];
    double J = 1.0;
    double h = 0.01;
    double beta = 1.0;

    long int numtrials = 1000000000;

    long int neigh[4];

    char filename[100];

    //get1Dcoord(2,2,N,out);
    //printf("*out = %d\n",*out);

    //get2Dcoord(8,N,outx,outy);
    //printf("*outx = %d, *outy = %d\n",*outx,*outy);

    //getNeighbors(2,neigh,N);
    //printf("%d, %d, %d, %d\n",neigh[0],neigh[1],neigh[2],neigh[3]);

    FILE *flip_file; FILE *prob_file; FILE *spins_file;
    FILE *write_file;
    flip_file = fopen("flips.bin","rb");
    prob_file = fopen("probs.bin","rb");
    spins_file = fopen("spins.bin","rb");
    

    long int *flips; flips = malloc(numtrials*sizeof(long int));
    double *probs; probs = malloc(numtrials*sizeof(double));
    int *spins; spins = malloc(Ns*sizeof(int));

    fread(flips , sizeof(long int), numtrials , flip_file); 

    fread(probs , sizeof(double), numtrials , prob_file); 

    fread(spins , sizeof(int), Ns , spins_file);

    long int i;

    double deltaE,boltz,acceptprob;

    int neighspins[4];

    int p;

    for(i=0;i<numtrials;i++)
    {
        
	    
        if (i%500000==0) 
	{
	    printf("i = %ld\n",i);
	    sprintf(filename,"data/%ld.bin",i);
	    write_file = fopen(filename,"wb");
	    fwrite(spins,sizeof(int),Ns,write_file);
	    fclose(write_file);
	}
	    
	p = flips[i];	
	getNeighbors(p,neigh,N);
	neighspins[0] = spins[neigh[0]];
	neighspins[1] = spins[neigh[1]];
	neighspins[2] = spins[neigh[2]];
	neighspins[3] = spins[neigh[3]];
	
        //printf("%d, %d, %d, %d\n",neighspins[0],neighspins[1],neighspins[2],neighspins[3]);
	//2.0*h*spins[p]+spins[p]*J*(neighspins[0]+neighspins[1]+neighspins[2]+neighspins[3]);
        	
	deltaE = 2.0*h*spins[p]+spins[p]*J*(neighspins[0]+neighspins[1]+neighspins[2]+neighspins[3]);
        if(i%500000==0) printf("deltaE = %.17g\n",deltaE);
	//printf("%.17g\n",deltaE);
	if (deltaE<0)
	{
	    spins[p] *= -1;
	    continue;
	}
	else
	{
	    if (beta*deltaE > 15.0) acceptprob = 0.0;
	    else acceptprob = exp(-beta*deltaE);

	    if (probs[i]<acceptprob) spins[p] *= -1;
	}
    }

    //for(i=0;i<Ns;i++) printf("%d\n",spins[i]);

    free(flips);free(probs);free(spins);


    fclose(flip_file);fclose(prob_file);fclose(spins_file);

    return 0;
}

int arg_min(double* numbers, int length) {
    double min = numbers[0];
    int argmin = 0;
    for (int j = 1; j < length; j++) 
    {
	if (numbers[j] < min) 
	{
	    min = numbers[j];
	    argmin = j;
	}
    }

    return argmin;
}

void get1Dcoord(long int i, long int j, int N, long int * out)
{
    *out = j+i*N;
}


void get2Dcoord(long int k, int N, long int *i, long int*j)
{
    *i = k/N;
    *j = k - (*i)*N;
}


void getNeighbors(long int p, long int *neigh, int N)
{
    long int ix[1], iy[1], out[1];
    get2Dcoord(p,N,ix,iy);
    long int neighbors[8];
    int i;
    neighbors[0] = *ix-1; neighbors[1]= *iy;
    neighbors[2] = *ix+1;neighbors[3]= *iy;
    neighbors[4] = *ix; neighbors[5]= *iy-1;
    neighbors[6] = *ix; neighbors[7]= *iy+1;
    for(i=0;i<8;i++)
    {
        if (neighbors[i] < 0) neighbors[i] = N-1;
	if (neighbors[i] > N-1) neighbors[i] = 0;
	//if (i%2==1) get1Dcoord(neighbors[i-1],neighbors[i],N,out);
	//printf("%d\n",neighbors[i]);
	//if (i%2==1) printf("%d\n\n",*out);
    }
    get1Dcoord(neighbors[0],neighbors[1],N,out);
    neigh[0] = *out;


    get1Dcoord(neighbors[2],neighbors[3],N,out);
    neigh[1] = *out;


    get1Dcoord(neighbors[4],neighbors[5],N,out);
    neigh[2] = *out;

    get1Dcoord(neighbors[6],neighbors[7],N,out);
    neigh[3] = *out;
}
