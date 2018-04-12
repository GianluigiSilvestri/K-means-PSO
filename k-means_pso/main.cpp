//
//  main.cpp
//  k-means_pso
//
//  Created by Gianluigi Silvestri on 05/12/16.
//  Copyright © 2016 Università degli Studi di Parma. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <random>
using  namespace  std;

#define	D_max 100  // Max number of dimensions of the search space
#define	S_max 3000 // Max swarm size

// Structures

struct velocity
{
    double v[D_max];
};

struct position
{
    double x[D_max];
    double f;
    int seed=-1;
};

struct cluster
{
    int size=0;
    double x[D_max];
    double sigma=0;
    struct position bp; // best particle
    bool best=false;
};

// Sub-programs
void identify_niches();
void k_means();
double alea( double a, double b );
double perf( int s, struct position H[S_max]); // Fitness evaluation
void update();

// Global variables
int best;
double pi; // Useful for some test functions
struct position P[S_max]; // Best positions found by each particle
int D; // Search space dimension
int S; // Swarm size
int K; // Number of seeds
double v[D_max]; // vector for distance
struct velocity V[S_max]; // Velocities
struct position X[S_max]; // Positions
struct cluster M[S_max]; // Definitive Seeds
struct cluster MT[S_max]; // Temporary Seeds
double xmin[D_max], xmax[D_max]; // Intervals defining the search space
int T;
double x;
double c; // acceleration
double w; // constriction factor


int main(int argc, const char * argv[]) {
    
    int d; // Current dimension
    int s; // Rank of the current particle
    int k; // current seed
    int c1; // intervals for identify niches
    
    ofstream myfile("example.txt");
    ofstream myfile1 ("example1.txt");
    ofstream seeds("seeds.txt");
    ofstream seeds1("seeds1.txt");
    
    pi = acos( -1 ); // for rastrigin function
    if (argc < 7)
    { // We expect 7 arguments: the program name, the source path and the destination path
        cerr << "Usage: D, S, K, x, T, c1" << endl;
        return 1;
    }
    else
    {
        D =atoi(argv[1]); // Search space dimension
        S=atoi(argv[2]);
        K=atoi(argv[3]);
        x=atof(argv[4]);
        T=atoi(argv[5]);
        c1=atoi(argv[6]);
    }
    
    w = 0.73; // defined in thesis passaro
    c=2.05;
    // D-cube data
    for ( d = 0; d < D; d++ )
    {
        xmin[d] = -x; xmax[d] = x;
    }
    
    //-----------------------INITIALIZATION
    srand( static_cast<unsigned int>(time(NULL)));
    
    for ( s = 0; s < S; s++ ) // create S particles
    {
        for ( d = 0; d < D; d++ )
        {
            X[s].x[d] = alea( xmin[d], xmax[d] );
            V[s].v[d] = (alea( xmin[d], xmax[d] ) - X[s].x[d])/2; // Non uniform
            // V[s].v[d] = ( xmin[d]-xmax[d] )*(0.5-alea(0,1)); //Uniform. 2006-02-24
        }
        X[s].f = perf( s,X);// calculate fitness
        P[s] = X[s]; // Best position = current one
    }
    
    for(s=0; s<S; s++)
    {
            if (myfile.is_open())
            {
                for(d=0;d<D;d++)
                {
                    myfile<<P[s].x[d];
                    myfile<< " ";
                }
                myfile<<P[s].f;
                myfile<<"\n";
            }
        
    }
    
    
    k_means(); //k-means algorithm
    

    //--------------------ITERATIONS
    for (int t=1; t<T; t++)
    {
        
        update();
        if (t % c1 ==0)
        {
            k_means();
            identify_niches();
        }
        
    }
    
    /////-----------------PRINTING RESULTS
    
    for(s=0; s<S; s++)
    {
        if (X[s].seed!=-1)
        {
            if (myfile1.is_open())
            {
                for(d=0;d<D;d++)
                {
                    myfile1<<P[s].x[d];
                    myfile1<< " ";
                }
                myfile1<<P[s].f;
                myfile1<<" ";
                myfile1<<X[s].seed*5;
                myfile1<< "\n";
            }
        }
        
    }
    for (k=0;k<K;k++)
    {
        if (seeds1.is_open())
        {
            for(d=0;d<D;d++)
            {
                seeds1<<M[k].x[d];
                seeds1<< " ";
            }
            
            seeds1<<"0";
            seeds1<<" ";
            seeds1<<k;
            seeds1<< "\n";
        }
    }

    return 0;
}


double alea( double a, double b )
{ // random number (uniform distribution) in [a b]
    double r;
    r=(double)rand(); r=r/RAND_MAX;
    return a + r * ( b - a );
}

double perf( int s, struct position H[S_max])
{ // Evaluate the fitness value for the particle of rank s
    int d;
    int k;
    double f, xd;
    struct position xs;
    
    xs = H[s];
    //rastriging
    k = 10;
    f = 0;
    for ( d = 0; d < D; d++ )
    {
        xd = xs.x[d];
        f += xd * xd - k * cos( 2 * pi * xd );
    }
    f += D * k;
    return f;
}

void k_means()
{
    int k, d, s;
    double k1, kt;
    bool change;
    bool insert;
    int seed=-1;
    double JT; // temporary J
    double J=0;
    
    for (s=0;s<S;s++) X[s].seed=-1;
    
    for(int rk=0; rk<10; rk++) // n times random seeds to find the best seeds
    {
        for (k=0; k<K; k++) //initialize seeds
        {
            for (d=0; d<D; d++)
            {
                MT[k].x[d]= alea( xmin[d], xmax[d] );
            }
        }
        
        do
        {
            change =false;
            for (k=0; k<K; k++)MT[k].size=0;
            for(s=0; s<S; s++) // for each particle i do
            {
                k1=0;
                insert=false; //doesn't belong to a cluster
                for (k=0; k<K; k++) // find the nearest seed mk
                {
                    for (d=0; d<D; d++)
                    {
                        v[d] = P[s].x[d]-MT[k].x[d];
                    }
                    kt=sqrt(inner_product(v, v+D, v, 0.0L)); // calculate distance p-m
                    if((insert==false ) || kt<k1 )
                        // if is the first evaluation or a smaller distance found
                    {
                        insert=true;
                        k1=kt; // set the smallest distance
                        seed=k;
                    }
                }
                // assign i to the cluster ck
                if(X[s].seed!=seed) // if found a nearer seed set it
                {
                    X[s].seed=seed;
                    change=true; // something has changed
                }
                MT[X[s].seed].size+=1;// increase the size of the cluster
                
            }
            for(k=0; k<K; k++) // for each cluster recalculate the new mean
            {
                if(MT[k].size>0)
                {
                    for(d=0; d<D; d++)
                    {
                        MT[k].x[d]=0; // set the position to 0 to calculate the new one
                        for (s=0; s<S; s++)
                        {
                            if (X[s].seed==k)MT[k].x[d]+=P[s].x[d];// for each particle in the cluster add the PB position
                        }
                        MT[k].x[d]=MT[k].x[d]/MT[k].size; // final new position
                    }
                }
            }
        }while(change==true);
        
        
        // now calculate J of each configuration, setting M to the best foud
        
        JT=0;
        for(k=0;k<K;k++)
        {
            MT[k].sigma=0;
            for(s=0; s<S; s++)
            {
                if (X[s].seed==k)
                {
                    for (d=0; d<D; d++)
                    {
                        v[d] = P[s].x[d]-MT[k].x[d];
                    }
                    MT[k].sigma+=inner_product(v, v+D, v, 0.0L); // distance (p-m)^2
                    MT[k].size+=1;
                }
            }
            MT[k].sigma=MT[k].sigma/(MT[k].size-1);
            JT+=M[k].sigma;
        }
        if(rk==0 || JT<J) // set the smallest J
        {
            J=JT;
            for (k=0; k<K; k++)
            {
                M[k]=MT[k]; // assign the definitive distribution to the best one
            }
        }
    }
    
    // reassign the particles to cluster
    for (k=0; k<K; k++)M[k].size=0;
    for(s=0; s<S; s++) // for each particle i do
    {
        k1=0;
        insert=false; //doesn't belong to a cluster
        for (k=0; k<K; k++) // find the nearest seed mk
        {
            for (d=0; d<D; d++)
            {
                v[d] = P[s].x[d]-M[k].x[d];
            }
            kt=sqrt(inner_product(v, v+D, v, 0.0L)); // calculate distance p-m
            if((insert==false ) || kt<k1 )
                // if is the first evaluation or a smaller distance found
            {
                insert=true;
                k1=kt; // set the smallest distance
                seed=k;
            }
        }
        // assign i to the cluster ck
        X[s].seed=seed;
        M[X[s].seed].size+=1;// increase the size of the cluster
        if( M[X[s].seed].best==false || P[s].f<M[X[s].seed].bp.f)
        {
            M[X[s].seed].bp=P[s];
            M[X[s].seed].best=true;
        }
    }
}

void update()
{
    int s, d;
    for (s=0; s<S; s++)
    {
        if(X[s].seed!=-1) // if belongs to a cluster
        {
            for ( d = 0; d < D; d++ )
            {
                V[s].v[d] = V[s].v[d] + alea( 0, c ) * ( P[s].x[d] - X[s].x[d] );
                V[s].v[d] = V[s].v[d] + alea( 0, c ) * ( M[X[s].seed].bp.x[d] - X[s].x[d] );
                V[s].v[d] = w*V[s].v[d];
                if(V[s].v[d]>2*(sqrt(M[X[s].seed].sigma))) V[s].v[d]=2*(sqrt(M[X[s].seed].sigma));
                if(V[s].v[d]<-2*(sqrt(M[X[s].seed].sigma))) V[s].v[d]=-2*(sqrt(M[X[s].seed].sigma));
                X[s].x[d] = X[s].x[d] + V[s].v[d];
            }
            X[s].f = perf( s,X);
            if (X[s].f<P[s].f)
            {
                P[s]=X[s];
                if(P[s].f<M[X[s].seed].bp.f)
                {
                    M[X[s].seed].bp=P[s];
                }
            }
        }
        else
        {
            for ( d = 0; d < D; d++ )
            {
                V[s].v[d] = V[s].v[d] + alea( 0, c ) * ( P[s].x[d] - X[s].x[d] ); // cognition only
                V[s].v[d] = w*V[s].v[d];
                X[s].x[d] = X[s].x[d] + V[s].v[d];
            }
            X[s].f =perf( s,X);
            if (X[s].f<P[s].f)
            {
                P[s]=X[s];
            }
        }
        
        // ... interval confinement (keep in the box)
        for ( d = 0; d < D; d++ )
        {
            if ( X[s].x[d] < xmin[d] )
            {
                X[s].x[d] = xmin[d]; V[s].v[d] = 0;
            }
            if ( X[s].x[d] > xmax[d] )
            {
                X[s].x[d] = xmax[d]; V[s].v[d] = 0;
            }
        }

    }

}

void identify_niches()
{
    int navg=0; // avarage number of particles per cluster
    int nu;
    double wf; //worst fitness
    double worst=-1; // worst particle
    int k, s, d;
    
    for(k=0; k<K; k++)
    {
        navg+=M[k].size;
    }
    navg=navg/K; // calculate average number of particles per cluster
    nu=0;
    for(k=0; k<K; k++)
    {
        if (M[k].size>navg)
        {
            for(int z=0; z<M[k].size-navg; z++)
            {
                wf=0;
                for(s=0;s<S;s++)
                {
                    if(X[s].seed==k)
                    {
                        if (X[s].f>wf)
                        {
                            wf=X[s].f;
                            worst=s;
                        }
                    }
                }
                for(s=worst;s<S;s++) // remove the nj-navg worst particles from cj
                {
                    X[s]=X[s+1];
                    P[s]=P[s+1];
                    V[s]=V[s+1];
                }
            }
            nu+=M[k].size-navg;
            M[k].size-=M[k].size-navg;
        }
    }
    for(s=S-nu;s<S;s++) // reinitialize the nu un-niched particles
    {
        for ( d = 0; d < D; d++ )
        {
            X[s].x[d] = alea( xmin[d], xmax[d] );
            V[s].v[d] = (alea( xmin[d], xmax[d] ) - X[s].x[d])/2; // Non uniform
            // V[s].v[d] = ( xmin[d]-xmax[d] )*(0.5-alea(0,1)); //Uniform. 2006-02-24
            //cout << X[s].x[0] << "\n";
        }
        X[s].f = perf( s,X);
        P[s] = X[s]; // Best position = current one
        X[s].seed=-1;
    }
}
