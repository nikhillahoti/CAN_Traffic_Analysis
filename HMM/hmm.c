//
// Hidden Markov Model program for written English
//
// Consistent with "Revealing Introduction" version dated January 12, 2017
//
// The program reads in data (English text) ignoring
// all punctuation, leaving only 26 letters and
// word spaces. Consequently, the A matrix is N x N and 
// the B matrix is N x M, where M = 27. The program
// begins with random (approximately uniform) A, B and pi,
// then attempts to climb to the optimal A and B.
// This program implements the work described in 
// Cave and Neuwirth's paper "Hidden Markov Models
// for English"
//
// Note: The data file must not have any formatting!
//
// To compile: gcc -o hmm hmm.c -O3
//

#include <string.h>
#include <stdlib.h>

#include "hmm.h"

void saveBToFile(double B[N][M]);

int main(int argc, const char *argv[])
{
    int maxIters,
        i,
        j,
        T,
        iter;
        
    int seed;

    double logProb,
           newLogProb;

    double pi[N],
           piBar[N],
           A[N][N],
           Abar[N][N],
           B[N][M],
           Bbar[N][M];
           
    char fname[80];
    int NumberOfModels = 1;
    
    struct stepStruct *step;

    if(argc != 5)
    {
oops:   fprintf(stderr, "\nUsage: %s filename maxIters SizeOfT NumberOfModels\n\n", argv[0]);
        fprintf(stderr, "where filename == input file,  maxIters = number of Iterations, SizeOfT is the length of the observations sequence and Number of Models\n");
        exit(0);
    }

//    sprintf(fname, argv[1]);
    strcpy(fname, argv[1]);
    maxIters = atoi(argv[2]);
    T = atoi(argv[3]);
    NumberOfModels = atoi(argv[4]);

    ////////////////////////
    // read the data file //
    ////////////////////////
    
    // determine number of observations
    printf("GetT... ");
    fflush(stdout);
    //T = 134200;
    
    printf("T = %d\n", T);

    // allocate memory
    printf("allocating %lu bytes of memory... ", (T + 1) * sizeof(struct stepStruct));
    fflush(stdout);
    if((step = calloc(T + 1, sizeof(struct stepStruct))) == NULL)
    {
        fprintf(stderr, "\nUnable to allocate alpha\n\n");
        exit(0);
    }
    printf("done\n");

    // read in the observations from file
    printf("GetObservations... ");
    fflush(stdout);
    T = GetObservations(fname, 
                        step, 
                        T);

    printf("T = %d\n", T);

    /////////////////////////
    // hidden markov model //
    /////////////////////////

    srandom(seed);

    int modelNumber = 0;

    while (modelNumber < NumberOfModels){
    
        // initialize pi[], A[][] and B[][]
        initMatrices(pi, A, B, seed);
        
        // print pi[], A[][] and B[][] transpose
        printf("\nN = %d, M = %d, T = %d\n", N, M, T);
        printf("initial pi =\n");
        printPi(pi);
        printf("initial A =\n");
        printA(A);
        printf("initial B^T =\n");
        printBT(B);

        // initialization
        iter = 0;
        logProb = -1.0;
        newLogProb = 0.0;

        // main loop
        while((iter < maxIters))
        {
            printf("\nbegin iteration = %d\n", iter);

            logProb = newLogProb;

            // alpha (or forward) pass
            printf("alpha pass... ");
            fflush(stdout);
            alphaPass(step, pi, A, B, T);
            printf("done\n");
            
            // beta (or backwards) pass
            printf("beta pass... ");
            fflush(stdout);
            betaPass(step, pi, A, B, T);
            printf("done\n");
            
            // compute gamma's and diGamma's
            printf("compute gamma's and diGamma's... ");
            fflush(stdout);
            computeGammas(step, pi, A, B, T);
            printf("done\n");
            
            // find piBar, reestimate of pi
            printf("reestimate pi... ");
            fflush(stdout);
            reestimatePi(step, piBar);
            printf("done\n");
            
            // find Abar, reestimate of A
            printf("reestimate A... ");
            fflush(stdout);
            reestimateA(step, Abar, T);
            printf("done\n");
            
            // find Bbar, reestimate of B
            printf("reestimate B... ");
            fflush(stdout);
            reestimateB(step, Bbar, T);
            printf("done\n");
            
    #ifdef PRINT_REESTIMATES
            printf("piBar =\n");
            printPi(piBar);
            printf("Abar =\n");
            printA(Abar);
            printf("Bbar^T = \n");
            printBT(Bbar);
    #endif // PRINT_REESTIMATES

            // assign pi, A and B corresponding "bar" values
            for(i = 0; i < N; ++i)
            {
                pi[i] = piBar[i];
            
                for(j = 0; j < N; ++j)
                {
                    A[i][j] = Abar[i][j];
                }

                for(j = 0; j < M; ++j)
                {
                    B[i][j] = Bbar[i][j];
                }
                
            }// next i

            // compute log [P(observations | lambda)], where lambda = (A,B,pi)
            newLogProb = 0.0;
            for(i = 0; i < T; ++i)
            {
                newLogProb += log(step[i].c);
            }
            newLogProb = -newLogProb;

            // a little trick so that no initial logProb is required
            if(iter == 0)
            {
                logProb = newLogProb - 1.0;
            }

            printf("completed iteration = %d, log [P(observation | lambda)] = %f\n", 
                    iter, newLogProb);

            ++iter;

        }// end while
        
        printf("\nT = %d, N = %d, M = %d, iterations = %d\n\n", T, N, M, iter);
        printf("final pi =\n");
        printPi(pi);
        printf("\nfinal A =\n");
        printA(A);
        printf("\nfinal B^T =\n");
        printBT(B);
        printf("\nlog [P(observations | lambda)] = %f\n\n", newLogProb);
        printf("\n\n Number of Iterations executed -> %d", iter);

        saveBToFile(B);

        modelNumber++;
    }
}// end hmm


void saveBToFile(double B[N][M]){
    FILE *fp = fopen("BMatrix3.txt", "a");

    for(int i = 0 ; i < M ; i++){
        for (int j = 0 ; j < N ; j++){
            fprintf(fp, "%f,", B[j][i]);
        }
        fprintf(fp, "%d\n", i);
    }
    fclose(fp);
}

//
// alpha pass (or forward pass) including scaling
//
void alphaPass(struct stepStruct *step,
               double pi[], 
               double A[][N],
               double B[][M],
               int T)
{
    int i,
        j,
        t;
        
    double ftemp;
    
    // compute alpha[0]'s
    ftemp = 0.0;
    for(i = 0; i < N; ++i)
    {
        step[0].alpha[i] = pi[i] * B[i][step[0].obs];
        ftemp += step[0].alpha[i];
    }
    step[0].c = 1.0 / ftemp;

    // scale alpha[0]'s
    for(i = 0; i < N; ++i)
    {
        step[0].alpha[i] /= ftemp;
    }

    // alpha pass
    for(t = 1; t < T; ++t)
    {
        ftemp = 0.0;
        for(i = 0; i < N; ++i)
        {
            step[t].alpha[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].alpha[i] += step[t - 1].alpha[j] * A[j][i];
            }
            step[t].alpha[i] *= B[i][step[t].obs];
            ftemp += step[t].alpha[i];
        }
        step[t].c = 1.0 / ftemp;
        
        // scale alpha's
        for(i = 0; i < N; ++i)
        {
            step[t].alpha[i] /= ftemp;
        }
    
    }// next t
    
}// end alphaPass


//
// beta pass (or backwards pass) including scaling
//
void betaPass(struct stepStruct *step,
              double pi[], 
              double A[][N],
              double B[][M],
              int T)
{
    int i,
        j,
        t;

    // compute scaled beta[T - 1]'s
    for(i = 0; i < N; ++i)
    {
        step[T - 1].beta[i] = 1.0 * step[T - 1].c;
    }

    // beta pass
    for(t = T - 2; t >= 0; --t)
    {
        for(i = 0; i < N; ++i)
        {
            step[t].beta[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].beta[i] += A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j];
            }
            
            // scale beta's (same scale factor as alpha's)
            step[t].beta[i] *= step[t].c;
        }

    }// next t
        
}// end betaPass


//
// compute gamma's and diGamma's including optional error checking
//
void computeGammas(struct stepStruct *step,
                   double pi[], 
                   double A[][N],
                   double B[][M],
                   int T)
{
    int i,
        j,
        t;
        
    double denom;

#ifdef CHECK_GAMMAS
    double ftemp,
           ftemp2;
#endif // CHECK_GAMMAS

    // compute gamma's and diGamma's
    for(t = 0; t < T - 1; ++t)// t = 0,1,2,...,T-2
    {
        
#ifdef CHECK_GAMMAS
        ftemp2 = 0.0;
#endif // CHECK_GAMMAS

        for(i = 0; i < N; ++i)
        {
            step[t].gamma[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].diGamma[i][j] = (step[t].alpha[i] * A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j]);
                step[t].gamma[i] += step[t].diGamma[i][j];
            }

#ifdef CHECK_GAMMAS
            // verify that gamma[i] == alpha[i]*beta[i] / sum(alpha[j]*beta[j])
            ftemp2 += step[t].gamma[i];
            ftemp = 0.0;
            for(j = 0; j < N; ++j)
            {
                ftemp += step[t].alpha[j] * step[t].beta[j];
            }
            ftemp = (step[t].alpha[i] * step[t].beta[i]) / ftemp;
            if(DABS(ftemp - step[t].gamma[i]) > EPSILON)
            {
                printf("gamma[%d] = %f (%f) ", i, step[t].gamma[i], ftemp);
                printf("********** Error !!!\n");
            }
#endif // CHECK_GAMMAS

        }// next i
            
#ifdef CHECK_GAMMAS
        if(DABS(1.0 - ftemp2) > EPSILON)
        {
            printf("sum of gamma's = %f (should sum to 1.0)\n", ftemp2);
        }
#endif // CHECK_GAMMAS
            
    }// next t
    
    // special case for t = T-1
    for(j = 0; j < N; ++j)
    {
        step[T-1].gamma[j] = step[T-1].alpha[j];
    }

}// end computeGammas


//
// reestimate pi, the initial distribution
//
void reestimatePi(struct stepStruct *step, 
                  double piBar[])
{
    int i;
    
    // reestimate pi[]        
    for(i = 0; i < N; ++i)
    {
        piBar[i] = step[0].gamma[i];
    }
        
}// end reestimatePi


//
// reestimate the A matrix
//
void reestimateA(struct stepStruct *step, 
                 double Abar[][N], 
                 int T)
{
    int i,
        j,
        t;
    
    double numer,
           denom;
           
    // reestimate A[][]
    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < N; ++j)
        {
            numer = denom = 0.0;

            // t = 0,1,2,...,T-2
            for(t = 0; t < T - 1; ++t)
            {
                numer += step[t].diGamma[i][j];
                denom += step[t].gamma[i];
                
            }// next t

            Abar[i][j] = numer / denom;
        
        }// next j
        
    }// next i
        
} // end reestimateA    


//
// reestimate the B matrix
//
void reestimateB(struct stepStruct *step, 
                 double Bbar[][M], 
                 int T)
{
    int i,
        j,
        t;
    
    double numer,
           denom;
           
    // reestimate B[][]
    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < M; ++j)
        {
            numer = denom = 0.0;

            // t = 0,1,2,...,T-1
            for(t = 0; t < T; ++t)
            {
                if(step[t].obs == j)
                {
                    numer += step[t].gamma[i];
                }
                denom += step[t].gamma[i];

            }// next t

            Bbar[i][j] = numer / denom;
        
        }// next j
        
    }// next i
        
}// end reestimateB


//
// initialize pi[], A[][] and B[][]
//
void initMatrices(double pi[], 
                  double A[][N], 
                  double B[][M],
                  int seed)
{
    int i,
        j;
        
    double prob,
           ftemp,
           ftemp2;
    
    // initialize pseudo-random number generator
    srandom(seed);

    // initialize pi
    prob = 1.0 / (double)N;
    ftemp = prob / 10.0;
    ftemp2 = 0.0;
    for(i = 0; i < N; ++i)
    {
        if((random() & 0x1) == 0)
        {
            pi[i] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
        }
        else
        {
            pi[i] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
        }
        ftemp2 += pi[i];
        
    }// next i
    
    for(i = 0; i < N; ++i)
    {
        pi[i] /= ftemp2;
    }

    // initialize A[][]
    prob = 1.0 / (double)N;
    ftemp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        ftemp2 = 0.0;
        for(j = 0; j < N; ++j)
        {
            if((random() & 0x1) == 0)
            {
                A[i][j] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
            }
            else
            {
                A[i][j] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
            }
            ftemp2 += A[i][j];
            
        }// next j
        
        for(j = 0; j < N; ++j)
        {
            A[i][j] /= ftemp2;
        }
        
    }// next i

//    A[0][0] = .275;
//    A[0][1] = .725;
//    A[1][0] = .780;
//    A[1][1] = .220;
    
    // initialize B[][]
    prob = 1.0 / (double)M;
    ftemp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        ftemp2 = 0.0;
        for(j = 0; j < M; ++j)
        {
            if((random() & 0x1) == 0)
            {
                B[i][j] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
            }
            else
            {
                B[i][j] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
            }
            ftemp2 += B[i][j];
            
        }// next j
        
        for(j = 0; j < M; ++j)
        {
            B[i][j] /= ftemp2;
        }
        
    }// next i
    
}// end initMatrices


//
// read and save observations
//
int GetObservations(char fname[], 
                    struct stepStruct *step,
                    int T)
{
    FILE *in;
            char temp[MAX_CHARS + 2];
    
    char space[1] = {" "};
    
    char alphabet[M] = ALPHABET;
    
    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }
    
    // read 'em in
    int num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL && num < T)
    {
        step[num++].obs = atoi(temp);
    }// end while
    
    fclose(in);

    return(T);

}// end GetObservations


//
// print pi[]
//
void printPi(double pi[])
{
    int i;
        
    double ftemp;

    ftemp = 0.0;
    for(i = 0; i < N; ++i)
    {
        printf("%8.5f ", pi[i]);
        ftemp += pi[i];
    }
    printf(",  sum = %f\n", ftemp);

}// end printPi


//
// print A[][]
//
void printA(double A[][N])
{
    int i,
        j;
        
    double ftemp;

    for(i = 0; i < N; ++i)
    {
        ftemp = 0.0;
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", A[i][j]);
            ftemp += A[i][j];
        }
        //printf(", Before sum = %f\n", ftemp);
        A[i][N-1] = 1.00 - ftemp;
        printf(", sum = %f\n", ftemp);

    }// next i

}// end printA


//
// print BT[][]
//
void printBT(double B[][M])
{
    int i,
        j;
        
    double ftemp;
    
    char alphabet[M] = ALPHABET;
    
    for(i = 0; i < M; ++i)
    {
        printf("%d ", i);
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", B[j][i]);
        }
        printf("\n");
    }
    for(i = 0; i < N; ++i)
    {
        ftemp = 0.0;
        for(j = 0; j < M; ++j)
        {
            ftemp += B[i][j];
        }
        //printf("Before sum[%d] = %f ", i, ftemp);
        B[i][M-1] = 1.00 - ftemp;
        printf("sum[%d] = %f ", i, ftemp);
    }
    printf("\n");

}// end printB
