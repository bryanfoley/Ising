/*------------------------------------------------------------------------------*
  * File Name: Ising Model.c                                                     *
  * Creation: ER 7/1/2001                                                        *
  * Purpose: Origin C file                                                       *
  * Copyright (c) OriginLab Corp. 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 *
  * All Rights Reserved                                                          *
  *                                                                              *
  * Modification Log:                                                            *
  *------------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////////
// you must include this header file for all Origin built-in functions and classes
#include <origin.h>
//
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// start your functions here
// 2D Ising Model simulation using Metropolis Monte Carlo algorithm
// The code for the main loop of this algorithm is from the website:
// http://www.npac.syr.edu/users/paulc/lectures/montecarlo/node31.html
// First prototype all functions that are called by the main function
     void  imInitAllDown(Matrix &matSpin, int isize);
     void  imInitAllUp(Matrix &matSpin, int isize);
     void  imInitWarm(Matrix &matSpin, int isize);
     void  imMetropolis(Matrix &matSpin, int size, double dBeta);
     void  imAvgSpin(Matrix &matSpin, int isize, Dataset &dsAvgSpin);
     void  imAvgEnergy(Matrix &matSpin, int isize, Dataset &dsAvgEnergy);
// Main function of the code follows
// A message variable is passed to this section, which is used to then call subsequent function
void IsingModel(int iMsg)
{
     // define datasets for spin lattice, temperature, average spin, and average energy
     Matrix matLattice();
     matLattice.Attach("Matrix1");
     Dataset dsTemp("data1_temp");
     Dataset dsAvgSpin("data1_avgspin");
     Dataset dsAvgEnergy("data1_avgenergy");
     // get size of lattice and check that it is square
     int iNRows = matLattice.GetNumRows();
     int iNCols = matLattice.GetNumCols();
     if (iNRows != iNCols)
     {
         printf("Not a square matrix!\n");
          return;
     }
     // if message = 0, call function to initialize all spins to -1
     if (iMsg == 0)
     {
         imInitAllDown(matLattice, iNRows);
          return;
     }
     // if message = 1, call function to initialize all spins to +1
     else if (iMsg == 1)
     {
         imInitAllUp(matLattice, iNRows);
          return;
     }
     // if message = 2, call function to initialize spins randomly to +1/-1
     else if (iMsg == 2)
     {
         imInitWarm(matLattice, iNRows);
          return;
     }
    // if message = 3, call the Metropolis Montecarlo function to loop over lattice to change s
    else if (iMsg == 3)
    {
        double dBeta = 1/dsTemp[0];              // beta = 1/T
        imMetropolis(matLattice, iNRows, dBeta);
        return;
    }
    // if message = 4, call function to compute average spin of lattice
    else if (iMsg == 4)
    {
        imAvgSpin(matLattice, iNRows, dsAvgSpin);
        return;
    }
    // if message = 5, call function to compute average energy
    else if (iMsg == 5)
    {
        imAvgEnergy(matLattice, iNRows, dsAvgEnergy);
        return;
    }
}
// end of main function
// function to intialize all spins to -1
// note that the lattice is passed by reference
void imInitAllDown(Matrix &matSpin, int isize)
{
    int ir, ic;
    double dRan;
    for (ir=0; ir<isize; ir++)
    {
        for (ic=0; ic<isize; ic++)
        {
            matSpin[ir][ic] = -1;
        }
    }
}
// function to intialize all spins to +1
// note that the lattice is passed by reference
void imInitAllUp(Matrix &matSpin, int isize)
{
    int ir, ic;
    double dRan;
    for (ir=0; ir<isize; ir++)
    {
        for (ic=0; ic<isize; ic++)
        {
            matSpin[ir][ic] = 1;
        }
    }
}
// function to intialize spins randomly
// note that the lattice is passed by reference
void imInitWarm(Matrix &matSpin, int isize)
{
    int ir, ic;
    double dRan;
    // use timer to seed random number generator
    uint wTick = GetTickCount();
    dRan = rnd(wTick);
    for (ir=0; ir<isize; ir++)
    {
        for (ic=0; ic<isize; ic++)
        {
            dRan = rnd();
            if (dRan <= 0.5)
            {
                matSpin[ir][ic] = 1;
            }
            else matSpin[ir][ic] = -1;
        }
    }
}
// function that performs the Metropolis Monte Carlo computation to change lattice spins
// note that the lattice is passed by reference
// the code here is almost the same as at the website:
// http://www.npac.syr.edu/users/paulc/lectures/montecarlo/node31.html
void imMetropolis(Matrix &matSpin, int size, double dBeta)
{
    int i,j,iPosRow,iPosCol,iNegRow,iNegCol;
    int iOldSpin,iNewSpin,iSpinSum;
    int iOldEnergy,iNewEnergy;
    double dEnergyDiff;
    // use timer to seed random number generator
    uint wTick =GetTickCount();
    double dRan = rnd(wTick);
    /* Loop over sites */
    for (i=0;i<size;i++)
    {
        for (j=0;j<size;j++)
        {
            /* periodic boundary conditions */
            iPosRow = (i+1) % size;
            iPosCol = (j+1) % size;
            iNegRow = (i+size-1) % size;
            iNegCol = (j+size-1) % size;
            iOldSpin = matSpin[i][j];
            iNewSpin = - iOldSpin;
            /* Sum of neighboring spins */
            iSpinSum = matSpin[i][iPosCol] + matSpin[iPosRow][j] + matSpin[i][iNegCol] + matSpi
            iOldEnergy = - iOldSpin * iSpinSum;
            iNewEnergy = - iNewSpin * iSpinSum;
            dEnergyDiff = dBeta * (iNewEnergy - iOldEnergy);   /* beta = J/kT */
            /* Do the update */
            if ( (dEnergyDiff <= 0.0) || (exp(-dEnergyDiff) > rnd() ) )
            {
                 /* Accept the change */
                matSpin[i][j] = iNewSpin;
            }
        }
    }  /* end of loop over sites */
}
// function that computes average spin of lattice
// note that the lattice and average spin datasets are passed by reference
void imAvgSpin(Matrix &matSpin, int isize, Dataset &dsAvgSpin)
{
    int ir, ic;
    double dRan;
    double dAvg;
    for (ir=0; ir<isize; ir++)
    {
        for (ic=0; ic<isize; ic++)
        {
            dAvg += matSpin[ir][ic];
        }
    }
    dAvg = dAvg/(isize*isize);
    int iLen = dsAvgSpin.GetSize();
    if(iLen == 50)
    {
        iLen = 0;
        dsAvgSpin.SetUpperBound(1);
    }
    iLen += 1;
    dsAvgSpin.SetSize(iLen);
    dsAvgSpin[iLen-1] = dAvg;
    dsAvgSpin.Update(FALSE, REDRAW_REALTIME_SCOPE);
}
// function that computes average energy
// note that the lattice and the average energy dataset are passed by reference
void imAvgEnergy(Matrix &matSpin, int isize, Dataset &dsAvgEnergy)
{
    int i, j, iPosRow, iPosCol, iNegRow, iNegCol, iSpinSum;
    double dEnergy, dAvg;
    dEnergy = 0.0;
    for (i=0; i<isize; i++)
    {
        for (j=0; j<isize; j++)
        {
            iPosRow = (i+1) % isize;
            iPosCol = (j+1) % isize;
            iNegRow = (i+isize-1) % isize;
            iNegCol = (j+isize-1) % isize;
            /* Sum of neighboring spins */
            iSpinSum = matSpin[i][iPosCol] + matSpin[iPosRow][j] + matSpin[i][iNegCol] + matSpi
            dEnergy += - matSpin[i][j] * iSpinSum;
        }
    }
    dAvg = dEnergy/(isize*isize);
    int iLen = dsAvgEnergy.GetSize();
    if(iLen == 50)
    {
        iLen = 0;
        dsAvgEnergy.SetUpperBound(1);
    }
    iLen += 1;
    dsAvgEnergy.SetSize(iLen);
    dsAvgEnergy[iLen-1] = dAvg;
    dsAvgEnergy.Update(FALSE, REDRAW_REALTIME_SCOPE);
}
// end of all functions
