#include<float.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>



//Function Headers
int nmax=10;
void main(void);
void calculate(double spin[nmax][nmax],int nmax,double t,int n_passes,double current_mag,double current_energy,double ave_mag,double ave_energy);
void spin(void);
void init_spin(double spin[nmax][nmax],int nmax,double current_m,double current_e);
void initialise(double spin[nmax][nmax],int nmax,double t_init,double t_final,int n_points,int n_passes,double current_m,double current_e);
void test_for_edge(int x, int nmax);
void test_for_flip(double spin[nmax][nmax],int j,int k,double t,double sum,double current_mag,double current_energy);
double neighbours(double spin[nmax][nmax],int nmax,int j,int k);

//Declare Arrays


//main
void main(void)
{
 int i,j,k;

 double mag[nmax];
 double temp[nmax];
 double energy[nmax];
 double spin[nmax][nmax];
 
 double t;
 double dt;
 double t_init;
 double t_final;
 int n_points;
 int n_passes;
 double current_m;
 double current_e;
 double ave_mag;
 double ave_energy;
 double sum;
 

 initialise(spin,nmax,t_init,t_final,n_points,n_passes,current_m,current_e);
 dt = (t_final - t_init) / (n_points - 1);
 t = t_init;

 init_spin(spin,nmax,current_m,current_e);

 for(i = 0;i<n_points;i++)
 {
      calculate(spin,nmax,t,n_passes,current_m,current_e,ave_mag,ave_energy);
      temp[i]= t;
      t = t + dt;
 }
}

void initialise(double spin[nmax][nmax],int nmax,double t_init,double t_final,int n_points,int n_passes,double current_m,double current_e)
{
 t_init = 1;
 t_final = 5;
 n_points = 9; 
 n_passes = 100;
 
 init_spin(spin,nmax,current_m,current_e);
}

void init_spin(double spin[nmax][nmax],int nmax,double current_m,double current_e)
{
 int i,j;
 for(i=0;i<nmax;i++)
   {
       for(j=0;j<nmax;j++)
          {
              spin[i][j] = 1;
          }
   }
   current_m = nmax*nmax;
   current_e = -2*(nmax*nmax);      
}

void calculate(double spin[nmax][nmax],int nmax,double t,int n_passes,double current_mag,double current_energy,double ave_mag,double ave_energy)
{
     double sum;
     int i,j,k;
     double m_tmp = 0;
     double e_tmp = 0;
     int n = 0;
     double m[nmax];
     double e[nmax];
     int ts[nmax];
     
     m[0] = current_mag;
     e[0] = current_energy;
     ts[0] = 0;
     
     for(i = 0; i< n_passes;i++)
           {
                for(j = 0; j < nmax;j++)
                      {
                           for(k = 0;k<nmax;k++)
                                 {
                                    sum = neighbours(spin,nmax,j,k);
                                    test_for_flip(spin,j,k,t,sum,current_mag,current_energy);
                                 }
                      }
                      m[i+1] = current_mag;
                      e[i+1] = current_energy;
                      ts[i+1] = i;
                      
                      if(i>(n_passes/2))
                         {
                            m_tmp = m_tmp + current_mag;
                            e_tmp = e_tmp + current_energy;
                            n = n+1;   
                         }
           }
           ave_mag = m_tmp/n;
           ave_energy = e_tmp/n;                     
}

double neighbours(double spin[nmax][nmax],int nmax,int j,int k)
{
   double sum;
   int x1,y1,x2,y2,x3,y3,x4,y4;
   x1 = j;
   y1 = k + 1;
   test_for_edge(y1,nmax);
   x2 = j;
   y2 = k - 1;
   test_for_edge(y2,nmax);
   x3 = j + 1;
   y3 = k;
   test_for_edge(x3,nmax);
   x4 = j - 1;
   y4 = k;
   test_for_edge(x4,nmax);
   sum = spin[x1][y1] + spin[x2][y2] + spin[x3][y3] + spin[x4][y4];
   return(sum);    
}

void test_for_edge(int x, int nmax)
{
   if( x < 1)
       x = nmax;
   if( x > nmax)
       x = 1;     
}

void test_for_flip(double spin[nmax][nmax],int j,int k,double t,double sum,double current_mag,double current_energy)
{
 double denergy;
 
 denergy = 2 * (spin[j][k])*sum;
 if(denergy < 0)
    {
       spin[j][k] = - spin[j][k];
       current_mag = current_mag + 2*spin[j][k];
       current_energy = current_energy + denergy;
    }
    else if(exp(-denergy/t) > rand()/RAND_MAX)
       {
          spin[j][k] = -spin[j][k];
          current_mag = current_mag + 2 * spin[j][k];
          current_energy = current_energy + denergy;
       }
}
