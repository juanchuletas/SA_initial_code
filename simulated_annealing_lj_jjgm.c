//--------------------------------------------------------------//
// Simulated Annealing Code for N particles into a Sphere	//
// Juan José García Miranda					//
// Maestría en Ciencias Químicas				//
// Universidad Autónomoa Metropolitana Unidad Iztapalapa	//
// Noviembre  2017						//
// -------------------------------------------------------------//	
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
struct Atom
{
	double rx,ry,rz;
};
void system_init(struct Atom atoms[],int n_atom)
{
	int i;
	for(i=0;i<n_atom;i++)
	{
		atoms[i].rx = 0.0;
		atoms[i].ry = 0.0;
		atoms[i].rz = 0.0;
	}
}
double random_g(double r_user)
{
	double target,up,down;
	up = 2.0*r_user;
	down = -1*2.0*r_user;

	target  = drand48()*(up-down) + down;

	return target;
}

/*Random asigments MODULE---------------------------------------------------------- */
void random_moves(struct Atom atoms[],int n_atom, struct Atom atomran[],double r)
{
	int i,j, flag;
  	flag = 0;
	for (i=0; i<n_atom;i++)
	{
	  atomran[i].rx = atoms[i].rx + random_g(r);
	  atomran[i].ry = atoms[i].ry + random_g(r);
	  atomran[i].rz = atoms[i].rz + random_g(r);
	}

}
/* SPHERE MODULE----------------------------------------------------------------------- */
int sphere_criterion(struct Atom atoms[],int n_atom,struct Atom atomran[],double r_user) 
{
	int i,flag;
	double r_xc,r_yc,r_zc,r_ic2,r_ic;
	flag = 0;
	for (i=0; i<n_atom; i++)
	{
		
		r_xc = atomran[i].rx;
		r_yc = atomran[i].ry;
		r_zc = atomran[i].rz;
		r_ic2 = (r_xc*r_xc + r_yc*r_yc + r_zc*r_zc);
		r_ic = sqrt(r_ic2);         /*Distance from the sphere center per particle*/
		if (r_ic > r_user)
		{
			flag = 1;
		}
	}
	return flag;
}
/* Energy  MODULE-------------------------------------------------------------------- */
double energy(struct Atom atomran[],int n_atom) 
{
  int i, j;
  double dx, dy, dz, r, r2, uij, sigma, sigma6, epsilon, ep_init;
  sigma = 0.3243;
  epsilon = 0.9518;
  ep_init = 0.f;
  for (i=0;i<n_atom-1; i++)
  {
    for (j=i+1;j<n_atom;j++)
    {
      dx=(atomran[i].rx-atomran[j].rx);
      dy=(atomran[i].ry-atomran[j].ry);//compute distances
      dz=(atomran[i].rz-atomran[j].rz);// xyz coordinates
      r2=(dx*dx + dy*dy + dz*dz);
      r=sqrt(r2);
      sigma6 = pow((sigma/r),(double) 6); //from Lennard-Jones equation
      uij= 4.f*epsilon*sigma6*(sigma6-1.f); //factorization of LJ eq.
      ep_init = ep_init + uij;
    }
  }
  return ep_init;
}
/*Cooling Schedule---------------------------------------------------------------------------------------------*/
double low_temperature(double t_c, int i)
{
        double alpha,Tf;
        alpha  = 0.3;

        //Tf = t_c*pow(alpha,i);
	Tf = (t_c)*exp(-alpha*i/3.0);
	return Tf;
}
/*Simulated Annealing & Metropolis MODULE---------------------------------------------------------------------*/
int metropolis(double epij,double T,double **in_energy)
{
	int flag;
	double fin_energy,delta_E;
	double beta,p;
	
	beta = 1.0/(T);
	flag = 0;
        fin_energy = epij;
	printf("Temperature in Metropolis Module %lf\n",T);


	delta_E = fin_energy - **in_energy;

	if(delta_E<0.f)
	{

                **in_energy = fin_energy;
		flag = 1;
	}
	else 
	{
		p = rand()/(double)RAND_MAX;
		if(p<exp((-beta*delta_E)))
		{
                	**in_energy = fin_energy;
			flag = 1;
		}
	}

	return flag;
}
/*Write Energy and Structures MODULE-----------------------------------------------------------*/
void write_structures(FILE *outputdata,struct Atom atomran[],int n_atom, double energy,double k,double r_user)
{
	int i;

	fprintf(outputdata,"%d\n",n_atom);
	fprintf(outputdata,"Energy of the system = %lf\n",energy);
	for(i=0; i<n_atom; i++)
		{
			fprintf(outputdata,"Na %lf	%lf	%lf\n",atomran[i].rx,atomran[i].ry,atomran[i].rz);
		}
}
/*----- Run the simulated annealing   ---------------------------------------------------------------------*/
void simulated_annealing(double r_user,int n_atoms,int n_steps,double temp, double *in_energy)
{
        FILE *outputdata;
        outputdata=fopen("Na_SA.dat","wb");
        int k,flag,i,cond;
        double target,band,aux,t_aux,t_trial;
        double tk,secs;
        struct Atom atoms[100],atomran[100];
	clock_t t_init,t_fin;
        t_init = clock();

	system_init(atoms,n_atoms);
	i=0;
	cond=0;
	do
	{
		printf("\n");
		tk = low_temperature(temp,i);
		printf("Initial Energy = %lf at T=%lf\n",*in_energy,tk);
		for(k=0;k<n_steps;k++)
		{
			random_moves(atoms,n_atoms,atomran,r_user);
			flag=sphere_criterion(atoms,n_atoms,atomran,r_user);

			if(flag==0)
			{
				target=energy(atomran,n_atoms);
				band=metropolis(target,tk,&in_energy);
				if(band==1)
				{

					write_structures(outputdata,atomran,n_atoms,target,k,r_user);
					printf("Accepted energy %lf in step %d\n",target,k+1);
					aux = target;
				}
			}
		}
		if(tk<=1.0e-4)
		{
			cond=1;
		}
		*in_energy = aux;
		i++;
		
	}while(cond==0);

	fclose(outputdata);
        printf("Na_SA.dat'\n");
	t_fin = clock();
        secs = (double)(t_fin - t_init) / CLOCKS_PER_SEC;
        printf("Done in %.16g ms\n", secs * 1000.0);

}
/* MAIN MODULE -------------------------------------------------------------------------------------------*/
int main  ()
{
	srand(time(NULL)); //random numbers
	int total_atoms,n_steps;
	double r_user,t_0,e_INIT=1.0E9;
	printf("Simulated Annealing program for LJ potential\n");
	printf("Restricted atoms inside a sphere \n");
	printf("Enter the sphere radius r:\n");
	scanf("%lf",&r_user);
	printf("Enter a high temperature\n");
	scanf("%lf",&t_0);
	printf("Enter the number of atoms:\n");
	scanf("%d",&total_atoms);
	printf("How many random moves?\n");
	scanf("%d",&n_steps);

	simulated_annealing(r_user,total_atoms,n_steps,t_0,&e_INIT);

	return 0;
}
/* ---------------------------------------------------------------------------------------------------*/
