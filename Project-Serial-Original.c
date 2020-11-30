#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#define PI 3.14159265

//variable declarations
FILE *fd;
int nx, ny, nt, auto_save;
double dx,dy,dt,sigma,nu,Lx,Ly;
double *x, *y;
double **u, **v, **un, **vn; // A 2D array of char arrays (a pointer to pointers to chars)
int char_counts[26];
char hostname[256];
double tstart, ttotal;

double myclock() 
{
	static time_t t_start = 0;  // Save and subtract off each time

	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	if( t_start == 0 ) t_start = ts.tv_sec;

	return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}

void * linspace(double a, double b, int n, double *array)
{
    double d;
    int i;
    
    /* Check if correct input*/
    if(n < 2 || array == 0)
        printf("Error in linspace\n");
    
    d = (b - a)/(n - 1);

    for(i = 0; i < n - 1; ++i)
        array[i] = a + i*d;
    array[n - 1] = b;
}


void * ones2D (double** array)
{
	int i,j;

	for ( i = 0; i < ny; i++) 
	{
		for ( j = 0; j < nx; j++ ) 
		{
			array[i][j] = 1;
		}
	}

}

void * create_arrays()
{


	x = (double *)malloc(nx*sizeof(double));
	y = (double *)malloc(ny*sizeof(double));

	// allocate an array of size ARRAY_SIZE pointers to chars
	u    = (double **)malloc(sizeof(double *)*ny);
	v    = (double **)malloc(sizeof(double *)*ny);
	un   = (double **)malloc(sizeof(double *)*ny);
	vn   = (double **)malloc(sizeof(double *)*ny);

	// for each row in char_array, we add to it an array of size STRING_SIZE
	// to make it a 2d array
	long i;
	for ( i=0; i < ny; i++) 
	{
		u[i]    = (double *)malloc(sizeof(double)*nx);
		v[i]    = (double *)malloc(sizeof(double)*nx);
		un[i]   = (double *)malloc(sizeof(double)*nx);
		vn[i]   = (double *)malloc(sizeof(double)*nx);
	}
}


void *write_output(int time)
{
	int i,j;

	char filename1[20], filename2[20];
	snprintf(filename1,20,"/homes/mcheikh/CIS_625/Project/Serial-Original/output/u-%d",time);
	snprintf(filename2,20,"/homes/mcheikh/CIS_625/Project/Serial-Original/output/v-%d",time);

	fd = fopen(filename1, "w" );
	for ( i = 0; i < ny; i++) 
	{
		for ( j = 0; j < nx; j++ ) 
		{
			fprintf(fd, "%lf",u[i][j]);
			if (j != nx-1) 	fprintf(fd, ",");
		}
		fprintf(fd,"\n");
	}
	fclose( fd );

	fd = fopen(filename2,"w" );
	for ( i = 0; i < ny; i++) 
	{
		for ( j = 0; j < nx; j++ ) 
		{
			fprintf(fd, "%lf",v[i][j]);
			if (j != nx-1) 	fprintf(fd, ",");
		}
		fprintf(fd,"\n");
	}
	fclose( fd );
}

void *write_boundary_data()
{
	char filename[20];
	snprintf(filename,20,"/homes/mcheikh/CIS_625/Project/Serial-Original/output/boundary");

	fd = fopen(filename,"w");
	fprintf(fd, "%lf,%lf,%d,%d,%d,%lf,%lf,%d\n",Lx,Ly,nx,ny,nt,sigma,nu,auto_save);
	fclose( fd );

}

void *init_condition()
{
	int i,j;
	for ( i = (int) (0.5/dy); i < (int) ((Ly-0.5)/dy); i++) 
	{
		for ( j = (int) (0.5/dx); j < (int) ((Lx-0.5)/dx); j++ ) 
		{
			u[i][j] = sin(dx*j)*cos(dy*i)+2;
			v[i][j] = sin(dx*j)*cos(dy*i)+2;
		}
	}
}

void * copy_arrays()
{
	int i,j;
	for ( i = 0; i < ny; i++) 
		for ( j = 0; j < nx; j++ ) 
		{
			un[i][j] = u[i][j];
			vn[i][j] = v[i][j];
		}
}

void * update_boundary()
{
	int i,j;
	for (j = 0; j < nx; j++)
	{
		u[0][j]    = u[ny-2][j];
		u[ny-1][j] = u[1][j];

		v[0][j]    = v[ny-2][j];
		v[ny-1][j] = v[1][j];
	}

	for (i = 0; i < ny; i++)
	{
		u[i][0]  	= u[i][nx-2];
		u[i][nx-1] 	= u[i][1];
		
		v[i][0] 	= v[i][nx-2];
		v[i][nx-1] 	= v[i][1];
	}
}

int main(int argc, char *argv[]) 
{
	tstart = myclock(); // Global Clock
	int n;
	int i,j;
	struct rusage ru;

	// Default Values:
	nx = 101; 
	ny = 101;
	nt = 6000;
	Lx = 6.0; // Domain length in x
	Ly = 6.0; // Domain length in y
	auto_save = 60; // Auto save 

	if (argc >= 6){
		nx = atol(argv[1]);
		ny = atol(argv[2]);
		nt = atol(argv[3]);
		Lx = atol(argv[4]);
		Ly = atol(argv[5]);
		auto_save = atol(argv[6]);
	}

	sigma = 0.0003;
	nu = 0.0001;

	dx = Lx / (nx - 1);
	dy = Ly / (ny - 1);
	dt = sigma * dx * dy / nu ;

	//write_boundary_data();

	create_arrays();

	linspace(0,Lx,nx, x);
	linspace(0,Ly,ny, y);

	ones2D(u);
	ones2D(v);
	ones2D(un);
	ones2D(vn);

	init_condition();
	
	//write_output(0);

	for (n=1;n<=nt;n++)
	{
		copy_arrays(); // copys u and v into un and vn

		for (i = 1; i < ny -1; i ++)
			for (j = 1; j < nx - 1; j++)
				{
					u[i][j] = (un[i][j] - dt/dx*un[i][j]*(un[i][j] - un[i][j-1]) 
					 	      - dt/dy*vn[i][j]*(un[i][j] - un[i-1][j]) 
					 	      + nu * dt / (dx*dx) * (un[i][j+1] - 2.0 * un[i][j] + un[i][j-1]) 
					 	      + nu * dt / (dy*dy) * (un[i+1][j] - 2.0 * un[i][j] + un[i-1][j]) );
				
					v[i][j] = (vn[i][j] - dt/dx*un[i][j]*(vn[i][j] - vn[i][j-1])
							  - dt/dy*vn[i][j]*(vn[i][j] - vn[i-1][j]) 
							  + nu * dt / (dx*dx) * (vn[i][j+1] - 2.0 * vn[i][j] + vn[i][j-1]) 
							  + nu * dt / (dy*dy) * (vn[i+1][j] - 2.0 * vn[i][j] + vn[i-1][j]) );
				}

		update_boundary();
		//if(n % auto_save == 0) write_output(n);
	}

	ttotal = myclock() - tstart;

	getrusage(RUSAGE_SELF, &ru);
    long MEMORY_USAGE = ru.ru_maxrss;   // Memory usage in Kb

	printf("Serial-Original: Simulation of nx:%d,ny:%d,nt:%d,Lx:%lf,Ly:%lf,auto_save:%d, took %lf seconds, with %ld Kb\n", nx,ny,nt,Lx,Ly,auto_save, ttotal, MEMORY_USAGE);	

	return 0;
}