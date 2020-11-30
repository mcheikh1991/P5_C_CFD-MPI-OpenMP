#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <mpi.h>

// Global Variable declarations:
//------------------------------
FILE *fd;
int nx, ny, nt, auto_save, nx_local;
double dx,dy,dt,sigma,nu,Lx,Ly;
double *x, *y, *x_local, *y_local;
double **u, **v, **un, **vn; // A 2D array original size

int char_counts[26];
char hostname[256];
double tstart, ttotal;
int MPIthreads = 1;
int i_start, i_end, j_start, j_end;
int neig_r, neig_l, neig_t, neig_b; // Neighbors MPI id
// Data Sent/Recieved Accross the cores
double *u_r, *u_l, *u_t, *u_b; 
double *v_r, *v_l, *v_t, *v_b; 

// Functions:
//-----------

#include "Functions.H"
/* This list include functions like:
------------------------------------
myclock
linspace
ones2D
create_arrays
create_matrices
write_boundary_data
write_domain_info
write_output
copy_arrays
*/

#include "Functions_MPI.H"
/* This list include functions like:
------------------------------------
init_condition
update_boundary
solve_equation
setup_MPI_Domain
send_boundary_data
recieve_boundary_date
*/

// Main Function:
//---------------

int main(int argc, char *argv[]) 
{
	tstart = myclock(); // Global Clock
	int n, i, rank;
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

	// Starting the MPI Process
	int rc;
    MPI_Request Request;
	MPI_Status Status;

	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS){
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

    MPI_Comm_size(MPI_COMM_WORLD,&MPIthreads); // Number of cores
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	   // rank of each core

	//write_boundary_data();

	create_arrays();

	linspace(0,Lx,nx, x);
	linspace(0,Ly,ny, y);

	setup_MPI_Domain(&rank); // This function divides the array according to the mpi-process

	//write_domain_info( &rank);

	create_matrices();
	
	ones2D(u);
	ones2D(v);
	ones2D(un);
	ones2D(vn);

	init_condition(&rank);
	
	copy_arrays(); // copys u and v into un and vn

	//write_output(0,&rank);

	for (n=2;n<=nt;n+=2)
	{
		// Inside foor loop

		solve_equation(1); // solve Bergurs equation for un and vn

		send_boundary_data(&Request, &rank, 2); // send data from un and vn

		recieve_boundary_date(&Status,&Request,&rank);

		update_boundary(&rank, 2); // update un and vn

		solve_equation(2); // solve Bergurs equation for u and v

		send_boundary_data(&Request, &rank, 1);

		recieve_boundary_date(&Status,&Request,&rank);

		update_boundary(&rank, 1); 

		//if(n % auto_save == 0) write_output(n,&rank);
	}

	ttotal = myclock() - tstart;

	getrusage(RUSAGE_SELF, &ru);
    long MEMORY_USAGE = ru.ru_maxrss;   // Memory usage in Kb

    if(rank ==0)
		printf("MPI: Simulation of nx:%d,ny:%d,nt:%d,Lx:%lf,Ly:%lf,auto_save:%d,nthreads:%d, took %lf seconds, with %ld Kb\n", nx,ny,nt,Lx,Ly,auto_save,MPIthreads, ttotal, MEMORY_USAGE);	


	MPI_Finalize();
	return 0;

}