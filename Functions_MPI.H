#ifndef Functions_MPI_H
#define Functions_MPI_H

void *init_condition( void *rank)
{
	int i,j;
	int myID =  *((int*) rank); // Real ID, will be used to move along x
	for ( i = (int) (0.5/dy); i < (int) ((Ly-0.5)/dy); i++) 
	{
		for ( j = (int) (0.5/dx); j < (int) ((Lx-0.5)/dx); j++ ) 
		{
			if(j>=j_start && j<=j_end)
			{
				//if (myID == 1) printf("%d,%d,%d,%d\n",i,j,i_start,j_start);
				u[i-i_start][j-j_start] = sin(dx*j)*cos(dy*i)+2;
				v[i-i_start][j-j_start] = sin(dx*j)*cos(dy*i)+2;
			}
		}
	}
}


void * update_boundary( void *rank, int type)
{
	int i,j;
	if (type == 1) // Update the boundaries of u & v
	{
		for (i = 0; i < ny; i++)
		{
			u[i][nx_local-1] 	= u_r[i] ;
			u[i][0] 			= u_l[i] ;
			v[i][nx_local-1] 	= v_r[i] ;
			v[i][0] 			= v_l[i] ;
		}

		for (j = 0; j < nx_local; j++)
		{
			u[0][j]    = u[ny-2][j];
			u[ny-1][j] = u[1][j];
			v[0][j]    = v[ny-2][j];
			v[ny-1][j] = v[1][j];
		}
	}
	else  // Update the boundaries of un & vn
	{
		for (i = 0; i < ny; i++)
		{
			un[i][nx_local-1] 	= u_r[i] ;
			un[i][0] 			= u_l[i];
			vn[i][nx_local-1] 	= v_r[i] ;
			vn[i][0] 			= v_l[i];
		}
		
		for (j = 0; j < nx_local; j++)
		{
			un[0][j]    	= un[ny-2][j];
			un[ny-1][j] 	= un[1][j];
			vn[0][j]    	= vn[ny-2][j];
			vn[ny-1][j] 	= vn[1][j];
		}
	}
}

void * solve_equation(int type)
{
	int i,j;
	if (type == 1) // Solve with u & v as input
	{
		for (i = 1; i < ny -1; i ++)
			for (j = 1; j < nx_local - 1; j++)
				{
					un[i][j] = (u[i][j] - dt/dx*u[i][j]*(u[i][j] - u[i][j-1]) 
					 	      - dt/dy*v[i][j]*(u[i][j] - u[i-1][j]) 
					 	      + nu * dt / (dx*dx) * (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) 
					 	      + nu * dt / (dy*dy) * (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) );
				
					vn[i][j] = (v[i][j] - dt/dx*u[i][j]*(v[i][j] - v[i][j-1])
							  - dt/dy*v[i][j]*(v[i][j] - v[i-1][j]) 
							  + nu * dt / (dx*dx) * (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]) 
							  + nu * dt / (dy*dy) * (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j]) );
				}
	}
	else // Solve with un & vn as input
	{
		for (i = 1; i < ny -1; i ++)
			for (j = 1; j < nx_local - 1; j++)
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
		
	}
}


void * setup_MPI_Domain(void *rank)
{
	int myID =  *((int*) rank); // Real ID, we be used to move along x
	int i, j;
	nx_local = (int) (nx/MPIthreads);

	if (myID ==0)
		{
			printf("Divide the x-direction into %d divisions.\n", MPIthreads);    
		}
	
	if(MPIthreads < 2)
		printf("YOOO ERROR!! Don't use this low number of cores.\n");

	// Diving the domain 
	j_start = nx_local*myID;
	j_end   = j_start+nx_local;
	neig_r = myID+1; // Neighbors to the right
	neig_l = myID-1; // Neighbors to the left

	//Correction by adding the ghost cells
	if(myID == 0) // First Rank
	{ 
		neig_l = MPIthreads - 1;
	}
	else if (myID == MPIthreads - 1) // Last Rank
	{
		neig_r = 0;
		j_start--;
	}
	else
	{
		j_start--;
	}

	nx_local = j_end - j_start + 1; // Correct number of points in each rank

	printf("Rank:%d,Neighbours r:%d,l:%d,nx:%d, x:%d,%d\n",myID,neig_r,neig_l,nx_local,j_start,j_end);

	x_local = (double *)malloc(nx_local*sizeof(double));
	for(j = j_start; j <=j_end ; j++)
		x_local[j-j_start] = x[j];

	// Creating arrays that will be sent/recieved by the MPI
	u_r = (double *)malloc(ny*sizeof(double));
	u_l = (double *)malloc(ny*sizeof(double));
	v_r = (double *)malloc(ny*sizeof(double));
	v_l = (double *)malloc(ny*sizeof(double));

}

void * send_boundary_data(void *Request, void * rank, int type)
{
	int myID =  *((int*) rank); // Real ID, we be used to move along x

	if(type == 1)
	{
		// Send Right and Left Boundary of u & v
		int i;
		for (i = 0; i < ny; i++)
		{
			u_r[i] = u[i][nx_local-2];
			u_l[i] = u[i][1];
			v_r[i] = v[i][nx_local-2];
			v_l[i] = v[i][1];
		}
		MPI_Isend(u_r, ny, MPI_DOUBLE, neig_r, 6660, MPI_COMM_WORLD, Request);
		MPI_Isend(v_r, ny, MPI_DOUBLE, neig_r, 6661, MPI_COMM_WORLD, Request);
		MPI_Isend(u_l, ny, MPI_DOUBLE, neig_l, 4440, MPI_COMM_WORLD, Request);
		MPI_Isend(v_l, ny, MPI_DOUBLE, neig_l, 4441, MPI_COMM_WORLD, Request);

	}
	else
	{
		// Send Right and Left Boundary of un & vn
		int i;
		for (i = 0; i < ny; i++)
		{
			u_r[i] = un[i][nx_local-2];
			u_l[i] = un[i][1];
			v_r[i] = vn[i][nx_local-2];
			v_l[i] = vn[i][1];
		}
		MPI_Isend(u_r, ny, MPI_DOUBLE, neig_r, 6660, MPI_COMM_WORLD, Request);
		MPI_Isend(v_r, ny, MPI_DOUBLE, neig_r, 6661, MPI_COMM_WORLD, Request);
		MPI_Isend(u_l, ny, MPI_DOUBLE, neig_l, 4440, MPI_COMM_WORLD, Request);
		MPI_Isend(v_l, ny, MPI_DOUBLE, neig_l, 4441, MPI_COMM_WORLD, Request);
	}
}

void * recieve_boundary_date(void *Status,void *Request, void * rank)
{

	int myID =  *((int*) rank); // Real ID, we be used to move along x

	int flag = 0; 
	// recieve from left
	while(!flag)
		MPI_Iprobe(neig_l, 6660, MPI_COMM_WORLD, &flag, Status );
	MPI_Irecv(u_l, ny, MPI_DOUBLE, neig_l, 6660, MPI_COMM_WORLD, Request);	

		flag = 0; 
	while(!flag)
		MPI_Iprobe(neig_l, 6661, MPI_COMM_WORLD, &flag, Status );
	MPI_Irecv(v_l, ny, MPI_DOUBLE, neig_l, 6661, MPI_COMM_WORLD, Request);	
	
	// recieve from right
	flag = 0; 
	while(!flag)
		MPI_Iprobe(neig_r, 4440, MPI_COMM_WORLD, &flag, Status );
	MPI_Irecv(u_r, ny, MPI_DOUBLE, neig_r, 4440, MPI_COMM_WORLD, Request);	
	
	flag = 0;  // ERROR HERE INVALID RANK
	
	while(!flag)
		MPI_Iprobe(neig_r, 4441, MPI_COMM_WORLD, &flag, Status );
	MPI_Irecv(v_r, ny, MPI_DOUBLE, neig_r, 4441, MPI_COMM_WORLD, Request);	

}
#endif