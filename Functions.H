#ifndef Functions_H
#define Functions_H

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
		for ( j = 0; j < nx_local; j++ ) 
		{
			array[i][j] = 1;
		}
	}

}

void * create_arrays()
{
	x = (double *)malloc(nx*sizeof(double));
	y = (double *)malloc(ny*sizeof(double));
}

void * create_matrices()
{
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
		u[i]    = (double *)malloc(sizeof(double)*nx_local);
		v[i]    = (double *)malloc(sizeof(double)*nx_local);
		un[i]   = (double *)malloc(sizeof(double)*nx_local);
		vn[i]   = (double *)malloc(sizeof(double)*nx_local);
	}
}



void *write_boundary_data()
{
	char filename[20];
	snprintf(filename,20,"output/boundary");

	fd = fopen(filename,"w");
	fprintf(fd, "%lf,%lf,%d,%d,%d,%lf,%lf,%d,%d,\n",Lx,Ly,nx,ny,nt,sigma,nu,auto_save, MPIthreads);
	fclose( fd );
}

/*
void *write_domain_info(void *rank)
{
	int myID =  *((int*) rank); 
	int i;
	char filename[20];

	snprintf(filename,20,"output/domain_x_%d",myID);
	fd = fopen(filename,"w");
	for ( i = 0; i < nx_local; i++) 
	{
		fprintf(fd, "%d,%lf\n",i,x_local[i]);
	}
	fclose( fd );

	snprintf(filename,20,"output/domain_y_%d",myID);
	fd = fopen(filename,"w");
	for ( i = 0; i < ny_local; i++) 
	{
		fprintf(fd, "%d,%lf\n",i,y_local[i]);
	}
	fclose( fd );

}*/

void *write_output(int time, void *rank)
{
	int i,j;
	int myID =  *((int*) rank); // Real ID, we be used to move along x

	char filename1[20], filename2[20];
	snprintf(filename1,20,"output/u-%d-%d",time,myID);
	snprintf(filename2,20,"output/v-%d-%d",time,myID);

	fd = fopen(filename1, "w" );
	for ( i = 0; i < ny; i++) 
	{
		for ( j = 0; j < nx_local; j++ ) 
		{
			fprintf(fd, "%lf",u[i][j]);
			if (j != nx_local-1) 	fprintf(fd, ",");
		}
		fprintf(fd,"\n");
	}
	fclose( fd );

	fd = fopen(filename2,"w" );
	for ( i = 0; i < ny; i++) 
	{
		for ( j = 0; j < nx_local; j++ ) 
		{
			fprintf(fd, "%lf",v[i][j]);
			if (j != nx_local-1) 	fprintf(fd, ",");
		}
		fprintf(fd,"\n");
	}
	fclose( fd );
}


void * copy_arrays()
{
	int i,j;
	for ( i = 0; i < ny; i++) 
		for ( j = 0; j < nx_local; j++ ) 
		{
			un[i][j] = u[i][j];
			vn[i][j] = v[i][j];
		}
}

#endif