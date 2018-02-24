#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "matrixCells.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    // Correctness variables
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );


    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    //
    //  allocate generic resources
    //
    //  Files (if asked) are generated by the rank = 0 process
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    // A MPI Datatype PARTICLE is defined containing the six values particle_t.
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );

    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];

    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    if( rank == 0 ) {
        init_particles( n, particles );
    }
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    matrixMapp::matrixCells* mesh = new matrixMapp::matrixCells(n, size, cutoff);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        //Correctness
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //

        // Insert the particles into the matrix
        mesh->clear();
        push2Mesh(n, particles, mesh);

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
            if( rank == 0 && fsave && (step % SAVEFREQ) == 0 )
                save( fsave, n, particles );

        //
        //  compute all forces
        //
        for (int i = partition_offsets[rank]; i < partition_offsets[rank + 1]; i++)
        {
            particles[i].ax = particles[i].ay = 0;

            // Only check the neighbors of the current particle: at most 8 cells
            matrixMapp::matrixCells::matrixIter adjIter;
            for (adjIter = mesh->AdjInitial(particles[i]); adjIter != mesh->AdjEnding(particles[i]); ++adjIter) {
                apply_force(particles[i], **adjIter, &dmin, &davg, &navg);
            }
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {

          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);


          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }


        //
        //  move particles
        //
        for (int i = partition_offsets[rank]; i < partition_offsets[rank + 1]; i++){
            int old_index = mesh->get_index(particles[i]);
            move( particles[i] );
            int new_index = mesh->get_index(particles[i]);
            if (old_index != new_index) {
                mesh->remove(particles[i], old_index);
                mesh->insert(particles[i]);
            }
        }

        // A distribution of the current work is performed.
        particle_t* local = &particles[partition_offsets[rank]];
        MPI_Allgatherv(local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );

    }
    simulation_time = read_timer( ) - simulation_time;

    if (rank == 0) {
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

      if( find_option( argc, argv, "-no" ) == -1 )
      {
        if (nabsavg) absavg /= nabsavg;
      //
      //  -The minimum distance absmin between 2 particles during the run of the simulation
      //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
      //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
      //
      //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
      //
      printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
      if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
      if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");

      //
      // Printing summary data
      //
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
    //
    //  release resources
    //
    delete mesh;
    if ( fsum )
        fclose( fsum );
    free( partition_offsets );
    free( partition_sizes );
    free( particles );

    if( fsave )
        fclose( fsave );

    MPI_Finalize( );

    return 0;
}