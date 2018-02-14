#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "common.h"
#include "matrixCells.h"
#include "omp.h"


#define min( i, j ) ( (i)<(j) ? (i): (j) )
#define max( i, j ) ( (i)>(j) ? (i): (j) )
//
//  benchmarking program
//
int main( int argc, char **argv )
{
    /*
    omp_lock_t write_lock;
    omp_init_lock(&write_lock);*/

    int navg,nabsavg=0,numthreads;
    double dmin, absmin=1.0,davg,absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    matrixMapp::matrixCells* mesh = new matrixMapp::matrixCells(n, size, cutoff);
    push2Mesh(n, particles, mesh);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    int max_locks = 320;
    int num_entries = min(mesh->get_cols() * mesh->get_rows(), max_locks);
    omp_lock_t* locks = (omp_lock_t*)malloc(sizeof(omp_lock_t) *  num_entries);
    #pragma omp parallel private(dmin)
    {
        #pragma omp for
        for (int i = 0; i < num_entries; i++) {
            omp_init_lock(&locks[i]);
        }
        numthreads = omp_get_num_threads();
    for( int step = 0; step < NSTEPS ; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;

            // Only check the neighbors of the current particle: at most 8 cells
            matrixMapp::matrixCells::matrixIter adjIter;
            matrixMapp::matrixCells::matrixIter adjIterEnd = mesh->AdjEnding(particles[i]);

            for (adjIter = mesh->AdjInitial(particles[i]); adjIter != adjIterEnd; ++adjIter) {
                apply_force(particles[i], **adjIter,&dmin,&davg,&navg);
            }
        }


        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ){
            int old_index = mesh->get_index(particles[i]);
            move( particles[i] );
            int new_index = mesh->get_index(particles[i]);
            if (old_index != new_index) {
                int lock_old_index = old_index * num_entries / (mesh->get_rows() * mesh->get_cols());
                omp_set_lock(&locks[lock_old_index]);
                mesh->remove(particles[i], old_index);
                omp_unset_lock(&locks[lock_old_index]);

                int lock_new_index = new_index * num_entries / (mesh->get_rows() * mesh->get_cols());
                omp_set_lock(&locks[lock_new_index]);
                mesh->insert(particles[i]);
                omp_unset_lock(&locks[lock_new_index]);
            }

            /* 
            fine grained locking?
            1) cache location of old particle
            2) move particle
            3) get new location
            4) if bucket_changed:
                // we can either used fine grained locking or overall locking
                5) lock old bucket, remove particle, unlock
                6) lock new bucket, add particle, unlock
            */
        }


        // Update grid hash set.
       /* 
        #pragma omp single
        {
            mesh->clear();
            push2Mesh(n, particles, mesh);
        }*/

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) {
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }

    }
        #pragma omp for
        for (int i = 0; i < num_entries; i++) {
            omp_destroy_lock(&locks[i]);
        }

}

        free(locks);
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
