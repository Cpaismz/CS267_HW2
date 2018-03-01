#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_set>
#include <unistd.h>
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
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );


    // invariant: at any given time, only the values within one mesh are valid

    // we have to broadcast particles because we want a consistent initial state
    if( rank == 0 ) {
        init_particles( n, particles );
    }
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    matrixMapp::matrixCells* mesh = new matrixMapp::matrixCells(n, size, cutoff, n_proc);
    // a set of pointers... each to be individually malloced
    // my code is disgusting
    std::unordered_set<particle_t*> owned;

    // filters a list of particles into a set to back the memory and adds that list to the grid
    push2Set(n, particles, mesh, owned, rank);


    int gridLen = mesh->get_cols();
    int* migrated_sizes = (int*)malloc(n_proc * sizeof(int));
    int* disp_sizes = (int*)malloc(n_proc * sizeof(int));
    //
    //  simulate a number of time steps
    //
    int size;
    double simulation_time = read_timer( );
    char a[30];
    sprintf(a, "outfiles/r%d.out", rank);
    FILE* fptr = fopen(a, "w");

    dup2(fileno(fptr), fileno(stdout));
    for( int step = 0; step < NSTEPS; step++ )
    {
        printf("a"); fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        size = owned.size();
        MPI_Allgather(&size, 1, MPI::INT, migrated_sizes, 1, MPI::INT, MPI_COMM_WORLD);
        int sum = 0;
        for (int i = 0; i < n_proc; i++) {
            sum += migrated_sizes[i];
        }
        if (sum != n) {
            std::cerr << "wrong size" << std::endl;
        }

        printf("b"); fflush(stdout);
        //Correctness
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
            if( rank == 0 && fsave && (step % SAVEFREQ) == 0 )
                save( fsave, n, particles );

        particle_t*  halo_buf = (particle_t*)malloc(sizeof(particle_t) * owned.size());
        // only send the number of necessary rows, as opposed to all
        int i = 0;
        for (auto & a : owned) {
            halo_buf[i] = *a;
            if (mesh->get_owner(halo_buf[i]) != rank) {
                
                printf("sn: %d ", mesh->get_owner(halo_buf[i])); fflush(stdout);
            }
            i++;
        }

        MPI_Request* trash = (MPI_Request*)malloc(n_proc * sizeof(MPI_Request));
        for (int i = 0; i < mesh->get_adj(); i += mesh->get_proc_rows()) {
            
            int left_addr = rank - i - 1;
            int right_addr = rank + i + 1;
            if (left_addr >= 0) {
                printf("send %d ", (int)owned.size());
                MPI_Isend(halo_buf, owned.size(), PARTICLE, left_addr, 0, MPI_COMM_WORLD, &(trash[left_addr]));
            }
            if (right_addr < n_proc) {
                printf("send %d ", (int)owned.size());
                MPI_Isend(halo_buf, owned.size(), PARTICLE, right_addr, 0, MPI_COMM_WORLD, &(trash[right_addr]));
            }
        }
        free(trash);

        printf("c"); fflush(stdout);
        particle_t*  rec_buf = (particle_t*)malloc(sizeof(particle_t) * n);
        int offset = 0;
        for (int i = 0; i < mesh->get_adj(); i += mesh->get_proc_rows()) {
            int left_addr = rank - i - 1;
            int right_addr = rank + i + 1;
            if (left_addr >= 0) {
                MPI_Status stat;
                MPI_Recv(rec_buf + offset, n, PARTICLE, left_addr, 0, MPI_COMM_WORLD, &stat);
                int num_rec;
                MPI_Get_count(&stat, PARTICLE, &num_rec);
                printf("recv %d ", num_rec);
                offset += num_rec;
            }
            if (right_addr < n_proc) {
                MPI_Status stat;
                MPI_Recv(rec_buf + offset, n, PARTICLE, right_addr, 0, MPI_COMM_WORLD, &stat);
                int num_rec;
                MPI_Get_count(&stat, PARTICLE, &num_rec);
                printf("recv %d ", num_rec);
                offset += num_rec;
            }
        }

        printf("%d %d", offset, n); fflush(stdout);
        push2Mesh(offset, rec_buf, mesh);
        //
        //  compute all forces
        //

        printf("d"); fflush(stdout);
        for (auto & part : owned) {
            part->ax = part->ay = 0;

            // Only check the neighbors of the current particle: at most 8 cells
            matrixMapp::matrixCells::matrixIter adjIter;
            for (adjIter = mesh->AdjInitial(*part); adjIter != mesh->AdjEnding(*part); ++adjIter) {
                // calling function by modifying non-const ref to part :(
                apply_force(*part, **adjIter, &dmin, &davg, &navg);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        mesh->clear_fringes(rank);
        free(rec_buf);
        free(halo_buf);

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
        std::vector<particle_t *> to_del;

        for (auto & part : owned) {
            int old_index = mesh->get_index(*part);
            move( *part );
            int new_index = mesh->get_index(*part);
            if (old_index != new_index) {
                mesh->remove(*part, old_index);
                if (mesh->owns_particle(*part, rank)) {
                    mesh->insert(*part);
                } else {
                    to_del.push_back(part);

                }
            }
        }


        std::vector<particle_t> migrated;
        for (auto a : to_del) {
            migrated.push_back(*a);
            owned.erase(a);
            free(a);
        }
        // TODO: this is still an all-to-all broadcast, but only of the particles that moved
        printf("e"); fflush(stdout);

        size = migrated.size();
        MPI_Allgather(&size, 1, MPI::INT, migrated_sizes, 1, MPI::INT, MPI_COMM_WORLD);

        printf("f"); fflush(stdout);
        disp_sizes[0] = 0;
        for (int i = 1; i < n_proc; i++) {
            disp_sizes[i] = disp_sizes[i-1] + migrated_sizes[i-1];
        }

        int tot_migrated = 0;
        for (int i = 0; i < n_proc; i++) {
            tot_migrated += migrated_sizes[i];
        }

        MPI_Allgatherv(&(migrated[0]), migrated.size(), PARTICLE, particles, migrated_sizes, disp_sizes, PARTICLE, MPI_COMM_WORLD );

        printf("g"); fflush(stdout);
        push2Set(tot_migrated, particles, mesh, owned, rank);
    }
    simulation_time = read_timer( ) - simulation_time;

    for (auto a : owned) {
        free(a);
    }

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
    free( particles );
    free(migrated_sizes);
    free(disp_sizes);

    if( fsave )
        fclose( fsave );

    MPI_Finalize( );

    return 0;
}
