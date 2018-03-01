// Class definition
#ifndef MATRIXCELLS
#define MATRIXCELLS

// vector and vector for pointer vectors
#include <vector>

// Common header
#include "common.h"

//New namespace definition
namespace matrixMapp {

    // Main class: Matrix containing the mesh for the particles
    class matrixCells {

    // Public elements
    public:
        // Utils for iterating through the matrix's cells: new helper class
        // New class is defined as a friend for accessing private data of the
        // main class
        friend class matrixIter;

        // "Sub-Class" definition: vector of pointers
        class matrixIter: std::iterator_traits<particle_t *> {

            // Main class is a friend
            friend class matrixCells;

            // Public components
            public:
                // Constructor
                matrixIter();
                matrixIter(matrixCells *, particle_t &, bool end = false);

                // Destructor
                ~matrixIter();

                // Flag operators for checking references
                matrixIter & operator++();
                pointer operator*() const;
                bool operator==(matrixIter const &);
                bool operator!=(matrixIter const &);

            // Private components
            private:
                // Pointer to mesh
                matrixCells * mmesh;

                // Boolean flag: TRUE if end of the vector is reached
                bool meshEnd;

                // Rows and columns
                int r, c;

                // Coordinates: upper left, bottom right
                int ulr, ulc, brr, brc;

                // Particle iterator (from vector)
                std::vector<particle_t *>::iterator pit;

                // Particle iterator: last element (from vector)
                std::vector<particle_t *>::iterator pit_last;

                // Next iterator function: updates the outer iteration
                void updateIter();

                // Particle from last cell
                matrixIter& lastCell();
        };

    // Matrix cells class constructor
    matrixCells(int n, double size, double cutoff_radius);
    matrixCells(int n, double size, double cutoff_radius, int grid_len_fac);

    // Destructor
    ~matrixCells();

    // Public Methods
    // Clears the matrix/mesh
    void clear();

    // Insert particles
    void insert(particle_t &);

    // Delete particles
    void remove(particle_t &, int);

    // Adjacent indexes for specific sequence containing the relevant particle
    std::vector<particle_t *>::iterator InitAdjPart(particle_t &);
    std::vector<std::vector<particle_t *>*>::iterator InitAdj();
    std::vector<std::vector<particle_t *>*>::iterator EndAdj();

    // Objects: initial and endint
    matrixIter AdjInitial(particle_t &);
    matrixIter AdjEnding(particle_t &);

    // Get methods: coordinates of a particle in rows and colums
    int getRow(particle_t &);
    int getCol(particle_t &);

    int get_index(particle_t &);
    int get_rows();
    int get_cols();
    int get_adj();
    
    void clear_fringes(int proc_n);
    int get_row_offset(int proc_n);
    int get_proc_rows();
    bool owns_particle(particle_t& p, int proc_n);

    // Private elements
    private:
        // Rows and columns
        int nrows;
        int ncols;
        int rows_per_proc;

        // Size of the matrix
        double msize;

        // Adjacent cells
        int Nadjacents;
        std::vector<std::vector<particle_t *>*>* cells;

        // Methods
        //
        void ULadj(int&, int&, particle_t&);
        void BRadj(int&, int&, particle_t&);

        // Get the index of a certain particle based on row/col or pointer
        int get_index(int row, int col);
    };
}

#endif
