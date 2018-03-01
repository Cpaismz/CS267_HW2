#include "matrixCells.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

// Definition of the newspace associated with the matrix/grid
namespace matrixMapp {

    ///////////////////////////////////
    //
    //          Main methods
    //
    ///////////////////////////////////


    // Constructor: a new matrix containing cells with the corresponding number of
    // particles and size is initialized
    matrixCells::matrixCells(int n, double size, double cutoff_radius) :
            // Main parameters are initialized
            // Number of rows of the matrix
            nrows(ceil(sqrt(n))),

            // Number of columns
            ncols(ceil(sqrt(n))),

            // Total size of the matrix
            msize(size / nrows),

            rows_per_proc(-1),

            // Number of adjacent/neighbor cells
            Nadjacents(ceil(cutoff_radius / msize)) {

                // New vector vector with references
                cells = new std::vector<std::vector<particle_t *>*>(nrows * ncols);

                // Iterator is initialized
                std::vector<std::vector<particle_t *>*>::iterator it;

                // Add references
                for (it = cells->begin(); it != cells->end(); ++it) {
                    *it = new std::vector<particle_t *>();
                }
            }


    matrixCells::matrixCells(int n, double size, double cutoff_radius, int grid_len_fac) {
        // for simplicity, we want the side of a grid to be a multiple of grid_len_fac (MPI)
        int grid_len = grid_len_fac;
        while (grid_len * grid_len < n) {
            grid_len += grid_len_fac;
        }

        // Number of rows of the matrix
        nrows = grid_len;

        // Number of columns
        ncols = grid_len;

        // Number of rows per processor
        rows_per_proc = nrows / grid_len_fac;

        // Total size of the matrix
        msize = size / nrows;

        // Number of adjacent/neighbor cells
        Nadjacents = ceil(cutoff_radius / msize);

        // New vector vector with references
        cells = new std::vector<std::vector<particle_t *>*>(nrows * ncols);

        // Iterator is initialized
        std::vector<std::vector<particle_t *>*>::iterator it;

        // Add references
        for (it = cells->begin(); it != cells->end(); ++it) {
            *it = new std::vector<particle_t *>();
        }
    }

    // Destructor method: delete the matrix's cells
    matrixCells::~matrixCells() {
        // Delete references (pointers to particles)
        std::vector<std::vector<particle_t *>*>::iterator it;
        for (it = cells->begin(); it != cells->end(); ++it) {
            delete *it;
        }

        // Delete the cells
        delete cells;
    }


    // Particles' pointers are eliminated from the matrix mesh
    void matrixCells::clear() {
         // Iterates over all the current pointers and delete them
        std::vector<std::vector<particle_t *>*>::iterator it;
        for (it = cells->begin(); it != cells->end(); ++it) {
            (*it)->clear();
        }
    }


    // Add a new particle via pointer inside the existing matrix via push_back function
    void matrixCells::insert(particle_t & p) {
        (*cells)[get_index(p)]->push_back(&p);
    }


    // Remove particle from old index bin
    void matrixCells::remove(particle_t & p, int index) {
        std::vector<particle_t *>* old_bin = (*cells)[index];
        particle_t * addr = &p;

        // maybe try swap too
        // https://stackoverflow.com/a/3385251
        old_bin->erase(std::remove(old_bin->begin(), old_bin->end(), addr), old_bin->end());
    }

    void matrixCells::clear_fringes(int proc_n) {
        for (int i = 0; i < Nadjacents; i++) {
            int left_row_index = proc_n * rows_per_proc - i - 1;
            int right_row_index = (proc_n + 1) * rows_per_proc + i;
            for (int j = 0; j < ncols; j++) {
                if (left_row_index >= 0) {
                    (*cells)[ncols * left_row_index + j]->clear();
                }
                if (right_row_index < nrows) {
                    (*cells)[ncols * right_row_index + j]->clear();
                }
            }
        }
    }

    // Adjacent indexes for specific sequence containing the relevant particle
    // Particle
    std::vector<particle_t *>::iterator matrixCells::InitAdjPart(particle_t & p) {
        return (*cells)[get_index(p)]->begin();
    }

    // Initial index
    std::vector<std::vector<particle_t *>*>::iterator matrixCells::InitAdj() {
        return cells->begin();
    }

    // Final index
    std::vector<std::vector<particle_t *>*>::iterator matrixCells::EndAdj() {
        return cells->end();
    }

    // Object: beginning
    matrixCells::matrixIter matrixCells::AdjInitial(particle_t & p) {
        return matrixIter(this, p);
    }

    // Object: ending
    matrixCells::matrixIter matrixCells::AdjEnding(particle_t & p) {
        return matrixIter(this, p).lastCell();
    }


    // Update the coordinates of adjacent cells: upper left to bottom right
    void matrixCells::ULadj(int& tlr, int& tlc, particle_t& p) {
        tlr = max(getRow(p) - Nadjacents, 0);
        tlc = max(getCol(p) - Nadjacents, 0);
    }

    void matrixCells::BRadj(int& brr, int& brc, particle_t& p) {
        brr = min(getRow(p) + Nadjacents, nrows - 1);
        brc = min(getCol(p) + Nadjacents, ncols - 1);
    }


    // Get the coordinate of the particle given row and column numbers
    int matrixCells::get_index(int row, int col) {
        return col + row * ncols;
    }


    // Get the coordinate of particle p inside the mesh
    int matrixCells::get_index(particle_t & p) {
        return get_index(getRow(p), getCol(p));
    }


    // Get the row value of particle p inside the mesh
    int matrixCells::getRow(particle_t & p) {
        int row = static_cast<int>(p.y / msize);
        assert(row < nrows);
        return row;
    }

    int matrixCells::get_adj() {
        return Nadjacents;
    }

    int matrixCells::get_rows() {
        return nrows;
    }

    int matrixCells::get_cols() {
        return ncols;
    }

    int matrixCells::get_row_offset(int proc_n) {
        return proc_n * rows_per_proc;
    }

    int matrixCells::get_proc_rows() {
        return rows_per_proc;
    }

    bool matrixCells::owns_particle(particle_t& p, int proc_n) {
        int p_row = getRow(p);
        return p_row >= get_row_offset(proc_n) && p_row < get_row_offset(proc_n) + rows_per_proc;
    }

    // Get the column value of particle p inside the mesh
    int matrixCells::getCol(particle_t & p) {
        int col = static_cast<int>(p.x / msize);
        assert(col < ncols);
        return col;
    }



    ///////////////////////////////////
    //
    //        Adjacent methods
    //
    ///////////////////////////////////

    // Empty constructor for avoiding initialization problems and completeness
    matrixCells::matrixIter::matrixIter() {}

    // Main constructor method: mesh bounds are calculated
    matrixCells::matrixIter::matrixIter(matrixCells* mmesh, particle_t& p, bool end) : mmesh(mmesh), meshEnd(end) {

        // Adjacent indices are updated
        mmesh->ULadj(ulr, ulc, p);
        mmesh->BRadj(brr, brc, p);

        // Rows and columns are initialized
        r = ulr;
        c = ulc;

        // vector vector is initialized
        pit = mmesh->cells->at(mmesh->get_index(r, c))->begin();
        pit_last = mmesh->cells->at(mmesh->get_index(r, c))->end();

        // Update the iteration value
        updateIter();
    }


    // Destructor
    matrixCells::matrixIter::~matrixIter() {}

    // References to vector: point to valid references
    matrixCells::matrixIter & matrixCells::matrixIter::operator++() {
        // vector reference is increased. Next iteration.
        ++pit;
        updateIter();
        return *this;
    }

    matrixCells::matrixIter::pointer matrixCells::matrixIter::operator*() const {
        // vector reference is dereferenced
        return *pit;
    }

    bool matrixCells::matrixIter::operator==(matrixIter const & other) {
        // End of the mesh or match between index
        return meshEnd == other.meshEnd && (meshEnd || (r == other.r && c == other.c && pit == other.pit));
    }

    bool matrixCells::matrixIter::operator!=(matrixIter const & other) {
        // Check if different
        return !((*this) == other);
    }


    // Push the iteration index to the last entry
    void matrixCells::matrixIter::updateIter() {
        // If pointer is useless: repeat
        while ((pit == pit_last || *pit == NULL) && !meshEnd) {
            // update the mesh: columns and rows
            // columns
            if (++c > brc) {
                c = ulc;
                ++r;
            }

            // rows
            if (r > brr) {
                meshEnd = true;
                return;
            }

            // Otherwise, indexes are updated for the vector
            else {
                pit = mmesh->cells->at(mmesh->get_index(r, c))->begin();
                pit_last = mmesh->cells->at(mmesh->get_index(r, c))->end();
            }
        }
    }


    // Particle from last cell
    matrixCells::matrixIter& matrixCells::matrixIter::lastCell() {
        // Next mesh
        r = brr + 1;
        c = ulc;

        // Reference from last particle
        meshEnd = true;
        return (*this);
    }

}

