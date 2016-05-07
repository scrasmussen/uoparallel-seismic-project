////////////////////////////////////////////////////////////////////////////////
// floatbox.h - 2016.05.06 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// build-tested using:
//   gcc 5.2
//   -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline
// report:
//   * all functions successfully inlined in get/put test
//   * boxindex() correctly computes indices when (x*y*z) > INT_MAX
// note:
//   For inlining purposes, these functions are fleshed out in this header file,
//   instead of using a separate compilation unit.
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct FLOATBOX {
    size_t sx, sy, sz; // array strides
    int nx, ny, nz; // dimensions
    float *flat; // [x][y][z] order
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
boxalloc (
    struct FLOATBOX *boxdst,
    int nx,
    int ny,
    int nz
)
// allocates boxdst->flat and sets strides and dimensions appropriately;
// on error: returns 0
// on success: returns non-zero
{
    size_t numbytes = (size_t)nx * ny * nz * sizeof(float);
    float *flat = malloc( numbytes );

    if( flat == NULL ) {
        // system couldn't allocate the memory
        return 0;
    }

    boxdst->sx = (size_t)ny * nz;
    boxdst->sy = nz;
    boxdst->sz = 1;

    boxdst->nx = nx;
    boxdst->ny = ny;
    boxdst->nz = nz;

    boxdst->flat = flat;

    return 1;
}


void
boxfree (
    struct FLOATBOX *boxdst
)
// frees boxdst->flat and sets it to NULL
{
    free( boxdst->flat );
    boxdst->flat = NULL;
}


inline extern
size_t
boxindex (
    struct FLOATBOX box,
    int x,
    int y,
    int z
)
// return an index into box.flat corresponding to the given coordinates
{
    return x * box.sx + y * box.sy + z * box.sz;
}


inline extern
float
boxget (
    struct FLOATBOX boxsrc,
    int x,
    int y,
    int z
)
// returns a single value from the given coordinates
{
    return boxsrc.flat[ boxindex( boxsrc, x, y, z ) ];
}


inline extern
void
boxput (
    struct FLOATBOX boxdst,
    int x,
    int y,
    int z,
    float val
)
// stores a single value at the given coordinates
{
    boxdst.flat[ boxindex( boxdst, x, y, z ) ] = val;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
