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

#include <stdio.h>
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
    struct FLOATBOX *box,
    int nx,
    int ny,
    int nz
)
// allocates box->flat and sets strides and dimensions appropriately;
// on error: returns 0
// on success: returns non-zero
{
    size_t numbytes = (size_t)nx * ny * nz * sizeof(float);
    float *flat = malloc( numbytes );

    if( flat == NULL ) {
        // system couldn't allocate the memory
        return 0;
    }

    box->sx = (size_t)ny * nz;
    box->sy = nz;
    box->sz = 1;

    box->nx = nx;
    box->ny = ny;
    box->nz = nz;

    box->flat = flat;

    return 1;
}


void
boxfree (
    struct FLOATBOX *box
)
// frees box->flat and sets it to NULL
{
    free( box->flat );
    box->flat = NULL;
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
    struct FLOATBOX box,
    int x,
    int y,
    int z
)
// returns a single value from the given coordinates
{
    return box.flat[ boxindex( box, x, y, z ) ];
}


inline extern
void
boxput (
    struct FLOATBOX box,
    int x,
    int y,
    int z,
    float val
)
// stores a single value at the given coordinates
{
    box.flat[ boxindex( box, x, y, z ) ] = val;
}


void
boxfprint (
    FILE *stream,
    const char *prefix,
    const char *indent,
    struct FLOATBOX box
)
// metadata friendly-printing
{
    if( stream == NULL ) stream = stdout;
    if( prefix == NULL ) prefix = "";
    if( indent == NULL ) indent = "  ";

    fprintf( stream, "%sFLOATBOX {\n", prefix );
    fprintf( stream, "%s%ssx: %zu, sy: %zu, sz %zu\n",
        prefix, indent, box.sx, box.sy, box.sz );
    fprintf( stream, "%s%snx: %d, ny: %d, nz: %d\n",
        prefix, indent, box.nx, box.ny, box.nz );
    fprintf( stream, "%s%sflat: %p\n",
        prefix, indent, (void*)box.flat );
    fprintf( stream, "%s}\n", prefix );
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
