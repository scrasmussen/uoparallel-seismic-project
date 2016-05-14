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
//
// Data:
//
//   struct FLOATBOX:
//      efficiently stores a 3D float volume in flat memory
//
//
// Functions:
//
//   int boxalloc( struct FLOATBOX*, int nx, int ny, int nz )
//   void boxsetall( struct FLOATBOX, float value )
//   void boxfree( struct FLOATBOX* )
//   float boxget( struct FLOATBOX, int x, int y, int z )
//   void boxput( struct FLOATBOX, int x, int y, int z, float )
//   void boxfprint( FILE*, char *prefix, char *indent, struct FLOATBOX )
//
//
// Example program:
//
//   #include "floatbox.h"
//   #include <stdio.h> /* for printf and stdout, not necessary otherwise */
//   int main() {
//     struct FLOATBOX box;
//     if( !boxalloc( &box, 241, 241, 51 ) ) { /* handle error */ }
//     boxfprint( stdout, "example: ", "\t", box );
//     boxput( box, 1, 2, 3, 4.567f );
//     printf( "value at (1,2,3): %g\n", boxget( box, 1, 2, 3 ) );
//     boxfree( &box );
//   }
//
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

void
boxinit (
    struct FLOATBOX *box
)
// initializes box to safe default values
{
    if( box == NULL ) return;

    box->sx = 0;
    box->sy = 0;
    box->sz = 0;

    box->nx = 0;
    box->ny = 0;
    box->nz = 0;

    box->flat = NULL;
}


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
// releases heap memory associated with box
{
    if( box == NULL ) return;
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
boxsetall (
    const struct FLOATBOX box,
    float val
)
// sets ALL values in the volume to the given value
{
    if( box.flat == NULL ) return;
    size_t i, last = boxindex( box, box.nx-1, box.ny-1, box.nz-1 );
    for( i = 0; i <= last; i++ ) box.flat[i] = val;
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
