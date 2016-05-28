////////////////////////////////////////////////////////////////////////////////
// velocitybox.h - 2016.05.26 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
// 
// A container for a regular-grid velocity volume using global coordinates.
//
// Requires:
//   floatbox.h
//
// Data:
//   struct VELOCITYBOX
//
// Functions:
//   vboxalloc
//   vboxfprint
//   vboxfree
//   vboxinit
// 
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "floatbox.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////
// structs 
////////////////////////////////////////////////////////////////////////////////

struct VELOCITYBOX {
    struct POINT3D min, max;
    struct FLOATBOX box; // contains dimensions and velocity data
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

void
vboxinit (
    struct VELOCITYBOX *vbox
)
// initializes vbox to safe default values
{
    if( vbox == NULL ) return;

    point3dset( &vbox->min, 0, 0, 0 );
    point3dset( &vbox->max, 0, 0 ,0 );

    boxinit( &vbox->box );
}


int
vboxalloc (
    struct VELOCITYBOX *vbox,
    const int ox, const int oy, const int oz,
    const int nx, const int ny, const int nz
)
// allocates and initializes a new VELOCITYBOX
// on error: returns 0 (failure to allocate memory)
// on success: returns non-zero
// note: you should use vboxfree(&vbox) when you're done with it.
{
    if( !boxalloc( &vbox->box, nx, ny, nz ) ) return 0;

    point3dset( &vbox->min, ox, oy, oz );
    point3dset( &vbox->max, ox + nx - 1, oy + ny - 1, oz + nz - 1 );

    return 1;
}


void
vboxfree (
    struct VELOCITYBOX *vbox
)
// releases heap memory associated with vbox
{
    if( vbox == NULL ) return;
    boxfree( &vbox->box );
}


void
vboxfprint (
    FILE *stream,
    const char *prefix,
    const char *indent,
    struct VELOCITYBOX vbox
)
// metadata friendly-printing
{
    if( stream == NULL ) stream = stdout;
    if( prefix == NULL ) prefix = "";
    if( indent == NULL ) indent = "  ";

    fprintf( stream, "%sVELOCITYBOX {\n", prefix );
    fprintf( stream, "%s%sminimum corner: (%d, %d, %d)\n",
        prefix, indent, vbox.min.x, vbox.min.y, vbox.min.z );
    fprintf( stream, "%s%smaximum corner: (%d, %d, %d)\n",
        prefix, indent, vbox.max.x, vbox.max.y, vbox.max.z );

    char *newprefix = malloc( strlen( prefix ) + strlen( indent ) + 1 );
    if( newprefix != NULL ) {
        strcpy( newprefix, prefix );
        strcpy( newprefix + strlen( prefix ), indent );
        boxfprint( stream, newprefix, indent, vbox.box );
        free( newprefix );
    } else {
        boxfprint( stream, prefix, indent, vbox.box );
    }

    fprintf( stream, "%s}\n", prefix );
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
