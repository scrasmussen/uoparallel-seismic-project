////////////////////////////////////////////////////////////////////////////////
// iovelocity.h - 2016.05.07 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// input/output for float-type, 3D, regular, velocity grid
//
// read two types of velocity files:
//   * verbose text velocity boxes
//   * dense binary velocity boxes
//
// write two types of velocity files:
//   * verbose text velocity boxes
//   * dense binary velocity boxes
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "../floatbox/floatbox.h"


#include <stdio.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct VELOCITYBOX {
    int ox, oy, oz; // offset of corner (0,0,0): used only for files
    struct FLOATBOX box; // contains dimensions and velocity data
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
vboxloadtext (
    struct VELOCITYBOX *vbox,
    const char *filename
)
// reads text-format velocity file of the format:
//   x_1,y_1,z_1,float_1
//   x_1,y_1,z_2,float_2
//   ...
//   x_nx,y_ny,z_nz,float_nxnynz
// example:
//   1,1,1,0.29762
//   ...
//   241,241,51,0.1771
// allocates vbox->box.flat and sets everything appropriately;
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vboxloadtext";

    FILE *infile = fopen( filename, "r" );
    if( infile == NULL ) {
        // there was a problem opening the given file
        fprintf( stderr, "%s: error opening file %s\n", fn, filename );
        return 0;
    }

    // read first line to get (x,y,z) of first velocity
    {
        int x, y, z;
        if( 3 != fscanf( infile, "%d,%d,%d", &x, &y, &z ) ) {
            // error
            fprintf( stderr, "%s: error reading first line from file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }
        vbox->ox = x;
        vbox->oy = y;
        vbox->oz = z;
    }

    // read last line to get (x,y,z) of last velocity;
    // this is used with the first coordinates to determine box size
    {
        int minlinelen = strlen( "1,1,1,0.0" );

        if( 0 != fseek( infile, -minlinelen, SEEK_END ) ) {
            // couldn't seek to estimated last line of file: maybe file is too short
            fprintf( stderr, "%s: error seeking to estimated last line in file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }

        // scan backward for a newline character
        for(;;) {
            int c = fgetc( infile );
            if( c == '\n' || c == '\r' ) {
                // found a newline: the file cursor is now at the start of the last line
                break;
            }
            if( c == EOF || 0 != fseek( infile, -2, SEEK_CUR ) ) {
                // some i/o error occurred
                fprintf( stderr, "%s: error scanning for last line in file %s\n",
                    fn, filename );
                fclose( infile );
                return 0;
            }
        }

        int x, y, z;
        if( 3 != fscanf( infile, "%d,%d,%d", &x, &y, &z ) ) {
            // error
            fprintf( stderr, "%s: error reading last line from file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }

        int nx, ny, nz;
        nx = x - vbox->ox + 1;
        ny = y - vbox->oy + 1;
        nz = z - vbox->oz + 1;

        if( nx <= 0 || ny <= 0 || nz <= 0 ) {
            // dimensions are nonsense
            fprintf( stderr, "%s: nonsense coordinates in file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }

        if( !boxalloc( &(vbox->box), nx, ny, nz ) ) {
            // couldn't allocate memory for floatbox
            fprintf( stderr, "%s: unable to allocate memory for a FLOATBOX with"
                "dimension: %d x %d x %d\n", fn, nx, ny, nz );
            fclose( infile );
            return 0;
        }
    }

    // starting from the beginning of the file, read and remember all velocity values
    {
        fseek( infile, 0, SEEK_SET );
        size_t numlines = (size_t)vbox->box.nx * vbox->box.ny * vbox->box.nz;
        size_t l;
        for( l = 0; l < numlines; l++ ) {
            int x, y, z;
            float vel;
            if( 4 != fscanf( infile, "%d,%d,%d,%f\n", &x, &y, &z, &vel ) ) {
                // problem parsing a line
                fprintf( stderr, "%s: I am confused by line %zu in %s\n",
                    fn, l+1, filename );
                boxfree( &(vbox->box) );
                fclose( infile );
                return 0;
            }
            vbox->box.flat[ l ] = vel;
        }
    }

    fclose( infile );

    // success
    return 1;
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
    fprintf( stream, "%s%sox: %d, oy: %d, oz: %d\n",
        prefix, indent, vbox.ox, vbox.oy, vbox.oz );

    char *newprefix = malloc( strlen( prefix ) + strlen( indent ) + 1 );
    strcpy( newprefix, prefix );
    strcpy( newprefix + strlen( prefix ), indent );
    boxfprint( stream, newprefix, indent, vbox.box );
    free( newprefix );

    fprintf( stream, "%s}\n", prefix );
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
