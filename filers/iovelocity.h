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


#include <stdint.h>
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
// unions
////////////////////////////////////////////////////////////////////////////////

union VBOX4BYTES {
    int8_t c4[4];
    int32_t i32;
    uint32_t u32;
    float f32; 
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

void
vboxfree (
    struct VELOCITYBOX *vbox
)
// releases heap memory associated with vbox
{
    boxfree( &(vbox->box) );
}


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


inline extern
union VBOX4BYTES
vboxreversebytes (
    const union VBOX4BYTES in
)
// reverses the byte-order of four bytes
{
    union VBOX4BYTES out;

    out.c4[0] = in.c4[3];
    out.c4[1] = in.c4[2];
    out.c4[2] = in.c4[1];
    out.c4[3] = in.c4[0];

    return out;
}


int
vboxwrite4bytes (
    const char *fn,
    FILE *outfile,
    const char *filename,
    union VBOX4BYTES fb,
    long *bytepos,
    uint32_t *checksum
)
// writes exactly 4 bytes to the given output FILE*;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    if( 4 != fwrite( &fb, 1, 4, outfile ) ) {
        fprintf( stderr, "%s: error writing to file %s at position %ld\n",
            fn, filename, *bytepos );
        return 0;
    }

    *bytepos += 4;
    *checksum += fb.u32;

    return 1;
}


int
vboxwrite4bytesreversed (
    const char *fn,
    FILE *outfile,
    const char *filename,
    union VBOX4BYTES fb,
    long *bytepos,
    uint32_t *checksum
)
// writes exactly 4 bytes to the given output FILE*, but in reverse byte order;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    union VBOX4BYTES rb = vboxreversebytes( fb );

    if( 4 != fwrite( &rb, 1, 4, outfile ) ) {
        fprintf( stderr, "%s: error writing to file %s at position %ld\n",
            fn, filename, *bytepos );
        return 0;
    }

    *bytepos += 4;
    *checksum += fb.u32;

    return 1;
}


int
vboxstorebinary (
    const char *filename,
    struct VELOCITYBOX vbox
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-0
{
    const char *fn = "vboxstorebinary";

    FILE *outfile = fopen( filename, "wb" );
    if( outfile == NULL ) {
        // there was a problem opening/creating the given file
        fprintf( stderr, "%s: error creating file %s\n", fn, filename );
        return 0;
    }

    // this union lets us do stuff like checksums and endianness conversion
    union VBOX4BYTES fb;

    // write header
    fb.c4[0] = (int8_t)'v';
    fb.c4[1] = (int8_t)'b';
    fb.c4[2] = (int8_t)'o';
    fb.c4[3] = (int8_t)'x';

    if( 4 != fwrite( &fb, 1, 4, outfile ) ) {
        // there was a problem writing to the file
        fprintf( stderr, "%s: error writing to file %s\n", fn, filename );
        fclose( outfile );
        return 0;
    }

    long bytepos = 4;
    uint32_t checksum = 0;

    // detect endianness of this machine: choose path accordingly
    if( fb.u32 == 0x786f6276 ) {
        // little-endian: write values directly

        int err = 0;

        // ox, oy, oz
        fb.i32 = vbox.ox;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.oy;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.oz;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum );

        // nx, ny, nz
        fb.i32 = vbox.box.nx;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.box.ny;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.box.nz;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        size_t numvals = (size_t)vbox.box.nx * vbox.box.ny * vbox.box.nz;
        size_t i;
        for( i = 0; i < numvals; i++ ) {
            fb.f32 = vbox.box.flat[i];
            if( !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum ) ) {
            fclose( outfile );
            return 0;
        }

    } else {
        // big-endian: must reverse byte order before writing

        int err = 0;

        // ox, oy, oz
        fb.i32 = vbox.ox;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.oy;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.oz;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum );

        // nx, ny, nz
        fb.i32 = vbox.box.nx;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.box.ny;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum );
        fb.i32 = vbox.box.nz;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        size_t numvals = (size_t)vbox.box.nx * vbox.box.ny * vbox.box.nz;
        size_t i;
        for( i = 0; i < numvals; i++ ) {
            fb.f32 = vbox.box.flat[i];
            if( !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !vboxwrite4bytesreversed( fn, outfile, filename, fb, &bytepos, &checksum ) ) {
            fclose( outfile );
            return 0;
        }
    }

    fclose( outfile );

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
