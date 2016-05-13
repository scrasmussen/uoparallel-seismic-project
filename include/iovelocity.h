////////////////////////////////////////////////////////////////////////////////
// iovelocity.h - 2016.05.07 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// input/output for float-type, 3D, regular, velocity box
//
// read two types of velocity files:
//   * text format velocity
//   * binary format velocity: see VBOXFORMAT.txt
//
// write one type of velocity file:
//   * binary velocity box: see VBOXFORMAT.txt
//
// features:
//   * checksum verification of file integrity for binary files
//   * friendly error messages to stderr
//   * friendly metadata printing: see vboxfprint(..)
//
//
////////////////////////////////////////////////////////////////////////////////
//
// Data:
//   struct VELOCITYBOX:
//      efficiently stores a 3D velocity volume in flat memory
//
//
// Functions:
//    void vboxfree( struct VELOCITYBOX* )
//    int vboxloadtext( struct VELOCITYBOX*, char *filename )
//    int vboxloadbinary( struct VELOCITYBOX*, char *filename )
//    int vboxstorebinary( char *filename, struct VELOCITYBOX )
//    void vboxfprint( FILE*, char *prefix, char *indent, struct VELOCITYBOX )
//
//
// Example program:
//
//    #include "iovelocity.h"
//    #include <stdio.h> /* for stdout, not necessary otherwise */
//    int main() {
//      struct VELOCITYBOX vbox;
//      if( !vboxloadtext( &vbox, "velocity-241-241-51.txt" ) ) { /* handle error */ }
//      if( !vboxstorebinary( "velocities.vbox", vbox ) ) { /* handle error */ }
//      vboxfree( &vbox );
//      if( !vboxloadbinary( &vbox, "velocities.vbox" ) ) { /* handle error */ }
//      vboxfprint( stdout, "example: ", "\t", vbox );
//      vboxfree( &vbox );
//    }
//
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
vboxinit (
    struct VELOCITYBOX *vbox
)
// initializes vbox to safe default values
{
    if( vbox == NULL ) return;

    vbox->ox = 0;
    vbox->oy = 0;
    vbox->oz = 0;
    
    boxinit( &(vbox->box) );
}

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

    vboxinit( vbox );

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
                vboxfree( vbox );
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


void
vboxchecksum (
    uint32_t *checksum,
    union VBOX4BYTES fb
)
// updates a checksum based on the bytes in fb;
// assumes little-endian byte order, so use this AFTER conversion, if necessary;
{
    *checksum += (uint32_t)fb.c4[0]
        + ((uint32_t)fb.c4[1] << 8)
        + ((uint32_t)fb.c4[2] << 16)
        + ((uint32_t)fb.c4[3] << 24);
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
    vboxchecksum( checksum, fb );

    // success
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
    vboxchecksum( checksum, rb );

    // success
    return 1;
}


int
vboxstorebinary (
    const char *filename,
    struct VELOCITYBOX vbox
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vboxstorebinary";

    if( vbox.box.flat == NULL ) {
        fprintf( stderr, "%s: provided vbox is empty\n", fn );
        return 0;
    }

    FILE *outfile = fopen( filename, "wb" );
    if( outfile == NULL ) {
        // there was a problem opening/creating the given file
        fprintf( stderr, "%s: error creating file %s\n", fn, filename );
        return 0;
    }

    // this union lets us do stuff like checksums and endianness conversion
    union VBOX4BYTES fb;

    // keep track of write position and checksum
    long bytepos = 0;
    uint32_t checksum = 0;

    // write signature
    fb.c4[0] = (int8_t)'v';
    fb.c4[1] = (int8_t)'b';
    fb.c4[2] = (int8_t)'o';
    fb.c4[3] = (int8_t)'x';
    if( !vboxwrite4bytes( fn, outfile, filename, fb, &bytepos, &checksum ) ) {
        fclose( outfile );
        return 0;
    }

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

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

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


int
vboxread4bytes (
    const char *fn,
    FILE *infile,
    const char *filename,
    union VBOX4BYTES *fb,
    long *bytepos,
    uint32_t *checksum
)
// reads exactly 4 bytes from the given input FILE*;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    if( 4 != fread( fb, 1, 4, infile ) ) {
        fprintf( stderr, "%s: error reading from file %s at position %ld\n",
            fn, filename, *bytepos );
        return 0;
    }

    *bytepos += 4;
    vboxchecksum( checksum, *fb );

    // success
    return 1;
}


int
vboxread4bytesreversed (
    const char *fn,
    FILE *infile,
    const char *filename,
    union VBOX4BYTES *fb,
    long *bytepos,
    uint32_t *checksum
)
// reads exactly 4 bytes from the given input FILE*, but in reverse byte order;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    union VBOX4BYTES rb;

    if( 4 != fread( &rb, 1, 4, infile ) ) {
        fprintf( stderr, "%s: error reading from file %s at position %ld\n",
            fn, filename, *bytepos );
        return 0;
    }

    *bytepos += 4;
    vboxchecksum( checksum, rb );

    *fb = vboxreversebytes( rb );

    return 1;
}


int
vboxloadbinary (
    struct VELOCITYBOX *vbox,
    const char *filename
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vboxloadbinary";

    vboxinit( vbox );

    FILE *infile = fopen( filename, "rb" );
    if( infile == NULL ) {
        // there was a problem opening the given file
        fprintf( stderr, "%s: error opening file %s\n", fn, filename );
        return 0;
    }

    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // keep track of read position and checksum
    long bytepos = 0;
    uint32_t checksum = 0;

    // read signature
    if( !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum ) ) {
        fclose( infile );
        return 0;
    }
    if( fb.c4[0] != (int8_t)'v' || fb.c4[1] != (int8_t)'b'
    || fb.c4[2] != (int8_t)'o' || fb.c4[3] != (int8_t)'x' ) {
        // signature doesn't match
        fprintf( stderr, "%s: input file %s is not a vbox binary file, or is corrupted\n",
            fn, filename );
        fclose( infile );
        return 0;
    }

    // detect endianness of this machine: choose path accordingly
    if( fb.u32 == 0x786f6276 ) {
        // little-endian: read values directly

        int err = 0;

        // ox, oy, oz
        err |= !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum );
        vbox->ox = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum );
        vbox->oy = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum );
        vbox->oz = fb.i32;

        // check for errors
        if( err ) {
            fclose( infile );
            return 0;
        }

        // nx, ny, nz
        int nx, ny, nz;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum );
        nx = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum );
        ny = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum );
        nz = fb.i32;

        // check for errors
        if( err ) {
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

        // flat array of velocities
        size_t numvals = (size_t)vbox->box.nx * vbox->box.ny * vbox->box.nz;
        size_t i;
        for( i = 0; i < numvals; i++ ) {
            if( !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum ) ) {
                fclose( infile );
                vboxfree( vbox );
                return 0;
            }
            vbox->box.flat[i] = fb.f32;
        }

        // read final checksum value and compare
        uint32_t mychecksum = checksum;
        if( !vboxread4bytes( fn, infile, filename, &fb, &bytepos, &checksum ) ) {
            fclose( infile );
            vboxfree( vbox );
            return 0;
        }
        if( mychecksum != fb.u32 ) {
            fprintf( stderr, "%s: checksum mismatch in input file %s: suspect corruption\n",
                fn, filename );
            vboxfree( vbox );
            return 0;
        }

    } else {
        // big-endian: reverse bytes after reading

        int err = 0;

        // ox, oy, oz
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum );
        vbox->ox = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum );
        vbox->oy = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum );
        vbox->oz = fb.i32;

        // check for errors
        if( err ) {
            fclose( infile );
            return 0;
        }

        // nx, ny, nz
        int nx, ny, nz;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum );
        nx = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum );
        ny = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum );
        nz = fb.i32;

        // check for errors
        if( err ) {
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

        // flat array of velocities
        size_t numvals = (size_t)vbox->box.nx * vbox->box.ny * vbox->box.nz;
        size_t i;
        for( i = 0; i < numvals; i++ ) {
            if( !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum ) ) {
                fclose( infile );
                vboxfree( vbox );
                return 0;
            }
            vbox->box.flat[i] = fb.f32;
        }

        // read final checksum value and compare
        uint32_t mychecksum = checksum;
        if( !vboxread4bytesreversed( fn, infile, filename, &fb, &bytepos, &checksum ) ) {
            fclose( infile );
            vboxfree( vbox );
            return 0;
        }
        if( mychecksum != fb.u32 ) {
            fprintf( stderr, "%s: checksum mismatch in input file %s: suspect corruption\n",
                fn, filename );
            vboxfree( vbox );
            return 0;
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
