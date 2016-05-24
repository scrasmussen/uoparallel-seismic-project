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

#include "floatbox.h"


#include <stdint.h>
#include <stdio.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct VELOCITYBOX {
    struct POINT3D min, max;
    struct FLOATBOX box; // contains dimensions and velocity data
};


struct VBOXOPENFILE {
    FILE *file;    // if non-NULL, then assume file is successfully open
    int is_little_endian;
    struct POINT3D min, dims;
    uint32_t checksum;
    const char *filename;
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

    point3dset( &vbox->min, 0, 0, 0 );
    point3dset( &vbox->max, 0, 0 ,0 );

    boxinit( &vbox->box );
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
        point3dset( &vbox->min, x, y, z );
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

        // read coordinates (assume they are correct)
        int x, y, z;
        if( 3 != fscanf( infile, "%d,%d,%d", &x, &y, &z ) ) {
            // error
            fprintf( stderr, "%s: error reading last line from file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }
        point3dset( &vbox->max, x, y, z );

        int nx, ny, nz;
        nx = vbox->max.x - vbox->min.x + 1;
        ny = vbox->max.y - vbox->min.y + 1;
        nz = vbox->max.z - vbox->min.z + 1;

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
        size_t numlines = (size_t)vbox->box.size.x * vbox->box.size.y * vbox->box.size.z;
        //size_t numlines = (size_t)vbox->box.nx * vbox->box.ny * vbox->box.nz;
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
    uint32_t *checksum
)
// writes exactly 4 bytes to the given output FILE*;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    if( 4 != fwrite( &fb, 1, 4, outfile ) ) {
        fprintf( stderr, "%s: error writing to file %s at position %ld\n",
            fn, filename, ftell( outfile ) );
        return 0;
    }

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
            fn, filename, ftell( outfile ) );
        return 0;
    }

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
    uint32_t checksum = 0;

    // write signature
    fb.c4[0] = (int8_t)'v';
    fb.c4[1] = (int8_t)'b';
    fb.c4[2] = (int8_t)'o';
    fb.c4[3] = (int8_t)'x';
    if( !vboxwrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
        fclose( outfile );
        return 0;
    }

    // detect endianness of this machine: choose path accordingly
    if( fb.u32 == 0x786f6276 ) {
        // little-endian: write values directly

        int err = 0;

        // ox, oy, oz
        fb.i32 = vbox.min.x;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.y;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.z;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // nx, ny, nz
        fb.i32 = vbox.box.size.x;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.y;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.z;
        err |= !vboxwrite4bytes( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        size_t numvals = boxvolume( vbox.box );
        //size_t numvals = (size_t)vbox.box.size.x * vbox.box.size.y * vbox.box.size.z;
        size_t i;
        for( i = 0; i < numvals; i++ ) {
            fb.f32 = vbox.box.flat[i];
            if( !vboxwrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !vboxwrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
            fclose( outfile );
            return 0;
        }

    } else {
        // big-endian: must reverse byte order before writing

        int err = 0;

        // ox, oy, oz
        fb.i32 = vbox.min.x;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.y;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.z;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum );

        // nx, ny, nz
        fb.i32 = vbox.box.size.x;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.y;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.z;
        err |= !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        size_t numvals = boxvolume( vbox.box );
        //size_t numvals = (size_t)vbox.box.nx * vbox.box.ny * vbox.box.nz;
        size_t i;
        for( i = 0; i < numvals; i++ ) {
            fb.f32 = vbox.box.flat[i];
            if( !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !vboxwrite4bytesreversed( fn, outfile, filename, fb, &checksum ) ) {
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
    uint32_t *checksum
)
// reads exactly 4 bytes from the given input FILE*;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    if( 4 != fread( fb, 1, 4, infile ) ) {
        fprintf( stderr, "%s: error reading from file %s at position %ld\n",
            fn, filename, ftell( infile ) );
        return 0;
    }

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
            fn, filename, ftell( infile ) );
        return 0;
    }

    vboxchecksum( checksum, rb );

    *fb = vboxreversebytes( rb );

    return 1;
}


int                                 // returns: 0 on error, else non-0
vboxopenbinary (
    struct VBOXOPENFILE *vbfile,    // out: stores metadata of open file
    const char *filename            // in: see VBOXFORMAT.txt for details
)
// note: when done with 'vbfile', be sure to run vboxclosebinary( &vbfile )
{
    const char *fn = "vboxopenbinary";

    vbfile->file = NULL;

    FILE *infile = fopen( filename, "rb" );
    if( infile == NULL ) {
        fprintf( stderr, "%s: error opening file %s\n", fn, filename );
        return 0;
    }

    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // keep track of read position and checksum
    uint32_t checksum = 0;

    // read signature
    if( !vboxread4bytes( fn, infile, filename, &fb, &checksum ) ) {
        fprintf( stderr, "%s: input file %s is not a vbox binary file, or is corrupted\n",
            fn, filename );
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

    // global minimum coordinates of velocity volume
    int ox, oy, oz;

    // dimensions of velocity volume
    int nx, ny, nz;

    // check after reading header: if err != 0, then there was an error
    int err = 0;

    // detect endianness of this machine: choose path accordingly
    if( (vbfile->is_little_endian = (fb.u32 == 0x786f6276)) ) { // 'vbox'
        // little-endian: read values directly

        // ox, oy, oz
        err |= !vboxread4bytes( fn, infile, filename, &fb, &checksum );
        ox = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &checksum );
        oy = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &checksum );
        oz = fb.i32;

        // nx, ny, nz
        err |= !vboxread4bytes( fn, infile, filename, &fb, &checksum );
        nx = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &checksum );
        ny = fb.i32;
        err |= !vboxread4bytes( fn, infile, filename, &fb, &checksum );
        nz = fb.i32;
    }
    else {
        // big-endian: reverse bytes after reading

        // ox, oy, oz
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &checksum );
        ox = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &checksum );
        oy = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &checksum );
        oz = fb.i32;

        // nx, ny, nz
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &checksum );
        nx = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &checksum );
        ny = fb.i32;
        err |= !vboxread4bytesreversed( fn, infile, filename, &fb, &checksum );
        nz = fb.i32;
    }

    // check for errors
    if( err ) {
        fprintf( stderr, "%s: error reading header in %s: suspect corruption\n",
            fn, filename );
        fclose( infile );
        return 0;
    }

    // as appropriate, store metadata from header
    vbfile->file = infile;
    point3dset( &vbfile->min, ox, oy, oz );
    point3dset( &vbfile->dims, nx, ny, nz );
    vbfile->checksum = checksum;
    vbfile->filename = filename;

    // success
    return 1;
}


void
vboxclosebinary (
    struct VBOXOPENFILE *vbfile
)
{
    if( vbfile == NULL || vbfile->file == NULL ) return;
    fclose( vbfile->file );
    vbfile->file = NULL;
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

    if( vbox == NULL ) return 0;

    // file metadata
    struct VBOXOPENFILE vbfile;

    // open file and read header
    if( !vboxopenbinary( &vbfile, filename ) ) return 0;

    // prepare FLOATBOX to store the contents of the file
    if( !boxalloc( &vbox->box, vbfile.dims.x, vbfile.dims.y, vbfile.dims.z ) ) {
        fprintf( stderr, "%s: unable to allocate memory for a FLOATBOX with"
            "dimension: %d x %d x %d\n", fn, vbfile.dims.x, vbfile.dims.y, vbfile.dims.z );
        vboxclosebinary( &vbfile );
        return 0;
    }

    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // read flat array of velocities from file and store in the FLOATBOX
    {
        size_t i, numvals = boxvolume( vbox->box );
        uint32_t stored_checksum;

        if( vbfile.is_little_endian ) {
            // read floats
            for( i = 0; i < numvals; i++ ) {
                if( !vboxread4bytes( fn, vbfile.file, filename, &fb, &vbfile.checksum ) ) {
                    fprintf( stderr, "%s: error reading value at byte position %zu"
                        " in %s\n", fn, ftell( vbfile.file ), filename );
                    vboxclosebinary( &vbfile );
                    vboxfree( vbox );
                    return 0;
                }
                vbox->box.flat[i] = fb.f32;
            }
            // read stored checksum value
            uint32_t dummy_checksum;
            if( !vboxread4bytes( fn, vbfile.file, filename, &fb, &dummy_checksum ) ) {
                fprintf( stderr, "%s: error reading stored checksum value from %s\n",
                    fn, filename );
                vboxclosebinary( &vbfile );
                vboxfree( vbox );
                return 0;
            }
            stored_checksum = fb.u32;
        }
        else { // big endian
            // read floats
            for( i = 0; i < numvals; i++ ) {
                if( !vboxread4bytesreversed( fn, vbfile.file, filename, &fb, &vbfile.checksum ) ) {
                    fprintf( stderr, "%s: error reading value at byte position %zu"
                        " in %s\n", fn, ftell( vbfile.file ), filename );
                    vboxclosebinary( &vbfile );
                    vboxfree( vbox );
                    return 0;
                }
                vbox->box.flat[i] = fb.f32;
            }
            // read stored checksum value
            uint32_t dummy_checksum;
            if( !vboxread4bytesreversed( fn, vbfile.file, filename, &fb, &dummy_checksum ) ) {
                fprintf( stderr, "%s: error reading stored checksum value from %s\n",
                    fn, filename );
                vboxclosebinary( &vbfile );
                vboxfree( vbox );
                return 0;
            }
            stored_checksum = fb.u32;
        }
        
        // done reading file
        vboxclosebinary( &vbfile );

        // verify that computed vbfile.checksum matches stored checksum:
        // else something is corrupt
        if( vbfile.checksum != stored_checksum ) {
            fprintf( stderr, "%s: checksum mismatch in input file %s: suspect corruption\n",
                fn, filename );
            vboxfree( vbox );
            return 0;
        }
    }

    // set volume metadata in vbox
    vbox->min = vbfile.min;
    point3dset(
        &vbox->max,
        vbfile.min.x + vbfile.dims.x - 1,
        vbfile.min.y + vbfile.dims.y - 1,
        vbfile.min.z + vbfile.dims.z - 1
    );
    
    // success
    return 1;
}

/*
int
vboxloadbinarysubset (
    struct VELOCITYBOX *vbox,
    const char *filename,

)
// see VBOXFORAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vboxloadbinarysubset";

    if( vbox == NULL ) return 0;

    // file metadata
    int little_endian;
    struct POINT3D boxmin, boxsize;
    uint32_t checksum;

    // open file and read header
    FILE *infile = vboxopenbinary (
        &little_endian,
        &boxmin,
        &boxsize,
        &checksum,
        filename
    );
    if( infile == NULL ) return 0;
}
*/

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
