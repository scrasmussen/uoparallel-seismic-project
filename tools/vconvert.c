// vconvert.c - Atlee Brink
// converts from old velocity file format to new (smaller) file format

#include "../include/iovelocity.h"

#include <stdio.h>

int
main (
    int argc,
    char *argv[]
)
{
    if( argc != 3 ) {
        printf( "vconvert: velocity file converter\n" );
        printf( "usage: %s <in:oldfile.txt> <out:newfile.vbox>\n", argv[0] );
        return 0;
    }

    char *infile = argv[1];
    char *outfile = argv[2];

    struct VELOCITYBOX vbox;

    printf( "reading old velocity model %s...", infile ); fflush( stdout );
    if( !vboxloadtext( &vbox, infile ) ) return 1;
    printf( " done.\n" ); fflush( stdout );

    printf( "writing new velocity model %s...", outfile ); fflush( stdout );
    if( !vboxstorebinary( outfile, vbox ) ) return 1;
    printf( " done.\n" ); fflush( stdout );

    vboxfree( &vbox );

    return 0;
}
