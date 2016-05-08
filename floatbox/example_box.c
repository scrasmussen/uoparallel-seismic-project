#include "floatbox.h"
#include <stdio.h> /* for printf and stdout, not necessary otherwise */
int main() {
    struct FLOATBOX box;
    if( !boxalloc( &box, 241, 241, 51 ) ) { /* handle error */ }
    boxfprint( stdout, "example: ", "\t", box );
    boxput( box, 1, 2, 3, 4.567f );
    printf( "value at (1,2,3): %g\n", boxget( box, 1, 2, 3 ) );
    boxfree( &box );
}
