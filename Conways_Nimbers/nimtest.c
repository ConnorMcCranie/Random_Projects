/*
** Make nim multiplication, inverse, and square root tables.
** David Eppstein / Columbia University / 19 Mar 1988
*/

#include "nimber.h"
#include <stdio.h>

#define TABHEIGHT 16			/* 16 rows */
#define TABWIDTH 16			/* of 16 columns in table */

main()
{
    int i, j;

    /* multiplication table */
    printf ("  * |");
    for (j = 0; j < TABWIDTH; j++) printf ("%3d", j);
    printf ("\n----+");
    for (j = 0; j < TABWIDTH; j++) printf ("---");
    for (i = 0; i < TABHEIGHT; i++) {
	printf ("\n%3d |", i);
	for (j = 0; j < TABWIDTH; j ++)
	    printf("%3d", (nimber(i)*nimber(j)).value());
    }

    /* unary functions: 1/x, sqrt(x) */
    printf ("\n\nx:      ");
    for (j = 0; j < TABWIDTH; j++) printf ("%3d", j);
    printf ("\n--------");
    for (j = 0; j < TABWIDTH; j++) printf ("---");
    printf ("\n1/x:    ");
    for (j = 0; j < TABWIDTH; j++) printf ("%3d", nimber(j).inverse().value());
    printf ("\nsqrt(x):");
    for (j = 0; j < TABWIDTH; j++) printf ("%3d", nimber(j).sqrt().value());
    printf ("\nx^2:    ");
    for (j = 0; j < TABWIDTH; j++) printf ("%3d", nimber(j).square().value());
    printf ("\n");
}
