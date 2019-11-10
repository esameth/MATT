/* Score.h -- Matt scoring function prototypes and alignment constants.
 *
 * Author: Matt Menke (2007)
 *
 * Matt is licensed under the GNU public license version 2.0. 
 *
 * If you would like to license Matt in an enviroment where the GNU public license is 
 * unacceptable (such as inclusion in a non-GPL software package) comercial Matt
 * licensing is available through the MIT and Tufts offices of Technology Transfer. 
 * Contact betawrap@csail.mit.edu or cowen@cs.tufts.edu for more information.
 */

#ifndef MAKE_RMSD_TABLES_H
#define MAKE_RMSD_TABLES_H

#include <stdio.h>
#include <math.h>
#include "Vector.h"

#define MIN_FRAGMENT_LEN 5
// Fragments may end up longer than this after extension phase.
#define MAX_FRAGMENT_LEN 5

#define MAX_EXTEND_DIST 5

extern const double RMSDTable[3001];

double RMSDScoreSimple(double rmsd, int length);

/* Saves some time using double angle identity.
 * Takes in cos(angle/2) as input.
 */
double AngleScoreHalfCos(double angleHalfCos);
double AngleScore(double angle);

double DisplacementScore(double displacement);

/* highest possible value of AngleScore() + DisplacementScore() */
#define MAX_BONUS (AngleScoreHalfCos(1)+DisplacementScore(0))
#define MAX_BONUS_TEST (AngleScore(0)+DisplacementScore(0))

/* highest possible value of DisplacementScore() */
#define MAX_POS_BONUS (DisplacementScore(0))

#endif
