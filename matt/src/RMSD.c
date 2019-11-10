/* RMSD.h -- RMSD calculation functions.
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

#include "RMSD.h"
#include <math.h>
#include <stdio.h>

void normalizeRow(Matrix *m, int row) {
	double len = sqrt(m->M[row][0] * m->M[row][0] + m->M[row][1]*m->M[row][1] + m->M[row][2]*m->M[row][2]);
	m->M[row][0] /= len;
	m->M[row][1] /= len;
	m->M[row][2] /= len;
}

void normalizeCol(Matrix *m, int col) {
	double len = sqrt(m->M[0][col] * m->M[0][col] + m->M[1][col]*m->M[1][col] + m->M[2][col]*m->M[2][col]);
	m->M[0][col] /= len;
	m->M[1][col] /= len;
	m->M[2][col] /= len;
}

double CalcRotation(Matrix *rotation, Matrix *A, Matrix *ATA2, double E0, double twoOverLength) {
	int i;
	int x, y, x1, x2, y1, y2;
	double sign = 1;
	Matrix vectors;
	Matrix u;
	Vector v;
	Matrix ATA = *ATA2;

	vectors = identity;
	for (i=0; i<80; i++) {
		double total = fabs(ATA.M[1][0]) + fabs(ATA.M[2][0]) + fabs(ATA.M[2][1]);
		double temp;
		if (total == 0) break;
		if (i >= 3) total = 0;
		else total *= 0.20;
		x = 0; y = 1; x1 = 0; x2 = 2; y1 = 1; y2 = 2;
		if (i > 3 && fabs(ATA.M[x][x]) + (temp = 100*fabs(ATA.M[y][x])) == fabs(ATA.M[x][x]) &&
			fabs(ATA.M[y][y])+temp == fabs(ATA.M[y][y])) {
		  ATA.M[y][x] = 0.0;
		} 
		else if (fabs(ATA.M[y][x]) > total) {
			double a = (ATA.M[y][y] - ATA.M[x][x]) / (2*ATA.M[y][x]);
			double t = (2*(a>=0)-1) / (fabs(a) + sqrt(a*a+1));
			double c = 1 / sqrt(t*t + 1);
			double s = c*t;
			double stau1 = 1.0 - s*s / (1.0+c);

			double vx = ATA.M[x2][x1];
			double vy = ATA.M[y2][y1];
			ATA.M[x2][x1] = vx * stau1 - s * vy;
			ATA.M[y2][y1] = vy * stau1 + s * vx;

			ATA.M[x][x] -= t*ATA.M[y][x];
			ATA.M[y][y] += t*ATA.M[y][x];
			ATA.M[y][x] = 0;

			vx = vectors.M[0][x];
			vy = vectors.M[0][y];
			vectors.M[0][x] = vx * stau1 - s * vy;
			vectors.M[0][y] = vy * stau1 + s * vx;
			vx = vectors.M[1][x];
			vy = vectors.M[1][y];
			vectors.M[1][x] = vx * stau1 - s * vy;
			vectors.M[1][y] = vy * stau1 + s * vx;
			vx = vectors.M[2][x];
			vy = vectors.M[2][y];
			vectors.M[2][x] = vx * stau1 - s * vy;
			vectors.M[2][y] = vy * stau1 + s * vx;
		}
		x = 0; y = 2; x1 = 0; x2 = 1; y1 = 1; y2 = 2;
		if (i > 3 && fabs(ATA.M[x][x]) + (temp = 100*fabs(ATA.M[y][x])) == fabs(ATA.M[x][x]) &&
			fabs(ATA.M[y][y])+temp == fabs(ATA.M[y][y])) {
		  ATA.M[y][x] = 0.0;
		} 
		else if (fabs(ATA.M[y][x]) > total) {
			double a = (ATA.M[y][y] - ATA.M[x][x]) / (2*ATA.M[y][x]);
			double t = (2*(a>=0)-1) / (fabs(a) + sqrt(a*a+1));
			double c = 1 / sqrt(t*t + 1);
			double s = c*t;
			double stau1 = 1.0 - s*s / (1.0+c);

			double vx = ATA.M[x2][x1];
			double vy = ATA.M[y2][y1];
			ATA.M[x2][x1] = vx * stau1 - s * vy;
			ATA.M[y2][y1] = vy * stau1 + s * vx;

			ATA.M[x][x] -= t*ATA.M[y][x];
			ATA.M[y][y] += t*ATA.M[y][x];
			ATA.M[y][x] = 0;

			vx = vectors.M[0][x];
			vy = vectors.M[0][y];
			vectors.M[0][x] = vx * stau1 - s * vy;
			vectors.M[0][y] = vy * stau1 + s * vx;
			vx = vectors.M[1][x];
			vy = vectors.M[1][y];
			vectors.M[1][x] = vx * stau1 - s * vy;
			vectors.M[1][y] = vy * stau1 + s * vx;
			vx = vectors.M[2][x];
			vy = vectors.M[2][y];
			vectors.M[2][x] = vx * stau1 - s * vy;
			vectors.M[2][y] = vy * stau1 + s * vx;
		}
		x = 1; y = 2; x1 = 0; x2 = 1; y1 = 0; y2 = 2;
		if (i > 3 && fabs(ATA.M[x][x]) + (temp = 100*fabs(ATA.M[y][x])) == fabs(ATA.M[x][x]) &&
			fabs(ATA.M[y][y])+temp == fabs(ATA.M[y][y])) {
			ATA.M[y][x] = 0.0;
		} 
		else if (fabs(ATA.M[y][x]) > total) {
			double a = (ATA.M[y][y] - ATA.M[x][x]) / (2*ATA.M[y][x]);
			double t = (2*(a>=0)-1) / (fabs(a) + sqrt(a*a+1));
			double c = 1 / sqrt(t*t + 1);
			double s = c*t;
			double stau1 = 1.0 - s*s / (1.0+c);

			double vx = ATA.M[x2][x1];
			double vy = ATA.M[y2][y1];
			ATA.M[x2][x1] = vx * stau1 - s * vy;
			ATA.M[y2][y1] = vy * stau1 + s * vx;

			ATA.M[x][x] -= t*ATA.M[y][x];
			ATA.M[y][y] += t*ATA.M[y][x];
			ATA.M[y][x] = 0;

			vx = vectors.M[0][x];
			vy = vectors.M[0][y];
			vectors.M[0][x] = vx * stau1 - s * vy;
			vectors.M[0][y] = vy * stau1 + s * vx;
			vx = vectors.M[1][x];
			vy = vectors.M[1][y];
			vectors.M[1][x] = vx * stau1 - s * vy;
			vectors.M[1][y] = vy * stau1 + s * vx;
			vx = vectors.M[2][x];
			vy = vectors.M[2][y];
			vectors.M[2][x] = vx * stau1 - s * vy;
			vectors.M[2][y] = vy * stau1 + s * vx;
		}
	}
	if (i == 50) {
		/* printf("Warning - RMSD calculation failed\n"); //*/
		return 100000;
	}
	if (ATA.M[0][0] < ATA.M[2][2]) {
		double v[3] = {vectors.M[0][0], vectors.M[1][0], vectors.M[2][0]};
		double temp = ATA.M[0][0];
		ATA.M[0][0] = ATA.M[2][2];
		ATA.M[2][2] = temp;

		vectors.M[0][0] = vectors.M[0][2];
		vectors.M[1][0] = vectors.M[1][2];
		vectors.M[2][0] = vectors.M[2][2];
		vectors.M[0][2] = v[0];
		vectors.M[1][2] = v[1];
		vectors.M[2][2] = v[2];
	}
	if (ATA.M[1][1] > ATA.M[0][0]) {
		double v[3] = {vectors.M[0][1], vectors.M[1][1], vectors.M[2][1]};
		double temp = ATA.M[0][0];
		ATA.M[0][0] = ATA.M[1][1];
		ATA.M[1][1] = temp;

		vectors.M[0][1] = vectors.M[0][0];
		vectors.M[1][1] = vectors.M[1][0];
		vectors.M[2][1] = vectors.M[2][0];
		vectors.M[0][0] = v[0];
		vectors.M[1][0] = v[1];
		vectors.M[2][0] = v[2];
	}
	else if (ATA.M[1][1] < ATA.M[2][2]) {
		double temp = ATA.M[2][2];
		ATA.M[2][2] = ATA.M[1][1];
		ATA.M[1][1] = temp;

		vectors.M[0][1] = vectors.M[0][2];
		vectors.M[1][1] = vectors.M[1][2];
		vectors.M[2][1] = vectors.M[2][2];
	}
	vectors.M[0][2] = vectors.M[1][0] * vectors.M[2][1] - vectors.M[1][1] * vectors.M[2][0];
	vectors.M[1][2] = vectors.M[2][0] * vectors.M[0][1] - vectors.M[2][1] * vectors.M[0][0];
	vectors.M[2][2] = vectors.M[0][0] * vectors.M[1][1] - vectors.M[0][1] * vectors.M[1][0];
	normalizeCol(&vectors, 0);
	normalizeCol(&vectors, 1);
	normalizeCol(&vectors, 2);

	mulMat3x3(&u, A, &vectors);
	v.x = u.M[1][0] * u.M[2][1] - u.M[1][1] * u.M[2][0];
	v.y = u.M[2][0] * u.M[0][1] - u.M[2][1] * u.M[0][0];
	v.z = u.M[0][0] * u.M[1][1] - u.M[0][1] * u.M[1][0];

	if (v.x * u.M[0][2] + v.y * u.M[1][2] + v.z * u.M[2][2] < 0) {
		sign = -1.0;
	}
	u.M[0][2] = v.x;
	u.M[1][2] = v.y;
	u.M[2][2] = v.z;
	normalizeCol(&u, 0);
	normalizeCol(&u, 1);
	normalizeCol(&u, 2);

	for (x=0; x<3; x++) {
		for (y=0; y<3; y++) {
			rotation->M[y][x] = u.M[y][0] * vectors.M[x][0] + u.M[y][1] * vectors.M[x][1] + u.M[y][2] * vectors.M[x][2];
		}
	}

	return fabs(twoOverLength * (E0 - sqrt(fabs(ATA.M[0][0])) - sqrt(fabs(ATA.M[1][1])) - sign*sqrt(fabs(ATA.M[2][2]))));
}

double CalcSqRMSDAndRotation(Matrix *rotation, Matrix *A, double E0, double twoOverLength) {
	Matrix ATA;
	int row;
	for (row = 0; row<3; row++) {
		ATA.M[row][0] =
			A->M[0][row] * A->M[0][0] +
			A->M[1][row] * A->M[1][0] +
			A->M[2][row] * A->M[2][0];
		ATA.M[row][1] =
			A->M[0][row] * A->M[0][1] +
			A->M[1][row] * A->M[1][1] +
			A->M[2][row] * A->M[2][1];
		ATA.M[row][2] =
			A->M[0][row] * A->M[0][2] +
			A->M[1][row] * A->M[1][2] +
			A->M[2][row] * A->M[2][2];
	}
	return CalcRotation(rotation, A, &ATA, E0, twoOverLength);
}


double CalcSqRMSDUnderestimate(Matrix *ATA, Matrix *A, double E0, double twoOverLength) {
	double angle;
	double val[3];
	double cubic[3];
	double cubic0_3;
	double p, q, t;
	int row;
	for (row = 0; row<3; row++) {
		ATA->M[row][0] =
			A->M[0][row] * A->M[0][0] +
			A->M[1][row] * A->M[1][0] +
			A->M[2][row] * A->M[2][0];
		ATA->M[row][1] =
			A->M[0][row] * A->M[0][1] +
			A->M[1][row] * A->M[1][1] +
			A->M[2][row] * A->M[2][1];
		ATA->M[row][2] =
			A->M[0][row] * A->M[0][2] +
			A->M[1][row] * A->M[1][2] +
			A->M[2][row] * A->M[2][2];
	}
	cubic[0] = 	-ATA->M[0][0]-ATA->M[1][1]-ATA->M[2][2];
	cubic[1] = ATA->M[0][0] * (ATA->M[1][1] + ATA->M[2][2]) + ATA->M[2][2]*ATA->M[1][1]
		      -ATA->M[1][2]*ATA->M[2][1] - ATA->M[2][0]*ATA->M[0][2] - ATA->M[1][0]*ATA->M[0][1];
	cubic[2] = ATA->M[0][0] * (ATA->M[2][1]*ATA->M[1][2] - ATA->M[1][1]*ATA->M[2][2]) +
			   ATA->M[1][0] * (ATA->M[0][1]*ATA->M[2][2] - ATA->M[2][1]*ATA->M[0][2]) +
			   ATA->M[2][0] * (ATA->M[1][1]*ATA->M[0][2] - ATA->M[1][2]*ATA->M[0][1]);

	cubic0_3 = cubic[0]*(1.0/3.0);
	t = cubic0_3*cubic0_3;

	p = t - cubic[1]*(1.0/3.0);
	q = (cubic[2] + cubic0_3*(2*t-cubic[1]))* 0.5;
	t = q * q - p * p * p;
	/* Don't think the latter ever happens, unless the first does as well. */
	if (t > -0.000001 || p < 0) {
		if (t < 0.000001) {
			/* Triple root. */
			val[0] = val[1] = val[2] = -cubic0_3;
		}
		else {
			/* Imaginary roots. */
			return 500000;
		}
	}
	else {
		double sqrtP = sqrt(p);
		double quadratic[2];
		double v, b2;

		angle = acos(q/(p * sqrtP)) * (1.0/3.0);

		val[0] = -2 * cos(angle) * sqrtP - cubic0_3;

		quadratic[0] = cubic[0] + val[0];
		quadratic[1] = cubic[1] + val[0]*quadratic[0];

		b2 = quadratic[0]*-0.5;
		v = b2*b2 - quadratic[1];
		if (v < 0) {
			if (v > -0.0001) {
				v = 0;
			}
			else {
				return 500000;
			}
		}
		v = sqrt(v);
		val[1] = b2+v;
		val[2] = b2-v;
	}
	return twoOverLength * (E0 - sqrt(fabs(val[0])) - sqrt(fabs(val[1])) - sqrt(fabs(val[2])));
}

double CalcPositionAllPairsRMSD(MultipleAlignment *ma, int *indices) {
  double sum = 0;
  int numSummed = 0;
  int i, j;
  for (i=0; i<ma->numChains; i++) {
    /* If location EXISTS (hasCoords) */
    if (indices[i] >= 0 && ma->chains[i]->res[indices[i]].hasCoords) {
      ResiduePosition *p1 = &ma->chains[i]->res[indices[i]];
      for (j=i+1; j<ma->numChains; j++) { /* Comp. to other residues in col. */
        if (indices[j] >= 0 && ma->chains[j]->res[indices[j]].hasCoords) {
          ResiduePosition *p2 = &ma->chains[j]->res[indices[j]];
          Vector v;
          sum += lengthSquaredVect(subVect(&v, &p1->coords, &p2->coords));
          numSummed++;
        }
      }
    }
  }
  if (!numSummed) return 0;
  return sqrt(sum/numSummed);
}


