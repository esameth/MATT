/* MultipleAlignment.c -- Contains functions for creating and manipulating multiple
 * alighments.
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

//#include <Windows.h>

//int time[6] = {0};


#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>

#include "MultipleAlignment.h"

#include "OctTree.h"
#include "RMSD.h"
#include "Score.h"
#include "Extend.h"

#define M_1_SQRT2PI 0.39894228040143270

struct SCORE_DATA {
	double minPVal;
	double maxSqRMSD[MAX_FRAGMENT_LEN+1];
	double angle;
	double angleInvHalfCos;
	double displacement;
	int firstStrandLocked;
	int lastStrandLocked;
	int selfAlignDistance;
	int pass;
};
typedef struct SCORE_DATA ScoreData;

ScoreData ScoreTables[3];

void InitMultipleAlignmentTables() {
	int i, len;
	for (i=0; i<3; i++) {
		ScoreData *s = ScoreTables+i;
		s->pass = i;
		if (i==0) {
			s->minPVal = 2.0;
			s->angle = 45 * M_PI / 180;
			s->displacement = 4;
		}
		else if (i==1) {
			s->minPVal = 1.6;
			s->angle = 45 * M_PI / 180;
			s->displacement = 5;
		}
		else if (i==2) {
			s->minPVal = 0.6;
			s->angle = 45 * M_PI / 180;
			s->displacement = 10;
		}
		s->angleInvHalfCos = cos(s->angle/2);
		for (len = MIN_FRAGMENT_LEN; len <= MAX_FRAGMENT_LEN; len++) {
			int pass;
			double min = 0, max = 5;
			for (pass = 0; pass < 50; pass++) {
				double mid = (min+max)*0.5;
				double score = RMSDScoreSimple(mid, len);
				if (score < s->minPVal) max = mid;
				else min = mid;
			}
			s->maxSqRMSD[len] = max*max;
		}
	}
}

/*double SequenceAlignmentScore(int p1, int p2, int count1, int count2, WeightedResiduePositions *r1, WeightedResiduePositions *r2, int len) {
	double score;
	int pos;
	for (pos=0; pos<len; pos++) {
		r1->res->residue;
	}
	return 0;
}//*/

typedef struct {
	Vector average;
	Vector com;
	double E0;
	double normalize;
} AlignmentCache;

/* Note: Doesn't handle non-existant residues. */
double AlignAlignmentBlocksCached(AlignmentCache *cache1, AlignmentCache *cache2, int len, Matrix *m, double maxSqRMSD, double twoOverLength) {
	Matrix A = {0,0,0,0,0,0,0,0,0,0,0,0};

	Matrix m2, m3;

	double E0;

	double sqRmsd;

	int i, j;

	for (i=0; i<len;  i++) {
		Vector comCoords1, comCoords2;


		subVect(&comCoords1, &cache1[i].average, &cache1->com);
		subVect(&comCoords2, &cache2[i].average, &cache2->com);

		for (j=0; j<3; j++) {
			A.M[j][0] += comCoords2.v[0] * comCoords1.v[j];
			A.M[j][1] += comCoords2.v[1] * comCoords1.v[j];
			A.M[j][2] += comCoords2.v[2] * comCoords1.v[j];
		}
	}

	E0 = cache1->E0 + cache2->E0;
	if (maxSqRMSD > 0) {
		Matrix ATA;
		sqRmsd = CalcSqRMSDUnderestimate(&ATA, &A, E0, twoOverLength);
		if (sqRmsd > maxSqRMSD) {
			/* Don't bother with rotation. */
			return sqRmsd;
		}
		sqRmsd = CalcRotation(&m2, &A, &ATA, E0, twoOverLength);
	}
	else {
		sqRmsd = CalcSqRMSDAndRotation(&m2, &A, E0, twoOverLength);
	}

	m2.M[0][3] = cache1->com.v[0];
	m2.M[1][3] = cache1->com.v[1];
	m2.M[2][3] = cache1->com.v[2];

	m3 = identity;

	m3.M[0][3] = -cache2->com.v[0];
	m3.M[1][3] = -cache2->com.v[1];
	m3.M[2][3] = -cache2->com.v[2];

	mulMat(m, &m2, &m3);

	return sqRmsd;
}

AlignmentCache *CalcAlignmentCache(int start, int count, WeightedResiduePositions *r, int len, int fragmentLen) {
	int i, j;
	double oneOverFragmentLen = 1/(double)fragmentLen;
	AlignmentCache *cache = (AlignmentCache*) malloc(sizeof(AlignmentCache) * len);
	for (i=0; i < len; i++) {
		double normalize = 0;
		Vector average = {0,0,0};
		for (j=0; j<count; j++) {
			int index = start+i;
			if (!r[j].res[index].exists) {
				/* Shouldn't ever happen, at the moment, so normalize should be the same everywhere. */
				continue;
			}
			normalize += r[j].weight;
			addMulVect(&average, &average, &r[j].res[index].coords, r[j].weight);
		}
		if (normalize) {
			cache[i].normalize = 1/normalize;
			mulVect(&average, &average, cache[i].normalize);
		}
		cache[i].average = average;
	}
	for (i=0; i <= len-fragmentLen; i++) {
		int k = 0;
		Vector com = {0,0,0};
		double E0 = 0;
		for (k=0; k < fragmentLen; k++) {
			addVect(&com, &com, &cache[i+k].average);
		}
		mulVect(&cache[i].com, &com, oneOverFragmentLen);
		for (k=0; k < fragmentLen; k++) {
			double E0delta = 0;
			for (j=0; j<count; j++) {
				int index = start+i+k;
				if (r[j].res[index].exists) {
					Vector temp1;
					subVect(&temp1, &r[j].res[index].coords, &cache[i].com);
					E0delta += r[j].weight * dotVect(&temp1, &temp1);
				}
			}
			E0 += E0delta * cache[i].normalize;
		}
		cache[i].E0 = E0 * 0.5;
	}
	return cache;
}

double AlignAlignmentBlocks(int p1, int p2, int count1, int count2, WeightedResiduePositions *r1, WeightedResiduePositions *r2, int len, Matrix *m) {
	AlignmentCache *cache1 = CalcAlignmentCache(p1, count1, r1, len, len);
	AlignmentCache *cache2 = CalcAlignmentCache(p2, count2, r2, len, len);
	double sqRmsd;
	
	sqRmsd = AlignAlignmentBlocksCached(cache1, cache2, len, m, -1, 2/(double)len);
	free(cache1);
	free(cache2);
	return sqRmsd;
}

double CalcAlignmentRMSD(MultipleAlignment *ma) {
	double total = 0;
	int i, j, k;
	if (ma->numChains <= 1) return ma->rmsd = 0;

	for (i=0; i<ma->numResidues; i++) {
		double weights = 0;
		double sum = 0;
		for (j=0; j<ma->numChains; j++) {
			if (ma->residues[j].res[i].hasCoords) {
				for (k=j+1; k<ma->numChains; k++) {
					if (ma->residues[k].res[i].hasCoords) {
						/* double w = residues[j].weight * residues[k].weight; */
						double w = 1;
						Vector v;
						sum += w * lengthSquaredVect(subVect(&v, &ma->residues[j].res[i].coords, &ma->residues[k].res[i].coords));
						weights += w;
					}
				}
			}
		}
		sum /= weights;
		total += sum;
	}
	ma->rmsd = sqrt(total/ma->numResidues);
	return ma->rmsd;
}

double GaussianIntegral(double x) {
	int swap = 0, steps = 1;
	double vals[2] = {-2,-2};
	if (x > 0) {
		x = -x;
		swap = 1;
	}


	do {
		double C, v = 0.5;
		int i;

		steps *= 2;
		C = fabs(x/steps) * (M_1_SQRT2PI/6);


		for (i=0; i<steps; i++) {
			double a = x*i/steps;
			double b = x*(i+1)/steps;
			double m = (a+b)/2;
			v -= C * (exp(a*a*-0.5) + 4 * exp(m*m*-0.5) + exp(b*b*-0.5));
		}
		vals[0] = vals[1];
		vals[1] = v;
	}
	while (fabs(vals[1] - vals[0]) > 0.00000001 && steps < 1000);

	/* Shouldn't happen, in general, but just in case... */
	if (vals[1] < 0) vals[1] = 0;
	if (swap) vals[1] = 1.0-vals[1];
	return vals[1];
}
double CalcAlignmentPvalue(MultipleAlignment *ma) {
	CalcAlignmentRMSD(ma);
	return ma->pvalue = GaussianIntegral((ma->rmsd - ma->numResidues * 0.155 + 1.6018)/2.13065);
}

double CalcStructalScore(MultipleAlignment *ma) {
	double score = 0;
	int i, j, k;
	for (i = 0; i<ma->numResidues; i++) {
		double weights = 0;
		double sum = 0;
		int gapCount = 0;
		if (i) {
			j = ma->numChains;
			while (j) {
				j--;
				if ((!ma->residues[j].res[i].exists && ma->residues[j].res[i-1].exists) || ma->residues[j].res[i].index-1 != ma->residues[j].res[i-1].index) {
					gapCount ++;
				}
			}
			/* Funky way of getting a number for multiple alignments.  Only currently used for pairwise, however. */
			if (ma->numChains == gapCount) gapCount --;
			score -= 10 * (gapCount+1) / ma->numChains;
		}

		for (j=0; j<ma->numChains; j++) {
			if (ma->residues[j].res[i].hasCoords) {
				for (k=j+1; k<ma->numChains; k++) {
					if (ma->residues[k].res[i].hasCoords) {
						/* double w = residues[j].weight * residues[k].weight; */
						double w = 1;
						Vector v;
						double d2 = lengthSquaredVect(subVect(&v, &ma->residues[j].res[i].coords, &ma->residues[k].res[i].coords));
						/* RMSD contribution to structal score. */
						sum += w * 20 / (1 + d2/(2.24*2.24));
						weights += w;
					}
				}
			}
		}
		score += sum / weights;
	}
	return score;
}

ResiduePositions *Distill(PDBChain *c, int id) {
	ResiduePositions *out = (ResiduePositions*) malloc(sizeof(ResiduePositions));
	int i;
	out->res = (ResiduePosition*) malloc(sizeof(ResiduePosition)*c->length);
	out->id = id;
	out->pdb = c;
	out->length = c->length;
	for (i=0; i<c->length; i++) {
		Atom *a = GetAtom(c, i, " CA ");
		if (a) {
			out->res[i].coords = a->pos;
			out->res[i].hasCoords = 1;
		}
		else {
			out->res[i].hasCoords = 0;
			out->res[i].coords.x = 0;
			out->res[i].coords.y = 0;
			out->res[i].coords.z = 0;
		}
		out->res[i].residue = ShortNames[c->seq->seq[i]];
		out->res[i].index = i;
	}

	return out;
}

void CalcWeights(MultipleAlignment *ma, int exclude) {
	int max = ma->numChains;
	int i;
	if (exclude<0) exclude = ma->numChains;
	else if (exclude < ma->numChains) max--;
	for (i=0; i<ma->numChains; i++) {
		if (i == exclude) ma->residues[i].weight = 0;
		else ma->residues[i].weight = 1.0/max;
	}
}

void RealignChain(MultipleAlignment *ma, int c, Matrix *matrices) {
	int i, j, k;
	WeightedResiduePositions tempResidues;

	CalcWeights(ma, c);

	tempResidues.weight = 1;
	tempResidues.res = ma->residues[c].res;
	for (i=0; i<ma->numResidues; i++) {
		tempResidues.res[i].coords = ma->chains[c]->res[tempResidues.res[i].index].coords;
	}
	for (j=0; j<ma->numBlocks; j++) {
		AlignAlignmentBlocks(ma->blocks[j].first, ma->blocks[j].first, ma->numChains, 1, ma->residues, &tempResidues, ma->blocks[j].last-ma->blocks[j].first+1, &matrices[j]);
		for (k=ma->blocks[j].first; k<=ma->blocks[j].last; k++) {
			transformVect(&tempResidues.res[k].coords, &matrices[j], &ma->chains[c]->res[tempResidues.res[k].index].coords);
		}
	}
}

Matrix ** Realign(MultipleAlignment *ma) {
	int block, k, chain;
	Matrix ** matrices = (Matrix **) malloc(sizeof(Matrix*) * ma->numChains);
	matrices[0] = (Matrix *) malloc(ma->numChains*ma->numBlocks * sizeof(Matrix));
	if (ma->numChains == 1) {
		for (chain=0; chain<ma->numBlocks; chain++) {
			matrices[0][chain] = identity;
		}
		ma->residues[0].weight = 1;
		return matrices;
	}
	for (chain=1; chain<ma->numChains; chain++) {
		matrices[chain] = matrices[chain-1] + ma->numBlocks;
	}

	RealignChain(ma, 0, matrices[0]);
	for (block=0; block<ma->numBlocks; block++) {
		Matrix m;
		invertMat(&m, &matrices[0][block]);
		for (k=ma->blocks[block].first; k<=ma->blocks[block].last; k++) {
			for (chain=0; chain<ma->numChains; chain++) {
				Vector v = ma->residues[chain].res[k].coords;
				transformVect(&ma->residues[chain].res[k].coords, &m, &v);
			}
		}
	}

	/* Recalculate alignments.of all blocks. */
	for (k=0; k<2; k++) {
		for (chain=ma->numChains-1; chain>=0; chain--) {
			RealignChain(ma, chain, matrices[chain]);
		}
	}

	/* Calculate new (And final) weights. */
	CalcWeights(ma, -1);
	return matrices;
}

void RealignExtendAndWeight(MultipleAlignment *ma) {
	/* Rotate alignment to coordinates of first chain. */
	int i, j, k;
	Matrix ** matrices = Realign(ma);
	Extend(ma, matrices, 1);

	free(matrices[0]);
	free(matrices);
	ma->conflictMap = (int**) malloc(sizeof(int*)*ma->numResidues);
	for (i=0; i<ma->numResidues; i++) {
		ma->conflictMap[i] = 0;
		for (j=0; j <= i; j++) {
			for (k = 0; k<ma->numChains; k++) {
				if (ma->residues[k].res[j].index >= ma->residues[k].res[i].index)
					break;
			}
			if (k != ma->numChains) {
				if (!ma->conflictMap[i]) {
					ma->conflictMap[i] = (int*) malloc(sizeof(int) * (1+i-j));
					ma->conflictMap[i][0] = j;
				}
				else {
					ma->conflictMap[i][j-ma->conflictMap[i][0]] = 1;
				}
			}
			else if (ma->conflictMap[i]) ma->conflictMap[i][j-ma->conflictMap[i][0]] = 0;
		}
	}
}












void TransformPDB(MultipleAlignment *ma) {
	/* Realign and rotate alignment to coordinates of first chain. */
	Matrix **matrices = Realign(ma);

	int i;
	/* Rotate residues of each chain. */
	for (i=0; i<ma->numChains; i++) {
		int block = 0;
		int pos = 0;
		Atom *atoms = ma->chains[i]->pdb->atoms;
		int numAtoms = ma->chains[i]->pdb->numAtoms;
		while (pos < numAtoms) {
			int newGroup = block+1;
			Vector v;
			if (newGroup<ma->numBlocks) {
				int temp = ma->blocks[newGroup].first;
				while (temp < ma->numResidues) {
					if (ma->residues[i].res[temp].exists) {
						while (newGroup<ma->numBlocks &&
							temp > ma->blocks[newGroup].last) {
							newGroup++;
						}
						break;
					}
					else {
						temp++;
					}
				}
				if (newGroup<ma->numBlocks && temp < ma->numResidues && ma->residues[i].res[temp].index <= atoms[pos].resIndex) {
					block = newGroup;
				}
			}
			v = atoms[pos].pos;
			transformVect(&atoms[pos].pos, &matrices[i][block], &v);
			pos ++;
		}
	}
	free(matrices[0]);
	free(matrices);
}

MultipleAlignment *DuplicateAlignment(const MultipleAlignment *srcMA, int first){
	/* Keep score and all length counts.
	 * Need to duplicate the assembly order and all arrays, except averages and
	 * conflictMap, which are calculated at the end, by another function. */
	WeightedResiduePositions tempResidues;
	ResiduePositions *tempChain;
	MultipleAlignment *ma = (MultipleAlignment*) memdup(srcMA, sizeof(MultipleAlignment));
	int i;

	ma->order = DuplicateAssemblyOrder(srcMA->order);
	ma->chains = (ResiduePositions **)memdup(srcMA->chains, srcMA->numChains * sizeof(ResiduePositions *));
	ma->blocks = (AlignedBlock *) memdup(srcMA->blocks, srcMA->numBlocks * sizeof(AlignedBlock));

	ma->conflictMap = 0;

	ma->residues = (WeightedResiduePositions*) malloc(srcMA->numChains*sizeof(WeightedResiduePositions));
	for (i=0; i<srcMA->numChains; i++) {
		ma->residues[i].res = (ResiduePosition*) memdup(srcMA->residues[i].res, srcMA->numResidues * sizeof(ResiduePosition));
		ma->residues[i].weight = srcMA->residues[i].weight;
	}

	/* Move the specified chain to the front. First one is the base alignment, which is (basically) unbent. */
	tempChain = ma->chains[0];
	ma->chains[0] = ma->chains[first];
	ma->chains[first] = tempChain;

	tempResidues = ma->residues[0];
	ma->residues[0] = ma->residues[first];
	ma->residues[first] = tempResidues;


	RealignExtendAndWeight(ma);
	return ma;
}


MultipleAlignment *CreateSingleStrandAlignment(ResiduePositions *chain) {
	MultipleAlignment *ma = (MultipleAlignment*) calloc(1, sizeof(MultipleAlignment));
	int i;
	ma->order = MakeAssemblyOrder(chain->id);
	ma->numChains = 1;
	ma->chains = (ResiduePositions **)memdup(&chain, sizeof(ResiduePositions *));
	ma->residues = (WeightedResiduePositions*) malloc(sizeof(WeightedResiduePositions));

	ma->numResidues = chain->length;
	ma->residues[0].res = (ResiduePosition*) malloc(ma->numResidues * sizeof(ResiduePosition));
	ma->numBlocks = 0;
	ma->blocks = 0;
	for (i=0; i<ma->numResidues; i++) {
		ma->residues[0].res[i] = chain->res[i];
		if (chain->res[i].exists) {
			if (!ma->numBlocks || ma->blocks[ma->numBlocks-1].last != i-1) {
				ma->blocks = (AlignedBlock *) realloc(ma->blocks, (++ma->numBlocks)*sizeof(AlignedBlock));
				ma->blocks[ma->numBlocks-1].first =
					ma->blocks[ma->numBlocks-1].last = i;
			}
			else {
				ma->blocks[ma->numBlocks-1].last = i;
			}
		}
	}
	ma->residues[0].weight = 1;

	return ma;
}



// Note:  Only use "average" in AlignmentCache, so don't need one for each length.
MatchedBlockPair *FindBestSet (MatchedBlockPair *matches, int numMatches, MultipleAlignment *a1, MultipleAlignment *a2, AlignmentCache *cache1, AlignmentCache *cache2, ScoreData *scoreData) {
	ScoreData s = *scoreData;
	MemoryManagementInfo info = {0,0};
	MatchedBlockPair **list = (MatchedBlockPair **) malloc((numMatches+MAX_OCT_TREE_VERTICES) * sizeof(MatchedBlockPair *));
	double displacementRange = s.displacement * 2;
	int best = 0;
	int i;
	int pos;
	OctTreeNode *trees = (OctTreeNode *) calloc(a2->numResidues, sizeof(OctTreeNode));
	double maxBonus = MAX_BONUS;
	double maxPosBonus = MAX_POS_BONUS;
	//int sum = 0;
	//int sum2 = 0;
	if (!numMatches) return 0;
	for (i=numMatches-1; i>=0; i--) {
		int pos = matches[i].p2 + matches[i].len-1;
		if (trees[pos].radius2 == 0) {
			// Flag value.  Will be overwritten later.
			trees[pos].radius2 = 1;
			trees[pos].corners[0] = trees[pos].corners[1] = matches[i].last;
		}
		else {
			if (matches[i].last.x < trees[pos].corners[0].x) {
				trees[pos].corners[0].x = matches[i].last.x;
			}
			else if (matches[i].last.x > trees[pos].corners[1].x) {
				trees[pos].corners[1].x = matches[i].last.x;
			}

			if (matches[i].last.y < trees[pos].corners[0].y) {
				trees[pos].corners[0].y = matches[i].last.y;
			}
			else if (matches[i].last.y > trees[pos].corners[1].y) {
				trees[pos].corners[1].y = matches[i].last.y;
			}

			if (matches[i].last.z < trees[pos].corners[0].z) {
				trees[pos].corners[0].z = matches[i].last.z;
			}
			else if (matches[i].last.z > trees[pos].corners[1].z) {
				trees[pos].corners[1].z = matches[i].last.z;
			}
		}
	}
	for (pos=0; pos<a2->numResidues; pos++) {
		InitOctTreeNode(&trees[pos], &trees[pos].corners[0], &trees[pos].corners[1]);
	}

	for (i=0; i<numMatches; i++) {
		Vector v, first;
		MatchedBlockPair *match = &matches[i];
		int *conflictMapA1 = a1->conflictMap[match->p1];
		int conflictMapA10_1 = conflictMapA1[0]+1;
		int *conflictMapA2 = a2->conflictMap[match->p2];
		int conflictMapA20 = conflictMapA2[0];
		match->bestPreviousMatch = -1;
		/* In this case, don't allow it to be added without a previous strand,
		 * unless it's the first strand.
		 */
		if (s.firstStrandLocked) {
			if (i) {
				match->bestAssembledScore = -5000;
			}
			else {
				/* Otherwise, keep previous assembled score.  Means I can
				 * just insert into old alignments without changing scores.
				 */
				match->bestAssembledScore -= match->score;
			}
		}
		else
			match->bestAssembledScore = 0;
		transformVect(&first, &match->m, &cache2[match->p2].average);
		for (pos=MIN_FRAGMENT_LEN-1; pos<match->p2; pos++) {
			int hits, j;
			// Skip empty trees.  Despite the fact that they just pass through all tests,
			// results in a fairly significant speedup (over 5%), at least on structures
			// that don't align well.
			if (!trees[pos].radius) continue;
			if (pos >= conflictMapA20) {
				// This test is only needed if pos is 0, which it isn't.
				// if (pos == conflictMapA20 ||
				if (conflictMapA2[pos - conflictMapA20])
					continue;
			}
			transformVect(&v, &match->m, &cache2[pos].average);

	//time[4] -= timeGetTime();
			hits = GetEntries2(list, (void**)(list+numMatches+MAX_OCT_TREE_VERTICES-1), &trees[pos], &v, displacementRange);
			//hits = GetEntries((void**)list, &trees[pos], &v, displacementRange);
	//time[4] += timeGetTime();

			//sum += hits;
			for (j=0; j < hits; j++) {
				Vector v1, v2;

				MatchedBlockPair *match2 = list[j];
				double score, invCosHalfTest, displacement;

				invCosHalfTest = invHalfCosQuat(&match->q, &match2->q);
				if (invCosHalfTest < s.angleInvHalfCos) {
					continue;
				}

				if (match2->p1 + match2->len >= conflictMapA10_1) {
					// This test is only needed if match2->p1 + match2->len is 0, which it isn't.
					// if (match2->p1 + match2->len == conflictMapA10_1 ||
					if (match2->p1 + match2->len > match->p1 ||
						conflictMapA1[match2->p1 + match2->len - conflictMapA10_1])
							continue;
				}

				score = match2->bestAssembledScore;

				if (score + maxBonus < match->bestAssembledScore) continue;

				/*
				angle = 0;
				if (invCosHalfTest<1) angle = 2*acos(invCosHalfTest);
				score += AngleScore(angle);
				//*/

				score += AngleScoreHalfCos(invCosHalfTest);

				if (score + maxPosBonus < match->bestAssembledScore) continue;

				subVect(&v1, &match2->last, &v);
				transformVect(&v2, &match2->m, &cache2[match->p2].average);
				subVect(&v2, &first,  &v2);

				displacement = (lengthVect(&v1) + lengthVect(&v2))/2;
				if (displacement > s.displacement) continue;
				score += DisplacementScore(displacement);

				if (score < match->bestAssembledScore) continue;

				match->bestAssembledScore = (float) score;
				match->bestPreviousMatch = (int)(match2 - matches);
			}
		}

		if (match->bestAssembledScore < 0) continue;

		match->bestAssembledScore += match->score;
		if (match->bestAssembledScore > matches[best].bestAssembledScore)
			best = i;
		pos = match->p2 + match->len-1;
		AddEntry(&trees[pos], match, &info, displacementRange);
	}

	FreeAllMemory(&info);
	free(trees);
	free(list);

	if (!s.lastStrandLocked)
		return &matches[best];
	else
		return &matches[numMatches-1];
}


/* One of the two main functions - Finds highest scoring set of consistent block pairs from a set,
 * using the given cutoffs.  The first and/or last block can be locked.
 */
int FindInterveningRegions (AlignmentCache **cache1, AlignmentCache **cache2, MatchedBlockPair **regions, int *numRegions, int index, MultipleAlignment *a1, MultipleAlignment *a2, ScoreData *s) {
	double ratios[MAX_FRAGMENT_LEN+1];
	double twoOverLength[MAX_FRAGMENT_LEN+1];
	double sqRmsd;
	int extend = 0;

	/* Note the extra one on the end - means don't have to reallocate for adding a locked strand
	 * at the end, if there is one.
	 */
	MatchedBlockPair *matches = (MatchedBlockPair *) malloc(10241 * sizeof(MatchedBlockPair));
	MatchedBlockPair *match;
	int maxMatches = 10240;
	int numMatches = 0;

	int i, block1, p1, block2, p2, len;
	//time[2] -= timeGetTime();

	for (i=1; i<=MAX_FRAGMENT_LEN; i++) {
		ratios[i] = ((i-1.0)/i)*((i-1.0)/i);
		twoOverLength[i] = 2.0/i;
	}

	s->firstStrandLocked = 0;
	s->lastStrandLocked = 0;

	if (index>=0) {
		s->firstStrandLocked = 1;
		matches[numMatches++] = regions[0][index];
	}
	if (index < *numRegions-1) {
		s->lastStrandLocked = 1;
	}

	for (block1 = 0; block1 < a1->numBlocks; block1++) {
		int maxGroup2 = a2->numBlocks;
		if (s->selfAlignDistance) maxGroup2 = block1+1;
		for (p1 = a1->blocks[block1].first; p1 <= a1->blocks[block1].last - MIN_FRAGMENT_LEN + 1; p1++) {
			if (s->firstStrandLocked && regions[0][index].p1 + regions[0][index].len > p1) {
				p1 = regions[0][index].p1 + regions[0][index].len;
				continue;
			}
			if (s->lastStrandLocked && regions[0][index+1].p1 < p1 + MIN_FRAGMENT_LEN) break;

			for (block2 = 0; block2 < maxGroup2; block2++) {
				int maxp2 = a2->blocks[block2].last - MIN_FRAGMENT_LEN + 1;
				if (s->selfAlignDistance && maxp2  + s->selfAlignDistance >= p1)
					maxp2 = p1 - s->selfAlignDistance;
				if (s->lastStrandLocked && regions[0][index+1].p2 < maxp2 + MIN_FRAGMENT_LEN) {
					maxp2 = regions[0][index+1].p2 - MIN_FRAGMENT_LEN;
				}

				for (p2 = a2->blocks[block2].first; p2 <= maxp2; p2++) {
					if (s->firstStrandLocked && regions[0][index].p2 + MIN_FRAGMENT_LEN > p2) {
						p2 = regions[0][index].p2 + MIN_FRAGMENT_LEN;
						continue;
					}
					sqRmsd = 0;
					for (len = MIN_FRAGMENT_LEN; len <= MAX_FRAGMENT_LEN; len++) {
						/* Underestimate of new RMSD.  Simple optimization. */
						sqRmsd = sqRmsd * ratios[len];
						if (sqRmsd < s->maxSqRMSD[len]) {
							Matrix m;
							double score;

							/* Make sure we aren't extending beyond a source block's limits. */
							if (p1 + len-1 > a1->blocks[block1].last ||
								p2 + len-1 > a2->blocks[block2].last) break;

							sqRmsd = AlignAlignmentBlocksCached(cache1[len]+p1, cache2[len]+p2, len, &m, s->maxSqRMSD[len], twoOverLength[len]);
							if (sqRmsd < s->maxSqRMSD[len]) {
								double rmsd = sqrt(sqRmsd);
								score = RMSDScoreSimple(rmsd, len);

								if (score >= s->minPVal) {
									if (numMatches == maxMatches) {
										maxMatches += 10240;
										matches = (MatchedBlockPair *) realloc(matches, (1+maxMatches) * sizeof(MatchedBlockPair));
									}
									//score += SequenceAlignmentScore(p1, p2, a1->numChains, a2->numChains, a1->residues, a2->residues, len);

									match = &matches[numMatches];
									match->m = m;
									matrixToQuaternion(&match->q, &m);
									match->score = (float)(score);
									match->p1 = p1;
									match->p2 = p2;
									match->len = len;
									transformVect(&match->last, &m, &cache2[len][p2+len-1].average);
									numMatches++;
								}
							}
						}
					}
				}
			}
		}
	}
	if (s->lastStrandLocked == 1) {
		matches[numMatches] = regions[0][index+1];
		numMatches++;
	}
	matches = (MatchedBlockPair *) realloc(matches, numMatches * sizeof(MatchedBlockPair));
	//time[2] += timeGetTime();

	//time[3] -= timeGetTime();
	match = FindBestSet(matches, numMatches, a1, a2, cache1[MIN_FRAGMENT_LEN], cache2[MIN_FRAGMENT_LEN], s);
	//time[3] += timeGetTime();
	if (!match) {
		free(matches);
		return 0;
	}
	if (s->lastStrandLocked == 1) {
		regions[0][index+1] = *match;
		if (match->bestPreviousMatch >= 0) {
			match = &matches[match->bestPreviousMatch];
		}
		else match = 0;
	}
	while (match) {
		if (match->bestPreviousMatch == -1 && s->firstStrandLocked) {
			regions[0][index] = *match;
			break;
		}
		else {
			extend++;

			regions[0] = (MatchedBlockPair*) realloc(regions[0], sizeof(MatchedBlockPair) * (*numRegions+1));
			memmove(regions[0]+index+2, regions[0]+index+1, sizeof(MatchedBlockPair) * (*numRegions-index-1));
			regions[0][index+1] = *match;
			numRegions[0]++;
		}
		if (match->bestPreviousMatch == -1) break;
		match = &matches[match->bestPreviousMatch];
	}
	free(matches);
	return extend;
}

void FindRegionsPass(AlignmentCache **c1, AlignmentCache **c2, MatchedBlockPair **regions, int *numRegions, MultipleAlignment *a1, MultipleAlignment *a2, ScoreData* s) {
	int i = -1;
	//time[1] -= timeGetTime();
	for (; i < *numRegions; i++) {
		int p = FindInterveningRegions(c1, c2, regions, numRegions, i, a1, a2, s);
		i+=p;
	}
	//time[1] += timeGetTime();
}


void FindRegions (MatchedBlockPair **regions, int *numRegions, MultipleAlignment *a1, MultipleAlignment *a2, int selfAlignDistance) {
	ScoreData s;
	int i;
	AlignmentCache *c1[MAX_FRAGMENT_LEN+1], *c2[MAX_FRAGMENT_LEN+1];
	for (i = MIN_FRAGMENT_LEN; i <= MAX_FRAGMENT_LEN; i++) {
		c1[i] = CalcAlignmentCache(0, a1->numChains, a1->residues, a1->numResidues, i);
		c2[i] = CalcAlignmentCache(0, a2->numChains, a2->residues, a2->numResidues, i);
	}
	memset(&s, 0, sizeof(s));

	*regions = 0;
	*numRegions = 0;
	s.selfAlignDistance = selfAlignDistance;

	for (i=0; i<sizeof(ScoreTables)/sizeof(ScoreData); i++) {
		s = ScoreTables[i];

		FindRegionsPass(c1, c2, regions, numRegions, a1, a2, &s);
	}
	for (i = MIN_FRAGMENT_LEN; i <= MAX_FRAGMENT_LEN; i++) {
		free(c1[i]);
		free(c2[i]);
	}
}

MultipleAlignment* AlignAlignments(MultipleAlignment *originalSrcMA1, MultipleAlignment *originalSrcMA2, int index1, int index2, int selfAlignDistance) {
	int i,j,k;

	/* New multiple alignment */
	MultipleAlignment *ma = (MultipleAlignment *)calloc(1, sizeof(MultipleAlignment));

	/* Duplicate source alignments. Moves the index1 chain to the front.
	 * Also realigns and extends the aligned blocks. */
	MultipleAlignment *srcMA1 = DuplicateAlignment(originalSrcMA1, index1);
	MultipleAlignment *srcMA2 = DuplicateAlignment(originalSrcMA2, index2);

	MatchedBlockPair *regions;
	int numRegions;

	ma->order = CombineAssemblyOrders(srcMA1->order, srcMA2->order);

	ma->numChains = srcMA1->numChains + srcMA2->numChains;
	ma->chains = (ResiduePositions **)malloc(ma->numChains * sizeof(ResiduePositions *));
	memcpy(ma->chains, srcMA1->chains, srcMA1->numChains * sizeof(ResiduePositions *));
	memcpy(ma->chains + srcMA1->numChains, srcMA2->chains, srcMA2->numChains * sizeof(ResiduePositions *));

	FindRegions(&regions, &numRegions, srcMA1, srcMA2, selfAlignDistance);

	for (i=0; i < numRegions; i++) {
		ma->numResidues += regions[i].len;
	}

	ma->score = -5000;
	ma->residues = (WeightedResiduePositions*) malloc(ma->numChains*sizeof(WeightedResiduePositions));
	for (i=0; i<ma->numChains; i++) {
		ma->residues[i].res = (ResiduePosition*) malloc(ma->numResidues*sizeof(ResiduePosition));
		ma->residues[i].weight = 0;
	}

	/* Alignment found (If no alignment, will return what we already have). */
	if (numRegions) {
		int pos = 0;
		ma->score = regions[numRegions-1].bestAssembledScore;

		ma->numBlocks = numRegions;
		ma->blocks = (AlignedBlock *) malloc(numRegions * sizeof(AlignedBlock));

		for (i=0; i < numRegions; i++) {
			ma->blocks[i].first = pos;
			ma->blocks[i].last = pos + regions[i].len-1;
			pos += regions[i].len;
			for (j=0; j<regions[i].len; j++) {
				for (k=0; k<srcMA1->numChains; k++) {
					ma->residues[k].res[ma->blocks[i].first+j] = srcMA1->residues[k].res[regions[i].p1 + j];
				}
				for (; k<ma->numChains; k++) {
					ma->residues[k].res[ma->blocks[i].first+j] = srcMA2->residues[k-srcMA1->numChains].res[regions[i].p2 + j];
					transformVect(&ma->residues[k].res[ma->blocks[i].first+j].coords, &regions[i].m, &srcMA2->residues[k-srcMA1->numChains].res[regions[i].p2 + j].coords);
				}
			}
		}
		free(regions);
	}
	CleanupAlignment(srcMA1);
	CleanupAlignment(srcMA2);
	return ma;
}

double AlignOrder (MultipleAlignment *ma, AssemblyNode *n, int *count, WeightedResiduePositions *r, int *indices, Matrix *matrices) {
	int i, j, k, chain;
	/* If only one chain, just find the index in ma->residues and copy the data over. */
	if (n->left == 0) {
		int index;
		for (index = 0; index < ma->numChains-1; index++) {
			if (ma->chains[index]->id == n->id) break;
		}
		*indices = index;
		r[0].res = (ResiduePosition*) memdup(ma->residues[index].res, sizeof(ResiduePosition)*ma->numResidues);
		r[0].weight = 1;
		*count = 1;
		matrices[index] = identity;
		return 0;
	}
	else {
		Matrix m;
		WeightedResiduePositions tempResidues;
		WeightedResiduePositions *r2;
		double test;
		int count1, count2;

		/* Get aligned sets of residues from both children. */
		AlignOrder(ma, n->left, &count1, r, indices, matrices);

		r2 = r+count1;
		AlignOrder(ma, n->right, &count2, r2, indices+count1, matrices);
		*count = count1+count2;
		for (i=0; i<count1; i++)
			r[i].weight = 1 / (double) count1;
		for (i=0; i<count2; i++)
			r2[i].weight = 1 / (double) count2;

		/* Align both sets of residues. */
		test = sqrt(AlignAlignmentBlocks(0, 0, count1, count2, r, r2, ma->numResidues, &m));
		for (chain = 0; chain<count2; chain++) {
			for (k=0; k<ma->numResidues; k++) {
				if (r2[chain].res[k].exists) {
					Vector v = r2[chain].res[k].coords;
					transformVect(&r2[chain].res[k].coords, &m, &v);
				}
			}
		}

		/* Realign all chains. */
		tempResidues.weight = 1;

		count1 = *count;
		for (j=0; j<2; j++) {
			for (chain = count1-1; chain>=0; chain--) {
				for (i=0; i<count1; i++)
					r[i].weight = 1 / (double) (count1-1);
				r[chain].weight = 0;

				tempResidues.res = r[chain].res;
				for (k=0; k<ma->numResidues; k++) {
					tempResidues.res[k].coords = ma->chains[indices[chain]]->res[tempResidues.res[k].index].coords;
				}
				AlignAlignmentBlocks(0, 0, count1, 1, r, &tempResidues, ma->numResidues, &m);
				matrices[indices[chain]] = m;
				for (k=0; k<ma->numResidues; k++) {
					transformVect(&tempResidues.res[k].coords, &m, &ma->chains[indices[chain]]->res[tempResidues.res[k].index].coords);
				}
			}
		}
		return test;
	}
}



int IsCoreResidue(MultipleAlignment *ma, int pos) {
  int chain;
  for (chain=0; chain<ma->numChains; chain++) {
    if (!ma->residues[chain].res[pos].hasCoords) return 0;
  }
  return 1;
}




void UnbentRMSDAlign(MultipleAlignment *ma) {
	WeightedResiduePositions *residuesTemp = (WeightedResiduePositions *)malloc(ma->numChains * sizeof(WeightedResiduePositions));
	int* indices = (int*) malloc(sizeof(int) * ma->numChains);
	Matrix *matrices = (Matrix*) malloc(sizeof(Matrix) * ma->numChains);
	int i, j;
	for (i=0; i<ma->numChains; i++) {
		for (j=0; j<ma->numResidues; j++) {
			ma->residues[i].res[j] = ma->chains[i]->res[ma->residues[i].res[j].index];
		}
	}
	AlignOrder(ma, ma->order->root, &i, residuesTemp, indices, matrices);

	for (i=0; i<ma->numChains; i++) {
		RotateChain(ma->chains[i]->pdb, &matrices[i]);
	}
	for (i=0; i<ma->numChains; i++) {
		for (j=0; j<ma->chains[i]->length; j++) {
			Atom *a = GetAtom(ma->chains[i]->pdb, j, " CA ");
			if (a)
				ma->chains[i]->res[j].coords = a->pos;
		}
	}
	for (i=0; i<ma->numChains; i++) {
		for (j=0; j<ma->numResidues; j++) {
			ma->residues[i].res[j] = ma->chains[i]->res[ma->residues[i].res[j].index];
		}
		free(residuesTemp[i].res);
	}
	free(indices);
	free(residuesTemp);
	free(matrices);
}

double CalcNonWeightedRMSD(MultipleAlignment *ma) {
	double total = 0;
	int i, j, k;
	if (ma->numChains <= 1) return ma->rmsd = 0;

	for (i=0; i<ma->numResidues; i++) {
		int weights = 0;
		double sum = 0;
		for (j=0; j<ma->numChains; j++) {
			if (ma->residues[j].res[i].hasCoords) {
				for (k=j+1; k<ma->numChains; k++) {
					if (ma->residues[k].res[i].hasCoords) {
						Vector v;
						sum += lengthSquaredVect(subVect(&v, &ma->residues[j].res[i].coords, &ma->residues[k].res[i].coords));
						weights += 1;
					}
				}
			}
		}
		sum /= weights;
		total += sum;
	}
	return ma->rmsd = sqrt(total/ma->numResidues);
}


void CleanupAlignment(MultipleAlignment *ma) {
	int i;

	free(ma->chains);
	free(ma->order);
	free(ma->blocks);
	if (ma->residues) {
		for (i=0; i<ma->numChains; i++) {
			free(ma->residues[i].res);
		}
		free(ma->residues);
	}
	if (ma->conflictMap) {
		for (i=0; i<ma->numResidues; i++) {
			free(ma->conflictMap[i]);
		}
		free(ma->conflictMap);
	}
	free(ma);
}






MultipleAlignment* Align(PDBChain **chains, int numChains, int displayStatus, int singleThreaded) {
	int i, j;

	int numAlignments = numChains;
	MultipleAlignment** alignments = (MultipleAlignment**) malloc(numAlignments * sizeof(MultipleAlignment*));
	int numAlignedAlignments = numChains * (numChains-1)/2;
	MultipleAlignment** alignedAlignments = (MultipleAlignment**) malloc(numAlignedAlignments * sizeof(MultipleAlignment*));
	double *pairScores;
	int numTotalAlignments = (numChains-1) * (numChains-1);
	int len = 1;
	int d = 10;
	MultipleAlignment *ma;

	/* Indices of chains used in pairwise alignments. Simplest way to divide up work
	 * for the first pass. */
	int *indices = (int*) malloc(sizeof(int) * 2 * numTotalAlignments);

	/* Used for output */
	char progressTemplate[64];
	int backspace2 = 0;
	int backspace = 0;
	/* Keeps track of which alignment we're on. For display only. */
	int count = 0;
	int pass = 1;

	for (i=0; i<numChains; i++) {
		alignments[i] = CreateSingleStrandAlignment(Distill(chains[i], i));
	}
	pairScores = (double*) malloc(sizeof(double) * numChains*numChains);


	while (numTotalAlignments >= d) {
		len ++;
		d*=10;
	}

	/* Structure pair list.  Simplest way of breaking up first pass for threading */
	for (i=0; i<numChains; i++) {
		for (j=0; j<i; j++) {
			indices[count*2] = i;
			indices[count*2+1] = j;
			count++;
		}
	}

	count = 0;

	//time[0] -= timeGetTime();

	if (displayStatus) {
		backspace2 = fprintf(stderr, "Pass %i, ", pass);
	}

	/* Single threaded version for OpenMP build. */
#ifdef _OPENMP
	if (singleThreaded)
	{
		if (displayStatus) {
			int threads = omp_get_num_threads();
			if (threads > 1) sprintf(progressTemplate, "%i threads [%%%ii/%i]", threads, len, numTotalAlignments);
			else sprintf(progressTemplate, "1 thread [%%%ii/%i]", len, numTotalAlignments);
			DisplayProgress(&backspace, progressTemplate, count);
		}
		for (i=0; i<numAlignedAlignments; i++) {
			alignedAlignments[i] = AlignAlignments(alignments[indices[i*2]], alignments[indices[i*2+1]], 0, 0, 0);
			pairScores[indices[i*2]*numChains + indices[i*2+1]] =
				pairScores[indices[i*2+1]*numChains + indices[i*2]] = alignedAlignments[i]->score;
			if (displayStatus) {
				DisplayProgress(&backspace, progressTemplate, ++count);
			}
		}
	}
	else
#endif

#ifndef _OPENMP	
	{
		if (displayStatus) {
			sprintf(progressTemplate, "Single threaded [%%%ii/%i]", len, numTotalAlignments);
			DisplayProgress(&backspace, progressTemplate, count);
		}
#else
#pragma omp parallel default(shared)
	{
		#pragma omp single
		if (displayStatus) {
			int threads = omp_get_num_threads();
			if (threads > 1) sprintf(progressTemplate, "%i threads [%%%ii/%i]", threads, len, numTotalAlignments);
			else sprintf(progressTemplate, "1 thread [%%%ii/%i]", len, numTotalAlignments);
			DisplayProgress(&backspace, progressTemplate, count);
		}
		#pragma omp for schedule(dynamic, 1)
#endif
		for (i=0; i<numAlignedAlignments; i++) {
			alignedAlignments[i] = AlignAlignments(alignments[indices[i*2]], alignments[indices[i*2+1]], 0, 0, 0);
			pairScores[indices[i*2]*numChains + indices[i*2+1]] =
				pairScores[indices[i*2+1]*numChains + indices[i*2]] = alignedAlignments[i]->score;
			if (displayStatus) {
#ifdef _OPENMP
#pragma omp critical
#endif
				DisplayProgress(&backspace, progressTemplate, ++count);
			}
		}
	}
	free(indices);

	pass++;

	while (numAlignedAlignments) {
		int best = 0;
		int index1, index2;
		for (i=1; i<numAlignedAlignments; i++) {
			if (alignedAlignments[i]->score > alignedAlignments[best]->score) {
				best = i;
			}
		}
		ma = alignedAlignments[best];
		alignedAlignments[best] = alignedAlignments[--numAlignedAlignments];

		/* Remove anything containing the elements of ma from both
		 * sets of alignments.
		 */
		index1 = ma->chains[0]->id;
		index2 = ma->chains[ma->numChains-1]->id;
		for (i=numAlignedAlignments-1; i >= 0; i--) {
			for (j=0; j<alignedAlignments[i]->numChains; j++) {
				if (alignedAlignments[i]->chains[j]->id == index1 ||
					alignedAlignments[i]->chains[j]->id == index2) {

					CleanupAlignment(alignedAlignments[i]);
					alignedAlignments[i] = alignedAlignments[--numAlignedAlignments];
					break;
				}
			}
		}

		for (i=numAlignments-1; i >= 0; i--) {
			for (j=0; j<alignments[i]->numChains; j++) {
				if (alignments[i]->chains[j]->id == index1 ||
					alignments[i]->chains[j]->id == index2) {

					CleanupAlignment(alignments[i]);
					alignments[i] = alignments[--numAlignments];
					break;
				}
			}
		}

		if (displayStatus) {
			backspace2 += backspace;
			backspace = 0;
			DisplayProgress(&backspace2, "Pass %i, ", pass++);
			DisplayProgress(&backspace, progressTemplate, count);
		}

#ifdef _OPENMP
		if (singleThreaded) {
			/* Run all the new alignments and add ma to the list of alignments I've kept. */
			for (i=numAlignments-1; i >= 0; i--) {
				double bestScore = 0;
				int best1=0, best2=0;

				/* Need private (thread) variables. */
				int j, k;
				for (j=0; j<alignments[i]->numChains; j++) {
					for (k=0; k<ma->numChains; k++) {
						double score = pairScores[ma->chains[k]->id * numChains + alignments[i]->chains[j]->id];
						if (score > bestScore) {
							bestScore = score;
							best1 = j;
							best2 = k;
						}
					}
				}

				alignedAlignments[numAlignedAlignments+i] = AlignAlignments(alignments[i], ma, best1, best2, 0);
				if (displayStatus) {
					DisplayProgress(&backspace, progressTemplate, ++count);
				}
			}
		}
		else
#endif
		{
#ifdef _OPENMP
			#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
			/* Run all the new alignments and add ma to the list of alignments I've kept. */
			for (i=numAlignments-1; i >= 0; i--) {
				double bestScore = 0;
				int best1=0, best2=0;

				/* Need private (thread) variables. */
				int j, k;
				for (j=0; j<alignments[i]->numChains; j++) {
					for (k=0; k<ma->numChains; k++) {
						double score = pairScores[ma->chains[k]->id * numChains + alignments[i]->chains[j]->id];
						if (score > bestScore) {
							bestScore = score;
							best1 = j;
							best2 = k;
						}
					}
				}

				alignedAlignments[numAlignedAlignments+i] = AlignAlignments(alignments[i], ma, best1, best2, 0);
				if (displayStatus) {
#ifdef _OPENMP
#pragma omp critical
#endif
					DisplayProgress(&backspace, progressTemplate, ++count);
				}
			}
		}
		numAlignedAlignments += numAlignments;
		alignments[numAlignments++] = ma;
	}
	if (displayStatus) {
		backspace += backspace2;
		DisplayProgress(&backspace, "", 0);
	}
	free(pairScores);

	/* The only alignment left. */
	ma = alignments[0];

	free(alignments);
	free(alignedAlignments);

	//time[0] += timeGetTime();
	/*for (i=0; i<sizeof(time)/sizeof(time[0]); i++) {
		printf("%i\t%i\n", i, time[i]);
	}//*/
	return ma;
}

