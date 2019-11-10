/* Matt.c -- Multiple alignment program (Multiple Alignment with Translations and Twists).
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

#ifdef _DEBUG
/* Memory leak detection library. */
/* #include <vld.h> */
#endif
#include <stdarg.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#define THREADING "OpenMP"
#else
#define THREADING "Single Threaded"
#endif

#include "Score.h"
#include "Protein.h"
#include "MultipleAlignment.h"
#include "MultipleAlignmentOutput.h"
#include "Extend.h"
//#include <Windows.h>

#define OUTPUT_PDB		 0x01
#define OUTPUT_TXT		 0x02
#define OUTPUT_FASTA	 0x04
#define OUTPUT_SPT		 0x08
#define OUTPUT_MSF		 0x10
#define OUTPUT_PIR		 0x20
#define OUTPUT_SPDBV	 0x40
#define OUTPUT_SSI		 0x80
#define OUTPUT_NULL		0x100
#define OUTPUT_RMSD     0x200     

const char formats[][6] = {"pdb", "txt", "fasta", "spt", "msf", "pir", "spdbv", "ssi", "null", "rmsd"};

#define VERSION "1.02"

void Usage() {
	printf("\nMatt version " VERSION " " THREADING " build.\n"
		"See http://matt.csail.mit.edu for more information.\n\n"
		"Usage: Matt -o outprefix [-c cutoff] [-t threads] [-[rslbdvp][01]]*\n"
		"            [-f <extension list>] [-L listfile]* [file[:<chain list>]]*\n"
	);
}

void Feedback(FILE *out, char *prefix, char *string, va_list arg, char *suffix) {
	fprintf(out, prefix);

	vfprintf(out, string, arg);

	fprintf(out, suffix);
}

/* Not necessarily accurate counts in OpenMP mode, but doesn't really matter.
 */
int errors = 0;
int warnings = 0;

void WarningLog(int console, FILE *f2, char *string, ...) {
	va_list arg;
	warnings ++;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, " * ", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, " * ", string, arg, "\n");
		va_end(arg);
	}
}

void StatusLog(int console, FILE *f2, char *string, ...) {
	va_list arg;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, "   ", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, "   ", string, arg, "\n");
		va_end(arg);
	}
}

void ErrorLog(int console, FILE *f2, char *string, ...) {
	va_list arg;
	errors++;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, "** ", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, "** ", string, arg, "\n");
		va_end(arg);
	}
}

void RootStatusLog(int console, FILE *f2, char *string, ...) {
	va_list arg;

	if (console) {
		va_start(arg, string);
		Feedback(stdout, "", string, arg, "\n");
		va_end(arg);
	}

	if (f2) {
		va_start(arg, string);
		Feedback(f2, "", string, arg, "\n");
		va_end(arg);
	}
}

struct SUMMARY {
	int coreSize;
	int lengths[2];
	double coreBentRMSD;
	double coreRMSD;
	// negative if more than 2 structures.
	double pvalue;
	double structalScores[3];
};
typedef struct SUMMARY Summary;

struct JOB_INFO {
	char *outName;
	char **files;
	int numFiles;
};
typedef struct JOB_INFO JobInfo;

void DisplaySummary(JobInfo *job, Summary * summary) {
	printf("%s:\t%i\t%0.3lf", job->outName, summary->coreSize, summary->coreRMSD);
	if (summary->pvalue >= -0.00001) {
		if (summary->pvalue >= 0.01)
			printf("\t%0.4lf", summary->pvalue);
		else
			printf("\t%0.4e", summary->pvalue);
		printf("\t%0.2lf", summary->structalScores[0]);
		printf("\t%0.2lf", summary->structalScores[1]);
		printf("\t%0.2lf", summary->structalScores[2]);
		printf("\t%i\t%i", summary->lengths[0], summary->lengths[1]);
	}
	printf("\n");
}

int OutputAllButPdb(MultipleAlignment *ma, char *outfile, char *pdb, int outFormats, int consoleOutput, FILE *log, int ref) {
	char *suffix = outfile + strlen(outfile);
	int problems = 0;

	strcpy(suffix, ".spt");
	if ((outFormats & OUTPUT_SPT) && !OutputAlignmentRasmol(ma, outfile)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}

	strcpy(suffix, ".spdbv");
	if ((outFormats & OUTPUT_SPDBV) && !OutputAlignmentSpdbv(ma, outfile)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}

	strcpy(suffix, ".txt");
	if ((outFormats & OUTPUT_TXT) && !OutputAlignmentText(ma, outfile, ref)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}

	strcpy(suffix, ".fasta");
	if ((outFormats & OUTPUT_FASTA) && !OutputAlignmentFasta(ma, outfile)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}

	strcpy(suffix, ".msf");
	if ((outFormats & OUTPUT_MSF) && !OutputAlignmentMsf(ma, outfile)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}

	strcpy(suffix, ".pir");
	if ((outFormats & OUTPUT_PIR) && !OutputAlignmentPir(ma, outfile)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}

	strcpy(suffix, ".ssi");
	if ((outFormats & OUTPUT_SSI) && !OutputAlignmentSsi(ma, outfile)) {
		problems++;
		WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
	}
    
    strcpy(suffix, ".rmsd");
    if ((outFormats & OUTPUT_RMSD) && !OutputAlignmentRmsd(ma, outfile)) {
      problems++;
      WarningLog(consoleOutput, log, "Can't open %s for writing.", outfile);
    }
    

	return problems;
}

int CountAlphaCarbons(PDBChain *c) {
	int i;
	int count = 0;
	for (i=0; i<c->length; i++) {
		Atom * atom = GetAtom(c, i, " CA ");
		if (atom) count++;
	}
	return count;
}

int RunJob(JobInfo *job, int seqres, int verbose, int renumber, int relabel, int displayStatus, int outputBent, double rmsdCutoff, int consoleOutput, int makeLog, int parallel, int outFormats, int partial, Summary *summary, int displaySummary, int singleThreaded) {
	/* Names of files to write to. */
	char *pdbFile = (char*)malloc(sizeof(char) * (strlen(job->outName) + 50));
	char *otherFile = (char*)malloc(sizeof(char) * (strlen(job->outName) + 50));
	int problems = 0;
	int noAlignment = 0;

	int i, j;
	int numChains = 0;
	PDBChain **chains = 0;
	FILE *log = 0;
	if (makeLog) {
		sprintf(otherFile, "%s.log", job->outName);
		log = fopen(otherFile, "wb");

		if (!log)
			problems++;
	}
	{
		int len=strlen(job->outName);
		char *commandLine, *end;
		int first = 1;
		for (i=0; i<job->numFiles; i++) {
			len += 2 + strlen(job->files[i]);
		}
		commandLine = malloc(200 + len);
		end = commandLine + sprintf(commandLine, "Matt -o %s -s%i -v%i -r%i -l%i -d%i -b%i -p%i -x%i -c %0.4lf -f ", job->outName, seqres, verbose, renumber, relabel, displayStatus, outputBent, partial, displaySummary, rmsdCutoff);
		for (i=0; i<sizeof(formats)/sizeof(formats[0]); i++) {
			if (outFormats & (1<<i)) {
				if (!first) {
					end[0] = ',';
					end++;
				}
				end += sprintf(end, "%s", formats[i]);
				first = 0;
			}
		}
		for (i=0; i<job->numFiles; i++) {
			end += sprintf(end, " %s", job->files[i]);
		}
		RootStatusLog(consoleOutput, log, commandLine);
		free(commandLine);
	}
	RootStatusLog(consoleOutput, log, "Loading files...");
	for (i = 0; i<job->numFiles; i++) {
		char *fileName = job->files[i];
		char *pos = strchr(fileName, ':');
		char *name;
		PDBData *pdb;

		/* Under Windows, check for absolute paths. */
#if defined(WIN32) || defined(__WIN32__)
		if (pos && pos - fileName == 1) {
			pos = strchr(pos + 1, ':');
		}
#endif
		if (pos) {
			*pos = 0;
			if (!pos[1]) pos = 0;
		}
		name = strrchr(fileName, DIR_DELIMITER);
		if (!name) name = fileName;
		pdb = LoadPDBFile(fileName, seqres, 2, 0, 1);
		if (!pdb || !pdb->numModels) {
			problems++;
			WarningLog(consoleOutput, log, "Error reading %s, skipping.", fileName);
			continue;
		}
		if (pos) {
			pos++;
			while(pos && pos[0]) {
				PDBModel *model = 0;
				char *nextColon = strchr(pos, ':');
				char *comma = strchr(pos, ',');
				char *next = 0;
				int badFormat = 0;
				int m = -1;
				if (comma) {
					if (!nextColon || comma < nextColon) {
						next = comma+1;
						nextColon = 0;
					}
				}
				if (nextColon) {
					char *end;
					nextColon++;
					if (nextColon[0] == 0 || nextColon[0] == ',') {
						next = nextColon+1;
					}
					else {
						m = strtol(nextColon, &end, 10);
						if (end == nextColon || m < 0) {
							WarningLog(consoleOutput, log, "Invalid model specified %s:%s.", fileName, pos);
							break;
						}
						next = end;
						if (next[0] == ',') next++;
						model = FindModel(pdb, m);
						if (!model) {
							WarningLog(consoleOutput, log, "Model %i doesn't exist in %s.", m, fileName);
							pos = next;
							continue;
						}
					}
				}
				if (!model) {
					model = pdb->models[0];
				}
				if (pos[0] == ':' || pos[0] == ',' || pos[0] == 0) {
					int i;
					chains = (PDBChain**) realloc(chains, sizeof(PDBChain*) * (numChains+model->numChains));
					for (i=0; i<model->numChains;i++) {
						PDBChain *c;
						if (model->numChains > 1 && model->chains[i]->chainName == ' ') continue;
						chains[numChains++] = c = DuplicateChain(model->chains[i]);
						if (m != -1) {
							char *temp = malloc(strlen(c->idString) + 50);
							StatusLog(consoleOutput, log, "Chain %s loaded, model %i (%i residues).", c->idString, model->number, CountAlphaCarbons(chains[numChains-1]));

							strcpy(temp, c->idString);
							free(c->idString);
							c->idString = temp;
							temp = strrchr(temp, ':');
							if (!temp || !temp[1] || temp[2]) {
								temp = strchr(temp, 0);
								temp[0] = ':';
								temp++;
							}
							else {
								temp+=2;
							}
							sprintf(temp, ":%i", model->number);
						}
						else
							StatusLog(consoleOutput, log, "Chain %s loaded (%i residues).", model->chains[i]->idString, CountAlphaCarbons(chains[numChains-1]));
						if (renumber == 1) Renumber(chains[numChains-1]);
					}
				}
				else {
					while (pos[0] && pos[0] != ':' && pos[0] != ',') {
						PDBChain *c = GetChain(model, pos[0]);
						char *end = pos+1;
						int p1=-1, p2;
						chains = (PDBChain**) realloc(chains, sizeof(PDBChain*) * (numChains+1));
						if (end[0] == '(') {
							int r1, r2;
							char i1[2]="\0", i2[2]="\0";
							char *middle = end+1;
							char *temp = middle;
							int happy = 0;
							end = strchr(middle, ')');

							r1 = strtol(temp, &middle, 10);
							if (!end) {
								WarningLog(consoleOutput, log, "Invalid range in %s:%s.", fileName, pos);
								break;
							}
							end++;
							while (1) {
								if (temp == middle || middle > end) {
									break;
								}
								while (middle[0] == ' ' || middle[0] == '\t') middle++;
								if (middle[0] != '-' && middle[0] != ')') {
									i1[0] = middle[0];
									middle++;
									while (middle[0] == ' ' || middle[0] == '\t') middle++;
								}

								if (middle[0] != '-' || middle > end) {
									break;
								}

								temp = middle + 1;
								r2 = strtol(temp, &middle, 10);
								if (temp == middle || middle > end) {
									break;
								}
								while (middle[0] == ' ' || middle[0] == '\t') middle++;
								if (middle[0] != ')') {
									i2[0] = middle[0];
									middle++;
									while (middle[0] == ' ' || middle[0] == '\t') middle++;
								}
								if (middle+1 != end) {
									break;
								}
								happy = 1;
								break;
							}
							if (!happy) {
								WarningLog(consoleOutput, log, "Invalid range in %s:%s.", fileName, pos);
								pos = end;
								continue;
							}
							if (c) {
								char *suffix;
								int atom = 0;
								int res = 0;
								int atom2;
								c = DuplicateChain(c);
								if (renumber == 1) {
									Renumber(c);
								}
								if (!GetResidueFirstIndex(c, &p1, r1, i1[0]) || !GetResidueFirstIndex(c, &p2, r2, i2[0]) || r2 < r1) {
									WarningLog(consoleOutput, log, "Invalid range in %s (%i%s-%i%s).", c->idString, r1,i1,r2,i2);
									pos = end;
									CleanupChain(c);
									continue;
								}
								for (atom2 = 0; atom2 < c->numAtoms; atom2++) {
									if (c->atoms[atom2].resIndex < p1 || c->atoms[atom2].resIndex > p2) continue;
									c->atoms[atom] = c->atoms[atom2];
									c->atoms[atom].res -= p1;
									atom++;
								}
								c->numAtoms = atom;
								atom = c->residues[p1].firstAtom;
								for (res = 0; res <= p2 - p1; res++) {
									c->residues[res] = c->residues[res+p1];
									if (c->residues[res].firstAtom == -1) continue;
									c->residues[res].firstAtom -= atom;
									c->residues[res].lastAtom -= atom;
								}
								c->length = p2-p1+1;
								memcpy(c->seq->seq, c->seq->seq+p1, c->length * sizeof(char));
								c->seq->length = c->length;
								pos = end;
								chains[numChains++] = c;
								suffix = malloc(150 + strlen(c->idString));
								strcpy(suffix, c->idString);
								temp = c->idString;
								c->idString = suffix;
								suffix = strrchr(c->idString, ':');
								if (!suffix || !suffix[1] || suffix[2]) {
									suffix = strchr(c->idString, 0);
									strcpy(suffix, ": ");
								}
								suffix += 2;
								sprintf(suffix, "(%i%s-%i%s)", r1,i1,r2,i2);
								if (m == -1) {
									StatusLog(consoleOutput, log, "Chain %s, residues %i%s-%i%s loaded (%i residues).", temp, r1,i1,r2,i2, CountAlphaCarbons(chains[numChains-1]));
								}
								else {
									StatusLog(consoleOutput, log, "Chain %s, residues %i%s-%i%s, model %i loaded (%i residues).", temp, r1,i1,r2,i2, model->number, CountAlphaCarbons(chains[numChains-1]));
									sprintf(strchr(suffix, 0), ":%i", model->number);
								}
								free(temp);
								continue;
							}
						}
						if (!c) {
							WarningLog(consoleOutput, log, "Couldn't find chain %c of %s.", pos[0], fileName);
							pos = end;
							continue;
						}
						c = DuplicateChain(c);
						if (renumber == 1) {
							Renumber(c);
						}
						chains[numChains++] = c;
						if (m == -1)
							StatusLog(consoleOutput, log, "Chain %s loaded (%i residues).", c->idString, CountAlphaCarbons(chains[numChains-1]));
						else {
							char *temp = malloc(strlen(c->idString) + 50);
							StatusLog(consoleOutput, log, "Chain %s loaded,  model %i (%i residues).", c->idString, model->number, CountAlphaCarbons(chains[numChains-1]));

							strcpy(temp, c->idString);
							free(c->idString);
							c->idString = temp;
							temp = strrchr(temp, ':');
							if (!temp || !temp[1] || temp[2]) {
								temp = strchr(temp, 0);
								temp[0] = ':';
								temp++;
							}
							else {
								temp+=2;
							}
							sprintf(temp, ":%i", model->number);
						}
						pos = end;
					}
				}
				pos = next;
			}
		}
		else {
			int failed = 1;
			chains = (PDBChain**) realloc(chains, sizeof(PDBChain*) * (numChains+pdb->models[0]->numChains));
			for (j=0; j<pdb->models[0]->numChains; j++) {
				if (pdb->models[0]->chains[j]->chainName && pdb->models[0]->chains[j]->chainName != ' ') {
					failed = 0;
					chains[numChains++] = PopChain(pdb->models[0], pdb->models[0]->chains[j]->chainName);
					if (renumber == 1) Renumber(chains[numChains-1]);
					if (verbose == 2) {
						StatusLog(consoleOutput, log, "Chain %s loaded (%i residues).", chains[numChains-1]->idString, CountAlphaCarbons(chains[numChains-1]));
					}
					j--;
				}
			}
			if (failed && pdb->models[0]->numChains == 1) {
				failed = 0;
				chains[numChains++] = PopChain(pdb->models[0], 0);
				if (renumber == 1) Renumber(chains[numChains-1]);
				if (verbose == 2) {
					StatusLog(consoleOutput, log, "Chain %s loaded (%i residues).", chains[numChains-1]->idString, CountAlphaCarbons(chains[numChains-1]));
				}
			}
			/* Should never happen. */
			if (failed) {
				problems++;
				WarningLog(consoleOutput, log, "Couldn't find chains in %s.", fileName);
			}
		}
		CleanupPDB(pdb);
	}
	for (i=0; i<numChains; i++) {
		int c = CountAlphaCarbons(chains[i]);
		if (c < 5) {
			break;
		}
	}
	if (i < numChains) {
		ErrorLog(consoleOutput, log, "Chain %s has less than 5 residues with alpha carbon coordinates.", chains[i]->idString);
		ErrorLog(consoleOutput, log, "Unable to run algorithm.");
	}
	else if (numChains < 2) {
		if (numChains == 0) {
			problems++;
			ErrorLog(consoleOutput, log, "No chains loaded, nothing to align.");
		}
		else {
			problems++;
			ErrorLog(consoleOutput, log, "Only one chain loaded, nothing to align.");
		}
	}
	else {
		MultipleAlignment *ma;
		RootStatusLog(consoleOutput, log, "Chains loaded, running alignment...\n");
		if (relabel) {
			RelabelChains(chains, numChains);
		}

		ma = Align(chains, numChains, displayStatus && !parallel, singleThreaded);
		if (ma->numResidues == 0) {
			problems++;
			noAlignment = 1;
			WarningLog(consoleOutput, log, "Unable to align structures, no output files generated.");
		}
		else {
			int *originalOrder = (int*) malloc(sizeof(int) * numChains * 2);
			int *newOrder = originalOrder + numChains;
			int ref, coreSize;

			/* Current order indicates reference structure. */
			for (i=0; i<ma->numChains; i++) {
				originalOrder[i] = i;
				newOrder[i] = ma->chains[i]->id;
			}
			ref = ma->chains[0]->id;

			if (outputBent) {
				/* Note - OutputAlignmentBentPDB aligns things relative to the first structure,
				 * and outputs them according to the ids in the passed in order.
				 * It needs to start with a bent alignment.
				 */
				sprintf(pdbFile, "%s_bent.pdb", job->outName);
				if ((outFormats & OUTPUT_PDB) && !OutputAlignmentBentPDB(ma, pdbFile, originalOrder)) {
					problems++;
					WarningLog(consoleOutput, log, "Can't open %s.", pdbFile);
				}
				CalcAlignmentRMSD(ma);

				ReorderAlignment(ma, originalOrder);

				sprintf(otherFile, "%s_bent", job->outName);
				problems += OutputAllButPdb(ma, otherFile, pdbFile, outFormats, consoleOutput, log, ref);

				ReorderAlignment(ma, newOrder);
			}

			if (ma->numChains == 2) {
				summary->structalScores[0] = CalcStructalScore(ma);
			}
			UnbentRMSDAlign(ma);
			if (ma->numChains == 2) {
				summary->structalScores[1] = CalcStructalScore(ma);
			}
			CalcAlignmentPvalue(ma);

			coreSize = ma->numResidues;
			while (1) {
				FinalExtend(ma, rmsdCutoff);
				if (coreSize == ma->numResidues) break;
				UnbentRMSDAlign(ma);
				coreSize = ma->numResidues;
				summary->coreBentRMSD = ma->rmsd;
			}
			summary->coreSize = coreSize;

			CalcAlignmentRMSD(ma);
			summary->coreRMSD = ma->rmsd;

			/* Makes partial algorithm a bit simpler, as it makes chain ids match chain order.
			 */
			ReorderAlignment(ma, originalOrder);

			summary->pvalue = -1;
			if (ma->numChains == 2) {
				summary->lengths[0] = ma->chains[0]->length;
				summary->lengths[1] = ma->chains[1]->length;
				summary->pvalue = ma->pvalue;
				summary->structalScores[2] = CalcStructalScore(ma);
			}
			else if (ma->numChains > 2) {
				if (partial) CalcPartial(ma, rmsdCutoff);
			}


			sprintf(pdbFile, "%s.pdb", job->outName);
			if ((outFormats & OUTPUT_PDB) && !OutputAlignmentPDB(ma, pdbFile, originalOrder)) {
				problems++;
				WarningLog(consoleOutput, log, "Can't open %s.", pdbFile);
			}

			strcpy(otherFile, job->outName);
			problems += OutputAllButPdb(ma, otherFile, pdbFile, outFormats, consoleOutput, log, ref);

			free(originalOrder);
		}
		for (i=0; i<numChains; i++) {
			free(ma->chains[i]->res);
			free(ma->chains[i]);
		}
		CleanupAlignment(ma);
	}

	for (i=0; i<numChains; i++) {
		CleanupChain(chains[i]);
	}
	free(chains);
	free(pdbFile);
	free(otherFile);
	if (log) fclose(log);
	if (noAlignment && problems == 1) problems = -1;

	return problems;
}

/* Get rid of a warning with MSVC and fastcall. */
int CDECL main(int realArgc, char **realArgv) {
	// int timeTest = -(int)timeGetTime();
	JobInfo *jobs = (JobInfo*) calloc(1, sizeof(JobInfo));
	int job;
	int numJobs = 1;
	int verbose = -1;
	int displayStatus = -1;
	int run = 0;
	int outFormats = -1;

	/* Used to figure out if there were any non-command line warnings on multiple jobs. */
	int commandWarnings = 0;

	/* Values indicate to use default. */
	double rmsdCutoff = -1;
	int relabel = -1;
	int renumber = -1;
	int seqres = -1;
	int outputBent = -1;
	int partial = -1;
	int displaySummary = -1;
	unsigned long numThreads = 0;
	/* Currently set according to how many output prefixes are given */
	int log = -1;
	int consoleOutput = -1;

	int i, term;

	if (realArgc == 1) {
		free(jobs);
		Usage();
		exit(1);
	}
	for (term = 0; term < 3; term++) {
		int argc = realArgc;
		char **argv = realArgv;
		int pos = 1;
		if (term > 0) {
			/* Perhaps a slightly odd way of handling defaults, it works and it's simple. */
			char *args2 = "-b0v2dr0slpcx 5.0 -f pdb,txt,fasta,spt";
			char *args;
			int maxPointers;
			int len;
			if (term == 1)
				args2 = getenv("MATT_PARAMS");
			if (!args2)
				continue;
			pos = 0;
			len = strlen(args2);
			maxPointers = 2+len/2;
			argv = (char**) malloc(sizeof(char) * (len+1) + sizeof(char*) * maxPointers);
			args = (char*) (argv+maxPointers);
			strcpy(args, args2);
			args = strtok(args, " \t");
			argc = 0;
			do {
				argv[argc++] = args;
			}
			while(args = strtok(0, " \t"));
		}

		while (pos < argc) {
			int nextpos = pos+1;
			int index = 1;
			if (argv[pos][0] != '-') {
				if (term) {
					ErrorLog(1, 0, "Can't specify pdb file in MATT_PARAMS (%s).", argv[pos]);
				}
				else {
					jobs[numJobs-1].numFiles++;
					jobs[numJobs-1].files = (char**) realloc(jobs[numJobs-1].files, sizeof(char*) * jobs[numJobs-1].numFiles);
					jobs[numJobs-1].files[jobs[numJobs-1].numFiles-1] = strdup(argv[pos]);
				}
				pos++;
				continue;
			}
			while (argv[pos][index]) {
				if (argv[pos][index] ==  'r' ||
					argv[pos][index] ==  's' ||
					argv[pos][index] ==  'd' ||
					argv[pos][index] ==  'v' ||
					argv[pos][index] ==  'l' ||
					argv[pos][index] ==  'b' ||
					argv[pos][index] ==  'x' ||
					argv[pos][index] ==  'p') {

						char arg = argv[pos][index];
						char val = argv[pos][index+1];
						if (val >= '0' && val <= '9') {
							val -= '0';
							index++;
						}
						else
							val = 1;
						if (arg ==  'v') {
							/* May add 2 back at some point. */
							if (val > 2) {
								ErrorLog(1, 0, "Invalid parameter value: -%c%i.", arg, val);
							}
							else if (verbose != -1) {
								if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
							}
							else verbose = val;
						}
						else {
							if (val > 1) {
								ErrorLog(1, 0, "Invalid parameter value: -%c%i.", arg, val);
							}
							else if (arg ==  'r') {
								if (renumber != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else renumber = val;
							}
							else if (arg ==  'p') {
								if (partial != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else partial = val;
							}
							else if (arg ==  'd') {
								if (displayStatus != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else displayStatus = val;
							}
							else if (arg ==  'x') {
								if (displaySummary != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else displaySummary = val;
							}
							else if (arg ==  's') {
								if (seqres != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else seqres = val;
							}
							else if (arg ==  'l') {
								if (relabel != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else relabel = val;
							}
							else if (arg ==  'b') {
								if (outputBent != -1) {
									if (!term) ErrorLog(1, 0, "Duplicate option: -%c.", arg);
								}
								else outputBent = val;
							}
						}
				}
				else if (argv[pos][index] == 'o' ||
						 argv[pos][index] == 'L' ||
						 argv[pos][index] == 'c' ||
						 argv[pos][index] == 't' ||
						 argv[pos][index] == 'f') {

					if (nextpos == argc) {
						ErrorLog(1, 0, "Argument for -%c missing.", argv[pos][index]);
					}
					else {
						if (argv[pos][index] == 'o') {
							if (jobs[0].outName) {
								if (!term) {
									numJobs++;
									jobs = (JobInfo*) realloc(jobs, sizeof(JobInfo) * numJobs);
									memset(jobs + (numJobs-1), 0, sizeof(JobInfo));
									jobs[numJobs-1].outName = strdup(argv[nextpos]);
								}
							}
							else {
								jobs[numJobs-1].outName = strdup(argv[nextpos]);
							}
						}
						else if (argv[pos][index] ==  'L') {
							FileReader *FileReader;
							char *line = 0;
							int len = 0;
							int count = 0;
							if (term) {
								ErrorLog(1, 0, "List files not allowed in MATT_PARAMS.");
							}
							else {
								FileReader = readerOpen(argv[nextpos]);
								if (!FileReader) {
									ErrorLog(1, 0, "Couldn't open list %s for reading.", argv[nextpos]);
								}
								else {
									while (readerGetLine(FileReader, &line, &len)) {
										int len;
										char *line2 = line;
										while (*line2 == ' ' || *line2 == '\t') line2++;
										len = (int)strlen(line);
										while (len && ((line[len-1] == ' ' && (len==1 || (line[len-2] != ',' && line[len-2] != ':'))) || line[len-1] == '\t')) len--;
										line[len] = 0;
										if (!len) continue;
										if (line2[0] == '-' && line2[1] == 'o' && (line2[2] == '\t' || line2[2] == ' ')) {
											line2 += 2;
											while (line2[0] == ' ' || line2[0] == '\t') line2++;
											if (!line2[0]) {
												ErrorLog(1, 0, "Argument for -o missing in list %s.", argv[nextpos]);
											}
											else {
												if (jobs[0].outName) {
													numJobs++;
													jobs = (JobInfo*) realloc(jobs, sizeof(JobInfo) * numJobs);
													memset(jobs + (numJobs-1), 0, sizeof(JobInfo));
												}
												jobs[numJobs-1].outName = strdup(line2);
											}
											count++;
										}
										else {
											jobs[numJobs-1].numFiles++;
											jobs[numJobs-1].files = (char**) realloc(jobs[numJobs-1].files, sizeof(char*) * jobs[numJobs-1].numFiles);
											jobs[numJobs-1].files[jobs[numJobs-1].numFiles-1] = strdup(line2);
											count ++;
										}
									}
									if (!count) {
										WarningLog(1, 0, "List file %s is empty.", argv[nextpos]);
									}
									free(line);
									readerClose(FileReader);
								}
							}
						}
						else if (argv[pos][index] ==  'c') {
							char *c;
							if (rmsdCutoff >= 0) {
								if (!term) ErrorLog(1, 0, "Only one -%c occurance allowed.", argv[pos][index]);
							}
							else {
								rmsdCutoff = strtod(argv[nextpos], &c);
								if (rmsdCutoff < -1 || c[0]) {
									ErrorLog(1, 0, "Invalid cutoff value.", argv[pos][index]);
								}
							}
						}
						else if (argv[pos][index] ==  't') {
							char *c;
							if (numThreads > 0) {
								if (!term) WarningLog(1, 0, "Only one -%c occurance allowed. Second occurance ignored.");
							}
							else {
								numThreads = strtoul(argv[nextpos], &c, 0);
								if (numThreads < 1 || numThreads > 256 || c[0]) {
									WarningLog(1, 0, "Thread count must be an integer between 1 and 256. Parameter ignored.");
								}
								else {
#ifdef _OPENMP
									omp_set_num_threads(numThreads);
#else
									if (!term) WarningLog(1, 0, "Not compiled with OpenMP support. -t ignored");
#endif
								}
							}
						}
						else if (argv[pos][index] ==  'f') {
							if (outFormats != -1) {
								if (!term) ErrorLog(1, 0, "Duplicate switch: -%c.", argv[pos][index]);
							}
							else {
								char *temp = strdup(argv[nextpos]);
								char *format = strtok(temp, ",");
								outFormats = 0;
								do {
									if (!strcmp("all", format)) {
										outFormats = 0xFFFF;
									}
									for (i=sizeof(formats)/sizeof(formats[0]); i>=0; i--) {
										if (!strcmp(formats[i], format)) break;
									}
									if (i < 0) {
										ErrorLog(1, 0, "Unsupported output format: %s.", format);
									}
									else if (outFormats & (1<<i)) {
										WarningLog(1, 0, "Duplicate output format ignored: %s.", format);
									}
									else {
										outFormats |= 1<<i;
									}
								}
								while (format = strtok(0, ","));
								free(temp);
							}
						}
						nextpos++;
					}
				}
				else {
					ErrorLog(1, 0, "Invalid option: %c.", argv[pos][index]);
				}
				index++;
			}
			pos = nextpos;
		}
		if (term > 0) {
			free(argv);
		}
		if (errors) break;
	}

	if (!errors) {
		InitMultipleAlignmentTables();
		if (verbose == 0) {
			log = consoleOutput = 0;
		}
		if (log == -1) log = (numJobs>1);
		if (consoleOutput == -1) consoleOutput = (numJobs<=1);
		if (!jobs[0].outName) {
			jobs[0].outName = strdup("MattAlignment");
			WarningLog(1, 0, "No output name specified. Using %s.", jobs[0].outName);
		}

		for (i=0; i<numJobs; i++) {
			if (jobs[i].numFiles) break;
		}

		if (i == numJobs) {
			ErrorLog(1, 0, "Nothing to do.");
		}
		else {
			commandWarnings = warnings;
			run = 1;
			if (numJobs == 1) {
				Summary summary;
				int problems = RunJob(jobs, seqres, verbose, renumber, relabel, displayStatus, outputBent, rmsdCutoff, consoleOutput, log, 0, outFormats, partial, &summary, displaySummary, 0);
				if (!problems && displaySummary) {
					DisplaySummary(jobs, &summary);
				}
			}
			else {
				char progressTemplate[64];
				int len = 1;
				int d = 10;
				int count = 0;
				int backspace = 0;
				while (numJobs >= d) {
					len ++;
					d*=10;
				}
#ifndef _OPENMP	
				{
					if (displayStatus) {
						sprintf(progressTemplate, "Single threaded [%%%ii/%i alignments]", len, numJobs);
						DisplayProgress(&backspace, progressTemplate, count);
					}
#else
#pragma omp parallel default(shared)
				{
#pragma omp single
					if (displayStatus) {
						int threads = omp_get_num_threads();
						if (threads > 1) sprintf(progressTemplate, "%i threads [%%%ii/%i alignments]", threads, len, numJobs);
						else sprintf(progressTemplate, "1 thread [%%%ii/%i alignments]", len, numJobs);
						DisplayProgress(&backspace, progressTemplate, count);
					}

#pragma omp for schedule(dynamic, 1)
#endif
					for (job=0; job<numJobs; job++) {
						Summary summary;
						int problems = RunJob(jobs + job, seqres, verbose, renumber, relabel, displayStatus, outputBent, rmsdCutoff, consoleOutput, log, 1, outFormats, partial, &summary, displaySummary, 1);

						if (problems | displaySummary | displayStatus) {
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (problems) {
									CleanupAndBackspace(backspace);
									backspace = 0;
									if (problems == -1) {
										printf("%s: Unable to align structures.\n", jobs[job].outName);
									}
									else {
										printf("%s: Problem with alignment.\n", jobs[job].outName);
									}
								}

								else if (displaySummary) {
									CleanupAndBackspace(backspace);
									backspace = 0;
									DisplaySummary(jobs+job, &summary);
								}

								if (displayStatus) {
									DisplayProgress(&backspace, progressTemplate, ++count);
								}
							}
						}
					}
				}
				if (displayStatus) {
					fprintf(stderr, "\n");
				}
			}
		}
	}
	if (run && (errors || commandWarnings != warnings) && numJobs > 1) {
		RootStatusLog(1, 0, "Done, but some problems logged.");
	}
	else if (errors) {
		RootStatusLog(1, 0, "No alignments run.");
		if (!run) {
			if (term == 1) {
				RootStatusLog(1, 0, "Check MATT_PARAMS.");
			}
			else Usage();
		}
	}
	else {
		RootStatusLog(1, 0, "Done.\n");
	}

	for (job=0; job<numJobs; job++) {
		for (i=0; i<jobs[job].numFiles; i++)
			free(jobs[job].files[i]);
		free(jobs[job].files);
		free(jobs[job].outName);
	}
	free(jobs);
	//timeTest += timeGetTime();
	//printf("Time: %0.3lf s\n", timeTest/1000.0);
	return 0;
}
