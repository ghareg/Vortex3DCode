#include <cmath>
#include "matrix.h"
#include "constants.h"
#include "operator.h"
#include "diag.h"
#include <thread>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <limits>

struct eigen {
	MatType value;
	const MatType *vector;
	
	bool operator<(eigen const &other) const {
		return value.real() < other.value.real();
	}
};

void GenerateMuMat(double*& muMat, Count& muSize);
void sortEValues(eigen* eval, const MatType* evalues, const MatType* evecs, Count bSize);
void ProbDensity(int st, double mu, double kz, const eigen* eval, Count bSize, const State* basis);
void free(SparseMat& baseHam, MatType* evalues, MatType* evecs);

int main(void)
{
	State* basis = nullptr;
	Count bSize;
	GenerateBasis(basis, bSize);
	
	double* muMat = nullptr;
	Count muSize;
	GenerateMuMat(muMat, muSize);

	SparseMat baseHam;
	bool zeroBased = false;
	MatType* evalues = new MatType[neigs];
	MatType* evecs = new MatType[neigs * bSize];
	eigen* eval = new eigen[neigs];
	GenMatProd op;
	bool matrixCons = false;
	int Ni = 0;

	Param pm;
	double kz = 0.0;

	FILE* evalFile;
	char str[100];
	sprintf (str, "Endepmu3DNx%dNy%dNz%d%sTIh%.2fxi%dlmb%dOTE%.2fFBESE%.2fPBCZ.dat", NPx, NPy, NPz, STI ? "S" : "W", hm, (int) xi, (int) lmb, DeltaOTE, Delta);


	evalFile = fopen(str, "w");
	muSize = Ni + 1;
	for (int muCount = Ni; muCount < muSize; ++muCount) {
		fprintf(evalFile, "%.6f", muMat[muCount]);
		pm.kz = kz;
		pm.mu =  muMat[muCount];
		pm.dmu = 0.0;
		pm.hm = hm;
		pm.DeltaS = Delta;
		pm.DeltaP = DeltaOTE;
		if (!matrixCons) {
			calcFullMat(baseHam, basis, bSize, pm, zeroBased);
			
			op.init(baseHam);
			matrixCons = true;
		}
		else {
			pm.dmu = muMat[muCount] - muMat[muCount - 1];
			changeFullUMat(baseHam, pm.dmu, zeroBased);
			fflush(evalFile);
			//op.restart(pm.dmu);
			op.restart(baseHam);
		}
		calcEValues(baseHam, op, evalues, evecs);
		sortEValues(eval, evalues, evecs, bSize);
/*
		ProbDensity(20, pm.mu, kz, eval, bSize, basis);	
		ProbDensity(21, pm.mu, kz, eval, bSize, basis);
		ProbDensity(22, pm.mu, kz, eval, bSize, basis);
		ProbDensity(23, pm.mu, kz, eval, bSize, basis);
		ProbDensity(24, pm.mu, kz, eval, bSize, basis);	
		ProbDensity(25, pm.mu, kz, eval, bSize, basis);
		ProbDensity(26, pm.mu, kz, eval, bSize, basis);
		ProbDensity(27, pm.mu, kz, eval, bSize, basis);
		ProbDensity(28, pm.mu, kz, eval, bSize, basis);
		ProbDensity(29, pm.mu, kz, eval, bSize, basis);
		ProbDensity(30, pm.mu, kz, eval, bSize, basis);
*/

		for (int i = 0; i < neigs; ++i) {
			fprintf(evalFile, "\t%.6f", eval[i].value.real());
		}
		fprintf(evalFile, "\n");	
	
		fflush(evalFile);
	}
	fclose(evalFile);

	free(baseHam, evalues, evecs);
	delete[] basis;
	delete[] muMat;
	delete[] eval;
}

void sortEValues(eigen* eval, const MatType* evalues, const MatType* evecs, Count bSize)
{
	for (int i = 0; i < neigs; ++i) {
		eval[i].value = evalues[i];
		eval[i].vector = &evecs[i * bSize];
	}
	std::sort(&eval[0], &eval[0] + neigs);
}

void GenerateMuMat(double*& muMat, Count& muSize)
{
	muSize = 0;
	int nmax = 600;
	muMat = new double[nmax + 1];
	double mu0 = 0.003;
	double murange = 4.5;
	for (int j = 0; j <= nmax; ++j) {
		muMat[muSize] = mu0 + murange * j / nmax;
		muSize++;
	}
}

void ProbDensity(int st, double mu, double kz, const eigen* eval, Count bSize, const State* basis)
{
	double* dens = new double[NPx * NPy * NPz];
	for (int i = 0; i < NPx * NPy * NPz; ++i) {
		dens[i] = 0.0;
	}
	for (Count i = 0; i < bSize; ++i) {
		dens[basis[i].mpos[0] * NPy * NPz + basis[i].mpos[1] * NPz + basis[i].mpos[2]] +=
			std::norm(eval[st].vector[i]);
	}		   
	char str[100];
//	sprintf (str, "ProbDens2DNx%dNy%dNz%dkz%smu%.2fh%.2fxi%dlmb%dOTE%.2fESE%.2fFBst%d.dat", NPx, NPy, NPz, kz > 0.1 ? "Pi" : "0.0", mu, hm, (int) xi, (int) lmb, DeltaOTE, Delta, st);
	sprintf (str, "ProbDens2DNx%dNy%dNz%dmu%.2fh%.2fxi%dlmb%dlp%dOTE%.2fESE%.2fVp%.2fFBst%d.dat", NPx, NPy, NPz, mu, hm, (int) xi, (int) lmb, (int) lp, DeltaOTE, Delta, Vp, st);
	
	FILE* evalFile = fopen(str, "w");
	for (int ix = 0; ix < NPx; ++ix) {
		for (int iy = 0; iy < NPy; ++iy) {
			for (int iz = 0; iz < NPz; ++iz) {
				fprintf(evalFile, "%d\t%d\t%d\t%.4f\n", ix, iy, iz, dens[ix * NPy * NPz + iy * NPz + iz]);
			}
		}
		fprintf(evalFile, "\n");
	}

	fclose(evalFile);
	delete[] dens;
}

void free(SparseMat& baseHam, MatType* evalues, MatType* evecs)
{
	delete[] baseHam.mat;
	delete[] baseHam.ia;
	delete[] baseHam.ja;

	baseHam.mat = nullptr;
	baseHam.ia = nullptr;
	baseHam.ja = nullptr;

	delete[] evalues;
	delete[] evecs;
}

