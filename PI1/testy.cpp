#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

#define M 2
#define N 1299827
#define THREADS 4

void countPrimes(bool isPrimeNumber[]) {
	int numberOfPrimes = 0;
	for (int i = 0; i <= N - M; i++) {
		// zlicz liczby pierwsze
		if (isPrimeNumber[i]) {
			numberOfPrimes++;
		}
	}
	printf("Usuwanie sekwencyjnie:\n");
	printf("W przedziale <%d, %d> jest %d liczb pierwszych\n", M, N, numberOfPrimes);
}

void byDivisionSeq(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {

	for (int currentNumber = M; currentNumber <= N; currentNumber++) {
		int maxDivisorValue = (int)sqrt(currentNumber); // maksymalna wartoœæ podzielnika dla danej liczby
		for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
			if (isPrimeNumber[currentDivisor - M]) {
				//sprawdzaj kolejne podzielniki
				if (currentNumber % currentDivisor == 0) {
					//je¿eli dzieli siê bez reszty to liczba jest z³o¿ona
					isPrimeNumber[currentNumber - M] = false;
					break;
				}
			}
		}
	}

	if (mode)
		countPrimes(isPrimeNumber);
}

void byDivisionParallel(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	omp_set_num_threads(4);
#pragma omp parallel
	{
#pragma omp for nowait schedule(static,8)
		for (int currentNumber = M; currentNumber <= N; currentNumber++) {
			int maxDivisorValue = (int)sqrt(currentNumber); // maksymalna wartoœæ podzielnika dla danej liczby

			for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
				//sprawdzaj kolejne podzielniki
				if (isPrimeNumber[currentDivisor - M]) {
					if (currentNumber % currentDivisor == 0) {
						//je¿eli dzieli siê bez reszty to liczba jest z³o¿ona
						isPrimeNumber[currentNumber - M] = false;
						break;
					}
				}
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}

void byDeletionSeq(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	for (long i = 2; i <= maxDivisorValue; i++) {
		if (isPrimeNumber[i - M]) {
			for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {
				long multiple = i * j;
				if (multiple >= M && multiple <= N)
					isPrimeNumber[multiple - M] = false;
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}



void byDeletionDomain(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
#pragma omp parallel 
	{
		int threadnum = omp_get_thread_num();
		int div = ((N - M + 1) / omp_get_num_threads());
		int threadM = !threadnum ? M : M + 1 + div * threadnum;
		int threadN = threadnum == omp_get_num_threads() - 1 ? N : M + div * (threadnum + 1);

		long threadMaxDivisorValue = (long)sqrt(threadN);

		for (long i = 2; i <= threadMaxDivisorValue; i++) {
			if (isPrimeNumber[i - M]) {
				for (long j = (long)max(2.0, ceil(threadM / i)); j <= (long)ceil(threadN / i); j++) {
					long multiple = i * j;
					if (multiple >= threadM && multiple <= threadN)
						isPrimeNumber[multiple - M] = false;
				}
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}


void byDeletionFunctional(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
#pragma omp parallel
	{
#pragma omp for nowait schedule(guided)
		for (long i = maxDivisorValue; i >= 2; i--) {
			if (isPrimeNumber[i - M]) {
				for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {
					long multiple = i * j;
					if (multiple >= M && multiple <= N)
						isPrimeNumber[multiple - M] = false;
				}
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}

bool isPrimeNumber[N + 1 - M];

int main() {
	omp_set_num_threads(THREADS);
	long maxDivisorValue = (long)sqrt(N);

	for (long i = 0; i <= N - M; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}


	byDivisionSeq(maxDivisorValue, isPrimeNumber);
	//byDivisionParallel(maxDivisorValue, isPrimeNumber, true);
	//byDeletionSeq(maxDivisorValue, isPrimeNumber, true);
	//byDeletionDomain(maxDivisorValue, isPrimeNumber, true);
	//byDeletionFunctional(maxDivisorValue, isPrimeNumber, true);

	return 0;
}