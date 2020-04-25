#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

#define M 2
#define N 500000000

void countPrimes(bool isPrimeNumber[], bool mode = false) {
	int numberOfPrimes = 0;
	for (int i = 0; i <= N - M; i++) {
		// zlicz liczby pierwsze
		if (isPrimeNumber[i]) {
			numberOfPrimes++;
			if (mode) {
				if (i<100 || i>N - M - 100)
					printf("%d. %d\n", numberOfPrimes, i + M);
			}
		}
	}
	printf("Usuwanie sekwencyjnie:\n");
	printf("W przedziale <%d, %d> jest %d liczb pierwszych\n", M, N, numberOfPrimes);
}

void byDivisionSequentially(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {

	for (int currentNumber = M; currentNumber <= N; currentNumber++) {
		int maxDivisorValue = (int)sqrt(currentNumber); // maksymalna wartość podzielnika dla danej liczby

		for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
			//sprawdzaj kolejne podzielniki
			if (currentNumber % currentDivisor == 0) {
				//jeżeli dzieli się bez reszty to liczba jest złożona
				isPrimeNumber[currentNumber - M] = false;
				break;
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}

void byDeletionSequentially(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	for (long i = 2; i <= maxDivisorValue; i++) {
		for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {
			long multiple = i * j;
			if (multiple >= M && multiple <= N)
				isPrimeNumber[multiple - M] = false;
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}

void byDivisionParallel(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	omp_set_num_threads(4);
#pragma omp parallel
	{
#pragma omp for nowait
		for (int currentNumber = M; currentNumber <= N; currentNumber++) {
			int maxDivisorValue = (int)sqrt(currentNumber); // maksymalna wartość podzielnika dla danej liczby

			for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
				//sprawdzaj kolejne podzielniki
				if (currentNumber % currentDivisor == 0) {
					//jeżeli dzieli się bez reszty to liczba jest złożona
					isPrimeNumber[currentNumber - M] = false;
					break;
				}
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}


void byDeletionDomain(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	omp_set_num_threads(4);
#pragma omp parallel 
	{
		int threadnum = omp_get_thread_num();
		int div = ((N - M + 1) / omp_get_num_threads());
		int threadM = !threadnum ? M : M + 1 + div * threadnum;
		int threadN = threadnum == omp_get_num_threads() - 1 ? N : M + div * (threadnum + 1);

		long threadMaxDivisorValue = (long)sqrt(threadN);

		for (long i = 2; i <= threadMaxDivisorValue; i++) {
			for (long j = (long)max(2.0, ceil(threadM / i)); j <= (long)ceil(threadN / i); j++) {
				long multiple = i * j;
				if (multiple >= threadM && multiple <= threadN)
					isPrimeNumber[multiple - M] = false;
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber);
}

void byDeletionDomainReduction(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	omp_set_num_threads(4);
#pragma omp parallel 
	{
		bool* threadIsPrimeNumber = new bool[N + 1 - M];
		memcpy(threadIsPrimeNumber, isPrimeNumber, N + 1 - M);

		int threadnum = omp_get_thread_num();
		int div = ((N - M + 1) / omp_get_num_threads());
		int threadM = !threadnum ? M : M + 1 + div * threadnum;
		int threadN = threadnum == omp_get_num_threads() - 1 ? N : M + div * (threadnum + 1);

		long threadMaxDivisorValue = (long)sqrt(threadN);

		for (long i = 2; i <= threadMaxDivisorValue; i++) {
			for (long j = (long)max(2.0, ceil(threadM / i)); j <= (long)ceil(threadN / i); j++) {
				long multiple = i * j;
				if (multiple >= threadM && multiple <= threadN)
					threadIsPrimeNumber[multiple - M] = false;
			}
		}


		for (int i = 0; i < N + 1 - M; i++) {
			isPrimeNumber[i] = isPrimeNumber[i] && threadIsPrimeNumber[i];
		}

	}

	if (mode)
		countPrimes(isPrimeNumber, true);
}

void byDeletionFunctional(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	omp_set_num_threads(4);
#pragma omp parallel
	{
#pragma omp for nowait schedule(guided)
		for (long i = maxDivisorValue; i >= 2; i--) {
			for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {
				long multiple = i * j;
				if (multiple >= M && multiple <= N)
					isPrimeNumber[multiple - M] = false;
			}
		}
	}
	if (mode)
		countPrimes(isPrimeNumber, true);
}

void byDeletionFunctionalReduction(long maxDivisorValue, bool isPrimeNumber[], bool mode = false) {
	omp_set_num_threads(4);
#pragma omp parallel
	{
		bool* threadIsPrimeNumber = new bool[N + 1 - M];
		memcpy(threadIsPrimeNumber, isPrimeNumber, N + 1 - M);

		int threadnum = omp_get_thread_num();

#pragma omp for nowait
		for (long i = 2; i <= maxDivisorValue; i++) {
			for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {
				long multiple = i * j;
				if (multiple >= M && multiple <= N)
					threadIsPrimeNumber[multiple - M] = false;
			}
		}


		for (int i = 0; i < N + 1 - M; i++) {
			if (threadIsPrimeNumber[i])
				isPrimeNumber[i] = isPrimeNumber[i] && threadIsPrimeNumber[i];
		}

	}
	if (mode)
		countPrimes(isPrimeNumber);
}


bool isPrimeNumber[N + 1 - M];

int main() {

	long maxDivisorValue = (long)sqrt(N); // maksymalna wartość podzielnika
	omp_set_num_threads(4);

	for (long i = 0; i <= N - M; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}

	//byDeletionDomain(maxDivisorValue, isPrimeNumber);
	//byDeletionFunctional(maxDivisorValue, isPrimeNumber);
	//byDeletionFunctional(maxDivisorValue, isPrimeNumber, true);
	//byDeletionFunctionalReduction(maxDivisorValue, isPrimeNumber, true);
	byDeletionSequentially(maxDivisorValue, isPrimeNumber);
	return 0;
}