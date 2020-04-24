#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

void countPrimes(long M, long N, bool isPrimeNumber[], bool mode=false) {
	int numberOfPrimes = 0;
	for (int i = 0; i <= N - M; i++) {
		// zlicz liczby pierwsze
		if (isPrimeNumber[i]) {
			numberOfPrimes++;
			if (mode) printf("%d\n", i + M);
		}
	}
	printf("Usuwanie sekwencyjnie:\n");
	printf("W przedziale <%d, %d> jest %d liczb pierwszych\n", M, N, numberOfPrimes);
}

void byDivisionSequentially(long M, long N, long maxDivisorValue, bool isPrimeNumber[]) {

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

	//countPrimes()
}

void byDeletionSequentially(long M, long N, long maxDivisorValue, bool isPrimeNumber[]) {
	for (long i = 2; i <= maxDivisorValue; i++) {
		for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {
			long multiple = i * j;
			if (multiple >= M && multiple <= N)
				isPrimeNumber[multiple - M] = false;
		}
	}
	countPrimes(M, N, isPrimeNumber);
}

void byDivisionParallel(long M, long N, long maxDivisorValue, bool isPrimeNumber[]) {
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
}


void byDeletionFunctional(long M, long N, long maxDivisorValue, bool isPrimeNumber[]) {
	omp_set_num_threads(4);
#pragma omp parallel
	{
#pragma omp for nowait
		for (long i = 2; i <= maxDivisorValue; i++) {
			for (long j = (long) max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {	
				long multiple = i * j;														
				if (multiple >= M && multiple <= N)											
					isPrimeNumber[multiple - M] = false;
			}
		}
	}
}

void byDeletionDomain(long M, long N, long maxDivisorValue, bool isPrimeNumber[]) {
	omp_set_num_threads(4);
#pragma omp parallel 
	{
		int threadnum = omp_get_thread_num();
		int div = ((N - M + 1) / omp_get_num_threads());
		int threadM = !threadnum ? M : M + 1 + div * threadnum;
		int threadN = threadnum == omp_get_num_threads() - 1 ? N : M + div * (threadnum + 1);


		long threadMaxDivisorValue = (long)sqrt(threadN);
		//printf("%d zakres: %d - %d, %d\n", omp_get_thread_num(), threadM, threadN, threadMaxDivisorValue);

		for (long i = 2; i <= threadMaxDivisorValue; i++) {
			for (long j = (long)max(2.0, ceil(threadM / i)); j <= (long)ceil(threadN / i); j++) {
				long multiple = i * j;
				if (multiple >= threadM && multiple <= threadN)
					isPrimeNumber[multiple - M] = false;
				printf("%d: %d * %d = %d\n", threadnum, i, j, multiple);
			}
		}
	}
	countPrimes(M, N, isPrimeNumber, true);
}


const long M = 2; // zakres od
const long N = 100; // zakres do

bool isPrimeNumber[N + 1 - M];

int main() {

	long maxDivisorValue = (long)sqrt(N); // maksymalna wartość podzielnika
	omp_set_num_threads(4);

	for (long i = 0; i <= N - M; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}

	//byDeletionForLoop(M, N, maxDivisorValue, isPrimeNumber);
	//byDeletionFunctional(M, N, maxDivisorValue, isPrimeNumber);
	//byDeletionDomain(M, N, maxDivisorValue, isPrimeNumber);
	//testdiv(M, N, maxDivisorValue, isPrimeNumber);

	byDeletionDomain(M, N, maxDivisorValue, isPrimeNumber);

	return 0;
}