#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

void byDivisionSequentially() {
	const int M = 100; // zakres od
	const int N = 997; // zakres do

	bool isPrimeNumber[N + 1];  //tu nie powinno być N+1-M?
								//chociaż w sumie N+1 jest wygodne i można by w przypadku byDeletion sobie po prostu stąd czytać te liczby pierwsze zamiast iść po kolei
								//i jest wygodniejsze w myśleniu, ale w sumie jakbyś dał przedział 1000-1100 to 1000 booli w dupę xD

	for (int i = M; i <= N; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}

	for (int currentNumber = M; currentNumber <= N; currentNumber++) {

		int maxDivisorValue = (int)floor(sqrt(currentNumber)); // maksymalna wartość podzielnika dla danej liczby

		//tu teoretycznie można by pokombinować z optymalizacją, że te które uznaliśmy za pierwsze zapisywać sobie do listy liczb pierwszych 
		//i wtedy w następnej iteracji jak currentDivisor dojedzie do pierwszej z tych liczb, to zamiast sprawdzać po kolei
		//można by brać liczby z tej listy, gdy ta lista się skończy można nadać currentDivisor wartość = ostatni element z listy + 1 i kontynować dalej,
		//z tym że zysk jakikolwiek pojawi się dopiero, gdy będziemy sprawdzać dzielniki większe od M (program będzie lepiej działał wtedy dla M=2, gorzej dla M bliższego N,
		//jeśli w ogóle będzie jakikolwiek zysk bo nie wiem czy operacja zapisu i odczytu z listy nie będzie bardziej zabierało więcej czasu niż po prostu dodawanie,
		//ale w sumie sporo liczb odpadnie i będzie mniej dzielenia

		//z treści zadania mi się wydaje, że takie rzeczy właśnie mamy też potestować ale nie jestem pewien

		for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
			//sprawdzaj kolejne podzielniki
			if (currentNumber % currentDivisor == 0) {
				//jeżeli dzieli się bez reszty to liczba jest złożona
				isPrimeNumber[currentNumber] = false;   //gdyby isPrimeNumber miało rozmiar N+1-M musiałoby być [currentNumber-M]
				break;
			}
		}
	}

	int numberOfPrimes = 0;
	for (int i = M; i <= N; i++) {
		// zlicz liczby pierwsze
		if (isPrimeNumber[i]) {
			numberOfPrimes++;
		}
	}
	printf("Dzielenie sekwencyjnie:\n");
	printf("W przedziale <%d, %d> jest %d liczb pierwszych\n", M, N, numberOfPrimes);
}

void byDeletionSequentially() {
	const int M = 100; // zakres od
	const int N = 997; // zakres do

	bool isPrimeNumber[N + 1];

	for (int i = M; i <= N; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}

	int maxDivisorValue = (int)floor(sqrt(N)); // maksymalna wartość podzielnika
	int currrentPrimeNumber = 2; // pierwsza liczba pierszwa
	int currentMultiplier = 2; // mnożnik do obliczania wielokrotności

	while (currrentPrimeNumber <= maxDivisorValue) {
		// dla każdej liczby pierwszej z zakresu <2, sqrt(N)>
		int currentMultiple = currrentPrimeNumber * currentMultiplier; // obecnie przetwarzana wielokrotność
		if (currentMultiple <= N) {
			isPrimeNumber[currentMultiple] = false; // wielokrotność ustawiamy na liczbę złożoną
			currentMultiplier++; // i zwiększamy mnożnik 
		}
		else {
			//jeśli wyjdziemy poza przedział z wielokrotnościami (currentMultiple > N) to przechodzimy do kolejnej liczby
			currrentPrimeNumber++;  //i tu tak samo - można by pomijać te liczby które już zostały usunięte
			currentMultiplier = 2;

			//currentMultiplier zamiast ustawiać na 2 ustawić na wartość ceil(M/currrentPrimeNumber)
			//przykład
			//M=100
			//currrentPrimeNumber=3
			//ceil(M/currentPrimeNumber) = 34
			// 3*34 = 102
			// 3*35 = 105
			//itd.
		}
	}

	int numberOfPrimes = 0;
	for (int i = M; i <= N; i++) {
		// zlicz liczby pierwsze
		if (isPrimeNumber[i]) {
			numberOfPrimes++;
		}
	}
	printf("Usuwanie sekwencyjnie:\n");
	printf("W przedziale <%d, %d> jest %d liczb pierwszych\n", M, N, numberOfPrimes);

}







//void byDivisionParallel(long M, long N, long maxDivisorValue, bool isPrimeNumber[]) {
//	omp_set_num_threads(4);
//#pragma omp parallel
//	{
//#pragma omp for nowait
//		for (int currentNumber = M; currentNumber <= N; currentNumber++) {
//			int maxDivisorValue = (int)floor(sqrt(currentNumber)); // maksymalna wartość podzielnika dla danej liczby
//
//			for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
//				//sprawdzaj kolejne podzielniki
//				if (currentNumber % currentDivisor == 0) {
//					//jeżeli dzieli się bez reszty to liczba jest złożona
//					isPrimeNumber[currentNumber - M] = false;
//					break;
//				}
//			}
//		}
//	}
//}





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
		for (long i = 2; i <= maxDivisorValue; i++) {
#pragma omp for nowait
			for (long j = (long)max(2.0, ceil(M / i)); j <= (long)ceil(N / i); j++) {			
				long multiple = i * j;														
				if (multiple >= M && multiple <= N)											
					isPrimeNumber[multiple - M] = false;
			}
		}
	}
}

const long M = 2; // zakres od
const long N = 120000000; // zakres do

bool isPrimeNumber[N + 1 - M];


int main() {

	long maxDivisorValue = (long)sqrt(N); // maksymalna wartość podzielnika
	omp_set_num_threads(4);

	for (long i = 0; i <= N - M; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}

	//byDeletionForLoop(M, N, maxDivisorValue, isPrimeNumber);
	//byDeletionFunctional(M, N, maxDivisorValue, isPrimeNumber);
	byDeletionDomain(M, N, maxDivisorValue, isPrimeNumber);

	return 0;
}