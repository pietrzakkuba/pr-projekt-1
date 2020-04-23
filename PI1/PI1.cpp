#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>

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

void byDivisionParallel() {
	const int M = 100; // zakres od
	const int N = 997; // zakres do

	bool isPrimeNumber[N + 1 - M];  //tu nie powinno być N+1-M?

	for (int i = 0; i < N - M; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}


	//#pragma omp for private(isPrimeNumber)  //isPrimeNumber może być prywatne, dzięki temu unikniemy unieważnień (chyba xD), a po zakończeniu obliczeń można mnożyć z globalnym te wartości,
											//które są false, to true*false da false i nie stanie się nic złego jak będzie false*false 
											//czyli nawet by nie trzeba tego robić atomowo tylko mnożyć na pałę

											//albo w sumie jednak nie, bo to po prostu te unieważnienia pojawią się przy tym mnożeniu i tak xD
#pragma omp parallel shared(isPrimeNumber)
	{

#pragma omp for nowait
		for (int currentNumber = M; currentNumber <= N; currentNumber++) {
			int maxDivisorValue = (int)floor(sqrt(currentNumber)); // maksymalna wartość podzielnika dla danej liczby

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

void byDeletionForLoop() {
	const int M = 100; // zakres od
	const int N = 120; // zakres do

	bool isPrimeNumber[N + 1 - M];

	for (int i = 0; i <= N - M; i++) {
		isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
	}

	int maxDivisorValue = (int)floor(sqrt(N)); // maksymalna wartość podzielnika

	//omp for domenowy
	for (int i = 2; i < maxDivisorValue; i++) {
		//omp for funkcyjny
		for (int j = ceil(M/i); j <= ceil(N/i); j++) {					//i*M/i daje nam M, i*N/i daje nam N, a więc mamy pokryty pełen zakres od M do N
			int multiple = i * j;										//ceil, bo dla M=100 daje i i=3 daje 33,(3), 3*33,(3) daje 99,(9), czyli liczbę poza zakresem
			if (multiple >= M && multiple <= N)							//ceil da nam 34*3 czyli 102, więc zaczynamy od pierwszej wielokrotności w zakresie
				isPrimeNumber[multiple - M] = false;
			cout << i << " * " << j << " = " << multiple << endl;
		}
	}

	int numberOfPrimes = 0;
	for (int i = 0; i < N+1-M; i++) {
		// zlicz liczby pierwsze
		if (isPrimeNumber[i]) {
			numberOfPrimes++;
			cout << i + M << endl;
		}
	}
	printf("Usuwanie sekwencyjnie:\n");
	printf("W przedziale <%d, %d> jest %d liczb pierwszych\n", M, N, numberOfPrimes);

}


int main() {
	//byDivisionSequentially();
	//byDeletionSequentially();
	byDeletionForLoop();

	return 0;
}