#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>

using namespace std;

void byDivisionSequentially() {
    const int M = 100; // zakres od
    const int N = 997; // zakres do

    bool isPrimeNumber[N + 1];

    for (int i = M; i <= N; i++) {
        isPrimeNumber[i] = true; // ustaw wszystkie liczby z zakresy jako pierwsze
    }

    for (int currentNumber = M; currentNumber <= N; currentNumber++) {
        
        int maxDivisorValue = (int)floor(sqrt(currentNumber)); // maksymalna wartość podzielnika dla danej liczby

        for (int currentDivisor = 2; currentDivisor <= maxDivisorValue; currentDivisor++) {
            //sprawdzaj kolejne podzielniki
            if (currentNumber % currentDivisor == 0) {
                //jeżeli dzieli się bez reszty to liczba jest złożona
                isPrimeNumber[currentNumber] = false;
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
            currrentPrimeNumber++;
            currentMultiplier = 2;
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


int main() {
    byDivisionSequentially();
    byDeletionSequentially();
	return 0;
}