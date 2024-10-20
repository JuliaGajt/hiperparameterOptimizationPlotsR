# Wizualizacja Optymalizacji Hiperparametrów

## Opis projektu

Projekt służy do analizy wyników optymalizacji hiperparametrów modeli uczenia maszynowego oraz metod walidacyjnych. Skrypt wykorzystuje różne techniki statystyczne, takie jak test Scott-Knott ESD i współczynnik zgodności Kendall W, aby porównywać efektywność metod optymalizacji i walidacji. Wyniki są przedstawiane w formie wykresów pudełkowych oraz pivot table, które ułatwiają wizualizację wyników dla różnych kombinacji metod.

## Zawartość Repozytorium

Repozytorium zawiera trzy główne pliki skryptów R, które analizują dane dotyczące optymalizacji hiperparametrów i walidacji modeli:

1. **KendallW.R**: Skrypt ten przeprowadza test Kendall W, mierzący zgodność wyników różnych metod regresji na podstawie błędów MAE. Skrypt generuje również wykres pudełkowy z wynkiami rang uzyskanych w teście SK-ESD zbiorczo dla wszystkich modeli i zbiorów danych. Zapisuje również tabele .csv ze średnimi rangami dla wszystkich kobinacji osobno względem zbioru danych i modelu.
2. **SK-ESD.R**: Skrypt wykonuje analizę i generuje wykresy pudełkowe wyników MAE za pomocą testu Scott-Knott ESD dla wszystkich kombinacji metod walidacji i optymalizacji i grupuje je według ich wydajności (osobno dla każdego zbioru danych i modelu uczenia maszynowego). Generuje również wykresy pudełkowe rang otrzymanych zbiorczo dla metod walidacji i technik optymalizacji. Tworzy również wykres pudełkowy wstępny rozkładów MAE względem metod optymalizacji i walidacji.
3. **SK-ESD-mean.R**: Ten skrypt rozszerza analizę Scott-Knott ESD, generuje rozkład rankingów oraz średnie dla różnych kombinacji metod uzyskanych za pomocą testu KS-ESD dla poszczególnych modeli i osobno dla poszczególnych zbiorów danych.

## Dane wejściowe

Skrypty przetwarzają pliki CSV zawierające wyniki eksperymentów optymalizacji hiperparametrów i walidacji modeli. Pliki CSV są wczytywane z lokalnych folderów, a następnie łączone w jedną tabelę danych. Ścieżki do folderów zawierających pliki CSV muszą być dostosowane w każdym skrypcie.

Uruchomienie analizy:
```r
source("KendallW.R")
source("SK-ESD.R")
source("SK-ESD-mean.R")
```