``` r
library(ggplot2)
library(dplyr)
library(tidyr)
```

Celem tego raportu jest przedstawienie metod wyceny opcji w modelu CRR w środowisku R. Rozważano europejskie i azjatyckie opcje call i dla każdego z przypadków zaimplementowano odpowiednią funkcję.

Wycena opcji europejskich
-------------------------

W tym przypadku wycena sprowadza się wyłącznie do wyznaczenia wartości oczekiwanej wypłaty z akcji, a następnie zdyskontowania na moment 0. Nie przedstawia to trudności, ponieważ wielkości wzrostu i spadku ceny oraz prawdopodobieństwa tych zdarzeń są określone przez parametry wejściowe, a wynik zależy tylko od ceny instrumentu bazowego w momencie wykupu. Wyceny opcji europejskich w modelu dwumianowym dla podanych przez użytkownika parametrów wejściowych dokonuje poniższa funkcja.

``` r
priceEuropean <- function(S_0, # cena instrumentu bazowego w chwili 0 
                          K, # cena wykonania
                          r, # intensywność oprocentowania
                          sigma, # zmienność
                          t, # czas do wykupu opcji
                          n) { # na ile odcinków jest dzielony przedział [0,t]
    
    delta_n <- t/n
    
    u_n <- exp(sigma * sqrt(delta_n))
    d_n <- exp(- sigma * sqrt(delta_n))
    
    p_n <- (exp(r * delta_n) - d_n)/(u_n - d_n)
    
    exp(-r * t) * sum(choose(n, 0:n) * p_n^(0:n) * (1 - p_n)^(n:0) 
                      * sapply(S_0 * u_n^(0:n) * d_n^(n:0) - K, max, 0))
}
```

Wycena opcji azjatyckich
------------------------

Ten przypadek jest nieco bardziej skomplikowany - wyznaczenie analitycznego wzoru na cenę opcji w tym modelu byłoby bardzo skomplikowane. Zamiast tego można zastosować podejście Monte Carlo - zasymulować pewną dużą liczbę trajektorii procesu ceny w oparciu o model dwumianowy, dla każdej z nich obliczyć wypłatę w momencie `t`, a następnie wszystkie wypłaty uśrednić i zdyskontować na moment 0.

``` r
priceAsian <- function(S_0, # cena instrumentu bazowego w chwili 0
                       K, # cena wykonania
                       r, # intensywność oprocentowania
                       sigma, # zmienność
                       t, # czas do wykupu opcji
                       n, # na ile odcinków dzielony jest przedział [0,t]
                       n_average, # z ilu punktów średnia w wypłacie z opcji azjatyckiej
                       n_MC, # liczba iteracji Monte Carlo
                       plotSimulations = TRUE) { # czy prezentować symulacje na wykresie
    
    delta_n <- t/n
    u_n <- exp(sigma * sqrt(delta_n))
    d_n <- exp(- sigma * sqrt(delta_n))
    p_n <- (exp(r * delta_n) - d_n)/(u_n - d_n)
    
    # Sprawdzanie czy model dwumianowy ma szansę zadziałać
    if(delta_n >= sigma^2/r^2) stop("Nieodpowiednie dane wejściowe!")
    
    plotTrajectories <- NULL
    price <- 0
    traject <- vector()
    
    for(i in 1:n_MC) {
        
        # Generowanie ruchów cen w górę
        traject[1] <- rbinom(n = 1, size = 1, prob = p_n)
        for(j in 2:n) traject[j] <- traject[j - 1] + rbinom(n = 1, size = 1, prob = p_n)
        
        # Wyznaczanie ceny w każdym momencie przy użyciu wygenerowanych ruchów
        traject <- S_0 * u_n ^ traject * d_n ^ (1:n - traject)
        
        # Odkładanie niektórych trajektorii do rysowania wykresów
        if(plotSimulations & i%%500 == 0) plotTrajectories <- cbind(plotTrajectories, traject)
        
        # Wypłata przy danej trajektorii jest uwzględniania w finalnej cenie
        price <- price + 1/n_MC * max(mean(traject[floor(n/n_average*1:n_average)]) - K, 0)
    }
    
    # Dyskontowanie na moment 0
    price <- exp(- r * t) * price
    
    if(plotSimulations) {
        plotTrajectories <- as.data.frame(plotTrajectories)
        colnames(plotTrajectories) <- 1:ncol(plotTrajectories)
        plotTrajectories <- plotTrajectories %>% 
            gather(key = Trajektoria, value = Cena) %>% 
            mutate(Czas = rep(t/n*1:n, ncol(plotTrajectories)))
        
        g <- ggplot(data = plotTrajectories, 
                    mapping = aes(x = Czas, y = Cena, colour = Trajektoria),
                    environment = environment())
        g <- g + geom_line() + geom_hline(yintercept = K) + theme(legend.position = "none") + 
            scale_color_brewer(palette = "Spectral")
        return(list(price = price, plot = g))
    }
    
    price
}
```

Przykładowe wykorzystanie powyższej funkcji:

``` r
set.seed(1234)
priceAsian(80, 90, .05, .1, 1, 300, 100, 3000, TRUE)
```

![](Model_dwumianowy_files/figure-markdown_github/example-1.png)

Metody Monte Carlo potrafią być czasochłonne, więc zbadano czas obliczeń w zależności od argumentu `n_MC`, przy ustalonych wszystkich pozostałych. Wybrano arbitralnie `n = 500` i `n_average = 30`.

``` r
set.seed(123)
n_MC <- c(1000, 2500, 5000, 7500, 10000, 25000, 50000, 100000, 150000)

time <- vector()

for(i in 1:length(n_MC)) {
    start <- Sys.time()
    priceAsian(100, 80, .05, .1, 1, 500, 30, n_MC[i], FALSE)
    end <- Sys.time()
    time[i] <- difftime(end, start, units = "mins")
}
```

``` r
library(ggplot2)

qplot(x = n_MC, y = time, xlab = "Liczba iteracji Monte Carlo", ylab = "Czas w minutach") + 
    geom_smooth(method = "lm")
```

![](Model_dwumianowy_files/figure-markdown_github/runtimePlot-1.png)

Czas wykonania wydaje się rosnąć liniowo wraz ze wzrostem iteracji Monte Carlo, co jest bardzo zachęcającym wynikiem.
