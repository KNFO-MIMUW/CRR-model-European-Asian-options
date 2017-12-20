
# COMPUTATIONAL FINANCE STUDENT ORGANIZATION, UNIVERSITY OF WARSAW
# JAN KOSCIALKOWSKI
# European and Asian Options Pricing Using CRR Model

# European option is priced using the derived exact formula

priceEuropean <- function(S_0, # underlying price at time 0
                          K, # strike price
                          r, # interest intensity
                          sigma, # volatility
                          t, # time until option expiration
                          n) { # number of subintervals of [0,t]
    
    delta_n <- t/n
    
    # Checking whether the CRR model will work
    if(delta_n >= sigma^2/r^2) stop("Incorrect input data!")
    
    u_n <- exp(sigma * sqrt(delta_n))
    d_n <- exp(- sigma * sqrt(delta_n))
    
    p_n <- (exp(r * delta_n) - d_n)/(u_n - d_n)
    
    exp(-r * t) * sum(choose(n, 0:n) * p_n^(0:n) * (1 - p_n)^(n:0) 
                      * sapply(S_0 * u_n^(0:n) * d_n^(n:0) - K, max, 0))
}


# Asian option is priced using a Monte Carlo simulation

priceAsian <- function(S_0, # underlying price
                       K, # strike price
                       r, # interest intensity
                       sigma, # volatility
                       t, # time until option expiration
                       n, # number of subintervals of [0,t]
                       n_average, # average over how many points in Asian option payoff
                       n_MC, # number of Monte Carlo iterations
                       plotSimulations = TRUE) { # whether to present simulations on a graph
    
    delta_n <- t/n
    u_n <- exp(sigma * sqrt(delta_n))
    d_n <- exp(- sigma * sqrt(delta_n))
    p_n <- (exp(r * delta_n) - d_n)/(u_n - d_n)
    
    # Checking whether the CRR model will work
    if(delta_n >= sigma^2/r^2) stop("Incorrect input data!")
    
    plotTrajectories <- NULL
    price <- 0
    traject <- vector()
    
    for(i in 1:n_MC) {
        
        # Generating price 'ups'
        traject[1] <- rbinom(n = 1, size = 1, prob = p_n)
        for(j in 2:n) traject[j] <- traject[j - 1] + rbinom(n = 1, size = 1, prob = p_n)
        
        # Calculating price for each moment using generated shifts
        traject <- S_0 * u_n ^ traject * d_n ^ (1:n - traject)
        
        # Saving some trajectories for plotting
        if(plotSimulations & i%%500 == 0) plotTrajectories <- cbind(plotTrajectories, traject)
        
        # Payoff for each trajectory is included in the final price
        price <- price + 1/n_MC * max(mean(traject[floor(n/n_average*1:n_average)]) - K, 0)
    }
    
    # Discounting to time 0
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
