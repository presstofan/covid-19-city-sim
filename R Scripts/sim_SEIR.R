library(tidyverse)

# Simulation --------------------------------------------------------------

area_to_id <- function(area, area_code_file = "tokyo_area_code.xlsx") {
  tokyo_area_code <- read.xlsx(area_code_file)
  tokyo_area_code$code <- str_pad(string = tokyo_area_code$code,
                                  width = 4, pad = "0", side = "left")
  tokyo_area_code_cell_id <- left_join(dfPop[, c("id", "code")], tokyo_area_code,  by = "code")
  
  id <- tokyo_area_code_cell_id %>%
    dplyr::filter(district == area) %>%
    dplyr::select(id) %>%
    sample_n(., 1) %>%
    as.numeric
  
  return(id)
}

area_to_id(area = "新宿区")

percent <- function(x, digits = 4, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

run_seir <- function(OD_full, OD_reduced = NULL, dfPop, seed = 324515, random_initiation = FALSE, 
                     max_sim = 600, zone0 = NULL, zone0_id = 1586, zone0_infected = 10,
                     beta = 0.75, beta_dist = "uniform", gamma = 0.3, sigma = 0.2,
                     reduced_beta_start = Inf, reduced_beta = 0.5, reduced_beta_stop = Inf,
                     reduced_OD_start = Inf, reduced_OD_stop = Inf){
  
  rowSums <- base::rowSums
  colSums <- base::colSums
  
  # checking parameters
  
  if (is.null(zone0_id) & is.null(zone0)) {
    stop("zone0_id and zone0 cannot both be NULL")
  }
  
  if (is.null(OD_reduced) & reduced_OD_start == Inf) {
    stop("OD_reduced is needed as you specified reduced_OD_start")
  }
  
  # browser()
  # debug
  
  # zone0_infected <- 100
  # random_initiation <- FALSE
  # max_sim <- 300
  # zone0 <- "新宿区" # 新宿区
  # zone0_id <- 821
  # beta <- 0.3 # the parameter controlling how often a susceptible-infected contact results in a new exposure
  # beta_dist <- "uniform"
  # gamma <- 0.1 # the rate an infected recovers and moves into the resistant phase (1/Recovery_Time)
  # sigma <- 0.2 # the rate at which an exposed person becomes infective (1/Incubation_Time)
  # OD <- OD_alltrip # OD_alltrip, OD_no_public, OD_no_school, OD_no_work, OD_no_work_school
  
  # parameters
  set.seed(seed) # set a random seed for reproducible results
  N_k <- dfPop$pop # population for each cell
  N_k_sum <- sum(N_k) # total population simulated
  locs_len <- nrow(dfPop) # number of cells simulated
  
  # randomly select a cell for the first infection if zone0_id not specified
  
  if (is.null(zone0_id)){
    zone0_id <- area_to_id(area = zone0)
  }
  
  R0 <- beta / gamma # calculate the reproduction rate
  
  # set transmission rate at each cell
  if (beta_dist == "uniform") {
    beta_vec <- rep(beta, locs_len)
  } else if (beta_dist == "gamma") {
    beta_vec <- rgamma(shape = beta / 2, scale = 2, n = locs_len)
    beta_vec[zone0_id] <- beta # this will make sure the disease will always transmit beyong the first cell
  } else {
    stop(print("Unknown distribution parameter for beta"))
  }
  
  sigma_vec <- rep(sigma, locs_len) # same incubation-infectious transition rate at each cell
  gamma_vec <- rep(gamma, locs_len) # same recovery rate at each cell
  
  # set up SEIR matrix
  SEIR <- matrix(nrow = locs_len, ncol = 4) # initiate an empty SEIR matrix
  colnames(SEIR) <- c("S", "E", "I", "R") # rename the vectors
  SEIR[, "S"] <- N_k # assign the number of successible people in each cell
  SEIR[, "E"] <- 0 # assign the number of exposed people in each cell
  SEIR[, "I"] <- 0 # assign the number of infected people in each cell
  SEIR[, "R"] <- 0 # assign the number of recovered people in each cell
  
  # first infection
  # assume no all are infectious and no one in incubation period at the beginning
  first_infections <- (dfPop$id == zone0_id) * zone0_infected
  
  if (random_initiation == TRUE){
    # radomply select cell for first infection, probability is weighted by the cell population
    zone0 <- sample(dfPop$id, zone0_infected, replace = TRUE, prob = dfPop$pop) %>%
      table %>%
      as.data.frame
    names(zone0) <- c("id", "Freq")
    zone0$id <- zone0$id %>% as.character %>% as.numeric
    first_infections <- data.frame(id = 1:locs_len)
    first_infections <- left_join(first_infections, zone0)
    first_infections <- first_infections$Freq
    first_infections <- replace_na(first_infections, 0)
  }
  
  SEIR[, "S"] <- SEIR[, "S"] - first_infections
  SEIR[, "I"] <- SEIR[, "I"] + first_infections
  
  # row normalize the SEIR matrix for keeping track of group proportions
  SEIR_n <- SEIR / rowSums(SEIR)
  SEIR_n[is.na(SEIR_n)] <- 0
  
  # make copy of the SEIR matrices 
  SEIR_sim <- SEIR
  SEIR_nsim <- SEIR_n
  
  ## -- RUN SIMULATION -- ##
  print(paste0("Checking if total population match...", sum(SEIR_sim) == sum(dfPop$pop)))
  
  # initialise vectors to store the SEIR percentage by time
  susceptible_pop_norm <- vector()
  exposed_pop_norm <- vector()
  infected_pop_norm <- vector()
  recovered_pop_norm <- vector()
  susceptible_pop_norm <- c(susceptible_pop_norm, sum(SEIR[, "S"]) / N_k_sum)
  exposed_pop_norm <- c(exposed_pop_norm, sum(SEIR[, "E"]) / N_k_sum)
  infected_pop_norm <- c(infected_pop_norm, sum(SEIR[, "I"]) / N_k_sum)
  recovered_pop_norm <- c(recovered_pop_norm, sum(SEIR[, "R"]) / N_k_sum)
  
  # store the numebr of infection by day for later use
  infected_by_day <- list()
  infected_by_day[[1]] <- SEIR_sim[, "I"]
  
  # initialise parameters for loop control
  day <- 0
  total_inflow <- 1
  total_new_exposed <- 1
  
  print("Starting simulation...")
  while (total_new_exposed > 0 & day < max_sim){
    day <- day + 1
    print(paste0("Day: ",day))
    
    # model
    # S(t) = S(t-1) - New E (from local) - New E (from foreign)
    # E(t) = E(t-1) + New E (from local) + New E (from foreign) - New I
    # I(t) = I(t-1) + New I - New R
    # R(t) = R(t-1) + New R
    
    # apply lockdown effect (reduce beta)
    if (reduced_beta_start != Inf) {
      if (day >= reduced_beta_start & day < reduced_beta_stop) {
        print(paste0("Reducing beta to ", reduced_beta))
        beta_vec <- rep(reduced_beta, locs_len)
      }
    }
    
    OD <- OD_full
    if (reduced_OD_start != Inf) {
      if (day >= reduced_OD_start & day < reduced_OD_stop) {
        print("Reducing traffic...")
        OD <- OD_reduced # use the raw OD flow number 
      }
    }
    
    # New E
    infected_mat <- replicate(locs_len, SEIR_nsim[, "I"])
    OD_infected <- round(OD * infected_mat) # use the raw OD flow number 
    inflow_infected <- colSums(OD_infected)
    total_inflow_infected <- sum(inflow_infected)
    print(paste0("Total infected inflow: ", total_inflow_infected))
    new_exposed <- 
      beta_vec * SEIR_sim[, "S"] * inflow_infected / (N_k + colSums(OD)) + # exposed by contacting with imported infectious cases 
      beta_vec * SEIR_sim[, "S"] * SEIR_sim[, "I"] / N_k # exposed by contacting with local infected cases
    new_exposed[is.na(new_exposed)] <- 0
    total_new_exposed <- round(sum(new_exposed))
    print(paste0("New exposed: ", total_new_exposed))
    new_exposed <- ifelse(new_exposed > SEIR_sim[, "S"], SEIR_sim[, "S"], new_exposed) # make sure the N exposed is not bigger than total susceptible
    
    # New I
    new_infected <- sigma_vec * SEIR_sim[, "E"]
    total_new_infected <- round(sum(new_infected, na.rm = T))
    print(paste0("New infected: ", total_new_infected))
    
    # New R
    new_recovered <- gamma_vec * SEIR_sim[, "I"]
    total_new_recovered <- round(sum(new_recovered, na.rm = T))
    print(paste0("New recovered: ", total_new_recovered))
    
    SEIR_sim[, "S"] <- SEIR_sim[, "S"] - new_exposed
    SEIR_sim[, "E"] <- SEIR_sim[, "E"] + new_exposed - new_infected
    SEIR_sim[, "I"] <- SEIR_sim[, "I"] + new_infected - new_recovered
    SEIR_sim[, "R"] <- SEIR_sim[, "R"] + new_recovered
    SEIR_sim <- ifelse(SEIR_sim < 0, 0, SEIR_sim)
    
    # recompute the normalized SEIR matrix
    SEIR_nsim <- SEIR_sim / rowSums(SEIR_sim)
    SEIR_nsim[is.na(SEIR_nsim)] <- 0
    S <- sum(SEIR_sim[, "S"]) / N_k_sum
    E <- sum(SEIR_sim[, "E"]) / N_k_sum
    I <- sum(SEIR_sim[, "I"]) / N_k_sum
    R <- sum(SEIR_sim[, "R"]) / N_k_sum
    print(paste("S:", percent(S), 
                "E:", percent(E),
                "I:", percent(I),
                "R:", percent(R),
                "Total:", N_k_sum, 
                sep = " "))
    
    susceptible_pop_norm <- c(susceptible_pop_norm, S)
    exposed_pop_norm <- c(exposed_pop_norm, E)
    infected_pop_norm <- c(infected_pop_norm, I)
    recovered_pop_norm <- c(recovered_pop_norm, R)
    
    # output the infected matrix
    infected_by_day[[day + 1]] <- round(SEIR_sim[, "I"])
    
  }
  
  df_pop_norm <- data.frame(day = 1:(day + 1),
                            s_p = susceptible_pop_norm,
                            e_p = exposed_pop_norm,
                            i_p = infected_pop_norm,
                            r_p = recovered_pop_norm)
  
  # total infected
  print(paste0("total infected: ", percent(1 - tail(df_pop_norm$s_p, 1))))
  
  # output
  sim_output <- list("OD_full" = OD_full,
                     "OD_reduced" = OD_reduced,
                     "reduced_beta_start" = reduced_beta_start, 
                     "reduced_beta_stop" = reduced_beta_stop,
                     "reduced_beta" = reduced_beta,
                     "reduced_OD_start" = reduced_OD_start,
                     "reduced_OD_stop" = reduced_OD_stop,
                     "beta" = beta,
                     "beta_dist" = beta_dist,
                     "seed" = seed,
                     "gamma" = gamma,
                     "sigma" = sigma,
                     "zone0_id" = zone0_id,
                     "zone0_infected" = zone0_infected,
                     "curves" = df_pop_norm,
                     "infected_by_day" = infected_by_day)
  
  print("Simulation completed")
  return(sim_output)
}

sim_all_seir <- run_seir(OD_full = OD_alltrip, dfPop = dfPop, seed = 342125)

sim_no_public_seir <- run_seir(OD_full = OD_alltrip, OD_reduced = OD_no_public, dfPop = dfPop, 
                               reduced_OD_start = 40, seed = 342125)

sim_only_essential_seir <- run_seir(OD_full = OD_alltrip, OD_reduced = OD_essential, dfPop = dfPop, 
                                    reduced_OD_start = 40,
                                    seed = 342125)

sim_lockdown_40_seir <- run_seir(OD_full = OD_alltrip, OD_reduced = OD_essential, dfPop = dfPop, 
                                 reduced_OD_start = 40, reduced_beta_start = 40, reduced_beta = 0.5,
                                 seed = 342125)

sim_lockdown_60_seir <- run_seir(OD_full = OD_alltrip, OD_reduced = OD_essential, dfPop = dfPop, 
                                 reduced_OD_start = 60, reduced_beta_start = 60, reduced_beta = 0.5,
                                 seed = 342125)

sim_lockdown_40_stop_seir <- run_seir(OD_full = OD_alltrip, OD_reduced = OD_essential, dfPop = dfPop, 
                                      reduced_OD_start = 40, reduced_OD_stop = 70,
                                      reduced_beta_start = 40, reduced_beta = 0.5, reduced_beta_stop = 70,
                                      seed = 342125)

saveRDS(sim_all_seir, "sim_all_seir.rds")
saveRDS(sim_no_public_seir, "sim_no_public_seir.rds")
saveRDS(sim_only_essential_seir, "sim_only_essential_seir.rds")
saveRDS(sim_lockdown_40_seir, "sim_lockdown_40_seir.rds")
saveRDS(sim_lockdown_60_seir, "sim_lockdown_60_seir.rds")
saveRDS(sim_lockdown_40_stop_seir, "sim_lockdown_40_stop_seir.rds")
