## simulate LD GWAS effect sizes
set.seed(123)

#specify parameters
n_loci = 2E3
prop_polymorphic <- 1
f_at_loci <- rbinom(n_loci, 1, prop_polymorphic) * rbeta(n_loci, 1, 9)
n_indiv = 7.5E2
prop_causal_loci = 0.1
n_generations = 1E2
recomb_prob = 0.1
mutation_rate_per_generation_per_locus = 1E-5
heritability <- 0.4

#causal loci wrt focal trait
n_causal_loci <- round(prop_causal_loci * n_loci)
causal_loci <- sample(1:n_loci, size = n_causal_loci, replace = F)
loci_effect_sizes <- rep(0, n_loci)
loci_effect_sizes[causal_loci] <- rnorm(n_causal_loci)

#selective advantage of loci to induce LD?
prop_fitness_loci <- 0.2
n_changes_to_fitness_loci <- 4
gens_when_fitness_loci_change <- sort(sample(1:n_generations, n_changes_to_fitness_loci))
n_fitness_loci <- round(prop_fitness_loci * n_loci)
fitnesss_loci <- sample(1:n_loci, size = n_fitness_loci, replace = F)
fitness_loci_relative_effect_sizes <- rep(1, n_loci)
fitness_loci_relative_effect_sizes[fitnesss_loci] <- exp(rnorm(n_fitness_loci, 0, 2))

#simulate initial population of individuals -- at linkage equilibrium, HW proportions
t0 <- lapply(1:n_indiv, function(indiv){
  cbind(chr1 = rbinom(n_loci, size = 1, prob = f_at_loci),
        chr2 = rbinom(n_loci, size = 1, prob = f_at_loci))
})

## simulate next generation (constant pop size, random mating, contant mutation, homogenous recombination probability throughout all loci)

single_indiv_mutation_recombination_randomMating <- function(n_indiv, current_pop, n_loci, 
                                                             mutation_rate_per_generation_per_locus, recomb_prob, current_pop_fitness){
  
  #simulate random mating
  indivs <- sample(1:n_indiv, 2, replace = F, prob = current_pop_fitness)
  i1 <- current_pop[[indivs[1]]]
  i2 <- current_pop[[indivs[2]]]
  
  #simulate mutation
  i1_mutation <- matrix(rbinom(n = n_loci*2, size = 1, prob = mutation_rate_per_generation_per_locus), ncol = 2)
  i1 <- (i1 + i1_mutation) %% 2
  i2_mutation <- matrix(rbinom(n = n_loci*2, size = 1, prob = mutation_rate_per_generation_per_locus), ncol = 2)
  i2 <- (i2 + i2_mutation) %% 2
  
  #simulate crossing over
  i1_xover_loc <- ifelse(rbinom(1, 1, recomb_prob) == 1, sample(1:n_loci, 1), NA)
  i2_xover_loc <- ifelse(rbinom(1, 1, recomb_prob) == 1, sample(1:n_loci, 1), NA)
  
  if(!is.na(i1_xover_loc)){
    i1_temp <- i1
    i1[i1_xover_loc:n_loci,1] <- i1_temp[i1_xover_loc:n_loci,2]
    i1[i1_xover_loc:n_loci,2] <- i1_temp[i1_xover_loc:n_loci,1]
    rm(i1_temp)
  }
  if(!is.na(i2_xover_loc)){
    i2_temp <- i2
    i2[i2_xover_loc:n_loci,1] <- i2_temp[i2_xover_loc:n_loci,2]
    i2[i2_xover_loc:n_loci,2] <- i2_temp[i2_xover_loc:n_loci,1]
    rm(i2_temp)
  }
  
  #generate offspring
  chromosomes_offspring <- sample(1:2, 2, T)
  offspring <- cbind(i1[,chromosomes_offspring[1]], i2[,chromosomes_offspring[2]])
  return(offspring)
  
}

#simulate population through time
current_pop <- t0
allele_freqs_through_time <- matrix(0, ncol = n_loci, nrow = n_generations) 
allele_freqs_through_time[1,] <- apply(do.call(cbind, current_pop), 1, mean)
fitness_effects <- matrix(0, ncol = n_loci, nrow = n_changes_to_fitness_loci + 1) 
fitness_effects[1,] <- fitness_loci_relative_effect_sizes
for(time in 2:n_generations){
  
  cat(paste0(time, " "))
  
  #change fitness landscape if it's time
  if(time %in% gens_when_fitness_loci_change){
    fitnesss_loci <- sample(1:n_loci, size = n_fitness_loci, replace = F)
    fitness_loci_relative_effect_sizes <- rep(1, n_loci)
    fitness_loci_relative_effect_sizes[fitnesss_loci] <- exp(rnorm(n_fitness_loci, 0, 2))
    fitness_effects[match(time, gens_when_fitness_loci_change)+1,] <- fitness_loci_relative_effect_sizes
  }
  
  #calculate fitness of current gen
  current_pop_fitness <- sapply(1:n_indiv, function(indiv) sum(apply(current_pop[[indiv]], 1, sum) * fitness_loci_relative_effect_sizes))
  current_pop_fitness <- current_pop_fitness / sum(current_pop_fitness)
  
  #sample next generation
  current_pop <- replicate(n_indiv, simplify = F,
                          single_indiv_mutation_recombination_randomMating(n_indiv, 
                                                                           current_pop = current_pop, 
                                                                           n_loci, 
                                                                           mutation_rate_per_generation_per_locus, 
                                                                           recomb_prob,
                                                                           current_pop_fitness))
  
  #track allele frequencies
  allele_freqs_through_time[time,] <- apply(do.call(cbind, current_pop), 1, mean)

  
}

#visualize allele frequencies through time
par(mfrow = c(1,1))
plot(NULL, ylim = c(0,1.15), xlim = c(1,n_generations), xlab = "generation", ylab = "allele frequency")
gen_bounds <- c(1, gens_when_fitness_loci_change, n_generations)
cols <- colorRampPalette(c("blue", "grey10", "red"))(100)
for(i in 1:n_loci){
  for(j in 1:(n_changes_to_fitness_loci+1)){
    lines(gen_bounds[j]:gen_bounds[j+1], allele_freqs_through_time[gen_bounds[j]:gen_bounds[j+1],i], lwd = 2,
          col = adjustcolor(cols[log(fitness_effects[j,i]) * (49 / max(abs(range(log(fitness_effects))))) + 50], 0.5))  
  }
}
abline(v = gens_when_fitness_loci_change, col = "darkgreen", lwd = 4, lty = 2)
legend(x = "topleft", legend = c("change in selection", "positive selection", "negative selection"), 
       col = c("darkgreen", "red", "blue"), lwd = 2, lty = c(2,1,1))

#simulate quantitative trait
allele_counts <- sapply(1:n_indiv, function(indiv) apply(current_pop[[indiv]], 1, sum))
genetic_values <- sapply(1:n_indiv, function(indiv) sum(allele_counts[,indiv] * loci_effect_sizes))
phenotypic_values <- genetic_values + rnorm(n_indiv, sd = sqrt(var(genetic_values) * (1-heritability) / heritability))

#fit model to quantitative trait & examine (marginally) estimated effect sizes
fits <- lapply(1:n_loci, function(locus) summary(lm(phenotypic_values ~ allele_counts[locus,]))$coefficients)
pvals <- sapply(1:n_loci, function(locus) ifelse(nrow(fits[[locus]]) == 2, fits[[locus]]["allele_counts[locus, ]","Pr(>|t|)"], NA))
pvals[is.na(pvals)] <- 1
par(mfrow = c(2,2))
barplot(-log(pvals), ylab = "-log(pval)", xlab = "pos on chr")
text(x = par('usr')[1], y = par('usr')[4], xpd = NA, pos = 4,
     labels = paste0("excess hits: ", round((1 - qvalue::pi0est(pvals)$pi0) - prop_causal_loci, 2)))
acf(log(pvals))
effect_size_estimates <- sapply(1:n_loci, function(locus) ifelse(nrow(fits[[locus]]) == 2, fits[[locus]]["allele_counts[locus, ]","Estimate"], NA))
effect_size_estimates[is.na(effect_size_estimates)] <- 0
plot(loci_effect_sizes, effect_size_estimates); abline(0, 1, col = "red", lty = 2)
legend(x = "topleft", legend = "1-to-1 line", col = "red", lty = 2)
acf(abs(effect_size_estimates))

#compute LD matrix
all_chr <- do.call(cbind, current_pop)
distance_mat <- as.matrix(dist(1:n_loci))
LDr2 <- cor(t(all_chr))
inds_to_plot <- sample(1:choose(n_loci, 2), 1E3)
plot(distance_mat[upper.tri(LDr2)][inds_to_plot], LDr2[upper.tri(LDr2)][inds_to_plot]^2, 
     xlab = "pairwise distance on chromosome", ylab = "pairwise squared pearson correlation")
allele_frequencies <- apply(all_chr, 1, mean)
pairwise_expectation <- t(t(allele_frequencies)) %*% t(allele_frequencies)
pairwise_realization <- outer(all_chr[,1], all_chr[,1], "+") == 2
for(chr in 2:ncol(all_chr)){cat(paste0(chr, " ")); pairwise_realization <- pairwise_realization + (outer(all_chr[,chr], all_chr[,chr], "+") == 2)}
pairwise_realization <- pairwise_realization / ncol(all_chr)
LDD <- pairwise_realization - pairwise_expectation 
# plot(allele_frequencies, diag(pairwise_realization)) #sanity check
plot(LDr2[upper.tri(LDr2)][inds_to_plot], LDD[upper.tri(LDD)][inds_to_plot], 
     xlab = "pairwise squared pearson correlation", ylab = "linkage disequilibrium coefficient (D)")
plot(distance_mat[upper.tri(LDr2)][inds_to_plot], LDD[upper.tri(LDD)][inds_to_plot],
     xlab = "pairwise distance on chromosome", ylab = "linkage disequilibrium coefficient (D)")

#simulate second quantitative trait
prop_causal_loci_2 <- 0.25
n_causal_loci_2 <- round(prop_causal_loci_2 * n_loci)
causal_loci_2 <- sample(1:n_loci, size = n_causal_loci_2, replace = F)
loci_effect_sizes_2 <- rep(0, n_loci)
loci_effect_sizes_2[causal_loci_2] <- rnorm(n_causal_loci_2)
genetic_values_2 <- sapply(1:n_indiv, function(indiv) sum(allele_counts[,indiv] * loci_effect_sizes_2))
coef_2traits <- 1.3
genetic_and_other_trait_values <- genetic_values_2 + phenotypic_values * coef_2traits
heritability_2 <- 0.6
residual_environmental_variance_trait2 <- ((var(genetic_values_2) + coef_2traits^2*var(genetic_values)) / heritability_2^2) - 
                                          (var(genetic_values_2) + coef_2traits^2*var(genetic_values)) - 
                                           coef_2traits^2*(var(genetic_values) * (1-heritability) / heritability)
phenotypic_values_2 <- genetic_and_other_trait_values + rnorm(n_indiv, sd = sqrt(residual_environmental_variance_trait2))
lm(phenotypic_values_2 ~ phenotypic_values)
fits_2 <- lapply(1:n_loci, function(locus) summary(lm(phenotypic_values_2 ~ allele_counts[locus,]))$coefficients)
pvals_2 <- sapply(1:n_loci, function(locus) ifelse(nrow(fits_2[[locus]]) == 2, fits_2[[locus]]["allele_counts[locus, ]","Pr(>|t|)"], NA))
pvals_2[is.na(pvals_2)] <- 1
effect_size_estimates_2 <- sapply(1:n_loci, function(locus) ifelse(nrow(fits_2[[locus]]) == 2, fits_2[[locus]]["allele_counts[locus, ]","Estimate"], NA))
effect_size_estimates_2[is.na(effect_size_estimates_2)] <- 0
plot(effect_size_estimates_2 ~ effect_size_estimates)
lm(effect_size_estimates_2 ~ effect_size_estimates)
plot(diff(effect_size_estimates_2) ~ diff(effect_size_estimates))
lm(diff(effect_size_estimates_2) ~ diff(effect_size_estimates))
