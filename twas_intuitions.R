#metaparameters
n_indiv <- 5000
n_loci <- 200

#parameters
snps <- t(replicate(n_loci, rbinom(n_indiv, 2, rbeta(1,1,1))))
b_snp_expr <- rnorm(n_loci)
b_snp_trait <- rnorm(n_loci)
b_expr_trait <- 0

#simulate expression
genetic_expression <- 0 + c(b_snp_expr %*% snps)
expression_noise <- rnorm(n_indiv, 0, sqrt(var(genetic_expression))) #hack to get heritability = 0.5
expression <- genetic_expression + expression_noise

#simulate trait
direct_genetic_trait <- 0 + c(b_snp_trait %*% snps)
genetic_trait <- direct_genetic_trait + genetic_expression * b_expr_trait
genetic_expression_trait <- genetic_trait + expression_noise * b_expr_trait
trait <- trait  + rnorm(n_indiv, 0, sqrt(var(trait))) #dunno what heritablity would be here exactly, but let's just add half again as much variance as noise

#estimate genetic component of expression 
expr_weights <- sapply(1:n_loci, function(locus_i) lm(expression ~ snps[locus_i,])$coefficients[2])
# plot(expr_weights, b_snp_expr); abline(0,1,col=2,lwd=2)

#perform twas
predicted_genetic_expression <- 0 + c(expr_weights %*% snps)
b_expr_trait
summary(lm(trait ~ predicted_genetic_expression))$coefficients
# plot(b_snp_trait, lm(trait ~ t(snps))$coefficients[-1]); abline(0,1,col=2,lwd=2)
cor.test(genetic_expression, genetic_expression_trait)
summary(lm(trait ~ genetic_expression))$coefficients


