#---------------------------------------------------------------------#
# Modeling Cyclic Patterns Using a Two-Stage Hybrid Bayesian Approach
#
# Real Data Analysis Example
#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#
#### Required Packages ####
#
# This R script requires the following packages:
#  1. `rblimp`
#  2. `remotes` (required to install below)
#  3. `mitml`
#
#---------------------------------------------------------------------#

## Install `rblimp` package from Github (Uncomment line below to install)
# remotes::install_github('blimp-stats/rblimp')

## Load required package
library(rblimp)
library(mitml)

#---------------------------------------------------------------------#
#### Reading Data ####
#
# Load data sets from CSV files for real data analysis. Next, clean
#  the data set to select only necessary variables
#
#---------------------------------------------------------------------#

# Read data in
original_data <- read.csv("real data example.csv", header = TRUE)

# Clean data
mydata <- data.frame(
    id = original_data$id,
    relationship = (original_data$relationship.status == 'in_relationship') |> 
        as.numeric(),
    desire = original_data$sexual.desire,
    time = original_data$t,
    omega = 2 * pi / (original_data$cycle.length - 1)
)

#---------------------------------------------------------------------#
#### Specify `rblimp` model ####
#
# Using Blimp's modeling syntax, we can specify the cyclic model.
#  Blimp allows for separating the model 'blocks' by names. Below,
#  specify the between-person (`level_1`) and within-person blocks
#  (`level_2`). Next, we create a `quantities` block to label the
#  generated quantities using the random effects.
#
#---------------------------------------------------------------------#

# Specify `rblimp` model
model <- list(
    # Person Level
    level_2 = c(
        'u ~ 1',
        'r1 ~ 1 relationship',
        'r2 ~ 1 relationship',
        'u r1 r2 ~~ u r1 r2'
    ),
    # Time Level
    level_1 = c(
        'desire ~ 1@u (cos(omega * time))@r1 (-sin(omega * time))@r2'
    ),
    # Calculate Transformed Quantities
    quantities = c(
        'Ri = sqrt( (r1)^2 + (r2)^2 )',
        'phii = switch(
            r1 > 0, atan(r2 / r1),
            r1 < 0, atan(r2 / r1) + _pi,
            r2 >= 0 and r1 == 0, _pi / 2,
            r2 <= 0 and r1 == 0, 3 * _pi / 2,
            nan # Should never happen
        )',
        'mui = u'
    )
)

#---------------------------------------------------------------------#
#### Run model via `rblimp`  ####
#
# rblimp is an interface to the Blimp computational engine. 
#  Below, we incorporate the model syntax above (`model`), the data
#  set above (`data`), specify latent variables (i.e., random
#  effects), and clustering variable. In addition, we incorporate
#  specifics for how long to run and how many imputed data sets to
#  save for the second stage. Finally, setting `add_save = FALSE`
#  will reduce some of the unused data in each imputed data set.
#---------------------------------------------------------------------#

# Run model
results <- rblimp(
    model,
    mydata,
    clusterid = 'id',
    latent = 'id = u r1 r2',
    nimps = 1000,
    burn = 10000,
    iter = 10000,
    seed = 197321,
    options = 'prior1',
    # Reduce size of imputed data sets
    add_save = FALSE
)

# Print out psr diagnostics
psr(results)

# Display trace plots
trace_plot(results)

# Print output for the run model
output(results)

# Display posterior plots
posterior_plot(results)

#---------------------------------------------------------------------#
#### Second stage estimation  ####
#
# Now that the imputations have been generated, we can analyze the
#  imputed latent variables specified in the `quanitites` block.
#  The imputations are saved in the `imputations` slot as a list.
#  First, Blimp outputs these as a level-1 data set, so we must
#  convert these to level-2 data sets.
#
#---------------------------------------------------------------------#

# Get imputations of random effects.
# Convert to level-2 units
imps <- results@imputations |> 
    lapply(aggregate, cbind(ri, mui, phii, relationship) ~ id, unique)

#---------------------------------------------------------------------#
#
# Next we can fit an OLS regression for each latent variable based
#  on the desire model. The `mitml` package can handle this for us,
#  automatically pooling the results.
#
#---------------------------------------------------------------------#

# Convert to mitml.list
imps <- as.mitml.list()

## Fit model via lm on each imputation and obtain estimates
mdl_ri <- with(imps, lm(ri ~ relationship)) |> testEstimates()
mdl_mui <- with(imps, lm(mui ~ 1)) |> testEstimates()
mdl_phii <- with(imps, lm(phii ~ 1)) |> testEstimates()

## Print Results and CI's for each model

mdl_ri
# Final parameter estimates and inferences obtained from 1000 imputed data sets.
# 
#               Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)      0.304     0.044     6.941  3887.817     0.000     1.028     0.507 
# relationship     0.014     0.062     0.222  5587.122     0.824     0.733     0.423 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

mdl_ri |> confint()
#                   2.5 %    97.5 %
# (Intercept)   0.2181164 0.3898584
# relationship -0.1078341 0.1353799

mdl_mui
# testEstimates(model = with(imps, lm(mui ~ 1)))
# 
# Final parameter estimates and inferences obtained from 1000 imputed data sets.
# 
#              Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)     2.202     0.108    20.443 2.931e+05     0.000     0.062     0.058 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

mdl_mui |> confint()
#                2.5 %   97.5 %
# (Intercept) 1.991008 2.413258

mdl_phii
# Final parameter estimates and inferences obtained from 1000 imputed data sets.
# 
#              Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)     1.584     0.390     4.060  3528.584     0.000     1.137     0.532 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

mdl_phii |> confint()
#                 2.5 %   97.5 %
# (Intercept) 0.8192295 2.349465

#---------------------------------------------------------------------#
#
# Alternatively, we can manually compute these within each
#  imputation using closed form solution formulas. Below we
#  illustrate how to do this and pool the results 'by hand'.
#
#---------------------------------------------------------------------#

## Compute manually
# Compute across 1000 imputed data sets
s <- sapply(imps, \(x) {
    # Compute Parameters
    Beta <- with(x, {
        sum((ri-mean(ri))*(relationship-mean(relationship))) /
            sum((relationship-mean(relationship))*(relationship-mean(relationship)))
    })
    R <- with(x, mean(ri)-Beta*mean(relationship))
    phi <- with(x, mean(phii))
    mu <- with(x, mean(mui))
    
    # Compute within imputation variance
    Rih <- with(x, R + Beta*relationship)
    rse <- with(x, sqrt(sum((ri - Rih)^2) / (NROW(ri) - 2)))
    w.Beta <- with(x, {
        ( rse / sqrt(sum((relationship - mean(relationship)) *
                             (relationship - mean(relationship)))))^2
    })
    w.R <- with(x, {
        (rse*sqrt(1/NROW(x)+mean(relationship)^2/
                      (sum((relationship-mean(relationship))*
                               (relationship - mean(relationship))))))^2
    })
    w.phi<- with(x, sd(phii)^2 / NROW(phii))
    w.mu <- with(x, sd(mui)^2 / NROW(mui))
    
    # Return
    c(
        # summaries
        Beta = Beta, R = R, phi = phi, mu = mu,
        # within variance 
        w.Beta = w.Beta, w.R = w.R, w.phi = w.phi, w.mu = w.mu
    )
}) 


### compute MI SE
vw.R <- mean(s['w.R', ])
vb.R <- var(s['R', ])
se.R <- sqrt(vw.R + vb.R + vb.R / NROW(s['R', ]))
CI.R <- c(mean(s['R', ]) - 1.96*se.R, mean(s['R', ]) + 1.96*se.R)

vw.Beta <- mean(s['w.Beta', ])
vb.Beta <- var(s['Beta', ])
se.Beta <- sqrt(vw.Beta + vb.Beta + vb.Beta / NROW(s['Beta', ]))
CI.Beta <- c(mean(s['Beta', ]) - 1.96*se.Beta, mean(s['Beta', ]) + 1.96*se.Beta)

vw.mu <- mean(s['w.mu', ])
vb.mu <- var(s['mu', ])
se.mu <- sqrt(vw.mu + vb.mu + vb.mu / NROW(s['mu', ]))
CI.mu <- c(mean(s['mu', ]) - 1.96*se.mu, mean(s['mu', ]) + 1.96*se.mu)

vw.phi <- mean(s['w.phi', ])
vb.phi <- var(s['phi', ])
se.phi <- sqrt(vw.phi + vb.phi + vb.phi / NROW(s['phi', ]))
CI.phi <- c(mean(s['phi', ]) - 1.96*se.phi, mean(s['phi', ]) + 1.96*se.phi)

# Create Results Table
rbind(
    R = c(mean = mean(s['R', ]), se = se.R, lower = CI.R[1], upper = CI.R[2]),
    Beta = c(mean = mean(s['Beta', ]), se = se.Beta, lower = CI.Beta[1], upper = CI.Beta[2]),
    mu = c(mean = mean(s['mu', ]), se = se.mu, lower = CI.mu[1], upper = CI.mu[2]),
    phi = c(mean = mean(s['phi', ]), se = se.phi, lower = CI.phi[1], upper = CI.phi[2])
)
#           mean         se      lower     upper
# R    0.3039874 0.04379889  0.2181416 0.3898332
# Beta 0.0137729 0.06203206 -0.1078099 0.1353557
# mu   2.2021330 0.10771852  1.9910047 2.4132613
# phi  1.5843472 0.39023942  0.8194779 2.3492164

