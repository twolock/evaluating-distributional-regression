# Analysis script for the probability distribution
# comparison in Wolock et al. (2021). Takes a model
# ID, `m_val`, identified in 'data/sa_model_list.csv',
# reads prepared data from 'data/', fits a model in
# BRMS, and saves a number of outputs

library(data.table)
library(brms)
library(rstan)
library(ggplot2)
source('R/helpers.R')

beta_max <- 150
result_dir <- '201013_single_age'

if (length(commandArgs(trailingOnly = T)) == 0) {
  m_val <- 2
} else {
  args <- commandArgs(trailingOnly = T)
  arg.l <- do.call(rbind, strsplit(args, '\\s*=\\s*'))
  formatted.args <- arg.l[,2]
  names(formatted.args) <- arg.l[,1]
  m_val <- as.integer(formatted.args['m_val'])
}

if (Sys.info()['sysname'] != 'Darwin') {
  output.file <- sprintf('logs/%s_brms.out', format(Sys.Date(), '%y%m%d'))
  sink(output.file)
}
set.seed(m_val)

# Get current model configuration
model.dt <- fread('data/sa_model_list.csv')
curr.config <- as.list(model.dt[model_id == m_val])

dep_var <- curr.config$dep_var
likelihood <- curr.config$likelihood
dataset <- curr.config$dataset

# Set dependent variable specification
use_log <- F
use_ratio <- F
use_diff <- F

switch (dep_var,
        'logratio' = {
          use_log <- T
          use_ratio <- T
        },
        'logage' = {
          use_log <- T
        },
        'diff' = {
          use_diff <- T
        },
        'linearage' = {},
        stop(sprintf('Independent variable %s is unknown', dep_var))
)

if (use_ratio & use_diff) {
  stop('Cannot model ratio and difference simultaneously')
}

# Set link and inverse link
if (use_log) {
  link <- log
  inv_link <- exp
} else {
  link <- identity
  inv_link <- identity
}

# Define function to make dependent variable
# and its inverse
make_dep_var <- function(p_age, age) {
  if (use_ratio) {
    return (link(p_age/age))
  } else if (use_diff) {
    return (p_age - age)
  } else {
    return (link(p_age))
  }
}
inv_dep_var <- function(y, age) {
  if (use_ratio) {
    return (inv_link(y)*age)
  } else if (use_diff) {
    return ((y + age))
  } else {
    return (inv_link(y))
  }
}

# Read prepared age-/sex-specific data
age.dt <- data.table(readRDS(sprintf('data/single age data %i.Rds', dataset)))

# Create stanvars object for BRMS
stanvars <- stanvar(scode = stan_funs, block = "functions")

# Restrict to target ages
age.dt <- age.dt[!is.na(sex) & age >= 15 & age < 65 & p_age >= 10]
age.dt[, age_diff := p_age-age]

# Change second argument of `sample(1:.N, .N)` to a desired sample size
# for subsampling
keep.i <- 1:nrow(age.dt)
# keep.i <- sample(1:nrow(age.dt), 2000)
subset.dt <- age.dt[keep.i]
subset.dt[, age := age - age %% 5]
subset.dt[, tmp_id := .I]

# Scale respondent age
age.mean <- subset.dt[, mean(age)]
age.sd <- subset.dt[, sd(age)]
subset.dt[, scaled_age := (age-age.mean)/age.sd]
subset.dt[, dep_var := make_dep_var(p_age, age)]

# Reflect and scale if using beta or gamma
if (likelihood == 'beta') {
  subset.dt[, dep_var := dep_var / beta_max]
}
if (likelihood == 'gamma') {
  subset.dt[, gamma_reflect := (sex == 'Male')*-2 +1]
  subset.dt[, dep_var := dep_var * gamma_reflect + beta_max]
}

# Create data set for prediction
pred.dt <- CJ(sex = subset.dt[, unique(sex)], age = min(subset.dt[, age]):max(subset.dt[, age]))
pred.dt[, scaled_age := (age-age.mean)/age.sd]
pred.dt[, Var2 := .I]

# Get BRMS family associated with configuration
# beta2 and gamma2 use mean-variance parametrisations
switch(likelihood,
       'norm' = {brm.fam <- 'gaussian'},
       'beta' = {brm.fam <- 'beta2'},
       'gamma' = {brm.fam <- 'gamma2'},
       'skewnorm' = {brm.fam <- 'skew_normal'},
       'sinh' = {brm.fam <- 'sinhasinh'},
       stop('Invalid distribution')
)
if (likelihood %in% c('gamma','sinh', 'beta')) {
  dpars <- get(brm.fam)$dpars
} else {
  dpars <- get(brm.fam)()$dpars
}
if (likelihood == 'norm') {
  dpars <- c('mu', 'sigma')
}

# Set priors
priors <- set_prior('normal(0, 5)', class = 'Intercept', dpar = c('', dpars[-1]))

# Generate baseline formula list
# (All constant in distribution comaprison)
f.list <- lapply(dpars, function(d) {
  formula(sprintf('%s ~ 1', d))
})
# Set outcome variable in first formula
f.list[[1]] <- formula(dep_var ~ 1)
f <- brmsformula(f.list[[1]], flist = f.list[-1])

# Create a junk variable for prediction
pred.dt[, dep_var := 0.5]

# Fit model
brm.fit <- brm(
  formula = f,
  prior = priors,
  data = data.frame(subset.dt),
  family = get(brm.fam),
  stanvars = stanvars,
  iter = 3000,
  warmup = 2500,
  cores = 4,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
)
# Get stan functions for posterior prediction
expose_functions(brm.fit, vectorize = T, show_compiler_warnings=F)

# Get posterior dpar predictions
draws <- extract_draws(restructure(brm.fit), newdata = pred.dt)
for (dp in names(draws$dpars)) {
  draws$dpars[[dp]] <- as.matrix(brms:::get_dpar(draws, dpar = dp))
}

# Get posterior predictive samples
post.samples <- dcast(data.table(melt(draws$dpars)), Var1 + Var2 ~ L1)
post.samples <- merge(post.samples, melt(predict(brm.fit, newdata = pred.dt, summary = F)))
setnames(post.samples, 'value', 'y_post')
post.samples <- merge(post.samples, pred.dt[, .(Var2, age, sex)], by=c('Var2'))
if (likelihood == 'beta') {
  post.samples[, y_post := y_post * beta_max]
}
if (likelihood == 'gamma') {
  post.samples[, gamma_reflect := (sex == 'Male')*-2 +1]
  post.samples[, y_post := (y_post - beta_max)/gamma_reflect]
}

post.samples[, p_age := inv_dep_var(y_post, age)]

# Combine everything
post.dt <- post.samples[, .(dep_var, likelihood, dataset, age, sex, sample_i=Var1, y_post, p_age)]
post.dt <- cbind(post.dt, post.samples[, mget(dpars)])
saveRDS(post.dt,
        file = sprintf('data/results/%s/post %i.Rds', result_dir, m_val))

# Calculate posterior summaries
par.summary <- melt(post.samples, id.vars = c('age', 'sex', 'Var1', 'Var2'))[,
  .(l = quantile(value, 0.025),
    m = median(value),
    u = quantile(value, 0.975)),
by=.(age, sex, variable)]


# Get posterior density and adjust as necessary
post.dens.dt <- melt(data.table(subset.dt[, .(age, sex, p_age, tmp_id)],
                                t(exp(log_lik(brm.fit)))),
                     id.vars = c('age','sex', 'p_age', 'tmp_id'))
if (use_log) {
  post.dens.dt[, value := value / p_age]
}
if (likelihood == 'beta') {
  post.dens.dt[, value := value / beta_max]
}
if(likelihood == 'gamma') {
  post.dens.dt[sex == 'Male', value := value * 1]
}

# Get WAIC (not used in the paper)
wide.waic.dt <- dcast(post.dens.dt, tmp_id + age + sex + p_age ~ variable, value.var = 'value')
waic.m <- as.matrix(wide.waic.dt[, -c(1,2,3,4)])
waic.m[waic.m == 0] <- .Machine$double.eps

lppd <- sum(log(rowMeans(waic.m)))
pwaic2 <- sum(apply(log(waic.m), 1, function(x) {
  sum((x - mean(x))^2) / (ncol(waic.m)-1)
}))

waic.dt <- data.table(likelihood = likelihood,
                      dep_var = dep_var,
                      dataset = dataset,
                      lppd = lppd,
                      pwaic = pwaic2,
                      waic = lppd - pwaic2,
                      n_par = length(names(brm.fit$fit))-1
)

# Get LOO-CV
# Needs to be done semi-manually to compare
# across different dependent variables
loo.cols <- 1:ncol(waic.m)
loo.m <- log(t(waic.m[, loo.cols]))

r_eff <- brms:::r_eff_log_lik(exp(loo.m), fit=brm.fit, FALSE)

loo.res <- loo(loo.m, r_eff = r_eff)

waic.dt[, (names(loo.res$estimates[, 'Estimate'])) := as.list(loo.res$estimates[, 'Estimate'])]
saveRDS(waic.dt, file = sprintf('data/results/%s/waic %i.Rds', result_dir, m_val))
saveRDS(loo.res, file = sprintf('data/results/%s/loo %i.Rds', result_dir, m_val))

print(brm.fit)

# Get QQ RMSE
quantiles <- seq(0, 1, 1/10)

obs.quant.dt <- subset.dt[, .(val_int = quantile(p_age, probs = quantiles),
                              val_jitter = quantile(p_age + runif(.N, -0.5, 0.5), probs = quantiles)),
                          by = .(age, sex)]
obs.quant.dt[, q_i := 1:.N, by=.(age, sex)]
obs.quant.dt[, q_val := quantiles[q_i]]

post.quant.dt <- post.dt[, .(val_est = quantile(p_age, probs = quantiles)),
                          by = .(age, sex)]
post.quant.dt[, q_i := 1:.N, by=.(age, sex)]
post.quant.dt[, q_val := quantiles[q_i]]
quant.dt <- merge(obs.quant.dt, post.quant.dt)
quant.dt[, dep_var := dep_var]
quant.dt[, likelihood := likelihood]
quant.dt[, dataset := dataset]

saveRDS(quant.dt, file = sprintf('data/results/%s/quantiles %i.Rds', result_dir, m_val))
