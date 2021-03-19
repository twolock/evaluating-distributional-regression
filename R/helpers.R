
rsinh <- function(n, mu, sigma, eps, delta) {
  mu + sigma * sinh((asinh(rnorm(n, 0, 1)) - eps) / delta)
}

paste2 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  gsub(paste0("(^",sep,"|",sep,"$)"),"",
       gsub(paste0(sep,sep),sep,
            do.call(paste,c(L,list(sep=sep)))))
}
sinhasinh <- custom_family(
  "sinhasinh", dpars = c("mu", "sigma", "eps", "delta"),
  links = c("identity", "log", "identity", "log"), lb = c(NA, NA, NA, NA),
  type = "real"
)

gamma2 <- custom_family(
  "gamma2", dpars = c("mu", "sigma"),
  links = c("log", "log"), lb = c(NA, NA),
  type = "real"
)

beta2 <- custom_family(
  "beta2", dpars = c("mu", "sigma"),
  links = c("logit", "logit"), lb = c(0, 0),
  ub = c(1, 0.5),
  type = "real"
)

sinhasinh2 <- custom_family(
  "sinhasinh2", dpars = c("mu", "sigma", "eps", "delta"),
  links = c("identity", "log", "identity", "log"), lb = c(NA, NA, NA, NA),
  type = "real"
)

log_lik_sinhasinh <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[, i]
  eps <- draws$dpars$eps[, i]
  delta <- draws$dpars$delta[, i]
  y <- draws$data$Y[i]
  sinhasinh_lpdf(y, mu, sigma, eps, delta)
}
predict_sinhasinh <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[i]
  eps <- draws$dpars$eps[i]
  delta <- draws$dpars$delta[i]
  sinhasinh_rng(mu, sigma, eps, delta)
}
posterior_predict_sinhasinh <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[,i]
  eps <- draws$dpars$eps[,i]
  delta <- draws$dpars$delta[,i]
  sinhasinh_rng(mu, sigma, eps, delta)
}

log_lik_gamma2 <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[, i]
  y <- draws$data$Y[i]
  
  gamma2_lpdf(y, mu, sigma)
}
predict_gamma2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[, i]
  gamma2_rng(mu, sigma)
}
posterior_predict_gamma2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[, i]
  gamma2_rng(mu, sigma)
}

log_lik_beta2 <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[, i]
  y <- draws$data$Y[i]
  beta2_lpdf(y, mu, sigma)
}
predict_beta2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[i]
  beta2_rng(mu, sigma)
}
posterior_predict_beta2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  sigma <- draws$dpars$sigma[,i]
  beta2_rng(mu, sigma)
}


stan_funs <- "
  real gamma2_lpdf(real y, real mu, real sigma) {
    real alpha;
    real beta;
    
    alpha = mu*mu/(sigma*sigma);
    beta = mu/(sigma*sigma);
    return gamma_lpdf(y | alpha, beta);
  }
  real gamma2_rng(real mu, real sigma) {
    real alpha;
    real beta;
    
    alpha = mu*mu/(sigma*sigma);
    beta = mu/(sigma*sigma);

    return gamma_rng(alpha, beta);
  }
  
  real beta2_lpdf(real y, real mu, real sigma) {
    real alpha;
    real beta;
    real adj_sigma;
    
    adj_sigma = sigma * 0.5;
    
    alpha = mu*mu*((1-mu)/(adj_sigma*adj_sigma) - 1/mu);
    beta = alpha * (1/mu - 1);
    return beta_lpdf(y | alpha, beta);
  }
  
  real beta2_rng(real mu, real sigma) {
    real alpha;
    real beta;
    real adj_sigma;
    
    adj_sigma = sigma * 0.5;

    alpha = mu*mu*((1-mu)/(adj_sigma*adj_sigma) - 1/mu);
    beta = alpha * (1/mu - 1);
    return beta_rng(alpha, beta);
  }

  real sinhasinh_lpdf(real y, real mu, real sigma, real eps, real delta) {
    real y_z;
    real sigma_star;
    real S_y;
    real S_y_2;
    real C_y;
    real nll;
    
    nll = 0;
    sigma_star = sigma * delta;
    y_z = (y - mu)/sigma_star;
    
    S_y = sinh(eps + delta * asinh(y_z));
    S_y_2 = S_y * S_y;
    C_y = sqrt(1 + S_y_2);
    nll += -0.5 * S_y_2 - log(sigma_star) - log(sqrt(2*pi()));
    nll += log(delta) + log(C_y) - log(sqrt(1 + y_z*y_z));
    return nll;
  }
  real sinhasinh_rng(real mu, real sigma, real eps, real delta) {
    return mu + sigma * delta * sinh((asinh(normal_rng(0, 1)) - eps)/delta);
  }
  
"
