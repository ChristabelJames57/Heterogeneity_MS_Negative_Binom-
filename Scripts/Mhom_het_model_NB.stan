  ////Mhom-het model stan
 functions {
  array[] real sir(
    real t,
    array[] real y,
    array[] real theta,
    array[] real x_r,
    array[] int x_i
  ) {
    real S = y[1];
    real I = y[2];
    real R = y[3];

    real N = x_r[1];        // EFFECTIVE population = pN * N
    real beta = theta[1];
    real v    = theta[2];

    real expo = 1 + v * v;  // v^2 ⇒ sign of v doesn’t matter

    real dS_dt = -beta * I * pow(S / N, expo);
    real dI_dt =  beta * I * pow(S / N, expo) - 0.2 * I;
    real dR_dt =  0.2 * I;

    return { dS_dt, dI_dt, dR_dt };
  }
}

data {
  int<lower=1> max_days;
  int<lower=1> n_states;
  real t0;
  array[n_states] int t_last;

  array[n_states] real N;                 //  EFFECTIVE POPULATION (pN * N)
  array[max_days, n_states] int cases;

  real<lower=0> spike_sd;
}

transformed data {
  array[0] int x_i;

  // 🔧 NEW: data-only ODE covariates (THIS FIXES THE ERROR)
  array[n_states, 1] real x_r;

  for (s in 1:n_states) {
    x_r[s, 1] = N[s];
  }
}

parameters {
  vector<lower=0>[n_states] beta; //states specific
  real<lower=0, upper=1> p_reported;

  vector<lower=0>[n_states] v; //states specific
  real<lower=0> slab_sd;
  real<lower=0, upper=1> theta_ss;

  vector<lower=0>[n_states] I0;
  real<lower=0> phi;// Negative-binomial dispersion
}

transformed parameters {
  array[n_states, 3] real y0;
  array[n_states, max_days, 3] real y;
  array[n_states, max_days - 1] real incidence;

  real MeanV;
  real meanbeta;
  real meanp_reported;
  real meanI0;

  for (s in 1:n_states) {
    y0[s, 1] = N[s] - I0[s];
    y0[s, 2] = I0[s];
    y0[s, 3] = 0;
  }

  for (s in 1:n_states) {
    array[2] real theta = { beta[s], v[s] };

    int tp = t_last[s];
    array[tp] real ts;
    for (i in 1:tp) ts[i] = i;


    // array[1] real x_r = { N[s] };

    // fIXED ODE CALL
    array[tp, 3] real y_tmp =
      integrate_ode_rk45(sir,y0[s],t0, ts,theta,x_r[s],   x_i);//  data-only

    for (i in 1:tp)
      y[s, i] = y_tmp[i];

    for (i in 1:(tp - 1)) {
      real deltaS = y[s, i, 1] - y[s, i + 1, 1];
      real mu = deltaS * p_reported;
      incidence[s, i] = mu + 1e-3;
    }
  }

  MeanV = mean(v);
  meanbeta = mean(beta);
  meanp_reported = p_reported;
  meanI0 = mean(I0);
}

//priors
model {
  beta ~ normal(0.5, 0.1);
  p_reported ~ beta(0.06, 8);

  slab_sd ~ gamma(1, 10);
  theta_ss ~ beta(1, 1);
  phi ~ exponential(1);  // NEW prior on overdispersion

  for (s in 1:n_states) {
    target += log_mix(
      theta_ss,
      normal_lpdf(v[s] | 0, spike_sd),
      normal_lpdf(v[s] | 0, slab_sd)
    );
  }

  for (s in 1:n_states) {
    I0[s] ~ normal(1 / p_reported, 1 / p_reported);
  }

//prior for 
  //for (s in 1:n_states) {
    //for (i in 1:(t_last[s] - 1)) {
      //cases[i, s] ~ normal(
       // incidence[s, i],
       // sqrt(incidence[s, i] * (1 - p_reported))
      //);
      // NEW: Negative binomial observation model prior
  for (s in 1:n_states) {
    for (i in 1:(t_last[s] - 1)) {
      cases[i, s] ~ neg_binomial_2(incidence[s, i], phi);
      
    }
  }
}

generated quantities {
  real R0;
  array[n_states, max_days - 1] real pred_cases;
  array[n_states, max_days - 1] real log_lik;

  R0 = mean(beta) / 0.2;
  real meanR0 = R0;


 // Likelihood
  for (s in 1:n_states) {
    for (i in 1:(max_days - 1)) {
      if (i < t_last[s]) {
        //pred_cases[s, i] =
          //normal_rng(incidence[s, i],
                    // sqrt(incidence[s, i] * (1 - p_reported)));
        //log_lik[s, i] =
         // normal_lpdf(cases[i, s] |
                      // incidence[s, i],
                      // sqrt(incidence[s, i] * (1 - p_reported)));
                       
       // NEW: NB predictive and log-likelihood
        pred_cases[s, i] =
          neg_binomial_2_rng(incidence[s, i], phi);
        log_lik[s, i] =
          neg_binomial_2_lpmf(cases[i, s] | incidence[s, i], phi);
          
      } else {
        pred_cases[s, i] = -1;
        log_lik[s, i] = 0;
      }
    }
  }
}
