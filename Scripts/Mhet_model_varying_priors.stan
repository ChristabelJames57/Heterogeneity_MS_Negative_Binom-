//The slab and spike code 2 code not v_raw Fixinf pr FOR REAL DATA model 1 ALL GLOBAL xpt beta n v
functions {
  array[] real sir(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];
    real beta = theta[1];
    real v    = theta[2];
    
    real expo = 1 + v * v;                 // v^2 ⇒ sign of v doesn’t matter
    real dS_dt = -beta * I * pow(S / N, expo);
    real dI_dt =  beta * I * pow(S / N, expo) - 0.2 * I;
    real dR_dt =  0.2 * I;
    
    return {dS_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> max_days;
  int<lower=1> n_states;
  real t0;
  array[n_states] int t_last;
  //array[n_states, 3] real y0;                // initial S,I,R per state not needed from R built in stan
  array[n_states] int N;                     // population per state
  array[max_days, n_states] int cases; //old      // [time, state]
  //array[max_days, n_states] real cases;//if u want cases as real new
  
  
  
  // Fix a very small spike scale to avoid identifiability problems with slab_sd
  real<lower=0> spike_sd;                    // e.g., 1e-3 or 1e-2 (set in R)
}

transformed data {
  array[0] real x_r;
  //vector[n_states] p_reported = rep_vector(0.015, n_states);// FIXED at 0.001 FOR REAL silence if not fixed
}

parameters {
  // epidemic parameters
  vector<lower=0>[n_states] beta;//states specific
  //real<lower=0> beta;//global
  //real<lower=0> phi_inv;//single phi
  real<lower=0, upper=1> p_reported;//glbal
  //vector<lower=0, upper=1>[n_states] p_reported;//silence since fixed and state specific
  
  // slab-and-spike on v (continuous mixture, no discrete indicators)
  vector<lower=0>[n_states] v;  //state specific      // unconstrained; used as v^2 in ODE
  //real<lower=0> v; //Global
  real<lower=0> slab_sd;                     // slab scale (>0)
  real<lower=0, upper=1> theta_ss;           // spike weight (prob of near-zero)
  vector<lower=0>[n_states] I0;      // initial infected (per-state), real-valued NEW
}

transformed parameters {
  array[n_states, 3] real y0; // initial S, I, R per state NEW   
  array[n_states, max_days, 3] real y;
  array[n_states, max_days - 1] real incidence;
  // real phi = 1.0 / phi_inv;//if single
  //vector[n_states] phi;  // per state dispersion
  real MeanV;
  real meanp_reported;//leave it
  real meanbeta;
  real meanI0; // add mean of I0 NEW
  // elementwise reciprocal
  // phi = 1.0 ./ phi_inv;//every state has its own
  
  // Construct y0 from I0 and N ( N is integer; y0 elements are real) NEW
  for (s in 1:n_states) {
    y0[s, 1] = N[s] - I0[s]; // S0 = N - I0
    y0[s, 2] = I0[s];        // I0
    y0[s, 3] = 0;            // R0 = 0
  }
  
  for (s in 1:n_states) {
    array[2] real theta = { beta[s], v[s] };//stste specific beta
    // array[2] real theta = { beta, v[s] };//GLobal beta
   //array[2] real theta = { beta, v };//GLobal v n Global beta
   // array[2] real theta = { beta[s], v };//GLobal v and beta state
    array[1] int x_i = { N[s] };
    
    int tp = t_last[s];
    array[tp] real ts;//data in days
    for (i in 1:tp) ts[i] = i;
    
    // (Optionally add tolerances: ..., 1e-6, 1e-6, 1e3)
    array[tp, 3] real y_tmp
    = integrate_ode_rk45(sir, y0[s], t0, ts, theta, x_r, x_i);
    
    for (i in 1:tp) y[s, i] = y_tmp[i];
    
    // Expected new cases; keep strictly positive for NB
    for (i in 1:(tp - 1)) {
      real deltaS = y[s, i, 1] - y[s, i + 1, 1];
      //real mu = deltaS * p_reported[s] + 1e-9;   // small offset prevents 0
      //incidence[s, i] = fmax(mu, 1e-9);
       //real mu = deltaS * p_reported[s]; //state specidic
      real mu = deltaS * p_reported;//global  
      incidence[s, i] = mu +1e-03;//added to avoid sd running to zero. incidence cant be 0
      
    }
  }
  
  MeanV = sum(v) / n_states;//states specific
  //meanp_reported = sum(p_reported) / n_states;//state specific
  meanp_reported = p_reported;//global
  meanbeta = sum(beta) / n_states;//stat specific
  //meanbeta = beta;
  meanI0 = sum(I0) / n_states;   //dded line
}


// varying one at a time
model {
  // Priors
//  varying one at a time priors now gamma alone now two at a time

  //beta ~ normal(0.5, 0.1);//epidemic good
   beta ~ normal(0.5, 0.08);// varyi or 0.45,0.1. vary mean but keep sd. or  0.45, 0.1
    //p_reported ~beta (0.06,8);//epidemic good
    p_reported ~beta (0.06,7.5);//vary worked
           
  slab_sd  ~ gamma(1,10);//good epidemic
  theta_ss ~ beta(1, 1); //  spike weight ~ Uniform(0,1)
  //slab_sd  ~ gamma(1,9);//varying
  

  // Spike-and-slab prior on v
  for (s in 1:n_states) {
    target += log_mix(theta_ss,
                      normal_lpdf(v[s] | 0, spike_sd),
                      normal_lpdf(v[s] | 0, slab_sd));
  }

  // Prior for initial infections I0[s]
  for (s in 1:n_states) {
 I0[s] ~ normal((1 / p_reported), (1 / p_reported));//used for epidemic  good
  }

  // Likelihood
  for (s in 1:n_states) {
    for (i in 1:(t_last[s] - 1)) {
      cases[i, s] ~ normal(incidence[s, i],
                           sqrt(incidence[s, i] * (1 - p_reported)));
    }
  }
}

generated quantities {
  //array[n_states] real R0;//stste specific R0
  real R0;//global R0
  //real meanR0;//identical
  array[n_states, max_days - 1] real pred_cases;
  array[n_states, max_days - 1] real log_lik;
  //real meanR0;//state specific then use
 // R0 = beta / 0.2;//global R0 alone
  R0 = mean(to_vector(beta)) / 0.2;  //  state specific beta but global R0 alone
  real meanR0 = R0; //for glabal and same with summary
  
  
  //for sate specifc R0
  for (s in 1:n_states) {
    //R0[s] = beta[s] / 0.2;//stste specific beta
    // R0[s] = beta / 0.2;//global beta but state specific R0 alone
    
    for (i in 1:(max_days - 1)) {
      if (i < t_last[s]) {
        //pred_cases[s, i] = neg_binomial_2_rng(incidence[s, i], phi);//single phi
        //log_lik[s, i] = neg_binomial_2_lpmf(cases[i, s] | incidence[s, i], phi);
        
        
        //state specific
        //pred_cases[s, i] = normal_rng( incidence[s, i],
                                           //sqrt( incidence[s, i] * (1 - p_reported[s]) ) );
       // log_lik[s, i] = normal_lpdf( cases[i, s] |
                                        // incidence[s, i],
                                      //sqrt( incidence[s, i] * (1 - p_reported[s]) ) );
        //global
        pred_cases[s, i] = normal_rng( incidence[s, i],
                                      sqrt( incidence[s, i] * (1 - p_reported ) ) );//global
        log_lik[s, i] = normal_lpdf( cases[i, s] |
                                       incidence[s, i],
                                    sqrt( incidence[s, i] * (1 - p_reported ) ) );
        
        
      } else {
        pred_cases[s, i] = -1;
        log_lik[s, i] = 0;
      }
    }
  }
  //meanR0 = mean(to_vector(R0));//for states specific
}
