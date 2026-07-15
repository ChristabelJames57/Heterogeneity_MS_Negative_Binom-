//Model2_Estimating PN_Mansucrip_Febru p_reported estimated
functions {
  array[] real sir(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];
    real beta = theta[1];
    real pN = theta[2]; // silence here if simulation bc not estimating pN

    // Use here if real data
    real dS_dt = -beta * I * (S / (pN * N));
    real dI_dt = beta * I * (S / (pN * N)) - 0.2 * I;
    real dR_dt = 0.2 * I;
    

    return {dS_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> max_days;
  int<lower=1> n_states;
  real t0;
  array[n_states] int t_last;
  array[n_states] int N;              
  array[max_days, n_states] int cases;
}

transformed data {
  array[0] real x_r;
   //vector[n_states] p_reported = rep_vector(0.001, n_states);  // FIXED at 0.001 FOR R
}

parameters {
  vector<lower=0, upper=1>[n_states] pN;      
  vector<lower=0>[n_states] beta;
  real<lower=0, upper=1> p_reported;           // still estimated globally
  vector<lower=0>[n_states] I0;      // initial infected (per-state), real-valued NEW
  //vector<lower=0>[n_states] I0;//
  real<lower=0> phi;// Negative-binomial dispersion
}

transformed parameters {
  array[n_states, 3] real y0; 
  array[n_states, max_days, 3] real y;
  array[n_states, max_days - 1] real incidence;
  real MeanpN;
  real meanp_reported;
  real meanbeta;
  real meanI0; // add mean of I0 NE not needed since estimate
  

  for (s in 1:n_states) {
    y0[s, 1] = pN[s]*N[s] - I0[s]; //NEW for manuscript 1
    y0[s, 2] = I0[s];//b c fixed
    y0[s, 3] = 0;
  }

  for (s in 1:n_states) {
    array[2] real theta = {beta[s], pN[s]};
    array[1] int x_i = {N[s]};
    int tp = t_last[s];
    array[tp] real ts;
    for (i in 1:tp) ts[i] = i;

    array[t_last[s], 3] real y_tmp;
    y_tmp = integrate_ode_rk45(sir, y0[s], t0, ts, theta, x_r, x_i);
    for (i in 1:t_last[s]) y[s, i] = y_tmp[i];

    for (i in 1:(t_last[s] - 1)) {
      incidence[s, i] = (y[s, i, 1] - y[s, i + 1, 1]) * p_reported + 1e-6;
    }
  }

  MeanpN = sum(pN) / n_states;
  meanp_reported = p_reported;
  meanbeta = sum(beta) / n_states;
  meanI0 = sum(I0) / n_states;  //now estimating I0 
}

model {
  beta ~ normal(0.5, 0.1);
  p_reported ~ beta(0.06, 8);

  for (s in 1:n_states) {
    pN[s] ~ beta(4, 4); 
     phi ~ exponential(1);//
  }
  
  
  for (s in 1:n_states) {
     I0[s] ~ normal((1 / p_reported), (1 / p_reported));//SD IS THE MEAN
  }

  for (s in 1:n_states) {
    for (i in 1:(t_last[s] - 1)) {
      //cases[i, s] ~ normal(incidence[s, i],
                           //sqrt(incidence[s, i] * (1 - p_reported)));//epidemic
      cases[i,s] ~ neg_binomial_2( incidence[s,i], phi); //new                    
    }
  }
}

generated quantities {

  real R0;

  real meanR0;

  array[n_states,max_days-1] real pred_cases;

  array[n_states,max_days-1] real log_lik;

  R0 = mean(to_vector(beta))/0.2;

  meanR0 = R0;

  for(s in 1:n_states){

    for(i in 1:(max_days-1)){

      if(i < t_last[s]){

        pred_cases[s,i] =neg_binomial_2_rng(incidence[s,i],phi );

        log_lik[s,i] = neg_binomial_2_lpmf( cases[i,s] |incidence[s,i],phi);

      } else {

        pred_cases[s,i] = -1;

        log_lik[s,i] = 0;

      }

    }

  }
}  
  
