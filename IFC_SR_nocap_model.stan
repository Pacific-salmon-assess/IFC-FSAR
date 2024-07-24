// Model with no prior on carrying capacity/srep
data {
  int<lower=0> N;           //number of years
  int<lower=0> C;           //number of CUs
  matrix[C,N] p_3;          //proportion of age 3 recruits
  matrix[C,N] p_4;          //proportion of age 3 recruits
  matrix[C,N] surv_3;       //logit SAS for age 3 recruits
  matrix[C,N] surv_4;       //logit SAS for age 4 recruits
  matrix[C,N] lrs;          //log recruits per spawner
  matrix[C,N] S;            //log recruits per spawner
  real mlogit_surv3;         //mean logit-transformed smolt-to-adult survival for age 3
  real mlogit_surv4;         //mean logit-transformed smolt-to-adult survival for age 4
}

parameters {
  vector[C] beta;           //population-specific Ricker beta parameter
  vector<lower=0>[C] alpha;          //population-specific Ricker alpha parameter
  vector[C] gamma;                   //population-specific survival index parameter
  vector<lower=0>[C] sigma;          //population-specific SD within the autocorrelated process
}

model {
    for(i in 1:N){
      lrs[1,i] ~ normal(p_3[1,i] * (beta[1] * S[1,i] + gamma[1] * surv_3[1,i] + alpha[1]) + p_4[1,i] * (beta[1] * S[1,i] + gamma[1] * surv_4[1,i] + alpha[1]), sigma[1]);
      lrs[2,i] ~ normal(p_3[2,i] * (beta[2] * S[2,i] + gamma[2] * surv_3[2,i] + alpha[2]) + p_4[2,i] * (beta[2] * S[2,i] + gamma[2] * surv_4[2,i] + alpha[2]), sigma[2]);
      lrs[3,i] ~ normal(p_3[3,i] * (beta[3] * S[3,i] + gamma[3] * surv_3[3,i] + alpha[3]) + p_4[3,i] * (beta[3] * S[3,i] + gamma[3] * surv_4[3,i] + alpha[3]), sigma[3]);
      lrs[4,i] ~ normal(p_3[4,i] * (beta[4] * S[4,i] + gamma[4] * surv_3[4,i] + alpha[4]) + p_4[4,i] * (beta[4] * S[4,i] + gamma[4] * surv_4[4,i] + alpha[4]), sigma[4]);
      lrs[5,i] ~ normal(p_3[5,i] * (beta[5] * S[5,i] + gamma[5] * surv_3[5,i] + alpha[5]) + p_4[5,i] * (beta[5] * S[5,i] + gamma[5] * surv_4[5,i] + alpha[5]), sigma[5]);
}
  beta ~ normal(0, 10);
  alpha ~ cauchy(0, 5);
  gamma ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
}

generated quantities {
  matrix[C,N] nu_Y;
  matrix[C,N] nu_rec;
  vector[C] srep;
  vector[C] smsy_80;
  vector[C] smsy;
  vector[C] umsy;
    for(i in 1:N){
      nu_Y[1,i] = normal_rng(p_3[1,i] * (beta[1] * S[1,i] + gamma[1] * surv_3[1,i] + alpha[1]) + p_4[1,i] * (beta[1] * S[1,i] + gamma[1] * surv_4[1,i] + alpha[1]), sigma[1]);
      nu_Y[2,i] = normal_rng(p_3[2,i] * (beta[2] * S[2,i] + gamma[2] * surv_3[2,i] + alpha[2]) + p_4[2,i] * (beta[2] * S[2,i] + gamma[2] * surv_4[2,i] + alpha[2]), sigma[2]);
      nu_Y[3,i] = normal_rng(p_3[3,i] * (beta[3] * S[3,i] + gamma[3] * surv_3[3,i] + alpha[3]) + p_4[3,i] * (beta[3] * S[3,i] + gamma[3] * surv_4[3,i] + alpha[3]), sigma[3]);
      nu_Y[4,i] = normal_rng(p_3[4,i] * (beta[4] * S[4,i] + gamma[4] * surv_3[4,i] + alpha[4]) + p_4[4,i] * (beta[4] * S[4,i] + gamma[4] * surv_4[4,i] + alpha[4]), sigma[4]);
      nu_Y[5,i] = normal_rng(p_3[5,i] * (beta[5] * S[5,i] + gamma[5] * surv_3[5,i] + alpha[5]) + p_4[5,i] * (beta[5] * S[5,i] + gamma[5] * surv_4[5,i] + alpha[5]), sigma[5]);
      nu_rec[1,i] = exp(nu_Y[1,i])*S[1,i]*1000;
      nu_rec[2,i] = exp(nu_Y[2,i])*S[2,i]*1000;
      nu_rec[3,i] = exp(nu_Y[3,i])*S[3,i]*1000;
      nu_rec[4,i] = exp(nu_Y[4,i])*S[4,i]*1000;
      nu_rec[5,i] = exp(nu_Y[5,i])*S[5,i]*1000;
  }
  for(c in 1:C){
    srep[c] = (alpha[c] + gamma[c]*mlogit_surv3)/(-1*beta[c]) * 1000;
    smsy_80[c] = 0.8*((1 - lambert_w0(exp(1 - (alpha[c] + gamma[c]*mlogit_surv3)))) / (-1*beta[c]))*1000;
    smsy[c] = ((1 - lambert_w0(exp(1 - (alpha[c] + gamma[c]*mlogit_surv3)))) / (-1*beta[c]))*1000;
    umsy[c] = 1 - lambert_w0(exp(1 - (alpha[c] + gamma[c]*mlogit_surv3)));
  }
}

