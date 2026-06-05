functions {
  real db_lambert(real x_b) {return 10 * log10(x_b * pi() / (10 ^ -6));}
  real db_spl(real x_l) {return 20 * log10(x_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real bright_to_loud(real x_bright, real alpha_l, real beta_l,
                        real beta_b, real w_1, real const_bl) {
    return db_spl(
       pow(w_1 * 1 / alpha_l * pow(db_inv_lambert(x_bright), beta_b) +
           inv(alpha_l) * const_bl, inv(beta_l))
    );
  }
}
data {
  int<lower=1>                     ntotal;
  int<lower=1>                      nx_bl;  // number of standards
  int<lower=1>                    nranges;
  vector<lower=1>[ntotal]            x_bl;
  vector<lower=1>[ntotal]            y_bl;
  array[ntotal] int<lower=1>     x_bl_idx;
  array[ntotal] int<lower=1>    range_idx;
  vector<lower=0>[ntotal]          sig_bl;
  real<lower=0>                     rho_b;  // fix light reference intensity
}
parameters {
  real<lower=0>          alpha_l;  // alpha_b = 1
  real<lower=0>           beta_l;
  real<lower=0>           beta_b;
  real<lower=0>              w_1;
  vector<lower=0>[nranges] rho_l;  // sound reference intensity
}
transformed parameters {
  vector[nranges] const_bl;
  const_bl = -w_1 * pow(db_inv_lambert(rho_b), beta_b) + 
                                alpha_l .* pow(rho_l, beta_l);
  vector[nx_bl] mu_bl;
  for (i in 1:ntotal) {
    mu_bl[x_bl_idx[i]] = bright_to_loud(x_bl[i], alpha_l, beta_l, beta_b, w_1, 
                                        const_bl[range_idx[i]]);
  }
}
model {
  target += normal_lpdf(w_1 | 1, .1) - normal_lccdf(0 | 1, .1);
  target += normal_lpdf(rho_l | 50, 15) - normal_lccdf(0 | 50, 15);
  target += normal_lpdf(alpha_l | 5, 3) - normal_lccdf(0 | 5, 3);
  target += beta_lpdf(beta_l | 6, 3);
  target += beta_lpdf(beta_b | 3, 6);
  target += normal_lpdf(y_bl | mu_bl[x_bl_idx], sig_bl);
}
generated quantities {
  vector[nranges] rho_l_db;
  for (r in 1:nranges) rho_l_db[r] = db_spl(rho_l[r]);
}

