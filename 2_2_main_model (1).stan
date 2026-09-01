
data {
  int<lower=1> T;                               // step number
  int<lower=1> n_step;                          // step for rel_sus, rel_inf
  int<lower=1> n_age;                           // age group
  array[T] int<lower=1, upper=n_step> step_id;  // t transport to step
  matrix[n_age, n_age] A;                       // normalized matirx
  array[T, n_age] int M;  
  matrix[T, n_age] Y;                           // observed infections = report_case / reporting_coverage
  matrix[T, n_age] s_mat;                       // suceptible population (daily by age、1 - cumulative infection）
  vector<lower=0>[n_age] n_pop;                 // population for age group
  vector<lower=0, upper=1>[T] reportcoverage;
  int<lower=1> g_len;                            // length of serial interval 
  vector[g_len] g_vec;                           // normalized discreted serial interval distribution
  int<lower=1, upper=T> fix_end;                 // last date when using observed infections
  
  //mask effect with two pathways　
  array[T] matrix[n_age, n_age]  mask_sus_mat; //suceptibility  (for wearer protection)
  array[T] matrix[n_age, n_age]  mask_inf_mat; //infectiousness  (for source contorol)
  array[T] matrix[n_age, n_age] mask_sus_mat_CF1; // suceptibility for high mask scenario
  array[T] matrix[n_age, n_age] mask_inf_mat_CF1; // infectiousness for high mask scenario
  array[T] matrix[n_age, n_age] mask_sus_mat_CF2; // suceptibility for low mask scenario
  array[T] matrix[n_age, n_age] mask_inf_mat_CF2; // infectiousness for low mask scenario
  
}

transformed data {
  array[T] row_vector[n_age] J_data;
  for (t in 1:T) {
    row_vector[n_age] tmp = rep_row_vector(0.0, n_age);
    int max_lag = (t - 1 < g_len) ? t - 1 : g_len;
    for (tau in 1:max_lag)
      tmp += Y[t - tau] * g_vec[tau];
    J_data[t] = tmp;
  }
}

parameters {
  matrix[n_step, n_age - 1] log_rel_sus; 
  matrix[n_step, n_age - 1] log_rel_inf;
  vector[T] log_Ct; 
}

transformed parameters {
  matrix[n_step, n_age - 1] rel_sus = exp(log_rel_sus);
  matrix[n_step, n_age - 1] rel_inf = exp(log_rel_inf);
  vector[T] Ct = exp(log_Ct);

  // ----convolution and susceptible propotion : using observed infections----
  array[T] matrix[n_age, n_age] K_out;
  for (t in 1:T) {
    int s = step_id[t];
    vector[n_age] rs = [rel_sus[s,1], rel_sus[s,2], rel_sus[s,3], 1.0]';
    vector[n_age] ri = [rel_inf[s,1], rel_inf[s,2], rel_inf[s,3], 1.0]';
    matrix[n_age, n_age] K = A;
    for (i in 1:n_age) K[i,] *= rs[i];
    for (j in 1:n_age) K[,j] *= ri[j];
    for (i in 1:n_age) K[i,] *= mask_sus_mat[t][i,i];
    for (j in 1:n_age) K[,j] *= mask_inf_mat[t][j,j];
    K_out[t] = K;
  }

  // ---- expected infecgtions and suceptible proportion calculating from (t-1)----
  array[T] vector[n_age] i_hat_tp;
  array[T] vector[n_age] M_hat_tp;

  {
    // ---- t = 1:fix_end using observed infections----
    for (t in 1:T) {
      vector[n_age] s_t = s_mat[t]';
      i_hat_tp[t] = diag_matrix(s_t) * (Ct[t] * K_out[t]) * J_data[t]';
      M_hat_tp[t] = i_hat_tp[t] * reportcoverage[t];
    }
  }
}

model {
  // ---- prior distribution ----
  to_vector(log_rel_sus) ~ normal(0, 1);
  to_vector(log_rel_inf) ~ normal(0, 1);
  log_Ct               ~ normal(0, 1);

  // ---- likelihood ----
  for (t in 1:T) {
    target += poisson_lpmf(M[t] | fmax(M_hat_tp[t], 1e-10));
  }
}


generated quantities {
  matrix[T, n_age] J;
  array[T] matrix[n_age, n_age] K_hat;
  array[T] vector[n_age] i_hat;       
  array[T] vector[n_age] s_hat;
  array[T] real Rt;

  array[T] real Ct_out;
  array[T] matrix[n_age, n_age] K_hat_out;
  array[T] matrix[n_age, n_age] K_hat_out_CF1;
  array[T] matrix[n_age, n_age] K_hat_out_CF2;
  array[T] vector[n_age] s_hat_out;
  array[T] vector[n_age] i_hat_out;  // true infections
  array[T] vector[n_age] M_hat_out;  // scaling to reported cases
  array[T] vector[n_age] j_hat;
  array[T] vector[n_age] y_rep;       // rmd scaling to relorted cases
  array[T] vector[n_age] M_hat_CF1_out;
  array[T] vector[n_age] M_hat_CF2_out;
  array[T] vector[n_age] y_rep_CF1;
  array[T] vector[n_age] y_rep_CF2;
  array[T] row_vector[n_age] J_out; 
  
  // --- couterfactual ---
  array[T] vector[n_age] s_hat_CF1;
  array[T] vector[n_age] s_hat_CF2;
  array[T] vector[n_age] i_hat_CF1;
  array[T] vector[n_age] i_hat_CF2;
  array[T] real Rt_CF1;
  array[T] real Rt_CF2;
  
  vector[n_age] cum_i_hat = rep_vector(0.0, n_age);
  vector[n_age] cum_i_hat_CF1 = rep_vector(0.0, n_age);
  vector[n_age] cum_i_hat_CF2 = rep_vector(0.0, n_age);
   
  
  // ---- initialize----
  for (t in 1:T) {
    i_hat[t]      = rep_vector(0, n_age);
    i_hat_CF1[t]  = rep_vector(0, n_age);
    i_hat_CF2[t]  = rep_vector(0, n_age);
    s_hat[t]      = rep_vector(0, n_age);
    s_hat_CF1[t]  = rep_vector(0, n_age);
    s_hat_CF2[t]  = rep_vector(0, n_age);
    j_hat[t]      = rep_vector(0, n_age);
    M_hat_out[t]  = rep_vector(0, n_age);
    J[t] = J_data[t];
    J_out[t] = J_data[t];
  }

  // ---- t=1:fix_end  : using observed infections ----
  for (t in 1:fix_end) {
    int s = step_id[t];
    
    // right：column scaling、left：row scaling
    vector[n_age] rel_sus_full = [ rel_sus[s,1], rel_sus[s,2], rel_sus[s,3], 1.0 ]';//suceptiblity maltiply to row
    vector[n_age] rel_inf_full = [ rel_inf[s,1], rel_inf[s,2], rel_inf[s,3], 1.0 ]';//infectiousness maltiply to column
    matrix[n_age, n_age] K = A;
    for (i in 1:n_age) K[i,] *= rel_sus_full[i];
    for (j in 1:n_age) K[,j] *= rel_inf_full[j];
    K_hat[t] = K;
    
    // mask effect by age
    K_hat_out[t] = K_hat[t];
    for (i in 1:n_age)
    K_hat_out[t][i,] = K_hat_out[t][i,] * mask_sus_mat[t][i,i];  
    for (j in 1:n_age)
    K_hat_out[t][,j] = K_hat_out[t][,j] * mask_inf_mat[t][j,j]; 

    
    // ---- t=1〜fix_end ----
    //if (t <= 2) {
      s_hat[t]      = s_mat[t]';
      s_hat_out[t]  = s_hat[t];
      s_hat_CF1[t]  = s_hat[t]; 
      s_hat_CF2[t]  = s_hat[t];
    
      //convolute to Y -> transformed data
    
      //actual
      i_hat[t]     = diag_matrix(s_hat[t]) * (Ct[t] * K_hat_out[t]) * J_data[t]';
       
      // counterfactual 
      K_hat_out_CF1[t] = mask_sus_mat_CF1[t] * K_hat[t] * mask_inf_mat_CF1[t];
      K_hat_out_CF2[t] = mask_sus_mat_CF2[t] * K_hat[t] * mask_inf_mat_CF2[t];
      i_hat_CF1[t] = diag_matrix(s_hat_CF1[t]) * (Ct[t] * K_hat_out_CF1[t]) * J_data[t]';
      i_hat_CF2[t] = diag_matrix(s_hat_CF2[t]) * (Ct[t] * K_hat_out_CF2[t]) * J_data[t]';
        
      // ----mu and poisson_rng ----
      vector[n_age] mu     = i_hat[t]     * reportcoverage[t];
      vector[n_age] mu_cf1 = i_hat_CF1[t] * reportcoverage[t];
      vector[n_age] mu_cf2 = i_hat_CF2[t] * reportcoverage[t];
       
      M_hat_out[t]     = mu;
      M_hat_CF1_out[t] = mu_cf1;
      M_hat_CF2_out[t] = mu_cf2;
       
      for (a in 1:n_age) {
         y_rep[t][a]     = poisson_rng(fmax(mu[a],     1e-10));
         y_rep_CF1[t][a] = poisson_rng(fmax(mu_cf1[a], 1e-10));
         y_rep_CF2[t][a] = poisson_rng(fmax(mu_cf2[a], 1e-10));
      }

    i_hat_out[t] = i_hat[t];
    j_hat[t] = J_data[t]';
    Ct_out[t] = Ct[t];

    // Rt 
    matrix[n_age, n_age] NG_matrix = (diag_matrix(s_hat[t]) * Ct[t] * K_hat_out[t]);
    Rt[t] = max(get_real(eigenvalues(NG_matrix)));

    matrix[n_age, n_age] NG_matrix_CF1 = diag_matrix(s_hat_CF1[t]) * (Ct[t]* K_hat_out_CF1[t]);
    matrix[n_age, n_age] NG_matrix_CF2 = diag_matrix(s_hat_CF2[t]) * (Ct[t] * K_hat_out_CF2[t]);
    Rt_CF1[t] = max(get_real(eigenvalues(NG_matrix_CF1)));
    Rt_CF2[t] = max(get_real(eigenvalues(NG_matrix_CF2)));
  }
  
  // ---- t=fix_end+1: using obeserved infection for last time  ----
if (fix_end < T) {
  
  int t = fix_end + 1;
  int s = step_id[t];

    vector[n_age] rel_sus_full = [ rel_sus[s,1], rel_sus[s,2], rel_sus[s,3], 1.0 ]';
    vector[n_age] rel_inf_full = [ rel_inf[s,1], rel_inf[s,2], rel_inf[s,3], 1.0 ]';

    K_hat[t] = A;
    for (i in 1:n_age) K_hat[t][i,] *= rel_sus_full[i];
    for (j in 1:n_age) K_hat[t][,j] *= rel_inf_full[j];

    K_hat_out[t] = K_hat[t];
    for (i in 1:n_age) K_hat_out[t][i,] *= mask_sus_mat[t][i,i];
    for (j in 1:n_age) K_hat_out[t][,j] *= mask_inf_mat[t][j,j];

    // counterfactual
    K_hat_out_CF1[t] = K_hat[t];
    for (i in 1:n_age) K_hat_out_CF1[t][i,] *= mask_sus_mat_CF1[t][i,i];
    for (j in 1:n_age) K_hat_out_CF1[t][,j] *= mask_inf_mat_CF1[t][j,j];

    K_hat_out_CF2[t] = K_hat[t];
    for (i in 1:n_age) K_hat_out_CF2[t][i,] *= mask_sus_mat_CF2[t][i,i];
    for (j in 1:n_age) K_hat_out_CF2[t][,j] *= mask_inf_mat_CF2[t][j,j];

    s_hat[t]     = s_mat[t]';
    s_hat_CF1[t] = s_hat[t];
    s_hat_CF2[t] = s_hat[t];
    
    s_hat_out[t] = s_hat[t]; 

    j_hat[t] = J_data[t]';               
    i_hat[t]     = diag_matrix(s_hat[t])     * (Ct[t] * K_hat_out[t])     * j_hat[t];
    i_hat_CF1[t] = diag_matrix(s_hat_CF1[t]) * (Ct[t] * K_hat_out_CF1[t]) * j_hat[t];
    i_hat_CF2[t] = diag_matrix(s_hat_CF2[t]) * (Ct[t] * K_hat_out_CF2[t]) * j_hat[t];
    
    i_hat_out[t] = i_hat[t];
    
    for (tt in 1:(fix_end + 1)) {
      cum_i_hat     += i_hat[tt];
      cum_i_hat_CF1 += i_hat_CF1[tt];
      cum_i_hat_CF2 += i_hat_CF2[tt];
      }
  
    M_hat_out[t]     = i_hat[t]     * reportcoverage[t];
    M_hat_CF1_out[t] = i_hat_CF1[t] * reportcoverage[t];
    M_hat_CF2_out[t] = i_hat_CF2[t] * reportcoverage[t];
    
    for (a in 1:n_age) {
      y_rep[t][a]     = poisson_rng(fmax(M_hat_out[t][a],     1e-10));
      y_rep_CF1[t][a] = poisson_rng(fmax(M_hat_CF1_out[t][a], 1e-10));
      y_rep_CF2[t][a] = poisson_rng(fmax(M_hat_CF2_out[t][a], 1e-10));
      }

    Ct_out[t] = Ct[t];

    matrix[n_age, n_age] NG = diag_matrix(s_hat[t]) * (Ct[t] * K_hat_out[t]);
    Rt[t] = max(get_real(eigenvalues(NG)));

    matrix[n_age, n_age] NG1 = diag_matrix(s_hat_CF1[t]) * (Ct[t] * K_hat_out_CF1[t]);
    matrix[n_age, n_age] NG2 = diag_matrix(s_hat_CF2[t]) * (Ct[t] * K_hat_out_CF2[t]);
    Rt_CF1[t] = max(get_real(eigenvalues(NG1)));
    Rt_CF2[t] = max(get_real(eigenvalues(NG2)));
  }

// ---- t=fix_end+2~: convolute with estimated ----
if (fix_end + 1 < T) {
   
  for (t in (fix_end + 2):T) {
  int max_lag = (t - 1 < g_len) ? t - 1 : g_len;
  j_hat[t] = rep_vector(0, n_age);
  vector[n_age] j_hat_CF1 = rep_vector(0, n_age);
  vector[n_age] j_hat_CF2 = rep_vector(0, n_age);

  for (tau in 1:max_lag) {
    j_hat[t]    += i_hat[t - tau]    * g_vec[tau];
    j_hat_CF1   += i_hat_CF1[t - tau] * g_vec[tau];
    j_hat_CF2   += i_hat_CF2[t - tau] * g_vec[tau];
  }
    int s = step_id[t];
 
    // A 
    vector[n_age] rel_sus_full = [ rel_sus[s,1], rel_sus[s,2], rel_sus[s,3], 1.0 ]';
    vector[n_age] rel_inf_full = [ rel_inf[s,1], rel_inf[s,2], rel_inf[s,3], 1.0 ]';

    K_hat[t] = A;
    for (i in 1:n_age) K_hat[t][i,] *= rel_sus_full[i];
    for (j in 1:n_age) K_hat[t][,j] *= rel_inf_full[j];

    K_hat_out[t] = K_hat[t];
    for (i in 1:n_age) K_hat_out[t][i,] *= mask_sus_mat[t][i,i];
    for (j in 1:n_age) K_hat_out[t][,j] *= mask_inf_mat[t][j,j];

    // counter factual
    K_hat_out_CF1[t] = K_hat[t];
    for (i in 1:n_age) K_hat_out_CF1[t][i,] *= mask_sus_mat_CF1[t][i,i];
    for (j in 1:n_age) K_hat_out_CF1[t][,j] *= mask_inf_mat_CF1[t][j,j];

    K_hat_out_CF2[t] = K_hat[t];
    for (i in 1:n_age) K_hat_out_CF2[t][i,] *= mask_sus_mat_CF2[t][i,i];
    for (j in 1:n_age) K_hat_out_CF2[t][,j] *= mask_inf_mat_CF2[t][j,j];
    
    s_hat[t] = fmax(0, 1 - cum_i_hat ./ n_pop);
    s_hat_out[t] = s_hat[t];
    s_hat_CF1[t]  = fmax(0, 1 - cum_i_hat_CF1 ./ n_pop);
    s_hat_CF2[t]  = fmax(0, 1 - cum_i_hat_CF2 ./ n_pop);

    // i_hat → mu → M_hat_out / y_rep (calculate prediction intervals)
    i_hat[t]     = diag_matrix(s_hat[t])     * (Ct[t]  * K_hat_out[t]) * j_hat[t];
    i_hat_CF1[t] = diag_matrix(s_hat_CF1[t]) * (Ct[t] * K_hat_out_CF1[t]) * j_hat_CF1;
    i_hat_CF2[t] = diag_matrix(s_hat_CF2[t]) * (Ct[t] * K_hat_out_CF2[t]) * j_hat_CF2;
    
    cum_i_hat     += i_hat[t];
    cum_i_hat_CF1 += i_hat_CF1[t];
    cum_i_hat_CF2 += i_hat_CF2[t];
    
    vector[n_age] mu     = i_hat[t]     * reportcoverage[t];
    vector[n_age] mu_cf1 = i_hat_CF1[t] * reportcoverage[t];
    vector[n_age] mu_cf2 = i_hat_CF2[t] * reportcoverage[t];
    
    M_hat_out[t]     = mu;
    M_hat_CF1_out[t] = mu_cf1;
    M_hat_CF2_out[t] = mu_cf2;
    
    for (a in 1:n_age) {
      y_rep[t][a]     = poisson_rng(fmax(mu[a],     1e-10));
      y_rep_CF1[t][a] = poisson_rng(fmax(mu_cf1[a], 1e-10));
      y_rep_CF2[t][a] = poisson_rng(fmax(mu_cf2[a], 1e-10));
      }

    i_hat_out[t] = i_hat[t];
    Ct_out[t] = Ct[t];

    // Rt
    matrix[n_age, n_age] NG_matrix     = diag_matrix(s_hat[t])     * (Ct[t]  * K_hat_out[t]);
    matrix[n_age, n_age] NG_matrix_CF1 = diag_matrix(s_hat_CF1[t]) * (Ct[t] * K_hat_out_CF1[t]);
    matrix[n_age, n_age] NG_matrix_CF2 = diag_matrix(s_hat_CF2[t]) * (Ct[t] * K_hat_out_CF2[t]);

    Rt[t] = max(get_real(eigenvalues(NG_matrix)));
    Rt_CF1[t] = max(get_real(eigenvalues(NG_matrix_CF1)));
    Rt_CF2[t] = max(get_real(eigenvalues(NG_matrix_CF2)));
    
  }
  
}

}


