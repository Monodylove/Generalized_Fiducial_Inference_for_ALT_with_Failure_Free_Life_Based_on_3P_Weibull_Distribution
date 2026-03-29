functions {
  real phi(real x) {
    return 1 / x; // Arrhenius函数
  }
  
  real term1(real t, real S, real delta, real b) {
    return t - exp(delta + b * phi(S));
  }
  
  real term2(real t, real S, real lambda, real b) {
    return exp(lambda + b * phi(S));
  }
  
  real f(real t, real S, real m, real delta, real lambda, real b) {
    real temp1 = term1(t, S, delta, b);
    real temp2 = term2(t, S, lambda, b);
    
    if (temp1 >= 0) {
      return m * pow(temp1, m - 1) * exp(-pow(temp1 / temp2, m)) / pow(temp2, m);
    } else {
      return 0.0;  // 当 temp1 < 0 时返回0
    }
  }
  
  real R(real t, real S, real m, real delta, real lambda, real b) {
    real temp1 = term1(t, S, delta, b);
    real temp2 = term2(t, S, lambda, b);
    
    if (temp1 >= 0) {
      return exp(-pow(temp1 / temp2, m));
    } else {
      return 1.0;  // 当 temp1 < 0 时返回1
    }
  }
  
  real L(matrix SAMPLE, real m, real delta, real lambda, real b) {
    int n_total = rows(SAMPLE);
    real log_likelihood = 0;
    
    for (i in 1:n_total) {
      real t = SAMPLE[i, 1];
      real status = SAMPLE[i, 2];
      real S = SAMPLE[i, 3];
      if (status == 1) {
        log_likelihood += log(f(t, S, m, delta, lambda, b));
      } else {
        log_likelihood += log(R(t, S, m, delta, lambda, b));
      }
    }
    
    return log_likelihood;
  }
}

data {
  int<lower=1> n; // TS的行数
  matrix[n, 2] TS; // 输入矩阵，包含T和S
  int<lower=1> n_sample; // SAMPLE的行数
  matrix[n_sample, 3] SAMPLE; // 输入的SAMPLE矩阵，包含T, status, S
}

parameters {
  real<lower=0> m; // 参数m
  real delta; // 参数delta
  real lambda; // 参数lambda
  real b; // 参数b                
}

model { 
  real log_likelihood = L(SAMPLE, m, delta, lambda, b);
  target += log_likelihood; 
}

