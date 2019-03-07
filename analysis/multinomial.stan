data{
    int<lower=0> N; // number of observations
    int<lower=0> Y[N, 3]; // vector of multinomial counts
    int<lower=0> K; // number of covariates
    int<lower=0> X[N, K]; // covariate matrix
    int<lower=0> C; // Number of communities
    int<lower=0> cidx[N]; // community index for each observation
}

parameters{
    matrix[K,2] beta; // betas for the two alternative options
    matrix[C,2] zeta; // random effects on community
    vector[2] sigma; // sd of random effects
    real<lower=-1, upper=1> rho; // 
}

transformed parameters{
    matrix[N, 3] pstar;
    vector[2] mu;
    matrix[2,2] Sigma;
    
    Sigma[1,1] = sigma[1]^2;
    Sigma[2,2] = sigma[2]^2;
    Sigma[1,2] = sigma[1]*sigma[2]*rho;
    Sigma[2,1] = sigma[1]*sigma[2]*rho;
    
    mu[1] = 0;
    mu[2] = 0;

    for(n in 1:N){
        pstar[n,1] = 0;
        pstar[n,2] = 0;
        pstar[n,3] = 0;
        for(j in 2:3){
            pstar[n,j] = zeta[cidx[n], j-1];
            for(k in 1:K){
                pstar[n,j] += X[n,k] * beta[k,j-1];
            }
        }
    }
}

model{
    // priors
    for(j in 1:2){
        for(k in 1:K){
            beta[k, j] ~ normal(0, 10);
        }
        sigma[j] ~ gamma(2,0.1);
    }
    
    for(c in 1:C){
        to_vector(zeta[c,]) ~ multi_normal(mu, Sigma);
    }
    
    (rho+1)/2 ~ beta(1,1);

    
    // data likelihood
    for(n in 1:N){
        Y[n,] ~ multinomial(softmax(to_vector(pstar[n,]))); // 1st is always 0
    }
}
