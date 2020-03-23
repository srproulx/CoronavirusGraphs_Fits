
    data {
    int total_datapoints ; // how much data is there?
    int<lower=0> maxdays[2]; // number of time points
    int total_cases[total_datapoints] ; // observed cases
    int days[total_datapoints] ; //day index within each dataset
    int dataset[total_datapoints] ; //index for which dataset the data are from
    real pstrength ; // strength of the prior
    }
    transformed data{
    int new_cases[total_datapoints] ; //really only need n-1 but keep the indexing the same for simplicity
    new_cases[1]=0;
    for(i in 2:maxdays[1]){
    new_cases[i]=total_cases[i]-total_cases[i-1];
    }
    new_cases[(maxdays[1]+1)]=0;
    for(i in (maxdays[1]+2):(maxdays[1]+maxdays[2])){
    new_cases[i]=total_cases[i]-total_cases[i-1];
    }
    }
    parameters {
    real <lower=0, upper=1> mean_lambda ; // the midpoint of the lambdas of the two datasets, weak uniform prior
    real lambda_diff ; // the log difference of the two lambdas
    }
    transformed parameters{
    real lambda[2]; // the two lambdas are given by log scale transformations using the parameter labmda_diff
    lambda[1] = mean_lambda * exp(-lambda_diff);
    lambda[2] = mean_lambda * exp(lambda_diff);
    }
    model {
    lambda_diff ~ normal(0,pstrength) ; // fairly weak prior on the difference between the two samples
    for(i in 2:maxdays[1]){
    new_cases[i] ~ poisson(total_cases[i-1]* lambda[1] ); // update for each time step is the sum of Poisson RVs with mean lambda
    } // update the likelihood for the first dataset
    for(i in (maxdays[1]+2):(maxdays[1]+maxdays[2])){
    new_cases[i] ~ poisson(total_cases[i-1]*lambda[2]); // update for each time step is the sum of Poisson RVs with mean lambda
    } // update for the second dataset
    }
    
