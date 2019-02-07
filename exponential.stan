functions{
    //returns the reciprocal of argument  
    real get_reciprocal(real a){
        	real b;
        	b = 1/a;
        	return b;   
    }
     //log likelihood function
    real dexp_lpdf(vector x, real lambda){
        vector[num_elements(x)] temp; //array that stores densities
        real log_prob;
        for(i in 1:num_elements(x)){  //loop over data and store each density
            temp[i] = lambda*exp(-lambda*x[i]); 
        }
        log_prob = sum(log(temp));   //return sum of log probabilities
        return log_prob;
    }
}


data{
 int LENGTH;
 vector[LENGTH] Y;
}

parameters{
   real<lower=0> lambda;
}


transformed parameters{
    real x_trans;
    x_trans = get_reciprocal(lambda);
}


model{
     real alpha;
     real beta;
     alpha = 1;
     beta = 1;
    lambda ~ gamma(alpha,beta);
    Y ~ dexp(lambda);
}

generated quantities{
    real pred;
    pred = exponential_rng(lambda);
}

