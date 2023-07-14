function [res] = sar_robust_est(y,W,X,tol,c1,c2,c3,max_iter,init)
%PURPOSE:   computes robust M-estimates of the SAR model based on the proposed iterative algorithm.
% ---------------------------------------------------
%  USAGE: res = sar_robust_est(y,W,X,tol,c1,c2,c3,max_iter,init);
%  where:   y = response vector
%           W = standardised spatial weight matrix
%           X = model matrix (with intercept term in first column if used)
%           tol = tolerance parameter for convergence of iterative algorithm
%           c1 = tuning parameter used for eta_{R,beta}
%           c2 = tuning parameter used for eta_{R,sigma}
%           c3 = tuning parameter used for eta_{R,rho}
%           max_iter = maximum number of iterations for the iterative algorithm
%           init = [init_beta', init_sigma, init_rho] where init_beta is a k x 1 vector of initial values for the
%           regression coefficient vector with k denoting the number of covariates (including intercept if used), init_sigma is the initial value for the standard deviation of error terms, 
%           and init_rho is the initial value for the spatial dependence parameter
% ---------------------------------------------------
%  RETURNS: a structure
%           res.beta = robust M-estimate for the regression coefficient vector   (k x 1) vector
%           res.sigma = robust M-estimate for the standard deviation of error terms
%           res.rho = robust M-estimate for the spatial dependence parameter
% --------------------------------------------------
%  NOTES: when the iterative algorithm doesn't converge, we recommend to try alternative initial values. 
%         The returned estimates are set to be NaN when the iterative algorithm doesn't converge.
% --------------------------------------------------  

%N = number of spatial units, k = number of covariates
[N,k] = size(X);

%Compute the smallest and largest eigenvalue for W
min_eig = min(eig(W)); 
max_eig = max(eig(W)); 

%Defining the Huber function with different tuning parameters
huber1=@(x) (x<-c1)*(-c1) + (x>c1)*c1 + x * (x<= c1 & x>= -c1);
huber2=@(x) (x<-c2)*(-c2) + (x>c2)*c2 + x * (x<= c2 & x>= -c2);
huber3=@(x) (x<-c3)*(-c3) + (x>c3)*c3 + x * (x<= c3 & x>= -c3);

%Defining htilde function as in the paper
htilde =@(c) 2*(c^2)*(1-normcdf(c)) - 2*c*normpdf(c)-1+2*normcdf(c);

%Initializing 
old_beta = (init(1:k))';
old_sigma = init(k+1);
old_rho = init(k+2);
beta_diff = 10;
sigma_diff = 10;
rho_diff = 10;
    
iter = 0;
%Iterative Algorithm    
    while beta_diff > tol || sigma_diff > tol || rho_diff > tol
        
        %If the algorithm doesn't converge, then set the resulting
        %estimates as NaN
        if iter > max_iter
            old_beta = NaN([k,1]);
            old_sigma = NaN(1);
            old_rho = NaN(1);
            new_beta = old_beta;
            new_sigma = old_sigma;
            new_rho = old_rho;
            
            break
        end
        
            
        iter = iter +1 ;
        
        %Beta step
        z = ((eye(N)-old_rho*W)*y - X*old_beta)/old_sigma;
        h1z = arrayfun(huber1,z);
        alpha = h1z ./ z;
        new_beta = (X' *diag(alpha) *X)\ (X' *diag(alpha) *(y-old_rho*W*y));
        
        %Sigma step
        z = ((eye(N)-old_rho*W)*y - X*new_beta)/old_sigma;
        h2z = arrayfun(huber2,z);
        new_sigma = ((((old_sigma)^2)/(N*htilde(c2))) * (h2z'*h2z))^(0.5);
           
        %rho step
        fun = @(x) ((1/new_sigma)*(W*inv(eye(N)-x*W)*X*new_beta)'*arrayfun(huber3,((eye(N)-x*W)*y - X*new_beta)/new_sigma) +(arrayfun(huber3,((eye(N)-x*W)*y - X*new_beta)/new_sigma))'*(W*inv(eye(N)-x*W))'*arrayfun(huber3,((eye(N)-x*W)*y - X*new_beta)/new_sigma)-trace(inv(eye(N)-x*W)*W)*htilde(c3))^2;
        new_rho = fminbnd(fun,(1/min_eig),(1/max_eig));
    
        %Compute the differences between current estimate and previous
        %loop's estimates to check for convergence
        beta_diff = max(abs(new_beta-old_beta));
        sigma_diff = abs(new_sigma - old_sigma);
        rho_diff = abs(new_rho - old_rho);
        
        old_beta = new_beta;
        old_sigma = new_sigma;
        old_rho = new_rho;
    end
    
    
    res.beta = new_beta;
    res.sigma = new_sigma;
    res.rho = new_rho;

    