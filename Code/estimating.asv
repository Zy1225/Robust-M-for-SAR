function F = estimating(param,y,W,X,c1,c2,c3)
%PURPOSE:   compute the estimating functions of the proposed robust M-estimator evaluated at given values of parameters.
% ---------------------------------------------------
%  USAGE: F = estimating(param,y,W,X,c1,c2,c3)
%  where:   param = [beta', sigma, rho] where beta is a k x 1 vector of regression coefficients, sigma is the standard deviation of error terms, and rho is the spatial dependence parameter 
%           y = response vector
%           W = standardised spatial weight matrix
%           X = model matrix (with intercept term in first column if used)
%           c1 = tuning parameter used for eta_{R,beta}
%           c2 = tuning parameter used for eta_{R,sigma}
%           c3 = tuning parameter used for eta_{R,rho}
% ---------------------------------------------------
%  RETURNS: a (k+2) x 1 vector denoted as [eta_{R,beta}; eta_{R,sigma}; eta_{R,rho}] consisting of the values of the estimating functions for
%  the regression coefficients (eta_{R,beta}), standard deviation of the error terms (eta_{R,sigma}), and spatial dependence parameter (eta_{R,rho}).
% --------------------------------------------------

%N = number of spatial units, k = number of covariates
[N,k] = size(X);

%Extract parameters from param vector
beta = (param(1:k))';
sigma = param(k+1);
rho = param(k+2);

G = W* inv(eye(N)-rho*W);

%Defining the Huber function with different tuning parameters
huber1=@(x) (x<-c1)*(-c1) + (x>c1)*c1 + x * (x<= c1 & x>= -c1);
huber2=@(x) (x<-c2)*(-c2) + (x>c2)*c2 + x * (x<= c2 & x>= -c2);
huber3=@(x) (x<-c3)*(-c3) + (x>c3)*c3 + x * (x<= c3 & x>= -c3);

%Defining htilde function as in the paper
htilde =@(c) 2*(c^2)*(1-normcdf(c)) - 2*c*normpdf(c)-1+2*normcdf(c);

F = [ X' * arrayfun(huber1,((eye(N)-rho*W)*y - X*beta)/sigma);
    (arrayfun(huber2,((eye(N)-rho*W)*y - X*beta)/sigma))'* (arrayfun(huber2,((eye(N)-rho*W)*y - X*beta)/sigma)) - N*htilde(c2);
    (1/sigma)*(G*X*beta)' *(arrayfun(huber3,((eye(N)-rho*W)*y - X*beta)/sigma)) + (arrayfun(huber3,((eye(N)-rho*W)*y - X*beta)/sigma))' * G' *(arrayfun(huber3,((eye(N)-rho*W)*y - X*beta)/sigma)) - (trace(G))*htilde(c3)  ];