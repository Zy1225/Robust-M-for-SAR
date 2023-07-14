function out = mle_var(beta,sigma,rho,X,W)
%PURPOSE:   computes the asymptotic covariance matrix estimate for MLE by transforming the expression in Lee (2004), as we are estimating
%           (beta',sigma,rho) instead of (beta,sigma^2,rho).
% ---------------------------------------------------
%  USAGE: out = mle_var_sigma(beta,sigma,rho,X,W)
%  where:   beta = MLE of the regression coefficient vector
%           sigma = MLE of the standard deviation of error term
%           rho = MLE of the spatial dependence parameter
%           X = model matrix (with intercept term in first column if used)
%           W = standardised spatial weight matrix
% ---------------------------------------------------
%  RETURNS: the asymptotic covariance matrix estimate for MLE of the SAR model.
% --------------------------------------------------

[N,k] = size(X);
G = W/(eye(N)-rho*W);

%Sigma_{beta,beta}, Sigma_{beta,sigma}, Sigma_{beta_rho}
beta_beta = (1/N) * (1/ sigma^2) * X' * X;
beta_sigma = zeros(k,1);
beta_rho =  (1/N) * (1/ sigma^2) * X' * G * X * beta;   

%Sigma_{sigma,beta}, Sigma_{sigma,sigma}, Sigma_{sigma_rho}
sigma_beta = zeros(1,k);
sigma_sigma = 2/(sigma^2);
sigma_rho = (2/N) * (1/ sigma) * trace(G);

%Sigma_{rho,beta}, Sigma_{rho,sigma}, Sigma_{rho_rho}
rho_beta = beta_rho'; 
rho_sigma = sigma_rho; 
rho_rho = (1/N) * (1/ sigma^2) * (G*X*beta)' * (G*X*beta)+ (1/N) * (3*diag_square(G)+cross_diag(G)+cross_off_diag(G)+off_diag_square(G) - trace(G) * trace(G)); 

Sigma = zeros(k+2);
Sigma(1:k,1:k) = beta_beta;
Sigma(1:k,k+1) = beta_sigma;
Sigma(1:k,k+2) = beta_rho;

Sigma(k+1,1:k) = sigma_beta;
Sigma(k+1,k+1) = sigma_sigma;
Sigma(k+1,k+2) = sigma_rho;

Sigma(k+2,1:k) = rho_beta;
Sigma(k+2,k+1) = rho_sigma;
Sigma(k+2,k+2) = rho_rho;

%Estimated covariance matrix is (1/N) * Sigma^{-1}
out = Sigma\eye(k+2) ;
out = (1/N) * out;
end

function out = diag_square(A)
%PURPOSE:   compute the sum of square of the diagonal elements of the input matrix
% ---------------------------------------------------
%  USAGE: out = diag_square(A)
%  where:   A is a square matrix
% ---------------------------------------------------
%  RETURNS: sum of square of the diagonal elements of the input matrix
% --------------------------------------------------
out = sum(diag(A).^2);
end 

function out = off_diag_square(A)
%PURPOSE:   compute the sum of square of all off-diagonal elements of the input matrix
% ---------------------------------------------------
%  USAGE: out = off_diag_square(A)
%  where:   A is a square matrix
% ---------------------------------------------------
%  RETURNS: sum of square of all off-diagonal elements of the input matrix
% --------------------------------------------------
[N,N] = size(A);
A(1:N+1:end) = zeros(N,1);
out = sum(sum(A.^2));
end

function out = cross_off_diag(A)
%PURPOSE:   compute the sum of a_{ij}*a_{ji} for all i!=j where a_{ij} denotes the (i,j)-th element of the input matrix A.
% ---------------------------------------------------
%  USAGE: out = cross_off_diag(A)
%  where:   A is a square matrix
% ---------------------------------------------------
%  RETURNS: sum of a_{ij}*a_{ji} for all i!=j
% --------------------------------------------------
[N,N] = size(A);
res = 0;
for i = 1:N
    for j = 1:i-1
        res = res + A(i,j) * A(j,i);
    end
end
out = 2*res;
end

function out = cross_diag(A)
%PURPOSE:   compute the sum of a_{ii}*a_{jj} for all i!=j where a_{ii} denotes the (i,i)-th element of the input matrix A.
% ---------------------------------------------------
%  USAGE: out = cross_diag(A)
%  where:   A is a square matrix
% ---------------------------------------------------
%  RETURNS: sum of a_{ii}*a_{jj} for all i!=j
% --------------------------------------------------
[N,N] = size(A);
res = 0;
for i = 1:N
    for j = 1:i-1
        res = res + A(i,i) * A(j,j);
    end
end
out = 2*res;
end