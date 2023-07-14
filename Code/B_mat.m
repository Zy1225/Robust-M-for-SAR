function out = B_mat(beta,sigma,rho,X,W,c1,c2,c3)
%PURPOSE:   computes B matrix defined in Section 3 of the paper and Section S.1 of the supplementary material, which is used as part of the
%computation for the estimated covariance matrix for the robust M-estimates.
% ---------------------------------------------------
%  USAGE: out = B_mat(beta,sigma,rho,X,W,c1,c2,c3)
%  where:   beta = robust M-estimate of the regression coefficient vector
%           sigma = robust M-estimate of the standard deviation of error term
%           rho = robust M-estimate of the spatial dependence parameter
%           X = model matrix (with intercept term in first column if used)
%           W = standardised spatial weight matrix
%           c1 = tuning parameter used for eta_{R,beta}
%           c2 = tuning parameter used for eta_{R,sigma}
%           c3 = tuning parameter used for eta_{R,rho}
% ---------------------------------------------------
%  RETURNS: B matrix as defined in Section 3 of the paper and Section S.1 of the supplementary material.
% --------------------------------------------------
[N,k] = size(X);
G = W /(eye(N)-rho*W);

%b_{beta,beta}, b_{beta,sigma}, b_{beta,rho}
beta_beta = (1/(N*sigma))* (2*normcdf(c1)-1)* X' * X;
beta_sigma = zeros(k,1);
beta_rho = (1/(N*sigma))*(2*normcdf(c1)-1)* X' * G * X * beta;

%b_{sigma,beta}, b_{sigma,sigma}, b_{sigma,rho}
sigma_beta = zeros(1,k);
sigma_sigma = (2/sigma) * (2*normcdf(c2)-1-2*(c2)*normpdf(c2));
sigma_rho = (2/N) * trace(G) * (2*normcdf(c2)-1-2*(c2)*normpdf(c2));

%b_{rho,beta}, b_{rhoa,sigma}, b_{rho,rho}
rho_beta = (1/(N * sigma^2)) * (2*normcdf(c3)-1) * (G*X*beta)' * X;
rho_sigma = (2/(N * sigma)) * trace(G) * (2*normcdf(c3)-1-2*(c3)*normpdf(c3));
rho_rho = (1/(N * sigma^2)) * (2*normcdf(c3)-1) * (G*X*beta)' * (G*X*beta) + (2/N)* (2*normcdf(c3)-1-2*(c3)*normpdf(c3)) * (diag_square(G)) + ((2*normcdf(c3)-1)^2) * (1/N) * (off_diag_square(G) + cross_off_diag(G));

out = zeros(k+2);
out(1:k,1:k) = beta_beta;
out(1:k,k+1) = beta_sigma;
out(1:k,k+2) = beta_rho;

out(k+1,1:k) = sigma_beta;
out(k+1,k+1) = sigma_sigma;
out(k+1,k+2) = sigma_rho;

out(k+2,1:k) = rho_beta;
out(k+2,k+1) = rho_sigma;
out(k+2,k+2) = rho_rho;
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