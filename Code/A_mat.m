function out = A_mat(beta,sigma,rho,X,W,c1,c2,c3)
%PURPOSE:   computes A matrix defined in Section 3 of the paper and Section S.1 of the supplementary material, which is used as part of the
%computation for the estimated covariance matrix for the robust M-estimates.
% ---------------------------------------------------
%  USAGE: out = A_mat(beta,sigma,rho,X,W,c1,c2,c3)
%  where:   beta = robust M-estimate of the regression coefficient vector
%           sigma = robust M-estimate of the standard deviation of error term
%           rho = robust M-estimate of the spatial dependence parameter
%           X = model matrix (with intercept term in first column if used)
%           W = standardised spatial weight matrix
%           c1 = tuning parameter used for eta_{R,beta}
%           c2 = tuning parameter used for eta_{R,sigma}
%           c3 = tuning parameter used for eta_{R,rho}
% ---------------------------------------------------
%  RETURNS: A matrix as defined in Section 3 of the paper and Section S.1 of the supplementary material.
% --------------------------------------------------

[N,k] = size(X);
max13 = max(c1,c3);
min13 = min(c1,c3);
max23 = max(c2,c3);
min23 = min(c2,c3);
htilde = @(c) 2*(c^2)*(1-normcdf(c)) - 2*c*normpdf(c) + 2 * normcdf(c) - 1;
G = W /(eye(N)-rho*W);

%a_{beta,beta}, a_{beta,sigma}, a_{beta,rho}
beta_beta = htilde(c1) * (1/N) * X' * X;
beta_sigma = zeros(k,1);
beta_rho = (1/(N * sigma)) * (2*c1*c3*(1-normcdf(max13))-2*(min13)*normpdf(max13) + 2*normcdf(min13)-1) * X' * G * X * beta;

%a_{sigma,beta}, a_{sigma,sigma}, a_{sigma,rho}
sigma_beta = zeros(1,k);
sigma_sigma = 2* ((c2)^4) * (1-normcdf(c2)) - 2* ((c2)^3) * normpdf(c2) + 3 * (2*normcdf(c2)-1-2*(c2)*normpdf(c2)) - (htilde(c2))^2 ;
sigma_rho = (1/N) *trace(G) * ( 2* ((c2)^2) * ((c3)^2) * (1-normcdf(max23)) - 2* ((min23)^3) * normpdf(min23) + 3 * (2*normcdf(min23) -1 - 2*(min23)*normpdf(min23)) + ((min23)^2)* (2* (min23) * normpdf(min23) - 2* (max23)*normpdf(max23) + 2 * normcdf(max23)-2*normcdf(min23)) - htilde(c2) * htilde(c3)        );

%a_{rho,beta}, a_{rho,sigma}, a_{rho,rho}
rho_beta = (beta_rho)' ; 
rho_sigma = sigma_rho ; 
rho_rho = (1/(N * (sigma)^2)) * htilde(c3) * (G*X*beta)' * (G*X*beta) + (1/N) * (2* ((c3)^4)* (1-normcdf(c3)) - 2* ((c3)^3)*normpdf(c3) + 3* (2*normcdf(c3)-1-2*(c3)*normpdf(c3)))* (diag_square(G)) + (1/N) * ((htilde(c3))^2) * (off_diag_square(G)+cross_off_diag(G)+cross_diag(G)) - (1/N) * (trace(G) * htilde(c3))^2;

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