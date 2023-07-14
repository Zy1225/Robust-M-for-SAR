function [result] = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter)
%PURPOSE:   perform simulations where data are generated with no
%           contamination in the error terms, to compare the performance of the
%           proposed robust M-estimator vs MLE. 
% ---------------------------------------------------
%  USAGE: result = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
%  where:   N = number of spatial units
%           W = standardised spatial weight matrix
%           X = model matrix (with intercept term in first column if used)
%           rho = true spatial dependence parameter
%           beta = true regression coefficient vector
%           sigma = true error standard deviation
%           B = number of simulations
%           tol = tolerance parameter for convergence of iterative algorithm
%           c1 = tuning parameter used for eta_{R,beta}
%           c2 = tuning parameter used for eta_{R,sigma}
%           c3 = tuning parameter used for eta_{R,rho}
%           max_iter = maximum number of iterations for the iterative algorithm
% ---------------------------------------------------
%  RETURNS: a structure
%           result.n = N
%           result.c = [c1,c2,c3]
%           result.B = B
%           result.flag = 1 indicating the use of MC lndet approximation to compute MLE
%           result.robust_array = a (B x 2(k+2)) matrix storing robust M-estimates and their biases for each of the B simulation runs where each row represent a simulation
%                                 run, k is the number of covariates, the first (k+2) columns represent the robust M-estimates, and the last (k+2) columns represent their biases.
%           result.robust_mean_vec = averaged (across B runs) robust M-estimates and their biases    (1 x 2(k+2)) vector
%           result.robust_esd = empirical standard deviation (across B runs) of robust M-estimates   (1 x (k+2)) vector
%           result.robust_RMSE = root mean squared error of robust M-estimates      (1 x (k+2)) vector
%           result.robust_average_se = average estimated standard error of robust M-estimates    (1 x (k+2)) vector
%           result.robust_coverage = empirical coverage probability of 95% confidence intervals based on robust M-estimates   (1 x (k+2)) vector
%           result.convergeratio = proportion of simulation runs with converged iterative algorithm
%           result.mle_array = a (B x 2(k+2)) matrix storing MLE and their biases for each of the B simulation runs where each row represent a simulation
%                                 run, k is the number of covariates, the first (k+2) columns represent the MLE, and the last (k+2) columns represent their biases.
%           result.mle_mean_vec = averaged (across B runs) MLE and their biases    (1 x 2(k+2)) vector
%           result.mle_esd = empirical standard deviation (across B runs) of MLE   (1 x (k+2)) vector
%           result.mle_RMSE = root mean squared error of MLE      (1 x (k+2)) vector
%           result.mle_average_se = average estimated standard error of MLE    (1 x (k+2)) vector
%           result.mle_coverage = empirical coverage probability of 95% confidence intervals based on MLE   (1 x (k+2)) vector

%Compute the smallest eigenvalue for W
 min_eig = min(eig(W));
 [N,k] = size(X);

%Array to store results related to the robust M estimator
robust_array = zeros([B,4+2*k]);

%Array to store results related to MLE
mle_array = zeros([B,4+2*k]);


%Array to store the estimated standard error of robust M estimates
se_robust_array = zeros([B,k+2]);

%Array to store the estimated standard error of MLE
se_mle_array = zeros([B,k+2]);


convergence_counter = 0;
non_convergence_counter = 0;
total_iter = 0;


while convergence_counter < B
    
    %Generate the error vector with no contamination
    eps = normrnd(0,sigma,[N,1]); 
  
    %Generate the y vector according to the SAR model
    y = inv(eye(N) - rho*W) * (X*beta) + inv(eye(N) - rho*W) * eps; 
   
    %Compute MLE using sar() from Spatial Econometrics Toolbox by James P. LeSage
    info.rmin = 1/min_eig;
    info.lflag = 1; %using MC lndet approximation
    mle = sar(y,X,W,info);
    
    %Compute robust M-estimate using the proposed iterative algorithm
    robust_est = sar_robust_est(y,W,X,tol,c1,c2,c3,max_iter,[1,1,0.5]);
    robust_beta = robust_est.beta;
    robust_sigma = robust_est.sigma;
    robust_rho = robust_est.rho;
    
  
    
    %Compute the biases for robust M estimator and MLE
    bias_rho = robust_rho - rho;
    bias_sigma = robust_sigma - sigma;
    bias_beta = robust_beta - beta;

    mle_bias_rho = mle.rho - rho;
    mle_bias_beta = mle.beta' - beta;
    mle_bias_sigma = sqrt(mle.sige) - sigma;
    
    
    if isnan(robust_rho)
        %If iterative algorithm failed to converge
        non_convergence_counter = non_convergence_counter +1;
        total_iter = total_iter + 1;
    else
        %If iterative algorithm converged
        convergence_counter = convergence_counter + 1;
        i = convergence_counter;
        if mod(convergence_counter,50)==0
            fprintf('Convergence No.%d\n', convergence_counter);
        end
        total_iter = total_iter + 1;
        
        %Compute and store standard error estimates for the robust M-estimates and MLE
        se_robust_array(i,:) =  sqrt(diag((1/N)* inv(B_mat(robust_beta,robust_sigma,robust_rho,X,W,c1,c2,c3))*A_mat(robust_beta,robust_sigma,robust_rho,X,W,c1,c2,c3)*(inv(B_mat(robust_beta,robust_sigma,robust_rho,X,W,c1,c2,c3)))'));
        se_mle_array(i,:) = sqrt(diag(mle_var(mle.beta',sqrt(mle.sige),mle.rho,X,W)));
        
    %Store robust M-estimates and MLE together with their biases
    robust_array(i,:) = [robust_beta',robust_sigma,robust_rho,bias_beta',bias_sigma,bias_rho];
    mle_array(i,:) = [mle.beta',sqrt(mle.sige),mle.rho,mle_bias_beta',mle_bias_sigma,mle_bias_rho];
    
    end
    
    
end

%Compute averaged (across B simulation runs) robust M-estimates and MLE, together with their biases
robust_mean_vec = mean(robust_array,'omitnan');
robust_esd = std(robust_array(:,1:k+2),'omitnan');

mle_mean_vec = mean(mle_array,'omitnan');
mle_esd = std(mle_array(:,1:k+2),'omitnan');


%Compute average estimated standard error for the robust M-estimate and MLE
robust_average_se = mean(se_robust_array,'omitnan');
mle_average_se = mean(se_mle_array,'omitnan');

%Compute empirical coverage probability for 95% CI based on the MLE
CI_lower_mle = zeros([N,k+2]);
CI_upper_mle = zeros([N,k+2]);
CI_lower_mle = mle_array(:,1:k+2) - norminv(0.975)* se_mle_array;
CI_upper_mle = mle_array(:,1:k+2) +norminv(0.975)* se_mle_array;
mle_coverage = sum([beta',sigma,rho] > CI_lower_mle & [beta',sigma,rho] < CI_upper_mle, 1);

%Compute empirical coverage probability for 95% CI based on the robust M-estimate
CI_lower_robust = zeros([N,k+2]);
CI_upper_robust = zeros([N,k+2]);
CI_lower_robust = robust_array(:,1:k+2) - norminv(0.975)* se_robust_array;   
CI_upper_robust = robust_array(:,1:k+2) + norminv(0.975)* se_robust_array;   
robust_coverage = sum([beta',sigma,rho] > CI_lower_robust & [beta',sigma,rho] < CI_upper_robust, 1);


%return results

%Function Input
result.n = N;
result.c = [c1,c2,c3];
result.B = B;
result.flag = info.lflag;

%Simulation Results related to robust M-estimate
result.robust_array = robust_array;
result.robust_mean_vec = robust_mean_vec;
result.robust_esd = robust_esd;
result.robust_RMSE = sqrt(robust_mean_vec(k+3:end).^2 + robust_esd.^2);
result.robust_average_se = robust_average_se;
result.robust_coverage = robust_coverage / B; 
result.convergeratio = convergence_counter / total_iter;

%Simulation Results related to MLE
result.mle_array = mle_array;
result.mle_mean_vec = mle_mean_vec;
result.mle_esd = mle_esd;
result.mle_RMSE = sqrt(mle_mean_vec(k+3:end).^2 + mle_esd.^2);
result.mle_average_se = mle_average_se;
result.mle_coverage = mle_coverage / B;
end

    