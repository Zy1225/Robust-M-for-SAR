function output = result_vector(result)
%PURPOSE:  Format the simulation results obtained from rob_sim() into a
%((k+2) x 8) matrix, similar to the tables in the paper
% ---------------------------------------------------
% USAGE: output = result_vector(result);
% where:    result = output from rob_sim()
% ---------------------------------------------------
%  RETURNS: a ((k+2) x 8) matrix where each row represent a different model
%  parameter, the first four columns represent (Bias,RMSE,ASE/ESD,CP) of
%  the robust M-estimates, and the last four columns represent
%  (Bias,RMSE,ASE/ESD,CP) of the MLE

k = length(result.mle_esd) - 2;
output = [result.robust_mean_vec(k+3:end)',result.robust_RMSE',result.robust_average_se' ./ result.robust_esd', result.robust_coverage', result.mle_mean_vec(k+3:end)',result.mle_RMSE',result.mle_average_se' ./ result.mle_esd',result.mle_coverage'];



end

