function results =  adf_new(x,p,l)
% PURPOSE: carry out DF tests on a time-series vector 
%---------------------------------------------------
% USAGE: results = adf(x,p,nlag)
% where:      x = a time-series vector
%             p = order of time polynomial in the null-hypothesis
%                 p = -1, no deterministic part
%                 p =  0, for constant term
%                 p =  1, for constant plus time-trend
%                 p >  1, for higher order polynomial
%             l = # of lagged changes of x included  
%---------------------------------------------------
% RETURNS: a results structure
%         results.meth  = 'adf'
%         results.alpha = estimate of the autoregressive parameter
%         results.adf   = ADF t-statistic
%         results.crit = (6 x 1) vector of critical values
%                        [1% 5% 10% 90% 95% 99%] quintiles    
%         results.nlag = nlag   
%---------------------------------------------------
% SEE ALSO: prt_coint()
%--------------------------------------------------- 
% References: Said and Dickey (1984) 'Testing for Unit Roots in
% Autoregressive Moving Average Models of Unknown Order', 
% Biometrika, Volume 71, pp. 599-607.

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

% Modeled after a similar Gauss routine by
% Sam Ouliaris, in a package called COINT
%
% fixed by
% Frauke Schleer, schleer@zew.de


% error checking on inputs
if (nargin ~= 3)
 error('Wrong # of arguments to adf');
end;

if (p < -1)
 error('p less than -1 in adf');
elseif (cols(x) > 1)
 error('adf cannot handle a matrix -- only vectors');
end;

nobs = rows(x);

if ((nobs - 2*l)+1 < 1)
    error('nlags too large in adf, negative dof');
end

ch = tdiff(x,1);
ch = trimr(ch,1,0);

% Gerard van den Hout suggested the fix below
% Erasmus University Rotterdam.
% The Netherlands.
k = 0 ;
z = [];
while (k < l);
    k = k+1;
    z = [z lag(ch,k)];
end;
z = trimr(z,k,0);

if (p > -1)
 z = [z ptrend(p,rows(z))];
end;

ylag=lag(x,1);
ylag=trimr(ylag,k+1,0);
dep=tdiff(x,1);
y2=trimr(dep,k+1,0);
regressor=[ylag,z];
results=ols(y2,regressor);

results.nlag = l;
results.alpha = results.beta(1,1);
results.adf = results.tstat(1,1);
results.crit = ztcrit(nobs,p); %do not coincide with EViews 7 critical values
results.meth = 'adf';

     
