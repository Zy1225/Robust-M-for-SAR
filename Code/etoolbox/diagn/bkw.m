function bkw(x,vnames,fmt)
% PURPOSE: computes and prints BKW collinearity diagnostics
%          variance-decomposition proportions matrix
%---------------------------------------------------
% USAGE: bkw(x,vnames,fmt)
% where:       x = independent variable matrix (from a regression model)
%         vnames = (optional) variable name vector    
%            fmt = (optional) format string, e.g., '%12.6f' or '%12d' 
%                  default = %10.2f   
%---------------------------------------------------
% NOTE: you can use either x-variable names or an ols
%       vnames argument containing x-variable names
% e.g.  vnames = strvcat('x1','x2') 
%---------------------------------------------------
% RETURNS:
%        nothing, just prints the table out
% --------------------------------------------------
% SEE ALSO: dfbeta, rdiag, diagnose
%---------------------------------------------------
% REFERENCES: Belsley, Kuh, Welsch, 1980 Regression Diagnostics
% ----------------------------------------------------

% written by:
% James P. LeSage, 
% Texas State University
% james.lesage@txstate.edu


[nobs, nvar] = size(x);
fid = 1;
% error checking on inputs
if nargin == 3
  Vname = vnames;  
  [nnames,junk] = size(Vname);
  if nnames ~= nvar
 fprintf(fid,'Wrong # of variable names in bkw -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 % make up some generic names
  Vname = 'var1';
  for i=2:nvar
     name = ['var' num2str(i)];
     Vname = strvcat(Vname,name);
  end
  end
  
elseif nargin == 2
  Vname = vnames;  
  [nnames,junk] = size(Vname);

  if nnames ~= nvar
 fprintf(fid,'Wrong # of variable names in bkw -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 % make up some generic names
  Vname = 'var1';
  for i=2:nvar
     name = ['var' num2str(i)];
     Vname = strvcat(Vname,name);
  end
  end
   fmt = '%10.2f';
elseif nargin == 1
% make up some generic names
  Vname = 'var1';
  for i=2:nvar
     name = ['var' num2str(i)];
     Vname = strvcat(Vname,name);
  end
   fmt = '%10.2f';
else
error('Wrong # of arguments to bkw');   
end

[u, d, v] = svd(x,0);

lamda = diag(d(1:nvar,1:nvar));
lamda2 = lamda.*lamda;
v = v.*v;

phi = zeros(nvar,nvar);
for i=1:nvar
phi(i,:) = v(i,:)./lamda2';
end


pi = zeros(nvar,nvar);
for i=1:nvar
phik = sum(phi(i,:));
pi(i,:) = phi(i,:)/phik;
end

% BUG fix suggested by 
% John P. Burkett <burkett@uriacc.uri.edu
lmax = lamda(1);
lmaxvec = lmax*ones(nvar,1);
lout = lmaxvec./lamda;


out = pi';


rnames = strvcat('K(x)',num2str(round(lout)));
in.fmt = fmt;
in.rnames = rnames;
in.cnames = Vname;
fprintf('\n Belsley, Kuh, Welsch Variance-decomposition \n');
mprint(out,in);


