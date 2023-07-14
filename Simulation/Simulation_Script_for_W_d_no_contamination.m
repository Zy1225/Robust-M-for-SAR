%% Add the Code folder to the search path
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
addpath(genpath(strcat(newdir,'\Code')))
%% Setting Up
%Setting true parameter value
beta = 1 ;
sigma = 1;
rho = 0.5;

%Setting tuning parameters and maximum number of iteration for the iterative algorithm
c1 = 1.4;
c2 = 2.4;
c3 = 1.65;
max_iter = 25;

%Number of simulation runs
B = 1000;

%Tolerance parameter to determine convergence of iterative algorithm
tol = 0.00001;



%% Simulation Script with n in {20,60,100,140,180,200,400,600} and no contamination
%Generate the spatial weight matrix W_d
N = 20;
W_ori = [];   
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);  
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(1);
core20_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core20_W_d.mat core20_W_d

%Generate the spatial weight matrix W_d
N = 60;
W_ori = [];   
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(2);
core60_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core60_W_d.mat core60_W_d

%Generate the spatial weight matrix W_d
N = 100;
W_ori = []; 
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(3);
core100_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core100_W_d.mat core100_W_d

%Generate the spatial weight matrix W_d
N = 140;
W_ori = [];  
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(4);
core140_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core140_W_d.mat core140_W_d

%Generate the spatial weight matrix W_d
N = 180;
W_ori = [];  
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(5);
core180_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core180_W_d.mat core180_W_d

%Generate the spatial weight matrix W_d
N = 200;
W_ori = [];   
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(6);
core200_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core200_W_d.mat core200_W_d

%Generate the spatial weight matrix W_d
N = 400;
W_ori = [];   
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(7);
core400_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core400_W_d.mat core400_W_d

%Generate the spatial weight matrix W_d
N = 600;
W_ori = [];   
for i=1:N
    for j=1:N
        if i ~=j
            W_ori(i,j) = 1/(abs(i-j));
        else
            W_ori(i,j) = 0;
        end
    end
end

%Row-standardise the spatial weight matrix
W = zeros(N);   
for i=1:N
    W(i,:) = W_ori(i,:) / sum(W_ori(i,:));
end
rng(1111);
X= normrnd(0,1,[N,1]) ; 

rng(8);
core600_W_d = rob_sim(N,W,X,rho,beta,sigma,B,tol,c1,c2,c3,max_iter);
save core600_W_d.mat core600_W_d

%% Simulation Results of 
load core20_W_d.mat
load core60_W_d.mat
load core100_W_d.mat
load core140_W_d.mat
load core180_W_d.mat
load core200_W_d.mat
load core400_W_d.mat
load core600_W_d.mat

%Get results from saved Matlab files
table_core = [result_vector(core20_W_d);result_vector(core60_W_d);result_vector(core100_W_d);result_vector(core140_W_d);result_vector(core180_W_d);result_vector(core200_W_d); result_vector(core400_W_d);result_vector(core600_W_d)];

%Formatting simulation results
result_table = array2table(table_core,...
    'VariableNames',{'Robust_Bias','Robust_RMSE','Robust_ASE/ESD', 'Robust_CP','ML_Bias','ML_RMSE','ML_ASE/ESD', 'ML_CP'});
table_label = table(repelem([20,60,100,140,180,200,400,600],3)',repelem({'beta','sigma','rho'},8)', ...
    'VariableNames',{'N' 'Parameter'});

%Output simulation results
horzcat(table_label,result_table)

