close all; clear; clc;

%% Load data
%D = csvread('./data/data_1.csv',1,0);
D = csvread('./data/data_2.csv',1,0);
t = D(:,1); y = D(:,2);

%% Single exponential model + t0 + noise floor
% Clearly not a very good model, but conveys the overall 
f = @(p) optimfun_no_irf(p,t,log10(y),[],false);
p0 = [3, 1, 1e4, 10];
options = optimset('Display','final','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[fval,m] = optimfun_no_irf(pf,t,log10(y),[],true); % Plot and get model function

fprintf('Single exponential fit results\n')
fprintf(' Rms residuals: %.4f log10(counts)\n',fval);
fprintf(' t0: %.2f ns\n',pf(1));
fprintf(' Noise floor: %.1f counts\n',pf(2));
fprintf(' Magnitude: %.1fk counts\n',pf(3)/1e3);
fprintf(' Time constant: %.2f ns\n',pf(4));

%% Single exponential model + t0 + noise floor
f = @(p) optimfun_no_irf(p,t,log10(y),[],false);
p0 = [3, 1, 1e4, 10, 1e4, 1];
options = optimset('Display','final','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[fval,m] = optimfun_no_irf(pf,t,log10(y),[],true);

fprintf('Double exponential fit results\n')
fprintf(' Rms residuals: %.4f log10(counts)\n',fval);
fprintf(' t0: %.2f ns\n',pf(1));
fprintf(' Noise floor: %.1f counts\n',pf(2));
fprintf(' 1st component\n')
fprintf('  Magnitude: %.1fk counts\n',pf(3)/1e3);
fprintf('  Time constant: %.2f ns\n',pf(4));
fprintf(' 2nd component\n')
fprintf('  Magnitude: %.1fk counts\n',pf(5)/1e3);
fprintf('  Time constant: %.2f ns\n',pf(6));
