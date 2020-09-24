close all; clear; clc;

%% Load data
D = csvread('./data/data_1.csv',1,0);
%D = csvread('./data/data_2.csv',1,0);
t = D(:,1); y = D(:,2);
y(y<=0) = 1; % Make sure there is at least one count in each bin
y = log10(y);

%% IRF model
% Gaussian lineshape
irf_fun = @(x,center,width,height) exp( -4.*log(2).*((x-center)./width).^2 ) .* height;

%% Make up weights vector
w = smooth(y,11); w(isnan(w)|isinf(w)) = 0;
w = w.^2; % Prop. to sqyare signal
w = w./max(w); % Normalise
w = 0.05+0.95*w; % Still care about noise floor somewhat

%% Single exponential + gaussian IRF model + t0 + noise floor
f = @(p) optimfun_fitted_irf(p,t,y,irf_fun,[],false);
p0 = [3, 1, 0.1, 1e4, 10];
pf = fminsearch(f,p0);
[~,m] = optimfun_fitted_irf(pf,t,y,irf_fun,[],true); % Plot and get the model function

%% Double exponential + gaussian IRF model + t0 + noise floor
f = @(p) optimfun_fitted_irf(p,t,y,irf_fun,w,false);
p0 = [3, 1, 0.1, 7e3, 4.5, 1e4, 0.3];
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_fitted_irf(pf,t,y,irf_fun,w,true);

%% Tripple exponential + gaussian IRF model + t0 + noise floor
f = @(p) optimfun_fitted_irf(p,t,y,irf_fun,w,false);
p0 = [3, 1, 0.3, 6e3, 4.5, 4.7e3, 1.7, 3.6e4, 0.2]; % <- Relatively good starting point
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[fval,m] = optimfun_fitted_irf(pf,t,y,irf_fun,w,true);

% Print fit result
fprintf('Tripple exponential fit results\n')
fprintf(' Rms residuals: %.4f log10(counts)\n',fval);
fprintf(' t0: %.2f ns\n',pf(1));
fprintf(' Noise floor: %.1f counts\n',pf(2));
fprintf(' IRF width: %.1f ps\n',pf(3)*1e3);
enum = {'1st','2nd','3rd'};
for j = 1:3
    offset = 2*(j-1);
    fprintf(' %s component\n',enum{j})
    fprintf('  Magnitude: %.0fk counts\n',pf(4+offset)/1e3);
    fprintf('  Time constant: %.2f ns\n',pf(5+offset));
end