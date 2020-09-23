close all; clear; clc;

%% Load data
D = csvread('data_1.csv',1,0);
%D = csvread('data_2.csv',1,0);
t = D(:,1); y1 = D(:,2);

D = csvread('irf.csv',1,0);
irf = D(:,2);

%% Plot
figure(1); clf;
plot(t,y1,t,irf);
a = gca; a.YScale = 'log';
grid on; box on; axis square;

%% Preprocess irf
% Noise floor ought to be removed from the irf, and then added back up
% after the convolution. Otherwise the parts of the signal which are
% ought to be flat (e.g. pre t0 counts), are not.
smooth = @(x,N) conv(x,(1-cos( 2*pi*linspace(0,1,N+3)))./(N+2),'same'); % Hann-window smoothing
irf2 = smooth(irf,0); % Smooth irf a little bit
irf_noise_floor = median(irf2); % Estimate the noise floor
irf2 = irf2 - irf_noise_floor;
irf2 = irf2./sum(irf2); % Normalise irf to unit 
irf_fun = @(k) interp1(t,irf2,k,'linear',0);
irf_cell = {irf_fun,irf_noise_floor}; % <- To be passed to the optimisation function

%% Make up weights vector
w = smooth(log10(y1),11); w(isnan(w)|isinf(w)) = 0;
w = w.^2; % Prop. to sqyare signal
w = w./max(w); % Normalise
w = 0.05+0.95*w; % Still care about noise floor somewhat

%% Fit single exponential with a noise floor and irf
f = @(p) optimfun_irf(p,t,log10(y1),irf_cell,[],false);
p0 = [0, 3, 1e4, 10];
pf = fminsearch(f,p0);
[~,m] = optimfun_irf(pf,t,log10(y1),irf_cell,[],true);

%% Fit double exponential with a noise floor and irf
f = @(p) optimfun_irf(p,t,log10(y1),irf_cell,w,false);
p0 = [0, 3, 7e3, 4.5, 1e4, 0.3];
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_irf(pf,t,log10(y1),irf_cell,w,true);

%% Fit tripple exponential with a noise floor and irf

f = @(p) optimfun_irf(p,t,log10(y1),irf_cell,w,false);
p0 = [0, 0, 6e3, 4.5, 4.7e3, 1.7, 3.6e4, 0.05]; % <- Good starting point
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
%options = optimoptions('patternsearch','Display','iter');
%pf = patternsearch(f,p0,[],[],[],[],[],[],[],options);
[~,m] = optimfun_irf(pf,t,log10(y1),irf_cell,w,true); % Plot and get

%% Same thing but with t0 and noise floor fixed at 0
% Also tripple exponential
f = @(p) optimfun_irf([0,0,p],t,log10(y1),irf_cell,w,false);
p0 = [3e3, 4.5, 2e4, 3, 2.5e4, 1];
p0 = [6e3 4.5 4.7e3 1.7 3.6e4, 0.05];
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_irf([0,0,pf],t,log10(y1),irf_cell,w,true); % Plot and get 