close all; clear; clc;

%% Load data
D = csvread('Decay_t0.csv',9,0);
t = D(:,1); y1 = D(:,2);
D = csvread('IRF_9nm_slit.csv',9,0);
irf = D(:,2);

%% Plot
figure(1); clf;
plot(t,y1,t,irf);
a = gca; a.YScale = 'log';
grid on; box on; axis square;

%% Preprocess irf
smooth = @(x,N) conv(x,(1-cos( 2*pi*linspace(0,1,N+3)))./(N+2),'same'); % Hann-window smoothing
irf2 = smooth(irf,0); % Smooth irf a little bit
irf_noise_floor = median(irf2);
irf2 = irf2 - irf_noise_floor;
irf2 = irf2./sum(irf2);
irf_fun = @(k) interp1(t,irf2,k,'linear',0);
irf_cell = {irf_fun,irf_noise_floor}; % <- To be passed to the optimisation functions

%% Make up weights vector
w = smooth(log10(y1),11); w(isnan(w)|isinf(w)) = 0;
w = w.^2; % Prop. to sqyare signal
w = w./max(w); % Normalise
w = 0.05+0.95*w; % Still care about noise floor somewhat

%% Fit single exponential with a noise floor
% Clearly not a very good model
f = @(p) optimfun_no_irf(p,t,log10(y1),w,false);
p0 = [3, 3, 1e4, 10];
pf = fminsearch(f,p0);
[fval,m] = optimfun_no_irf(pf,t,log10(y1),w,true); % Plot and get 

fprintf('Single exponential fit results\n')
fprintf(' t0: %.2f ps\n',pf(1));
fprintf(' Noise floor: %.1f counts\n',pf(2));
fprintf(' Exponential magnitude: %.0f counts\n',pf(3));
fprintf(' Time constant: %.2f ps\n',pf(4));
fprintf(' Rms residuals');

%% Fit double exponential with a noise floor
f = @(p) optimfun_no_irf(p,t,log10(y1),[],false);
p0 = [3, 3, 7e3, 4.5, 1e4, 0.3];
pf = fminsearch(f,p0);
[~,m] = optimfun_no_irf(pf,t,log10(y1),[],true);

%% Fit single exponential with a noise floor and irf
f = @(p) optimfun_irf(p,t,log10(y1),irf_cell,[],false);
p0 = [0, 3, 1e4, 10];
pf = fminsearch(f,p0);
[~,m] = optimfun_irf(pf,t,log10(y1),irf_cell,[],true);

%% Fit double exponential with a noise floor and irf
f = @(p) optimfun_irf(p,t,log10(y1),irf_cell,w,false);
p0 = [0, 3, 7e3, 4.5, 1e4, 0.3];
options = optimset('Display','iter','MaxFunEvals',1e6);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_irf(pf,t,log10(y1),irf_cell,w,true);

%% Fit tripple exponential with a noise floor and irf

w = smooth(log10(y1),11); w(isnan(w)|isinf(w)) = 0;
w = w./max(w); % Normalise
w = 0.1+0.9*w; % Still care about noise floor somewhat

f = @(p) optimfun_irf(p,t,log10(y1),irf_cell,w,false);
p0 = [0, 0, 3e3, 4.5, 2e4, 1, 2.5e4, 0.05]; % <- Good starting point
%options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e6);
%pf = fminsearch(f,p0,options);
options = optimoptions('patternsearch','Display','iter');
pf = patternsearch(f,p0,[],[],[],[],[],[],[],options);
[~,m] = optimfun_irf(pf,t,log10(y1),irf_cell,w,true); % Plot and get

%% Same thing but with t0 and noise floor fixed at 0
% Also tripple exponential
f = @(p) optimfun_irf([0,0,p],t,log10(y1),irf_cell,w,false);
p0 = [3e3, 4.5, 2e4, 3, 2.5e4, 1];
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e6);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_irf([0,0,pf],t,log10(y1),irf_cell,w,true); % Plot and get 
