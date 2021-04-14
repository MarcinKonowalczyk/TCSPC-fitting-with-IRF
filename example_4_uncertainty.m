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
smooth = @(x,N) conv(x,(1-cos( 2*pi*linspace(0,1,N+3)))./(N+2),'same'); % Hann-window smoothing

%% Make up weights vector
w = smooth(y,11); w(isnan(w)|isinf(w)) = 0;
w = w.^2; % Prop. to square signal
w = w./max(w); % Normalise
w = 0.05+0.95*w; % Still care about noise floor somewhat

%% Double exponential + gaussian IRF model + t0 + noise floor
f = @(p) optimfun_fitted_irf(p,t,y,irf_fun,w,false);
p0 = [3, 1, 0.1, 7e3, 4.5, 1e4, 0.3];
options = optimset('Display','none','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_fitted_irf(pf,t,y,irf_fun,w,false);

%%
fit = log10(m(t,pf));
r = (y-fit).*w;
nr = r(randsample(numel(r),numel(r),true));
figure(1); clf; hold on;
plot(t,fit+r./w,t,log10(round(10.^max(fit+nr./w,0))));