close all; clear; clc;

%% Load data
D = csvread('./data/data_1.csv',1,0);
%D = csvread('./data/data_2.csv',1,0);
t = D(:,1); y = D(:,2);
y(y<=0) = 1; % Make sure there is at least one count in each bin
y = log10(y);

D = csvread('./data/irf.csv',1,0);
irf = D(:,2);

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

% Determine the width of the irf
gaussian = @(x,center,width,height) exp( -4.*log(2).*((x-center)./width).^2 ) .* height;
[max_irf,I] = max(irf);
p0 = [t(I),0.1];
f = @(p) mean( (irf-gaussian(t,p(1),p(2),max_irf)).^2 );
pf = fminsearch(f,p0);
irf_width = pf(2);

% To be passed to the optimisation function
irf_cell = {irf_fun,irf_noise_floor,irf_width};

%% Make up weights vector
w = smooth(y,11); w(isnan(w)|isinf(w)) = 0;
w = w.^2; % Prop. to square signal
w = w./max(w); % Normalise
w = 0.05+0.95*w; % Still care about noise floor somewhat

%% Single exponential + IRF + t0 + noise floor
% Don't use weights for this fit, since it only gives an overall idea of the signal anyway
f = @(p) optimfun_measured_irf(p,t,y,irf_cell,[],false);
p0 = [0, 0, 1e4, 10];
pf = fminsearch(f,p0);
[~,m] = optimfun_measured_irf(pf,t,y,irf_cell,[],true); % Plot and get the model function

%% Double exponential + IRF + t0 + noise floor
f = @(p) optimfun_measured_irf(p,t,y,irf_cell,w,false);
p0 = [0, 3, 7e3, 4.5, 1e4, 0.3];
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_measured_irf(pf,t,y,irf_cell,w,true);

%% Tripple exponential + IRF + t0 + noise floor
f = @(p) optimfun_measured_irf(p,t,y,irf_cell,w,false);
p0 = [0, 0, 6e3, 4.5, 4.7e3, 1.7, 3.6e4, 0.2]; % <- Relatively good starting point
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',1e4);
pf = fminsearch(f,p0,options);
[fval,m] = optimfun_measured_irf(pf,t,y,irf_cell,w,true);

% Print fit result
fprintf('Tripple exponential fit results\n')
fprintf(' Rms residuals: %.4f log10(counts)\n',fval);
fprintf(' t0: %.2f ns\n',pf(1));
fprintf(' Noise floor: %.1f counts\n',pf(2));
enum = {'1st','2nd','3rd'};
for j = 1:3
    offset = 2*(j-1);
    fprintf(' %s component\n',enum{j})
    fprintf('  Magnitude: %.1fk counts\n',pf(3+offset)/1e3);
    fprintf('  Time constant: %.2f ns\n',pf(4+offset));
end

%% Tripple exponential + IRF (t0 and noise floor fixed at 0)
% Since t0 and noise floor are expressed by the IRF, they ought to be 0.
% Hence, might as well reduce model complexity and remove them.
f = @(p) optimfun_measured_irf([0,0,p],t,y,irf_cell,w,false);
p0 = [6e3 4.5 4.7e3 1.7 3.6e4, 0.2];
options = optimset('Display','iter','MaxFunEvals',1e6,'MaxIter',5e4);
pf = fminsearch(f,p0,options);
[~,m] = optimfun_measured_irf([0,0,pf],t,y,irf_cell,w,true);

% Print fit result
fprintf('Tripple exponential fit results\n')
fprintf(' Rms residuals: %.4f log10(counts)\n',fval);
enum = {'1st','2nd','3rd'};
for j = 1:3
    offset = 2*(j-1);
    fprintf(' %s component\n',enum{j})
    fprintf('  Magnitude: %.0fk counts\n',pf(1+offset)/1e3);
    fprintf('  Time constant: %.2f ns\n',pf(2+offset));
end
