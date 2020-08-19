close all; clear; clc;

%% Load data

D = csvread('Decay_t0.csv',9,0);
t = D(:,1); y1 = D(:,2);

D = csvread('IRF_9nm_slit.csv',9,0);
y2 = D(:,2);

%%
figure(1); clf;
plot(t,y1,t,y2);
a = gca; a.YScale = 'log';
grid on; box on; axis square;

%% Fit single exponential with a noise floor.
f = @(p) optimfun_no_irf(p,t,log10(y1),[],false);
p0 = [3, 3, 1e4, 10];
pf = fminsearch(f,p0);
[~,m] = optimfun_no_irf(pf,t,log10(y1),[],true);

%% Fit double exponential with a noise floor
f = @(p) optimfun_no_irf(p,t,log10(y1),[],false);

p0 = [3, 3, 7e3, 4.5, 1e4, 0.3];
pf = fminsearch(f,p0);
[~,m] = optimfun_no_irf(pf,t,log10(y1),[],true);

%% Fit tripple exponential with a noise floor
% Fit function
plt = false; % Whether to plot over the course of fitting
f = @(p) optimfun_no_irf(p,t,log10(y1),[],plt);

% Parameter guess and constraints
p0  = [3,  3, 5e3, 4.6, 4e3, 0.2, 5e3, 0.01];
plb = [0,  0, 1e2,   0, 1e2,   0, 1e2,    0];
pub = [5, 10, 1e5,  10, 1e5,  10, 1e5,   10];

% Linear constraint tau1 > tau2 > tau 3
A = [0 0 0 -1 0 1 0 0; 0 0 0  0 0 -1 0 1];
b = [0 0]';

options = optimoptions('fmincon','Display','Iter');
pf = fmincon(f,p0,A,b,[],[],plb,pub,[],options);
[~,m] = optimfun_no_irf(pf,t,log10(y1),[],true); % Plot and get the model function

%%