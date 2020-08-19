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

f = @(p) optimfun_1(p,t,log10(y1),false);

p0 = [10,1e4,10,5];
pf = fminsearch(f,p0);
[~,m] = optimfun_1(pf,t,log10(y1),true);
% figure(1); clf;
% plot(t,y1,t,m(t,p0),t,m(t,pf));
% a = gca; a.YScale = 'log';
% grid on; box on; axis square;

%% Fit double exponential with a noise floor

f = @(p) optimfun_2(p,t,log10(y1),false);

p0 = [3,1e4,4,1e4,0.1,3];
pf = fminsearch(f,p0);
[~,m] = optimfun_2(pf,t,log10(y1),true);
subplot(4,1,1:3);
plot(t,log10(m(t,p0)));
% grid on; box on; axis square;