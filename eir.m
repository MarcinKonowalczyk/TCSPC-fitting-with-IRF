function y = eir(t,t0,h,tau)
%% y = eir(t,t0,h,tau)
% Exponential impulse response model
%
% t   - time axis
% t0  - time of the impulse
% h   - height of the response
% tau - time constant of the resulting exponential
%

y = h*exp(-t/tau).*(t>t0);