function y = eir(t,t0,h,tau)
%% y = eir(t,t0,h,tau)
% Exponential impulse response funciton
%
% t   - time axis
% t0  - time of the impulse
% h   - height of the response
% tau - time constant of the resulting exponential

y = zeros(size(t));
region = (t>t0);
y(region) = h*exp(-(t(region)-t0)./tau);
end