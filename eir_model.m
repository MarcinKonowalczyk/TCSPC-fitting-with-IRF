function y = eir_model(t,p)
%% y = eir_model(t,p)
% Exponential impulse response model - Sum of n expoentials
n = (numel(p)-2)/2; % Number of exp. components
y = zeros(size(t));
for j = 1:n
    offset = 2*(j-1);
    y = y + eir(t,p(1),p(3+offset),p(4+offset));
end
end