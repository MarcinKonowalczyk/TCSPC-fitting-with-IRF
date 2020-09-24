function y = eir_sum(t,t0,ep)
%% y = eir_sum(t,t0,ep)
% Sum of n expoenntial impulse responses
%
% t       - time axis
% t0      - time of the impulse
% ep(1:2) - height and time constant of the 1st impulse
% ep(3:4) - -//- of the 2nd impulse
% ep(5:6) - ...

y = zeros(size(t));
for j = 1:(numel(ep)/2)
    offset = 2*(j-1);
    y = y + eir(t,t0,ep(1+offset),ep(2+offset));
end
end