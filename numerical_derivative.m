function [dydx,dx] = numerical_derivative(x,y,dx)
%% [dydx,dx] = numerical_derivative(x,y,dx)
% Take a numerical derivative `dydx` of the 1D vector specified by `x` and 
% `y`. The x-axis is assumed to be monotonic.
%
% dydx - numerical derivative vector
% dx - x spacing
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

%% Parse the input
x = x(:); y = y(:); % Make sure x and y are collumn vectors
n = numel(x);
if n~=numel(y)
   error('Input vectors must have the same size');
end

if nargin < 3
    dx = mean(diff(x));
end

%% Handle the trivial cases
if n == 0
    dydx = [];
    dx = [];
    return
elseif n == 1
    dydx = 0;
    dx = [];
    return
end

%% Calculate the derivative for unit spacing
dy = zeros(size(y));
% Forward and backward differences at the edges
dy(1) = y(2) - y(1);
dy(n) = y(n) - y(n-1);

% Center difference everywhere else
dy(2:n-1) = ( y(3:n)-y(1:n-2) )/2;

%% Correct for x spacing
dydx = dy/dx;

end