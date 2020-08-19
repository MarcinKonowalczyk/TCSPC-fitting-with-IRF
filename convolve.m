function [x,yc,bm] = convolve(fun,range,bm,h,N,method)
%% [x,yc,bm] = convolve(fun,range,bm,h,N)
% Convolve a funciton specified by `fun` over the `range`
% 
% OUTPUT
%  x  - x axis 
%  yc - convolved signal
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

%% Parse input
if nargin < 5 || isempty(N)
    N = 1001;
else
    if ~mod(N,2) % Make sure N is odd to alwyas have the middle sample
        N = N + 1;
    end
end

if nargin < 6 || isempty(method)
    method = 'conv';
end

bm = abs(bm); % Kernel width must be positive

%% Get the x axis
range2 = [min(range) max(range)];
x = linspace(range2(1),range2(2),N);
%dx = mean(diff(x));
dx = x(2)-x(1); % x spacing

% Special case for bm == 0
if bm == 0
    yc = zeros(size(x));
    return
end

% Pad x on either side to allow for edges where the convolution is invaid to be cut off
pad_N = ceil(bm/dx);
lower_pad = linspace(range2(1)-pad_N*dx,range2(1)-dx,pad_N);
upper_pad = linspace(range2(2)+dx,range2(2)+pad_N*dx,pad_N);
x2 = [lower_pad x upper_pad];

%% Do the convolution
% Evaluate fun at x
y = fun(x2);
[dy,~] = numerical_derivative(x2,y,dx);
%dy = gradient(y)/dx;

% Get the convolution kernel. Coersce size to the closest n of points
K_N = 2*pad_N+1;
bm = pad_N*dx; % Actuall bm used
error('Work In Prrogress');
K = 1;%

% Convolve
switch method
    case 'conv'
        yc = conv(dy,K,'same');
        yc = yc(pad_N+1:end-pad_N);
    case 'fft'
        yc = real(ifft(fft(dy).*fft(K',numel(dy))));
        yc = yc((2*pad_N+1):end);
end

% Interpolate onto the range if needed
if numel(range) > 2
    F = griddedInterpolant(x(:),yc(:),'linear');
    yc = F(range);
    x = range;
end