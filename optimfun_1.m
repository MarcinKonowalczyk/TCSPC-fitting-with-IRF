function [diff,m] = optimfun_1(p,t,data,plt)

t0 = abs(p(1));
h = abs(p(2));
tau = abs(p(3));
nf = abs(p(4));

eir = @(t,t0,h,tau) h*exp(-(t-t0)/tau).*(t>t0);
m = @(t,p) eir(t,t0,h,tau) + nf*ones(size(t));

fit = log10(m(t,p));

if plt % Plot
    figure(1); clf;
    s = subplot(4,1,1:3);
    plot(t,data,t,fit);
    grid on; box on;
    s = subplot(4,1,4);
    plot(t,data-fit,'k');
    s.YLim = [-1 1]*max(s.YLim);
    grid on; box on;
    drawnow;
end

diff = (data-fit).^2;
% Make sure no infs (-ve data) or nans (incomplete data)
diff(isinf(diff) | isnan(diff)) = [];
diff = mean(diff);
end