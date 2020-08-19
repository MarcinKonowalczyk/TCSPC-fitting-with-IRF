function [diff,m] = optimfun_2(p,t,data,weights,plt)

if any(p<0)
    diff = Inf;
    return
end
    

eir = @(t,t0,h,tau) h*exp(-(t-t0)/tau).*(t>t0);
m = @(t,p)...
    eir(t,p(1),p(2),p(3)) + ...
    eir(t,p(1),p(4),p(5)) + ...
    p(6)*ones(size(t));

fit = log10(m(t,p));

diff = (data-fit).*weights;

if plt % Plot
    figure(1); clf;
    s = subplot(4,1,1:3); hold on;
    plot(t,data,t,fit,'k--');
    plot(t,log10(eir(t,p(1),p(2),p(3)) + p(6)*ones(size(t))));
    plot(t,log10(eir(t,p(1),p(4),p(5)) + p(6)*ones(size(t))));
    grid on; box on;
    s = subplot(4,1,4);
    plot(t,diff,'k');
    s.YLim = [-1 1]*max(s.YLim);
    grid on; box on;
    drawnow;
end

diff = diff.^2;
% Make sure no infs (-ve data) or nans (incomplete data)
diff(isinf(diff) | isnan(diff)) = [];
diff = mean(diff);
end