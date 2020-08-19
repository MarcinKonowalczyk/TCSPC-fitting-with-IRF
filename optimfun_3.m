function [diff,m] = optimfun_3(p,t,data,weights,plt)

if any(p<0)
    diff = Inf;
    return
end

e1 = @(t,p) eir(t,p(1),p(2),p(3));
e2 = @(t,p) eir(t,p(1),p(4),p(5));
e3 = @(t,p) eir(t,p(1),p(6),p(7));
nf = @(t,p) p(8)*ones(size(t));

% Cumulative model
m = @(t,p) e1(t,p) + e2(t,p) + e3(t,p) + nf(t,p);

% Data is in log10(counts), so take log10 of the model too
fit = log10(m(t,p));

diff = (data-fit).*weights;

if plt % Plot
    figure(1); clf;
    s = subplot(4,1,1:3); hold on;
    plot(t,data,t,fit,'k--');
    plot(t,log10(e1(t,p)+nf(t,p)));
    plot(t,log10(e2(t,p)+nf(t,p)));
    plot(t,log10(e3(t,p)+nf(t,p)));
    grid on; box on;
    xlim([min(t),max(t)]);
    s.XTickLabel = [];
    
    s = subplot(4,1,4);
    plot(t,diff,'k');
    xlim([min(t),max(t)]);
    s.YLim = [-1 1]*max(s.YLim);
    grid on; box on;
    drawnow;
end

diff = diff.^2;
% Make sure no infs (-ve data) or nans (incomplete data)
diff(isinf(diff) | isnan(diff)) = [];
diff = mean(diff);
end