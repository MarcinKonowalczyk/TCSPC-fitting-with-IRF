function [delta,m] = optimfun_no_irf(p,t,data,weights,plot_flag)
%% [delta,m] = optimfun_no_irf(p,time,data,weights,plot_flag)
% Sum of n expnential impulse response functions with a noise floor.
% 
% 'p' if the vector of model parameters
% p(1) - t0
% p(2) - noise floor
% p(3:4) - amplitude and time constant of 1st exp. component
% p(5:6) - -//- of 2nd exp. component
% p(7:8) - etc...

%%
if nargin < 5 || isempty(plot_flag), plot_flag = false; end
if nargin < 4 || isempty(weights), weights = ones(size(t)); end

% Make sure parameters are all positive
if any(p<0), delta = Inf; return; end
% Make sure the time constants are not shorter than time axis spacing
if any(p(4:2:end)<2*(t(2)-t(1))), delta = 1e6; return; end

fit_noise_floor = p(2);

%% Model function
m = @(t,p) eir_model(t,p) + fit_noise_floor; % Alias the function so it can be returned
fit = log10(m(t,p)); % Fit
delta = (data-fit).*weights; % Difference between data and fit

%% Plot
if plot_flag
    figure(1); clf;
    s = subplot(4,1,1:3); hold on;
    plot(t,data,t,fit,'k--');
    % Plot all the exp. impulse response functions
    ne = (numel(p)-2)/2; % Number of eir components
    for j = 1:ne
        offset = 2*(j-1);
        y = eir(t,p(1),p(3+offset),p(4+offset)) +  p(2)*ones(size(t));
        plot(t,log10(y))
    end
    grid on; box on;
    title(sprintf('TCSPC data fit, no IRF, %d eir functions',ne))
    xlim([min(t),max(t)]); s.YLim = [0 s.YLim(2)];
    s.XTickLabel = []; ylabel('log_{10}( Intensity )');
    
    s = subplot(4,1,4);
    plot(t,delta,'k');
    xlim([min(t),max(t)]); s.YLim = [-1 1]*max(abs(s.YLim));
    grid on; box on; title('Residuals');
    xlabel('time / ns'); ylabel('\Delta log_{10}( I )');
    drawnow;
end

%% Root mean square residuals
% Make sure no infs (-ve data) or nans (incomplete data)
delta(isinf(delta) | isnan(delta)) = [];
delta = sqrt(mean(delta.^2));
end