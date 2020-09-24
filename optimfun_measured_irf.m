function [delta,m] = optimfun_measured_irf(p,t,data,irf,weights,plot_flag)
%% [delta,m] = optimfun_measured_irf(p,t,data,irf,weights,plot_flag)
% 
% 'p' if the vector of model parameters
% p(1) - t0
% p(2) - Noise floor
% p(3:4) - amplitude and time constant of 1st exp. component
% p(5:6) - -//- of 2nd exp. component
% p(7:8) - etc...

%% 
if nargin < 6 || isempty(plot_flag), plot_flag = false; end
if nargin < 5 || isempty(weights), weights = ones(size(t)); end

% Make sure parameters are all positive
if any(p<0), delta = Inf; return; end

irf_fun = irf{1};
irf_noise_floor = irf{2};
irf_width = irf{3};

% Make sure the time constants are not shorter than irf_width/5
if any(p(4:2:end)<irf_width/5), delta = Inf; return; end

%% Fit
% Convolve eir_sum with irf, and add the noise floor back up
m = @(t,p) irf_conv(@(k)eir_sum(k,p(1),p(3:end)),irf_fun,t) + irf_noise_floor + p(2);
fit = log10(m(t,p)); % Log10 for fitting
delta = (data-fit).*weights; % Difference between data and fit

%% Root mean square residuals
% Make sure no infs (-ve data) or nans (incomplete data)
delta(isinf(delta) | isnan(delta)) = [];
delta = sqrt(mean(delta.^2));

%% Plot
if plot_flag
    figure(1); clf;
    s = subplot(4,1,1:3); hold on;
    plot(t,data,t,fit,'k--');
    % Plot all the exp. impulse response functions
    t0 = p(1); fit_noise_floor = p(2); exp_params = p(3:end);
    for j = 1:(numel(exp_params)/2)
        offset = 2*(j-1);
        y = @(t) eir(t,t0,exp_params(1+offset),exp_params(2+offset));
        y = irf_conv(y,irf_fun,t) + irf_noise_floor + fit_noise_floor;
        plot(t,log10(y))
    end
    grid on; box on;
    title(sprintf('TCSPC data fit, with IRF, %d eir functions',numel(exp_params)/2))
    xlim([min(t),max(t)]); s.YLim = [0 s.YLim(2)];
    s.XTickLabel = []; ylabel('log_{10}( Intensity )');
    
    s = subplot(4,1,4);
    plot(t,(data-fit).*weights,'k');
    xlim([min(t),max(t)]);
    s.YLim = [-1 1]*max(abs(s.YLim));
    grid on; box on;
    title('Residuals');
    xlabel('time / ns'); ylabel('\Delta log_{10}( I )');
    drawnow;
end
end