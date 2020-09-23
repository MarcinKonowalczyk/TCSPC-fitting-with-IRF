function [diff,m] = optimfun_irf(p,t,data,irf,weights,plot_flag)
%% [diff,m] = optimfun_irf(p,t,data,irf,weights,plot_flag)
% Sum of n expnential impulse responses with a noise floor.
% 
% 'p' if the vector of model parameters
% p(1) - t0. This is given by the IRF and 
% p(2) - Noise floor
% p(3:4) - amplitude and time constant of 1st exp. component
% p(5:6) - -//- of 2nd exp. component
% p(7:8) - etc...

%% 
if nargin < 6 || isempty(plot_flag), plot_flag = false; end
if nargin < 5 || isempty(weights), weights = ones(size(t)); end

% Make sure parameters are all positive
if any(p<0), diff = Inf; return; end

% Unpack input
fit_noise_floor = p(2);
irf_fun = irf{1};
irf_noise_floor = irf{2};

%% Fit
fit = irf_conv(@(k)eir_model(k,p),irf_fun,t); % Convolve eir_model with irf
fit = fit + irf_noise_floor + fit_noise_floor; % Add the noise floor back up
fit = log10(fit); % Log10 for fitting
diff = (data-fit).*weights; % Difference between data and fit

%% Plot
if plot_flag
    figure(1); clf;
    s = subplot(4,1,1:3); hold on;
    plot(t,data,t,fit,'k--');
    
    % Plot all the exp. impulse response functions
    ne = (numel(p)-2)/2; % Number of eir components
    for j = 1:ne
        offset = 2*(j-1);
        y = @(t) eir(t,p(1),p(3+offset),p(4+offset));
        plot(t,log10(irf_conv(y,irf_fun,t) + irf_noise_floor + fit_noise_floor))
    end
    grid on; box on;
    title(sprintf('TCSPC data fit, with IRF, %d eir functions',ne))
    xlim([min(t),max(t)]); ylim([-1 5]);
    ylabel('log10( Intensity )');
    s.XTickLabel = [];
    
    s = subplot(4,1,4);
    plot(t,diff,'k');
    xlim([min(t),max(t)]);
    s.YLim = [-1 1]*max(abs(s.YLim));
    grid on; box on;
    title('Residuals');
    xlabel('time / ps'); ylabel('\Delta log10( I )');
    drawnow;
end

%% Root mean square residuals
% Make sure no infs (-ve data) or nans (incomplete data)
diff(isinf(diff) | isnan(diff)) = [];
diff = sqrt(mean(diff.^2));

% Alias the model function so it can be an output argument
if nargout > 1
    m = @(t,p) irf_conv(@(k)eir_model(k,p),irf_fun,t) + irf_noise_floor + fit_noise_floor;
end
end

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