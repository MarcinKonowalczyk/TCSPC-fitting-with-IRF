function [diff,m] = optimfun_no_irf(p,time,data,weights,plot_flag)
%% [diff,m] = optimfun_no_irf(p,time,data,weights,plot_flag)
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
if nargin < 4 || isempty(weights), weights = ones(size(time)); end

% Make sure parameters are all positive
if any(p<0), diff = Inf; return; end

fit_noise_floor = p(2);

%% Model function
m = @(t,p) eir_model(t,p) + fit_noise_floor; % Alias the function so it can be returned
fit = log10(m(time,p)); % Fit
diff = (data-fit).*weights; % Difference between data and fit

%% Plot
if plot_flag
    figure(1); clf;
    s = subplot(4,1,1:3); hold on;
    plot(time,data,time,fit,'k--');
    % Plot all the exp. impulse response functions
    ne = (numel(p)-2)/2; % Number of eir components
    for j = 1:ne
        offset = 2*(j-1);
        y = eir(time,p(1),p(3+offset),p(4+offset)) +  p(2)*ones(size(time));
        plot(time,log10(y))
    end
    grid on; box on;
    title(sprintf('TCSPC data fit, no IRF, %d eir functions',ne))
    xlim([min(time),max(time)]); s.YLim = [0 s.YLim(2)];
    s.XTickLabel = []; ylabel('log10( Intensity )');
    
    s = subplot(4,1,4);
    plot(time,diff,'k');
    xlim([min(time),max(time)]); s.YLim = [-1 1]*max(abs(s.YLim));
    grid on; box on; title('Residuals');
    xlabel('time / ns'); ylabel('\Delta log10( I )');
    drawnow;
end

%% Root mean square residuals
% Make sure no infs (-ve data) or nans (incomplete data)
diff(isinf(diff) | isnan(diff)) = [];
diff = sqrt(mean(diff.^2));
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