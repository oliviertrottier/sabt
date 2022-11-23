function [df,df_err,f] = boxdim_fit(R,N_boxes,varargin)
% Function to find the box dimension (Haussdorff dimension) by fitting the 
% number of boxes of size R needed to cover a shape.
%
% Input
% R = radius of each set of boxes used to cover the shape.
% N_boxes = number of boxes of radius R needed for the cover.
%
% Output
% df = box fractal dimension
% df_err = error of box fractal dimension calculated from non-linear fit
% f = figure handle of fit, if requested
%% Parse input
N_R = numel(R);
%% Parse optional parameters
p = inputParser;
addParameter(p, 'R_min', R(floor(1/3*N_R))); % Minimum box size for fit.
addParameter(p, 'R_max', R(ceil(2/3*N_R))); % Maximum box size for fit.
addParameter(p, 'Plot', nargout==0); % Plot the fit.
parse(p, varargin{:});
options = p.Results;
%%
% Define the power law to be fitted.
modelfun = @(beta, x) beta(2).*x.^(beta(1));

% Define the minimal number of datapoints to perform the non-linear
% regression of N_boxes vs rad.
df_min_datapoints = 5;

Fitted_ind = find(R >= options.R_min & R <= options.R_max);
if length(Fitted_ind) >= df_min_datapoints
    R_fit = R(Fitted_ind);
    N_boxes_fit = N_boxes(Fitted_ind);
    
    fit_model = fitlm(log(R_fit),log(N_boxes_fit));
    df = -fit_model.Coefficients.Estimate(2);
    df_err = fit_model.Coefficients.SE(2);
    beta = [fit_model.Coefficients.Estimate(2), exp(fit_model.Coefficients.Estimate(1))];
else
    df = nan;
    df_err = nan;
end
%% Plot the fit.
if options.Plot && ~isnan(df)
    f = figure;
    loglog(R,N_boxes,'ko');
    hold on
    xlabel('Box size r (\mum)');
    ylabel('Number of boxes N(r)');
    
    fith = plot(R_fit,modelfun(beta,R_fit),'k-','DisplayName',sprintf('N ~ r^{-%.2f \\pm %.2f}',df(1),df_err(1)));
    legend(fith,'Location','NorthEast');
else
    f = [];
end
end