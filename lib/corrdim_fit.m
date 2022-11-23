function [df,df_err,f] = corrdim_fit(R,C,varargin)
% Function to fit the correlation function and calculate the associated
% correlation dimension.
% 
% Input
% R = Radii (\mum) at which the correlation integral (C) is calculated.
% C = Correlation integral calculated at each radius.
%
% Output
% df = correlation dimension
% df_err = error of the correlation dimension inferred from fit.
% f = handle of figure if plot is requested.
%
% Reference:
% https://en.wikipedia.org/wiki/Correlation_dimension
N_R = numel(R);
%% Parse optional parameters
p = inputParser;
addParameter(p, 'R_min', R(floor(1/3*N_R))); % Minimum radius for correlation fit.
addParameter(p, 'R_max', R(ceil(2/3*N_R))); % Maximum radius for correlation fit.
addParameter(p, 'Plot', nargout==0); % Plot the fit.
addParameter(p, 'Piecewise_fit', false); % Perform a 2-pieces piecewise fit of C vs R.
parse(p, varargin{:});
options = p.Results;
%%
% Define the power law to be fitted.
modelfun = @(beta, x) beta(2).*x.^(beta(1));

% Define the minimal number of datapoints needed to perform the non-linear
% regression of C vs R.
df_min_datapoints = 5;

fit_ind = find(R > options.R_min & R < options.R_max);
N_fitted_datapoints = numel(fit_ind);
has_enough_points = N_fitted_datapoints >= df_min_datapoints;
if has_enough_points
    R_fit = R(fit_ind);
    C_fit = C(fit_ind);
    
    if options.Piecewise_fit
        % SLM fit.
        %N_segments = 2;
        %SLM_model = slmengine(log(R),log(C),'degree',1,'plot','on','knots',N_segments+1);
        %SLM_model = slmengine(log(R_fit),log(C_fit),'degree',1,'plot','on','knots',N_segments+1);
        
        % ARES fit.
        %params = aresparams2('cubic', false, 'useMinSpan',1,'maxfuncs',100,'maxFinalFuncs',N_segments+1);
        %model = aresbuild(log(R_fit),log(C_fit), params, [], [], [], [], 0);
        
        % Construct piecewise fit.
        Piecewise_model = @(df1,b1,c1,df2,x) (x<c1).*(df1.*x + b1) + (x>=c1).*(df2.*(x-c1) + df1.*c1 + b1);
        Piecewise_fittype = fittype(Piecewise_model);
        Fit_options = fitoptions(Piecewise_fittype);
        fit_options = fitoptions(fit_options,'Lower',[0 -Inf min(log(R_fit)) -Inf],'StartPoint',[1 1 1 mean(log(R_fit))]);
        [Piecewise_fit,Piecewise_gof] = fit(log(R_fit),log(C_fit),Piecewise_fittype,fit_options);
        
        % Calculate the standard errors of the coefficients.
        alpha = 0.95;
        ConfInt = confint(Piecewise_fit, alpha);
        t = tinv((1+alpha)/2, Piecewise_gof.dfe);
        SE = diff(ConfInt)/(2*t);
        df = [Piecewise_fit.df1 Piecewise_fit.df2];
        df_err = SE([1 4]);
    else
        % Fit the logged data with a linear fit.
        fit_model = fitlm(log(R_fit),log(C_fit));
        
        % Send warning if the fit is poor.
        if fit_model.Rsquared.Ordinary < 0.95
            warning('The R^2 value (%.3f) for the correlation dimension fit is below 0.95.',fit_model.Rsquared.Ordinary);
        end
        
        df = fit_model.Coefficients.Estimate(2);
        df_err = fit_model.Coefficients.SE(2);
        beta = [fit_model.Coefficients.Estimate(2), exp(fit_model.Coefficients.Estimate(1))];
        
        % Fit the data with a non-linear fit.
%         beta0 = [fit_model.Coefficients.Estimate(2), exp(fit_model.Coefficients.Estimate(1))];
%         [beta, ~, ~, CovB] = nlinfit(R_fit, C_fit, modelfun, beta0);
%         
%         % Record the fractal dimension and its error.
%         df = beta(1);
%         df_err = sqrt(CovB(1, 1));
    end
else
    warning('The correlation dimension did not have enough data points (%d < %d)',N_fitted_datapoints,df_min_datapoints);
    df = nan;
    df_err = nan;
end

% Plot the fit.
if options.Plot && has_enough_points
    f = figure;
    loglog(R,C,'ko');
    hold on
    xlabel('Radius r (\mum)');
    ylabel('Correlation integral C(r)');
    if options.Piecewise_fit
        R1_fitted = logspace(log10(R_fit(1)),log10(exp(Piecewise_fit.c1)));        
        C1_fitted = exp(Piecewise_fit(log(R1_fitted)));
        
        R2_fitted = logspace(log10(exp(Piecewise_fit.c1)),log10(R_fit(end)));        
        C2_fitted = exp(Piecewise_fit(log(R2_fitted)));
        
        Fit_h(1) = plot(R1_fitted,C1_fitted,'r-','DisplayName',sprintf('C ~ r^{%.2f}',df(1)));
        Fit_h(2) = plot(R2_fitted,C2_fitted,'k-','DisplayName',sprintf('C ~ r^{%.2f}',df(2)));
    else
        Fit_h = plot(R_fit,modelfun(beta,R_fit),'k-','DisplayName',sprintf('C ~ r^{%.2f \\pm %.2f}',df(1),df_err(1)));
    end
    legend(Fit_h,'Location','NorthWest');
else
    f = [];
end
end