function [Lp,Lp_err,f] = calculate_PersistenceLength(Tree,varargin)
% Function to calculate the persistence length (Lp) from the paths of
% branches in a tree structure.
%
% Input
% Tree = tree structure
%
% Output
% Lp = Persistence length
% Lp_err = Error of the persistence length
% f = handle of figure if plot is requested
%% Parse optional parameters
p = inputParser;
addParameter(p, 'Lengthscale', 1); % Multiply the tree nodes position by a given scale.
addParameter(p, 'Plot', nargout==0); % Plot the fit
addParameter(p, 'Max_fit_length', 'mean_length' ,@(x) isnumeric(x) || ismember(x,{'mean_length','2x_mean_length'})); % Maximum path length used in fit.
addParameter(p, 'Method', 'tangent_vec'); % Method used to fit the branches.
addParameter(p, 'PathLength_binwidth', []); % Binwidth used to bin the path lengths. The default is the minimum stepsize found.
parse(p, varargin{:});
options = p.Results;
%%
% Terminate early if the tree is empty.
N_branches = numel(Tree);
if N_branches == 0
    Lp = nan;
    Lp_err = nan;
    f = [];
    return;
end

% Calculate useful data about the branches.
for i=1:N_branches
    % Scale the tree nodes positions.
    Tree(i).PointsPos = Tree(i).PointsPos * options.Lengthscale;
    Tree(i).Length = size(Tree(i).PointsPos,1) - 1;
    
    % Define the tangent vectors that make the branch.
    Tree(i).TangentVec = diff(Tree(i).PointsPos);
    
    % Normalize the tangent vectors.
    Nodes_dist = sqrt(sum(Tree(i).TangentVec.^2,2));
    Tree(i).Nodes_dist = Nodes_dist;
    Tree(i).TangentVec = bsxfun(@rdivide,Tree(i).TangentVec,Nodes_dist);
    Tree(i).PathLength = [0;cumsum(Nodes_dist)];
    
    % Calculate the angle of each segment and the change in angle.
    Tree(i).Theta = atan2(Tree(i).TangentVec(:,2),Tree(i).TangentVec(:,1));
    Tree(i).Delta_Theta = mod(diff(Tree(i).Theta) + pi,2*pi) - pi;
    
    % Define the squared distance fromt the origin.
    Tree(i).SquaredDist = sum(bsxfun(@minus,Tree(i).PointsPos,Tree(i).PointsPos(1,:)).^2,2);
end
PathLength_max = max(cell2mat({Tree.PathLength}'));
Mean_length = mean(cellfun(@(x) x(end),{Tree.PathLength}));

% Determine the maximum path length to fit.
switch options.Max_fit_length
    case 'mean_length'
        PathLengthMax_fit = Mean_length;
    case '2x_mean_length'
        PathLengthMax_fit = 2*Mean_length;
    otherwise
        PathLengthMax_fit = options.Max_fit_length;
end

Nodes_dist_all = cell2mat({Tree.Nodes_dist}');
Nodes_dist_avg = mean(Nodes_dist_all);
Nodes_dist_std = std(Nodes_dist_all);
%%
switch options.Method
    case 'tangent_vec'
        % Fit the correlation function of the tangent vectors.
        
        % Initialize cell array to save the tangent vector inner products.
        if isempty(options.PathLength_binwidth)
            PathLengths_binwidth = min(Nodes_dist_all);
        else
            PathLengths_binwidth = options.PathLength_binwidth;
        end
        PathLengths_BinCenters = 0:PathLengths_binwidth:PathLength_max;
        PathLengths_BinEdges = [PathLengths_BinCenters - PathLengths_binwidth/2,PathLengths_BinCenters(end) + PathLengths_binwidth/2];
        N_bins = numel(PathLengths_BinCenters);
        
        % Remove the branches that are too short to calculate the
        % tangential angle.
        Is_branch_short = cellfun(@(x) x < 2,{Tree.Length});
        Tree = Tree(~Is_branch_short);
        N_branches = numel(Tree);
        
        for i=1:N_branches
            % Calculate cos(\Delta \theta) and the path length separation
            % for all pairs possible. The last point is omitted since it
            % has no definition of \theta.
            Branch_N_nodes = Tree(i).Length;
            [ind_upper,ind_lower] = triind2sub(Branch_N_nodes,(1:Branch_N_nodes*(Branch_N_nodes-1)/2)');
            Pair_indices = [ind_lower,ind_upper];
            if numel(ind_upper) == 1
                ds = diff(reshape(Tree(i).PathLength(Pair_indices),1,2),1,2);
                dtheta = diff(reshape(Tree(i).Theta(Pair_indices),1,2),1,2);
            else
                ds = diff(Tree(i).PathLength(Pair_indices),1,2);
                dtheta = diff(Tree(i).Theta(Pair_indices),1,2);
            end
            
            % Bin the path length separation.
            [~,Pair_stepsize_binind] = histc(ds,PathLengths_BinEdges);
            
            % Calculate the average cos(\Delta \theta) at each path length
            % difference.
            Tree(i).cos_theta_avg = accumarray(Pair_stepsize_binind(:),cos(dtheta),[N_bins,1],@mean);
        end
        N_branches_per_bin = histc(cell2mat({Tree.PathLength}'),PathLengths_BinEdges);
        N_branches_per_bin = N_branches_per_bin(1:end-1);
        
        cos_theta_all = cell2mat({Tree.cos_theta_avg});
        cos_theta_sum = sum(cos_theta_all,2);
        InnerProducts_avg = cos_theta_sum./N_branches_per_bin;
        InnerProducts_avg(1) = 1;
        
        % Fit the average of the tangent vector inner products to an exponential
        % function. The length scale of the exponential will correspond to the
        % persistence length.
        AutoCorrModel = @(b,x) exp(-(x)./(b(1)));
        PathLengthMaxInd = find(PathLengths_BinEdges >= PathLengthMax_fit,1);
        Percent20Ind = round(0.2*numel(PathLengths_BinEdges));
        if isempty(PathLengthMaxInd)
            FitMaxInd = Percent20Ind;
        else
            FitMaxInd = min(Percent20Ind,PathLengthMaxInd);
        end
        
        PathLengths_Fit = PathLengths_BinCenters(1:FitMaxInd);
        InnerProducts_avg_Fit = InnerProducts_avg(1:FitMaxInd)';
        
        % Find a good guess for beta by fitting the logged data.
        Fit_approx = fitlm(PathLengths_Fit,log(InnerProducts_avg_Fit));
        beta_guess = -1/Fit_approx.Coefficients.Estimate(2);
        
        % Fit with a non-linear least squares.
        [Lp,~,~,Lp_err] = nlinfit(PathLengths_Fit,InnerProducts_avg_Fit,AutoCorrModel,beta_guess);
        
    case 'squared_dist'
        % Fit the end-to-end distance squared averaged over all branches.
        % Use all points on a branch when assembling end-to-end distance
        % info
        
        % Calculate the end-to-end distance of each branch.
        EndToEndDistance = cell(N_branches,1);
        PathLengths = cell(N_branches,1);
        for i=1:N_branches
            EndToEndDistance{i} = Tree(i).SquaredDist;
            PathLengths{i} = Tree(i).PathLength;
        end
        PathLengths_all = cell2mat(PathLengths);
        EndToEndDistance_all = cell2mat(EndToEndDistance);
        
        % Bin the path lengths and calculate the average end-to-end
        % distance in each bin.
        PathLengths_BinEdges = (0:.1:ceil(max(PathLengths_all)/0.1 + 1)*0.1)';
        PathLengths_BinCenters = PathLengths_BinEdges(1:end-1) + diff(PathLengths_BinEdges(1:2))/2;
        [BinCounts,PathLengthsBinInd] = histc(PathLengths_all,PathLengths_BinEdges);
        
        % Calculate the end-to-end distance at each path length averaged over all branches.
        EndToEndDistance_avg = accumarray(PathLengthsBinInd,EndToEndDistance_all,[],@mean);
        EndToEndDistance_se = accumarray(PathLengthsBinInd,EndToEndDistance_all,[],@(x) std(x)/numel(x));
        
        % Remove empty bins.
        PathLengths_avg = PathLengths_BinCenters(BinCounts>0);
        EndToEndDistance_avg = EndToEndDistance_avg(BinCounts>0);
        EndToEndDistance_se = EndToEndDistance_se(BinCounts>0);
        
        % Fit the average end-to-end distance.
        FitMaxInd = find(PathLengths_avg >= min(max(PathLengths_avg),PathLengthMax_fit),1);
        
        EndToEndDistance_FitFunc = @(b,x) double(b>0)*(2*b*x - 2*b^2*(1-exp(-x/b)));
        beta_guess = 75; %um
        PathLengths_Fit = PathLengths_avg(1:FitMaxInd);
        EndToEndDistance_avg_Fit = EndToEndDistance_avg(1:FitMaxInd);
        
        Optim_options = optimoptions('lsqcurvefit','Display','off');
        [Lp,~,residual,~,~,~,J] = lsqcurvefit(EndToEndDistance_FitFunc,beta_guess,PathLengths_Fit,EndToEndDistance_avg_Fit,0,inf,Optim_options);
        Lp_cov = mean(residual.^2)*full(inv(J'*J));
        Lp_err = sqrt(Lp_cov);
    
    case 'squared_dist_normed'
        % Fit the normalized end-to-end distance squared for all branches.
        % Use only the endpoints of each branch when assembling data.
        EndToEndDistance = zeros(N_branches,1);
        PathLengths = zeros(N_branches,1);
        for i=1:N_branches
            EndToEndDistance(i) = Tree(i).SquaredDist(end);
            PathLengths(i) = Tree(i).PathLength(end);
        end
        EndToEndDistance_normalized = EndToEndDistance./PathLengths.^2;
        plot(1./PathLengths,EndToEndDistance_normalized,'o')
        beta_guess = 75; %um
        FitMaxInd = numel(PathLengths);
        PathLengths_Fit = PathLengths(1:FitMaxInd);
        EndToEndDistance_Fit = EndToEndDistance(1:FitMaxInd);
        EndToEndDistance_FitFunc = @(b,x) 2*b*x - 2*b^2*(1-exp(-x/b));
        
        [Lp,~,~,Lp_cov] = nlinfit(PathLengths_Fit,EndToEndDistance_Fit,EndToEndDistance_FitFunc,beta_guess);
        Lp_err = sqrt(Lp_cov);
    
    case 'delta_theta'
        %% Fit the change in the tangent vector angle.
        if Nodes_dist_std/Nodes_dist_avg > 0.1
            warning('The spread in the Nodes_dist is high. The persistence length may be inaccurate.');
        end
        
        % Get the distribution of all angular changes.
        Delta_thetas = cell2mat({Tree.Delta_Theta}');
        
        % Calculate the Mahalanobis distance
        d_mahal = mahal(Delta_thetas,Delta_thetas);
        d_mahal_max = 20;
        is_Delta_theta_outlier = d_mahal > d_mahal_max;
        is_Delta_theta_outlier = is_Delta_theta_outlier | abs(Delta_thetas) <= 1e-10;
        Delta_thetas_pruned = Delta_thetas(~is_Delta_theta_outlier);
        
        % Fit the distribution with a gaussian.
        fitted_dist = fitdist(Delta_thetas_pruned,'normal');
        Lp = 2*Nodes_dist_avg/fitted_dist.sigma^2;
        Lp_err = nan;
    case 'curvature'
        %% Curvature fit
        % Reference:
        % Wisanpitayakorn 2018, "Measurement of the Persistence Length of Cytoskeletal Filaments using Curvature Distributions"
        N_decimations = 30;
        delta_s = Nodes_dist_avg*(1:N_decimations)';
        mu = zeros(N_decimations,1);
        
        % Calculate the curvature of each filament at each level of
        % decimation.
        N_branches = numel(Tree);
        Tree_decimated = repmat(Tree(:),[1 N_decimations]);
        for i=1:N_branches
            for j=1:N_decimations
                % Define the tangent vectors.
                Nodes_Pos = Tree(i).PointsPos(1:j:end,:);
                if size(Nodes_Pos,1) > 1
                    Tree_decimated(i,j).TangentVec = diff(Nodes_Pos);
                    
                    % Normalize the tangent vectors.
                    Nodes_dist = sqrt(sum(Tree_decimated(i,j).TangentVec.^2,2));
                    Tree_decimated(i,j).Nodes_dist = Nodes_dist;
                    Tree_decimated(i,j).TangentVec = bsxfun(@rdivide,Tree_decimated(i,j).TangentVec,Nodes_dist);
                    
                    % Calculate the angle of each segment and the change in angle.
                    Tree_decimated(i,j).Theta = atan2(Tree_decimated(i,j).TangentVec(:,2),Tree_decimated(i,j).TangentVec(:,1));
                    Tree_decimated(i,j).Delta_Theta = mod(diff(Tree_decimated(i,j).Theta) + pi,2*pi) - pi;
                    Tree_decimated(i,j).kappa = abs(Tree_decimated(i,j).Delta_Theta)./Tree_decimated(i,j).Nodes_dist(1:end-1,:);
                else
                    Tree_decimated(i,j).kappa = [];
                end
            end
        end
        
        % Define the scale parameters \mu.
        for n=1:N_decimations
            mu(n) = n/((n*2 - 1)*(1 - 1/n)/3 + 1); % n = 2^m
        end
        Fitted_indices = 1:N_decimations;
        delta_s_fit = delta_s(Fitted_indices).*mu(Fitted_indices);
        
        % Fit the branch curvature distribution on bootstrapped datasets.
        N_bootstraps = 100;
        Lambda = zeros(N_decimations,N_bootstraps);
        Bootstrap_branches_indices = randi(N_branches,N_branches,N_bootstraps);
        Lp_fit_bootstrap = zeros(N_bootstraps,1);
        
        kappas_cell = reshape({Tree_decimated(:).kappa},[N_branches,N_decimations]);
        for k = 1:N_bootstraps
            kappas_sample = kappas_cell(Bootstrap_branches_indices(:,k),:);
            for n = 1:N_decimations
                % Fit the curvatures.
                kappas_decimated_sample = cat(1,kappas_sample{:,n});
                Lambda(n,k) = fit_curvature(kappas_decimated_sample,false);
            end
            
            % Fit the curvature inverse variance to find the persistence
            % length.
            Lambda_fit = Lambda(Fitted_indices,k);
            % Calculate the slode of the linear fit.
            Lambda_fit_centered = Lambda_fit - mean(Lambda_fit);
            delta_s_fit_centered = delta_s_fit - mean(delta_s_fit);
            Lp_fit_bootstrap(k) = sum(Lambda_fit_centered.*delta_s_fit_centered)/sum(delta_s_fit_centered.^2);
        end
        
        % Fit the bootstrapped distribution of persistence length to a lognormal
        % distribution and recover the mean of the distribution.
        Lp_lognormal_fit = fitdist(Lp_fit_bootstrap,'Lognormal');
        Lp = Lp_lognormal_fit.mean;
        Lp_err = Lp_lognormal_fit.std;
    otherwise
        error('Method %s is not supported.',options.Method);
end

% Send warning if the fit is poor.
Lp_relative_err = Lp_err/Lp;
if Lp_relative_err > 0.01
    warning('The persistence length fit is poor (Delta L_p/L_p=%.2f%%)',Lp_relative_err*100)
end
%% Plot the empirical end-to-end distance along with the fit.
if options.Plot
    f = figure;
    hold on
    f.Position(3:4)=[800 800];
    a = f.CurrentAxes;
    a.FontSize = 24;
    
    switch options.Method
        case 'tangent_vec'
            % Plot the average inner products.
            plot(PathLengths_BinCenters,InnerProducts_avg,'k.')
            ylim([0,1]);
            
            % Plot the fit of the average inner products.
            Fit_h = plot(PathLengths_Fit,AutoCorrModel(Lp,PathLengths_Fit),'r-','LineWidth',1);
            
            a.YScale = 'log';
            title({'Persistence length estimation using the inner products of tangent vectors.',['N_{branches}=',num2str(N_branches)]})
            xlabel('Path length \Deltas (\mum)')
            ylabel('$\langle t(\Delta s) \cdot t(0) \rangle_{\textrm{Branches}}$','interpreter','latex')
            
            Fit_leg=['$e^{-\frac{\Delta s}{L_p}}$ ($L_p$=',num2str(Lp(1),'%.1f'),'$\mu m$)'];
            xlim([0, 2*PathLengthMax_fit]);
            ylim([0.1,1]);
            
        case 'squared_dist'
            % Plot the calculated squared end-to-end distance.
            plot(PathLengths_avg,EndToEndDistance_avg,'k.')
            
            % Plot the fit to the squared end-to-end distance.
            Fit_h = plot(PathLengths_Fit,EndToEndDistance_FitFunc(Lp,PathLengths_Fit),'r-','LineWidth',1);
            
            title(['Persistence length estimation using the end-to-end distance. N_{branches}=',num2str(N_branches)])
            xlabel('Path length s [\mum]')
            ylabel('$\langle R(s)^2 \rangle $  $[\mu m^2]$','interpreter','latex')
            
            Fit_leg = sprintf('%s\n%s','$\langle R(s)^2 \rangle = 2 L_p (s - L_p (1 - e^{-\frac{s}{L_p}}))$',['$L_{p,Fit} = ',num2str(Lp),' \mu m $']);
        case 'squared_dist_normed'
            plot(PathLengths,EndToEndDistance./(PathLengths_Fit.^2),'.')
            PathLengths_fit = linspace(min(PathLengths),max(PathLengths));
            
            xlabel('Path Length s [\mum]')
            ylabel('$\frac{\langle R(s)^2 \rangle}{s^2} $','interpreter','latex')
            Fit_leg = sprintf('%s\n%s','$\langle R(s)^2 \rangle = 2 L_p (s - L_p (1 - e^{-\frac{s}{L_p}}))$',['$\L_{p,Fit} = ',num2str(Lp),' \mu m $']);
        case 'delta_theta'
            % Plot the distribution of angles.
            histogram(Delta_thetas_pruned,'BinWidth',0.01,'Normalization','pdf');
            
            % Plot the disitribution
            x_dist = linspace(min(Delta_thetas_pruned),max(Delta_thetas_pruned),1000);
            y_dist = pdf(fitted_dist,x_dist);
            Fit_h = plot(x_dist,y_dist);
            
            xlabel('Angle step \Delta\theta (rad)');
            ylabel('Probability Density');
            Fit_leg = ['Gaussian fit, $L_{p,Fit} = ',num2str(Lp),' \mu m $'];
        case 'curvature'
            f.Name = 'Inverse Curvature Variance Fit';
            plot(delta_s_fit,mean(Lambda,2),'o');
            
            x_fit = linspace(min(delta_s_fit),max(delta_s_fit),100)';
            y_fit = Lp*x_fit;
            Fit_h = plot(x_fit,y_fit);
            
            xlim([0,inf]);
            ylim([0,inf]);
            
            xlabel('\mu \Deltas (\mum)');
            ylabel('Curvature inverse variance \Lambda (\mum^2)')
            Fit_leg = ['$L_{p,Fit} = ',num2str(Lp),' \mu m $'];
    end
    hold off
    h = legend(Fit_h,{Fit_leg},'Interpreter','latex','location','southwest');
    h.FontSize = 24;
    drawnow
    if exist('tightax','file')
        tightax;
    end
    
    % Plot extra plots.
    if strcmp(options.Method,'curvature') && N_bootstraps > 1
        figure('Name','Persistence Length Bootstrap Distribution');
        histogram(Lp_fit_bootstrap,'Normalization','pdf','HandleVisibility','off')
        hold on
        beta_lognormal_fit_x = linspace(min(Lp_fit_bootstrap),max(Lp_fit_bootstrap),100);
        beta_lognormal_fit_y = pdf(Lp_lognormal_fit,beta_lognormal_fit_x);
        plot(beta_lognormal_fit_x,beta_lognormal_fit_y,'-','DisplayName',sprintf('Lognormal fit (mean = %.2f \\mum)',exp(Lp_lognormal_fit.mu + Lp_lognormal_fit.sigma^2/2)))
        xlabel('Persistence length (\mum)');
        ylabel('Probability density (\mum)');
        legend;
        
        % Plot an example of the curvature distribution.
        kappas = cell2mat({Tree_decimated(Bootstrap_branches_indices(:,1),1).kappa}');
        fit_curvature(kappas,true);
    end
end
end
%% Function to fit the curvation when using Method='curvature'.
function Lambda = fit_curvature(kappa,Plot)
% Input
% kappa = sample of curvatures to be fitted
%
% Output
% Lambda = inverse-variance of the curvature sample.

if nargin == 1
    Plot = false;
end

% Fit all curvatures
if Plot
    Gaussian_fit = prob.HalfNormalDistribution.fit(kappa(:));
    curvature_var = Gaussian_fit.sigma^2;
else
    curvature_var = sum(kappa.^2)/numel(kappa);
end
Lambda = 1/curvature_var;

% Plot for check
if Plot
    figure('Name','Curvature Distribution');
    histogram(kappa,'Normalization','Pdf','HandleVisibility','off')
    hold on;
    
    pdf_x = linspace(min(kappa),max(kappa),100);
    pdf_y = pdf(Gaussian_fit,pdf_x);
    plot(pdf_x,pdf_y,'DisplayName',sprintf('Gaussian Fit:\\mu = %.2f, \\sigma = %.2f',Gaussian_fit.mu,Gaussian_fit.sigma))
    
    legend;
    xlabel('Curvature \kappa (rad/\mum)');
    ylabel('Probability density P(\kappa)');
end
end