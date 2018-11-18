function productivity_analysis(catalog,res)
% This function is built to analize mainshock productivity in light of
% source parameters

%% 'standardize catalogs'

switch nargin 
    case 0
        [cat,res] = make_catalog('analysis_catalog.mat');
    case 2
        cat = catalog; % this catalog needs to be in the right format (as saved from make_cat)
    otherwise 
        disp('Something is wrong here')
end

%% multivariate analysis
% prep data:
    % 1) choose relevant variables  
    % 2) log transform data
    % 3) seperate data into predictors (X) and an observable (Y)
    
    
    % remove rows with earhtquakes that don't have aftershocks, nan entries
    % or inf
    
    [cleanCat] = cat(~any(isinf(cat.Variables) | isnan(cat.Variables),2) & ~isinf(res), :); 
    res        = res(~any(isinf(cat.Variables) | isnan(cat.Variables),2) & ~isinf(res));
    
% create covariance scatter matrix
%plot_scattermat(cleanCat,fms);
plot_scattermat(cleanCat,res);

% analyse
% pca_analysis(X);

end

function varargout = plot_scattermat(inputTable, Y)

h           = figure;
varNames    = inputTable.Properties.VariableNames;
nVarNames   = length(varNames);
X           = inputTable.Variables;

category = zeros(length(Y),1);
category(Y>=prctile(Y,10) & Y<(prctile(Y,90))) = 1;
category(Y>=prctile(Y,90)) = 2;


figure
gplotmatrix(X,[],category,['r' 'g' 'b'],[],[],false)
text(linspace(0,1, nVarNames), repmat(-.1,1,nVarNames), varNames, 'FontSize',8);
text(repmat(-.12,1,nVarNames), fliplr(linspace(0,1, nVarNames)), varNames, 'FontSize',8, 'Rotation',90);

figure
parallelcoords(X,'group',category,'standardize','on', 'labels',varNames,'quantile',.25)
%parallelcoords(X,'group',category,'standardize','on', 'labels',varNames)

varargout{1} = h;
end

function [Zx,PC] = pca_analysis(X)
% perform a principle component analysis of the input matrix X as a
% predictor of Y
dimX = size(X);
Z = zeros(dimX(1),dimX(2));

    for n = 1:dimX(2)
        centeredX       = X(:,n) - mean(X(:,n));    % center data
        standardizedX   = centeredX/std(centeredX); % standardize data
        Z(:,n)          = standardizedX;            % store
    end
    
ZTZ = Z'*Z;
[P,lambda] = eig(ZTZ);
[PC,I] = sort(diag(lambda),'descend');
Px = P(:,I);
Zx = Z*Px;

figure
plot(PC/sum(PC))
hold on
plot(cumsum(PC)/sum(PC))

    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    