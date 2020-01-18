function [validationPredictions, validationRMSE,regressionSVM] = train_SVM(predictors,response)
% Train an SVM model for aftershock forecasting. 

% Why use a random seed? - From matlab documentation:
% Kernel scale parameter, specified as the comma-separated pair consisting
% of 'KernelScale' and 'auto' or a positive scalar. The software divides
% all elements of the predictor matrix X by the value of KernelScale. Then,
% the software applies the appropriate kernel norm to compute the Gram
% matrix.

%If you specify 'auto', then the software selects an appropriate scale
%factor using a heuristic procedure. This heuristic procedure uses
%subsampling, so estimates can vary from one call to another. Therefore, to
%reproduce results, set a random number seed using rng before training.
rng(2)

% Perform cross-validation
KFolds = length(response);
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;

for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    responseScale = iqr(trainingResponse);
    if ~isfinite(responseScale) || responseScale == 0.0
        responseScale = 1.0;
    end
    boxConstraint = responseScale/1.349;
    epsilon = responseScale/13.49;
    regressionSVM = fitrsvm(...
        trainingPredictors, ...
        trainingResponse, ...
        'KernelFunction', 'polynomial', ...
        'PolynomialOrder',2, ...
        'KernelScale', 'auto', ...
        'BoxConstraint', boxConstraint, ...
        'Epsilon', epsilon, ...
        'Standardize', true);
    
    % Create the result struct with predict function
    svmPredictFcn = @(x) predict(regressionSVM, x);
    validationPredictFcn = @(x) svmPredictFcn(x);
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));


end
