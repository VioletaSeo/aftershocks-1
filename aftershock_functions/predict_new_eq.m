function [res,N] = predict_new_eq(M,logAnorm, Dip, Age)

load('SVMdata.mat');
rng(1)

trainingPredictors = pred;
trainingResponse = resp;
    
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

res = predict(regressionSVM,[Age,Dip,logAnorm]);

N = 10^res*6.7338e-7*10^(1.0412*M);
    
end