load analysis_catalog.mat
% Basic variables
% All variables are presented in SI units if not otherwise indicated

M                   = CAT.MSmag;
Mo                  = 10.^((M+10.7)*3/2)/10^7;
depth               = CAT.MSdepth*1000;
productivity        = CAT.MSprod;
residualProductivity= CAT.MSres;
MeBB                = CAT.MeBB_appended_cat1;
MeHF                = CAT.MeHF_appended_cat1;
meanSlip            = CAT.av_slip_appended_cat2;
faultWidth          = CAT.rwid_appended_cat2*1000;
faultLength         = CAT.rlen_appended_cat2*1000;
faultPerimeter      = CAT.Perimeter_appended_cat2;
stressDrop          = CAT.StressDrop_appended_cat2*10^6;
heterogeneity       = CAT.EllipseErrorNorm_appended_cat2;

% normalized/linearized parameters:

logSD               = log10(stressDrop);
MeBBNorm            = MeBB./M;
MeHFNorm            = MeHF./M;
depthNorm           = abs(depth-mean(depth))/mean(depth);
Ar                  = faultWidth./faultLength;
Ar(Ar<1)            = 1./Ar(Ar<1);
logAr               = log10(Ar);
PNorm               = faultPerimeter./(2*pi*(Mo*7/16./stressDrop).^(1/3));
H                   = heterogeneity;

% store variables into X
X = [logSD, MeBBNorm, MeHFNorm, depthNorm, logAr, PNorm, H];
Y = residualProductivity;
badI = any([isnan(X),isinf(X),isnan(Y),isinf(Y)],2);
X = X(~badI,:);
Y = Y(~badI);
SX= size(X);
X = X - repmat(mean(X,1),SX(1),1); % remove average value
Y = Y - mean(Y);

BETA = (X'*X)\X'*Y;

% remove outlier (rinse and repeat)
var = (X*BETA - Y).^2;
outlierInd = var > 3*mean(var);
X = X(~outlierInd,:);
Y = Y(~outlierInd);
SX= size(X);
X = X - repmat(mean(X,1),SX(1),1); % remove average value
Y = Y - mean(Y);


BETA = (X'*X)\X'*Y;

% comput R^2
R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2)
R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2)




