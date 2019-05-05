load('IRIS_DMC_with_FMS_and_energy.mat')
CAT = iris_dmc_cat_with_fms_and_energy;

inc     = 2;
years   = 1980:inc:2019;
nYr     = length(years);
bAry    = zeros(nYr,1);

for n = 1:(nYr-1)
    dateNumRange = datenum(years(n:n+1),[1 1], [1 1]);
    I = CAT.time > dateNumRange(1) & CAT.time < dateNumRange(2);
    yrCat = CAT.M(I);
    bAry(n) = calc_McMaxCurvature(yrCat);
end

figure
plot(years,bAry);

function [fMc] = calc_McMaxCurvature(mCatalog)
  % Get maximum and minimum magnitudes of the catalog
  fMaxMagnitude = max(mCatalog);
  fMinMagnitude = min(mCatalog);
  if fMinMagnitude > 0
    fMinMagnitude = 0;
  end
  
  % Number of magnitudes units
  nNumberMagnitudes = (fMaxMagnitude*10) + 1;
  
  % Create a histogram over magnitudes
  [vHist, vMagBins] = hist(mCatalog, fMinMagnitude:0.1:fMaxMagnitude);
  
  % Get the points with highest number of events -> maximum curvature  
  fMc = vMagBins(max(find(vHist == max(vHist))));
  if isempty(fMc)
    fMc = nan;
  end

end