load('IRIS_DMC_with_FMS_and_energy.mat')
CAT = iris_dmc_cat_with_fms_and_energy;
xlat = -20;
lonR = [-80,-50];
N = 200;
lons = linspace(lonR(1),lonR(2),N);

figure; plot(lons, -get_litho(repmat(xlat,1,N),lons));
hold on
I = CAT.lat>(xlat-0.5) & CAT.lat<(xlat+0.5);
scatter(CAT.lon(I),-CAT.depth(I),'.','MarkerEdgeAlpha',0.05)
xlim(lonR)