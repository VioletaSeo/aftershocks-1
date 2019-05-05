load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;
minMag  = 6.2;
% completensess = 4.3;
% magRange = [completeness,10];
maxDepth = 55;
[ASinfo,k,alpha] = aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'MinMainshockMag',minMag, ...
    'DepthRange',[0,maxDepth], ...
    'ReturnCatalog', 'yes', ...
    'SaveCatalog', 'no', ...
    'PlotYN','no', ...
    'Completeness',5);
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;

%%
wgs84   = wgs84Ellipsoid('kilometer');
loc     = geodetic2ecef(wgs84,MSCat.lat,MSCat.lon,-MSCat.depth);
res     = MSCat.MSres;

%% hold data out for validation
indices = crossvalind('Kfold',res,5);

figure; hold on
for i = 1:10
    test = (indices == i); 
    train = ~test;
    resHat = knn_regress(loc(train,:),res(train),loc(test,:),15);
    scatter(res(test),resHat,'k','filled','MarkerFaceAlpha',0.1);
end
axis equal
xlim([-1.5,1.5])
ylim([-1.5,1.5])
mm = minmax(res(~isinf(res))');
plot(mm,mm,'--','Color',[0 0 0 0.3],'LineWidth', 3);

%% temporal hold out

%% 
function Yq = knn_regress(Xt,Yt,Xq,k)
Mdl = KDTreeSearcher(Xt);
I  = knnsearch(Mdl,Xq,'K',k);
nQ = length(Xq);
Yq = zeros(nQ,1);
for n = 1:nQ
    Yq(n) = median(Yt(I(n,:)));
end

end