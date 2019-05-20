function xval_knn_prediction(lat,lon,depth,res)

%%
wgs84   = wgs84Ellipsoid('kilometer');
loc     = geodetic2ecef(wgs84,lat,lon,-depth);

%% hold data out for validation
indices = crossvalind('Kfold',res,5);
hold on
for i = 1:5
    test = (indices == i); 
    train = ~test;
    resHat = knn_regress(loc(train,:),res(train),loc(test,:),5);
    scatter(res(test),resHat,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.5);
end

mm = minmax(res(~isinf(res))');
plot(mm,mm,'--','Color','k','LineWidth', 2);
xlabel('Relative productivity')
ylabel('Prediction')

end 

function Yq = knn_regress(Xt,Yt,Xq,k)
Mdl = KDTreeSearcher(Xt);
I  = knnsearch(Mdl,Xq,'K',k);
nQ = length(Xq);
Yq = zeros(nQ,1);
for n = 1:nQ
    Yq(n) = median(Yt(I(n,:)));
end

end