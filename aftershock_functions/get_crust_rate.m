function CRUSTRATE = get_crust_rate(lat,lon,depth)


load('rate.mat'); 
rateAry = rate_3_6(~isnan(rate_3_6(:,3)),:);
wgs84 = wgs84Ellipsoid('kilometers'); 
[Xa,Ya,Za] = geodetic2ecef(wgs84,rateAry(:,2),rateAry(:,1),zeros(length(rateAry),1));
[Xe,Ye,Ze] = geodetic2ecef(wgs84,lat,lon,-depth);
age = rateAry(:,3);



%%
nMS = length(lat);
CRUSTRATE = zeros(nMS,1);
Dcutoff = 30; %km

for n = 1:nMS
    D = sqrt((Xe(n)-Xa).^2 + (Ye(n)-Ya).^2 + (Ze(n)-Za).^2);
    minD = min(D);
    if minD < (Dcutoff)
        I = D == minD;
        CRUSTRATE(n) = min(age(I));
    else
        CRUSTRATE(n) = nan;
    end
end

end
