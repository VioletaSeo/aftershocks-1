function CRUSTAGE = get_crust_age(lat,lon,depth)


load('age.mat'); 
ageAry = age_3_6(~isnan(age_3_6(:,3)),:);
wgs84 = wgs84Ellipsoid('kilometers'); 
[Xa,Ya,Za] = geodetic2ecef(wgs84,ageAry(:,2),ageAry(:,1),zeros(length(ageAry),1));
[Xe,Ye,Ze] = geodetic2ecef(wgs84,lat,lon,-depth);
age = ageAry(:,3);



%%
nMS = length(lat);
CRUSTAGE = zeros(nMS,1);
Dcutoff = 30; %km

for n = 1:nMS
    D = sqrt((Xe(n)-Xa).^2 + (Ye(n)-Ya).^2 + (Ze(n)-Za).^2);
    minD = min(D);
    if minD < (Dcutoff)
        I = D == minD;
        CRUSTAGE(n) = min(age(I));
    else
        CRUSTAGE(n) = nan;
    end
end

end
