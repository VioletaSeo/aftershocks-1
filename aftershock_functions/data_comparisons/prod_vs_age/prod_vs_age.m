tempCatFN= 'temp_earthquake_cat_merge.mat';
minMag = 6.5;
delete_temp         = @() system(sprintf('rm %s',tempCatFN));
load('IRIS_DMC_with_FMS.mat')
aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'SaveCatalog',tempCatFN, ...
    'MinMainshockMag',minMag, ...
    'DepthRange',[0,55], ...
    'PlotYN', 'no');
clear CAT
CAT = load(tempCatFN); delete_temp();
MSage = get_crust_age(CAT.MSlat,CAT.MSlon,CAT.MSdepth);

%% seafloor age

figure;
subplot(2,3,1:2); title(sprintf('Mainshocks MW%g +',minMag)); hold on
numBin      = 30;
edges       = linspace(0,max(MSage),numBin);
yyaxis right; ylim([0,30]); ylabel('No aftershocks (N)')
I = isinf(CAT.MSres);
histogram(MSage(I),edges,'EdgeColor','none');


% H           = cell(3,1);
% yyaxis right
% 
% colors = {[     0    0.4470    0.7410],...
%           [0.6350    0.0780    0.1840],...
%           [     0    1         0     ]};
% for n = 1:3
%     I = CAT.MSfms == n;
%     H{n} = histogram(MSage(I),edges,'FaceColor',colors{n});
% end
% legend([H{:}],{'Strike-slip','Normal','Reverse'})
% set(gca,'YScale','log')

yyaxis left; ylim([-3,2]);
scatter(MSage,CAT.MSres,25,'filled','MarkerFaceAlpha',0.5);
xlabel('Sea floor age')
ylabel('Relative Productivity')
[Mage,Mres,Merr] = mov_mean(MSage,CAT.MSres,150,30);
plot(Mage,Mres,'-','LineWidth',3,'Color',[1 0 0 0.5])
plot(Mage,Merr(:,1),'--','LineWidth',1,'Color',[1 0 0 0.5])
plot(Mage,Merr(:,2),'--','LineWidth',1,'Color',[1 0 0 0.5])
grid on

% sqrt age
subplot(2,3,3)
scatter(sqrt(MSage),CAT.MSres,30,'filled','MarkerFaceAlpha',0.3);
xlabel('Square root of Sea floor age')
ylabel('Relative Productivity')
hold on
[sqrtMage,Mres,Merr] = mov_mean(sqrt(MSage),CAT.MSres,200,sqrt(10));
plot(sqrtMage,Mres,'LineWidth',3,'Color',[1 0 0 0.5] )
plot(sqrtMage,Merr(:,1),'--','LineWidth',1,'Color',[1 0 0 0.5])
plot(sqrtMage,Merr(:,2),'--','LineWidth',1,'Color',[1 0 0 0.5])
grid on


% by focal mechanism
titleStr = {'Strike-slip','normal','reverse'};
for n = 1:3
subplot(2,3,n+3);hold on;
title(titleStr{n})

yyaxis left; ylim([-3,2]);
I = CAT.MSfms ==  n;
SSage = MSage(I);
SSres = CAT.MSres(I);
scatter(SSage, SSres,'Filled','MarkerFaceAlpha',0.5);
[SSMage,SSMres,SSMerr] = mov_mean(SSage,SSres,150,30);
plot(SSMage,SSMres,'-','LineWidth',3,'Color',[1 0 0 0.5])
plot(SSMage,SSMerr(:,1),'--','LineWidth',1,'Color',[1 0 0 0.5])
plot(SSMage,SSMerr(:,2),'--','LineWidth',1,'Color',[1 0 0 0.5])
xlabel('Sea floor age')
ylabel('N^*')

numBin      = 30;
edges       = linspace(0,max(MSage),numBin);
yyaxis right; ylim([0,30]); ylabel('No aftershocks (N)')
Iinf = isinf(CAT.MSres(I));
histogram(MSage(Iinf),edges,'EdgeColor','none');
grid on

end





%%
function CRUSTAGE = get_crust_age(lat,lon,depth)


load('age.3.6.xyz'); 
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
