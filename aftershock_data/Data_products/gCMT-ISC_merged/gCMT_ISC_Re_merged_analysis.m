%% radiated energy
load('ISC-gCTM_mainshock_catalog.mat')
    MSformat = {'MSt','MSlat','MSlon','MSmag'};
load('EQ_Energy_cat.mat');                  
    REformat        = {'t','lat','lon','mag'};    
load('agu_2018_source_attribute_data.mat');

[mergedCat,~] = merge_eq_catalog(CAT,MSformat, ...
                                 cat, REformat);
X = mergedCat.MeBB_appended_cat1./mergedCat.mag_appended_cat1;

                                                    
% X = [log10(mergedCat.EBB_appended_cat1), ...
%      log10(mergedCat.EHF_appended_cat1), ...
%      log10(mergedCat.TBB_appended_cat1), ...
%      mergedCat.MSdepth, ...
%      mergedCat.MeBB_appended_cat1, ...
%      mergedCat.MeHF_appended_cat1];


Y = mergedCat.MSres_appended_cat1;
SourceParameters = {'Scaled Energy Moment'};
%SourceParameters = {'ER_{BB}', 'ER_{HF}','T_{BB}','Depth','Me_{BB}','ME_{HF}'};                             
Rval = r2stem(X,Y,SourceParameters,10000);

%% finite fault inversions

load('this_is_the_catalog_I_am_using_god_dammit2.mat')
c = aguCAT;
% X = [c.Width, ...
%      c.Length,...
%      log10(c.AspectRatio),....
%      c.RuptureVelocity,...
%      log10(c.RuptureTime),...
%      log10(c.Heterogeneity),...
%      c.Young,...
%      c.Poisson,...
%      log10(c.StressDrop), ...
%      c.MSdepth_appended_cat1];
 X = [sqrt(c.RuptureTime)./c.Mo, ...
     sqrt(c.Width)./c.Mo, ...
     sqrt(c.Length)./c.Mo, ...
     c.AspectRatio,....
     c.RuptureVelocity,...
     c.Heterogeneity,...
     (c.StressDrop), ...
     abs(c.Depth-mean(c.Depth))];

Y = c.MSres_appended_cat1;
SourceParameters = {'t*','W*','L*','A_r','V_r','H','\Delta\sigma','\Delta d'};
r2stem(X,Y,SourceParameters,10000);

%% just depth
load('ISC-gCTM_mainshock_catalog.mat')
X = CAT.MSdepth;
Y = CAT.MSres_appended_cat1;
SourceParameters = {'depth'};
r2stem(X,Y,SourceParameters,10000);


%% 
load('EQ_Energy_cat.mat');                  
    REformat        = {'t','lat','lon','mag'};    
load('agu_2018_source_attribute_data.mat');
c = aguCAT;
    SDformat        = {'Time','Lat','Lon','Mw'};

[mergedCat,~] = merge_eq_catalog(c, SDformat,cat, REformat);

%% Where shit gets done ...
function [Rvalues] = r2stem(X,Y,SourceParameters,mit)    

% store variables into X
badI = any([isnan(X),isinf(X),isnan(Y),isinf(Y)],2);
X = X(~badI,:);
Y = Y(~badI);
SX= size(X);
X = X-mean(X,1);
X = X./std(X);
Y = Y - mean(Y);
Y = Y./std(Y); % ?

%% whole modelRrand
BETA = (X'*X)\X'*Y;
R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2);
R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2);
% remove outlier (rinse and repeat)
var = (X*BETA - Y).^2;
outlierInd = var > 4*mean(var);
X = X(~outlierInd,:);
Y = Y(~outlierInd);
SX= size(X);
X = X - repmat(mean(X,1),SX(1),1); % remove average value
Y = Y - mean(Y);
BETA = (X'*X)\X'*Y;
R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2);
R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2);

%% Random whole model
 for m = 1:mit
        Yshfl  = Y(randperm(length(Y)));
        BETArand= (X'*X)\X'*Yshfl;
        Rrand(m)= 1  -  sum((X*BETArand - Yshfl).^2)/sum(Y.^2);
 end
 

 %% individual models
for n = 1:SX(2)
    Xcol    = X(:,n);
    BETA    = (Xcol'*Xcol)\Xcol'*Y;
    R(n)    = 1  -  sum((Xcol*BETA - Y).^2)/sum(Y.^2);
end
Rvalues = R; %table(SourceParameters',R2');

%% random shuffled data
maxRshuffle = zeros(mit,1);
for m = 1:mit
    Rshuffle = zeros(mit,1);
    for n = 1:SX(2)
        Xcol    = X(:,n);
        Yshuffle= Y(randperm(length(Y)));
        BETAshuffle= (Xcol'*Xcol)\Xcol'*Yshuffle;
        Rshuffle(n)= 1  -  sum((Xcol*BETAshuffle - Yshuffle).^2)/sum(Yshuffle.^2);
    end
    maxRshuffle(m) = max(Rshuffle);
end

%% plots

figure 
histogram(maxRshuffle,'FaceColor', [0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on


xlabel('Variance Reduction')
ylabel('N')
set(gca,'Xlim',[0,0.3])
legend({'V.R. for shuffled data', 'Source parameters'},'Location','northeast')
ah = gca;
yLim = ah.YLim;
recX = linspace(0,max(maxRshuffle),20);
recX = recX(2:end);
for iDim = recX
    rectangle('Position',[0,0,iDim,yLim(2)],'FaceColor',[1 0 0 0.005],'EdgeColor','none')
end
stem(R,3*sqrt(mit)*ones(length(R), 1), ...
    'Color', [0.4660    0.6740    0.1880], ...
    'LineWidth',2)
texth = text(R, 3.3*sqrt(mit)*ones(length(R),1), SourceParameters);
set(texth, 'Rotation', 45)
set(gca,'YLim',yLim)
set(findall(gcf,'-property','FontSize'),'FontSize',12)

end





































