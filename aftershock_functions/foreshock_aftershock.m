maxDepth = 55; % km

figure
[prefactor,alpha, FSprefactor,FSalpha] = ...
    aftershock_productivity_kernel( 'DepthRange',           [0,maxDepth],...
                                    'ForeshockTimeWindow',  450);                % days

MScat       = load('mainshock_info.mat');
MScat       = struct2table(MScat);
allAS       = @(x,y,prefactor,exponent) log10(y) - log10(prefactor*10.^(exponent*x));
MScat.ASres = allAS(MScat.MSmag, MScat.MSprod, prefactor,     alpha);
MScat.FSres = allAS(MScat.MSmag, MScat.FSprod, FSprefactor,   FSalpha);

MTH       	= 7.5;
MScatTH     = MScat(MScat.MSmag > MTH, :);

figure
scatter(MScatTH.FSres, MScatTH.ASres, ...
       (MScatTH.MSmag-4).^3.5,MScatTH.MSfms, ...
        'Filled','MarkerFaceAlpha',0.5)
colormap jet
xlabel('Foreshock Abundance log(N/{k_010^{\alpha M}})')
ylabel('Aftershock Abundance')

hold on

% best fit
data    = [MScatTH.FSres, MScatTH.ASres];
data    = data(~isinf(prod(data,2)),:); % remove pesky -Inf
coef = polyfit(data(:,1), data(:,2), 1);
plot(minmax(data(:,1)),coef(1)*minmax(data(:,1))+coef(2), ...
    'LineWidth', 1, ...
    'Color', 'k');


pVal    = calcSignificance(data,coef);
disp(sprintf('P-Value = %g', pVal))


%% aggregate data
perWin = 30;



% old way better correlation but probably wrong
% % get foreshocks
% figure
% aftershock_productivity_kernel('DepthRange', [0,maxDepth],'TimeWindow',-900/2,'TimeSelectionWindow',-900/2)
% MScat           = load('mainshock_info.mat');
% MScat           = struct2table(MScat);
% MScat.MSres     = getMSres(MScat.MSmag,MScat.MSprod);
% foreshockCat    = MScat;
% 
% % get aftershocks
% figure
% aftershock_productivity_kernel('DepthRange', [0,maxDepth],'TimeWindow',900,'TimeSelectionWindow',900)
% MScat           = load('mainshock_info.mat');
% MScat           = struct2table(MScat);
% MScat.MSres     = getMSres(MScat.MSmag,MScat.MSprod);
% aftershockCat   = MScat;
% 
% MSformat        = {'MSt','MSlat','MSlon','MSdepth','MSmag'};
% 
% [mergedCat,I]   = merge_eq_catalog(aftershockCat, MSformat, ...
%                                    foreshockCat,  MSformat);
%                                
% MTH             = 6;
% mergedCatTH       = mergedCat(mergedCat.MSmag>MTH,:);
% 
% figure
% scatter(mergedCatTH.MSres_appended_cat1, mergedCatTH.MSres, ...
%         (mergedCatTH.MSmag-4).^3.5,mergedCatTH.MSfms, ...
%         'Filled','MarkerFaceAlpha',0.5)
% colormap jet
% xlabel('Foreshock Abundance log(N/{k_010^{\alpha M}})')
% ylabel('Aftershock Abundance')
%                              
% function MSRES = getMSres(MSMAG,MSPROD)
% h = findobj(gca,'Type','scatter');
% xData=get(h,'Xdata');
% yData=get(h,'Ydata');
% 
% if iscell(xData)
%     MAG     = [xData{:}];
%     NUMEQ 	= [yData{:}];
% else
%     MAG     = xData;
%     NUMEQ   = yData;
% end
% 
% [prefactor,exponent] = getproductivity(MAG,NUMEQ);
% 
% allAS       = @(x,y) log10(y) - log10(prefactor*10.^(exponent*x));
% MSRES       = allAS(MSMAG,MSPROD);
% 
% end
