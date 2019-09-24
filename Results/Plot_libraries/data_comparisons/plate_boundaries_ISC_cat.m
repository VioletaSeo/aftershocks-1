% compare all tectonic evironments using just ISC (no focal mechanism)

load('ISC_1976-2018.mat'); % loads quakes 
tempCatFN           = 'temp_earthquake_cat_merge.mat';
minM = 6;

%% 1) compute productivity catalogs
defaultInpStr = {'ForeshockTimeWindow', 10      , ...
                 'SaveCatalog',tempCatFN,...
                 'DepthRange',          [0,55]  , ...
                 'ShowOverviewYN',      'no', ...
                 'PlotYN', 'no', ...
                 'MinMainshockMag',     6.5}; % default input
PBClassArray  = {'all','OSR','OTF','OCB','CRB','CTF','CCB','SUB'};
CAT	          = [];

delete_temp         = @() system(sprintf('rm %s',tempCatFN));
for iPB = PBClassArray
    PB = iPB{:};
    aftershock_productivity_kernel(quakes.otime, ...
        quakes.lat                  , ...
        quakes.lon                  , ...
        quakes.depth                , ...
        quakes.mag                  , ...
        defaultInpStr{:}            , ...
        'PlateBoundaryClass',PB     , ...
        'PlateBoundaryDist',200      );
    CAT.(PB)  = load(tempCatFN);
    delete_temp;
end

%% 2) plot each output
PBClassSubset = PBClassArray;
plot_prod(CAT,PBClassSubset);

%%
plot_res_prod_comparison(CAT,PBClassSubset,minM+1);

%% *accessory functions
function plot_prod(CATSTRUCT,FIELDS)
figure; hold on;
for iField = FIELDS
    scatter(CATSTRUCT.(iField{:}).MSmag, CATSTRUCT.(iField{:}).MSprod,'filled','MarkerFaceAlpha',0.5);
end
legend(FIELDS)
set(gca,'YScale','log')    
end

function plot_res_prod_comparison(CATSTRUCT,FIELDS,Mc)

figure; hold on;
colors     = get(gca, 'ColorOrder');
sz         = 50;
N          = 0;
numNegInf     = zeros(length(FIELDS),1);
tickName      = cell(length(FIELDS),1);

for iField = FIELDS
    N   = N+1;
    mag = CATSTRUCT.(iField{:}).MSmag;
    Im  = mag>Mc;
    res = CATSTRUCT.(iField{:}).MSres(Im);
    numNegInf(N) = sum(isinf(res))/length(res);
    scatter(res,        N*ones(size(res)),  sz*2,'filled','MarkerFaceColor',colors(mod(N,7)+1,:), 'MarkerFaceAlpha',0.01)
    scatter(median(res),N,                 sz,  'filled','MarkerFaceColor',colors(mod(N,7)+1,:))
    plot(prctile(res,[40,60]),[N,N],'LineWidth',2,'Color',colors(mod(N,7)+1,:))
    tickName{N} = sprintf('%s, N0 = %g', iField{:}, numNegInf(N));
end

set(gca,'YTick',1:N,'YTickLabel',tickName','YLim',[0.1,N + 0.1])
xlabel('Residual productivity')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

end