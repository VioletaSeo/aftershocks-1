function CAT = plate_boundaries(varargin)

% compare all tectonic evironments:

%% 1) compute productivity catalogs
defaultInpStr = {'ForeshockTimeWindow', 10      , ...
                 'DepthRange',          [0,55]  , ...
                 'ShowOverviewYN',      'no'}; % default input
PBClassArray  = {'all','OSR','OTF','OCB','CRB','CTF','CCB','SUB'};
FMSArray      = {'all',2    ,1    ,3    ,2    ,1    ,3    ,3    };
CAT	          = [];

for iPB = PBClassArray
    PB = iPB{:};
    CAT.(PB) = aftershock_productivity_GCMT_ISC(defaultInpStr{:},  ...
                                                'PlateBoundaryClass',PB, ...
                                                'PlateBoundaryDist',200);
end

%% 2) plot each output
if strcmp(varargin{1},'yes')
PBClassSubset = PBClassArray;
plot_prob(CAT,PBClassSubset);

%%
plot_res_prod_comparison(CAT,PBClassSubset,6,FMSArray);
end
end

%% *accessory functions
function plot_prob(CATSTRUCT,FIELDS)
figure; hold on;
for iField = FIELDS
    scatter(CATSTRUCT.(iField{:}).MSmag, CATSTRUCT.(iField{:}).MSprod_appended_cat1,'filled','MarkerFaceAlpha',0.5);
end
legend(FIELDS)
set(gca,'YScale','log')    
end

function plot_res_prod_comparison(CATSTRUCT,FIELDS,Mc,fmsChoice)
if length(fmsChoice) == 1; fmsChoice = repmat(fmsChoice,1,length(fieldnames)); end

figure; hold on;
colors     = get(gca, 'ColorOrder');
sz         = 50;
N          = 0;
numNegInf     = zeros(length(FIELDS),1);
tickName      = cell(length(FIELDS),1);

for iField = FIELDS
    N   = N+1;
    mag = CATSTRUCT.(iField{:}).MSmag;
    fms = CATSTRUCT.(iField{:}).MSfms;
    if strcmp(fmsChoice(N),'all'); Ifms = ones(size(fms)); else; Ifms = fms == fmsChoice{N}; end
    Im  = mag>Mc;
    mag = mag(Im & Ifms);
    res = CATSTRUCT.(iField{:}).MSres_appended_cat1(Im & Ifms);
    prod= CATSTRUCT.(iField{:}).MSprod_appended_cat1(Im & Ifms);
    [kMo,alpha] = getproductivity(mag,prod);
    
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