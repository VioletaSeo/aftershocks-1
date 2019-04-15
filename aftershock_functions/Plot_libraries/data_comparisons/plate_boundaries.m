function PBCAT = plate_boundaries(varargin)

% compare all tectonic evironments:
load('IRIS_DMC_with_FMS_and_energy.mat');
CAT     = iris_dmc_cat_with_fms_and_energy;
minMag = 6.5;
maxDepth = 55;
%% 1) compute productivity catalogs
defaultInpStr = {CAT.time, ...
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
    'Completeness',4.5}; % default input
PBClassArray  = {'all','OSR','OTF','OCB','CRB','CTF','CCB','SUB'};
FMSArray      = {'all',2    ,1    ,3    ,2    ,1    ,3    ,3    };
DISTArray     = [nan, 40    ,100   ,300  ,300  ,30   ,400  ,400];

PBCAT	          = [];

for iPB = 1:length(PBClassArray)
    PB = PBClassArray{iPB};
    ASinfo = aftershock_productivity_kernel(defaultInpStr{:},  ...
        'PlateBoundaryClass',PB, ...
        'PlateBoundaryDist',DISTArray(iPB));
    
    PBCAT.(PB) = CAT(ASinfo.ID,:);
    PBCAT.(PB).MSres = ASinfo.MSres;
    PBCAT.(PB).MSprod= ASinfo.MSprod;
    
end

PB = 'none';
ASinfo = aftershock_productivity_kernel(defaultInpStr{:},  ...
        'PlateBoundary','none', ...
        'PlateBoundaryDist',500);
PBCAT.(PB) = CAT(ASinfo.ID,:);
PBCAT.(PB).MSres = ASinfo.MSres;
PBCAT.(PB).MSprod= ASinfo.MSprod;

%% 2) plot each output
PBClassSubset = PBClassArray;
% plot_prob(CAT,PBClassSubset);

%%

plot_res_prod_comparison(PBCAT,[PBClassSubset,{'none'}],6,[FMSArray,{'all'}]);

end

%% *accessory functions
function plot_prob(CATSTRUCT,FIELDS)
figure; hold on;
for iField = FIELDS
    scatter(CATSTRUCT.(iField{:}).M, CATSTRUCT.(iField{:}).MSprod,'filled','MarkerFaceAlpha',0.5);
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
    mag = CATSTRUCT.(iField{:}).M;
    fms = CATSTRUCT.(iField{:}).fms;
    if strcmp(fmsChoice(N),'all'); Ifms = ones(size(fms)); else; Ifms = fms == fmsChoice{N}; end
    Im  = mag>Mc;
    res = CATSTRUCT.(iField{:}).MSres(Im & Ifms);  
    numNegInf(N) = sum(isinf(res))/length(res);
    scatter(res,        N*ones(size(res)),  sz*2,'filled','MarkerFaceColor',colors(mod(N,7)+1,:), 'MarkerFaceAlpha',0.1)
    scatter(median(res),N,                  sz,  'filled','MarkerFaceColor',colors(mod(N,7)+1,:)/1.4)
    plot([prctile(res,33),prctile(res,66)],[N,N],'LineWidth',2,'Color',colors(mod(N,7)+1,:)/1.4)
    tickName{N} = sprintf('%s', iField{:});
end

set(gca,'YTick',1:N,'YTickLabel',tickName','YLim',[0.1,N + 0.1])
xlabel('Residual productivity')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

end