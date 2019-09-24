% summary plot ridgecrest

p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);
itsp = @(x) {x{1:2},x{3}+1};
subplotNum = {2,2,1};
%% a)

% GLOBAL catalog
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;


eqOfInterest = [7.1,22;6.4,6];
CAT     = CAT(CAT.time > datenum(1990,01,01),:);

inputCell = {CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',6.2, ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'Completeness', 4.5,...
            'TimeSelectionWindow',34/24};
[M,prod,cAry] = get_cat(CAT,inputCell);
resAll = get_res([M;eqOfInterest(:,1)],[prod,;eqOfInterest(:,2)]);
disp(['The mainshock was in the ',num2str(sum(resAll < resAll(end-1))/length(resAll)*100),' percentile of strike slip earthquakes'])
disp(['The mainshock was in the ',num2str(sum(resAll < resAll(end))/length(resAll)*100),' percentile of strike slip earthquakes'])

plot_productivity_law(M,prod,cAry,eqOfInterest,subplotNum,'a)');
subplotNum = itsp(subplotNum);
ylabel({'Global (M_c=4.5):','Number of Aftershocks'})
title('34-Hour Time Window:')


%% b)

% Observed counts
eqOfInterest = [7.1,28];

inputCell = {CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',6.2, ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'Completeness', 4.5,...
            'TimeSelectionWindow',28};

[M,prod,cAry] = get_cat(CAT,inputCell);
plot_productivity_law(M,prod,cAry,eqOfInterest,subplotNum,'b)');
subplotNum = itsp(subplotNum);
title('28-Day Time Window:')


%% c) Regional 34 
% SCEDC
load('SCEDC_CAT_all.mat')
CAT = table(t',lat',lon',depth',M',fms',...
    'VariableNames',{'time','lat','lon','depth','M','fms'});
eqOfInterest = [7.1,762;6.4,224];

inputCell = {CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',4.5, ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'Completeness', 2.5,...
            'TimeSelectionWindow',34/24};


[M,prod,cAry] = get_cat(CAT,inputCell);
plot_productivity_law(M,prod,cAry,eqOfInterest,subplotNum,'c)');
subplotNum = itsp(subplotNum); 
ylabel({'Regional (M_c=2.5):','Number of Aftershocks'})
xlabel('Mainshock Magnitude (M_w)')

        
%% d)Regional 28

eqOfInterest = [7.1,1355];

inputCell = {CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',4.5, ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'Completeness', 2.5,...
            'TimeSelectionWindow',28};


[M,prod,cAry] = get_cat(CAT,inputCell);
plot_productivity_law(M,prod,cAry,eqOfInterest,subplotNum,'d)');   
xlabel('Mainshock Magnitude (M_w)')

%%
ftsz(gcf,10)
setsize(gcf,6,6)
%%
savefigure(gcf,'ridgecrest_summary','yes')

%% comparison to global strike-slip:

% GLOBAL catalog
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;

eqOfInterest = [7.1,22;6.4,6];
CAT     = CAT(CAT.time > datenum(1990,01,01),:);

inputCell = {CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',6.2, ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'Completeness', 4.5,...
            'TimeSelectionWindow',34/24};
[M,prod,cAry,fms,~,PB] = get_cat(CAT,inputCell);
I = fms == 1;
resAll = get_res([M(I);eqOfInterest(:,1)],[prod(I),;eqOfInterest(:,2)]);
disp(['The mainshock was in the ',num2str(sum(resAll < resAll(end-1))/length(resAll)*100),' percentile of strike slip earthquakes'])
disp(['The mainshock was in the ',num2str(sum(resAll < resAll(end))/length(resAll)*100),' percentile of strike slip earthquakes'])

subplotNum = {1,1,1};
figure
plot_productivity_law(M(I),prod(I),cAry(I,:),eqOfInterest,subplotNum,'');

ylabel({'Global (M_c=4.5):','Number of Aftershocks'})
xlabel('Mainshock Magnitude (M_w)')
title({'34-Hour Time Window:','Strike-Slip Earthquakes'})


ftsz(gcf,10)
setsize(gcf,4,3.5)
savefigure(gcf,'ridgecrest_strike_slip','yes')
%% CTF's 34j

I = strcmp(PB,'CTF');
resAll = get_res([M(I);eqOfInterest(:,1)],[prod(I),;eqOfInterest(:,2)]);
disp(['The mainshock was in the ',num2str(sum(resAll < resAll(end-1))/length(resAll)*100),' percentile of strike slip earthquakes'])
disp(['The mainshock was in the ',num2str(sum(resAll < resAll(end))/length(resAll)*100),' percentile of strike slip earthquakes'])

subplotNum = {1,1,1};
figure
plot_productivity_law(M(I),prod(I),cAry(I,:),eqOfInterest,subplotNum,'');

ylabel({'Global (M_c=4.5):','Number of Aftershocks'})
xlabel('Mainshock Magnitude (M_w)')
title({'34-Hour Time Window:','Continental Transform Earthquakes'})


ftsz(gcf,10)
setsize(gcf,4,3.5)
savefigure(gcf,'ridgecrest_CTF','yes')

%%


function [M,prod,cAry,fms,res,PB] = get_cat(CAT,inputCell)% aftershock statistics
[ASinfo,~,~] = aftershock_productivity_kernel(inputCell{:});
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;
MSCat = MSCat(MSCat.fms ~= 0,:);


% shorter variable names
M       = MSCat.M;
fms     = MSCat.fms;
lat     = MSCat.lat;
lon     = MSCat.lon;

% aftershock stats
prod    = MSCat.MSprod;
res     = MSCat.MSres;

% plotting colors (be focal mechanism)
colors = {[0      0       0], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
cAry = assign_color(fms+1,colors);
distance2pb = 50;
PB = assign_PB_class(lat,lon, distance2pb,'yes',fms);
end

%%plots:
function plot_productivity_law(M,prod,c,eqOfInterest,subplotnum,abcd)

alpha = 0.2;
mRange= [min(M)-0.2,max(M+0.2)];
tpos = {-0.25,1};

% plot a)
subplot(subplotnum{:})
text(tpos{:},abcd,'Units','normalized')
hold on

[magArray,prodArray,K,A] = productivity_law(M,prod); 

scatter(M,prod,40,c,'filled','MarkerFaceAlpha',alpha);
set(gca,'Yscale','log');
squareH = scatter(magArray,prodArray,'s','k');
blankH  = scatter([],[],'k');
lineH   = plot(magArray,K*10.^(A*magArray),'--k','LineWidth',2);
set(gca,'xlim', mRange)


eqS = scatter(eqOfInterest(:,1),eqOfInterest(:,2),80,'rd','filled') ;

kexp = log10(K);
legend([lineH], ...
    {sprintf('10^{%0.1f M %0.1f}',A,kexp)}, ...
    'Location','southeast', ...
    'Box','off');

end


%%

%% other functions

function RES = get_res(M,N)
    [~,~,K,A]   = productivity_law(M,N);
    allAS    	= @(x,y) log10(y) - log10(K*10.^(A*x));
    RES         = allAS(M,N);
end

function C = assign_color(X,c)
nX = length(X);
C = zeros(nX,3);
for n = 1:nX
    C(n,:) = c{X(n)};
end
end

function [magArray,prodArray,K,A] = productivity_law(MAGNITUDE,AFTERSHOCKS)

increment = 0.1; 
minmaxMag = minmax(MAGNITUDE');

magArray    = (minmaxMag(1)):increment:minmaxMag(2);
numMag      = length(magArray);
prodArray   = zeros(1,numMag);

for iM = 1:(numMag)
    % fit a normal distribution parameters mu and sigma to the
    % distribution of earthquake productivities:
    asSub    = AFTERSHOCKS(MAGNITUDE>=magArray(iM) & MAGNITUDE<(magArray(iM)+increment));
    prodArray(iM) = median(asSub);
end

[K,A] = getproductivity(magArray,prodArray);

end

function fade_plot(X,yPos,c)
sz     = 100;
scatter(X,repmat(yPos,1,length(X)),sz,'filled','MarkerFaceAlpha',0.03,'MarkerFaceColor',c)
plot(prctile(X,[25,75]),[yPos,yPos],'LineWidth',3,'Color',c/1.4)  
scatter(median(X),yPos,sz,'filled','MarkerFaceColor',c/1.4)
end

function savefigure(fh,figureName,plotYN)
if strcmp(plotYN,'yes')
saveas(fh,[figureName]);
print(fh,[figureName,'.png'],'-dpng','-r300');
end
end

function OUT = it(start)
    persistent out
    if isempty(out)
        out = 1;
    end
    
    if nargin == 1
        out = start;
    end
    
    out = out+1;
    OUT = out;
end

function Yq = knn_regress(Xt,Yt,Xq,k)
Mdl = KDTreeSearcher(Xt);
I  = knnsearch(Mdl,Xq,'K',k);
nQ = length(Xq);
Yq = zeros(nQ,1);
for n = 1:nQ
    Yq(n) = median(Yt(I(n,2:end))); % skips itself
end

end

function X = get_coord(CAT)
[x,y,z] = geodetic2ecef(wgs84Ellipsoid('kilometers'),CAT.lat,CAT.lon,-10*CAT.depth);
X = [x,y,z];
end

function varargout = rem_row(I,varargin)
    for n = 1:length(varargin)
       arg = varargin{n};
       arg(I,:) = [];
       varargout{n} = arg;
    end
end

