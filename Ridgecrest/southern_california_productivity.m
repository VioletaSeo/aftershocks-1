

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROSS catalog
% rCAT = load('ross_cat.mat');
% rCAT = CAT.ross_cat;
% CAT = [];
% CAT.time= datenum(CAT.YEAR,CAT.MONTH,CAT.DAY,CAT.HOUR,CAT.MINUTE,CAT.SECOND);
% CAT.lat = rCAT.LATITUDE;
% CAT.lon = rCAT.LONGITUDE;
% CAT.lat = rCAT.DEPTH;
% CAT.fms = zeros(height(rCAT),1);

% Southern California catalog


% GLOBAL catalog
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;

% SCEDC
% load('SCEDC_CAT_all.mat')
% CAT = table(t',lat',lon',depth',M',fms',...
%     'VariableNames',{'time','lat','lon','depth','M','fms'});

% MAGIC NUMBERS                                                 Justification
% -------------                                                 -----------
minMag      = 6.3; % Minimum magnituce of mainshocks considered (Mc + 2)
maxDepth    = 55;  % Maximum mainshock depth considered         (No deep focus earthquakes)
completeness= 4.5; % Completeness of the catalog                  (from KS test for Mid ocean rideges)
startDate   = 1980;% start of the catalog                       big drop in MC around start of 1990s 
timeSelectionWindow = 34/24;
distance2pb = 60; % Maximum distance to plate                  ~ maxDepth/tand(10) 

% eqOfInterest = {7.1,1355};
% eqOfInterest = {7.1,762};
% eqOfInterest = {6.4,224};
eqOfInterest = {7.1,22};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);


CAT     = CAT(CAT.time > datenum(startDate,01,01),:);

inputCell = {CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',minMag, ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'Completeness', completeness,...
            'TimeSelectionWindow',timeSelectionWindow};

% aftershock statistics
[ASinfo,k,alpha] = aftershock_productivity_kernel(inputCell{:});
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;
MSCat = MSCat(MSCat.fms ~= 0,:);


% shorter variable names
M       = MSCat.M;
lat     = MSCat.lat;
lon     = MSCat.lon;
depth   = MSCat.depth;
t       = MSCat.time;
fms     = MSCat.fms;

% aftershock stats
prod    = MSCat.MSprod;
res     = get_res(M,prod);

% plotting colors (be focal mechanism)
colors = {[0      0       0], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
cAry = assign_color(fms+1,colors);

PB = assign_PB_class(lat,lon, distance2pb,'yes',fms);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) productivity law 
 I = strcmp(PB,'CTF');
% I = fms == 1;
% I = ones(size(fms));
plot_productivity_law(M(I),prod(I),res(I),cAry(I,:),eqOfInterest)
ftsz(gcf,8)
setsize(gcf,4.5,3.5)
savefigure(gcf,'ridgecrest_Mw71_ctf','no')

%%plots:
function plot_productivity_law(M,prod,res,c,eqOfInterest)

alpha = 0.2;
mRange= [min(M)-0.2,max(M+0.2)];
tpos = [-0.15,1];
figure

% plot a)
subplot(3,3,[1:2,4:5])
t = title('a)'); set(t,'Position',tpos,'Units','normalized')
hold on

[magArray,prodArray,K,A] = productivity_law(M,prod);
scatter(M,prod,40,c,'filled','MarkerFaceAlpha',alpha);
set(gca,'Yscale','log');
squareH = scatter(magArray,prodArray,'s','k');
blankH  = scatter([],[],'k');
lineH   = plot(magArray,K*10.^(A*magArray),'--k','LineWidth',2);
ylabel('Number of Aftershocks, N')
set(gca,'xlim', mRange,'XTickLabel',[])


eqS = scatter(eqOfInterest{:},80,'rd','filled') ;

kexp = log10(K);
lg = legend([blankH,squareH,lineH, eqS], ...
    {'Mainshock', ...
    'Median',...
    sprintf('10^{%0.1f M %0.1f}',A,kexp), ...
    sprintf('M_w%0.1f - %g',eqOfInterest{:})}, ...
    'Location','southeast');

newPosition = [0.48 0.42 0.1 0.1];
newUnits = 'normalized';
set(lg,'Position', newPosition,'Units', newUnits);

% plot b)
subplot(3,3,7:8)
t = title('b)'); set(t,'Position',tpos,'Units','normalized')

hold on
scatter(M,res,40,c,'filled','MarkerFaceAlpha',alpha);
set(gca,'xlim',mRange)
xlabel('M_w')
ylabel({'Relative', 'productivity','\Delta log(N)'})

% plot c)
subplot(3,3,9); hold on
t = title('c)'); set(t,'Position',tpos,'Units','normalized')
histogram(res,15,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
view(90, -90)
set(gca,'XTickLabel',[])
ylabel('Number of Mainshocks')


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

