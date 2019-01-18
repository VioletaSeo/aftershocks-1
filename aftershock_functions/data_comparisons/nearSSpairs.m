% robustly compare whether co located earthquakes have differences in
% productivity based on focal mechanism

%% load catalog
defaultInpStr = {'ForeshockTimeWindow', 10      , ...
                 'DepthRange',          [0,55]  , ...
                 'ShowOverviewYN',      'no'}; % default input

CAT      = aftershock_productivity_GCMT_ISC(defaultInpStr{:});
OTFCAT   = aftershock_productivity_GCMT_ISC(defaultInpStr{:}, ...
                                           'PlateBoundaryClass','OTF', ...
                                           'PlateBoundaryDist',150);                                       
CTFCAT  = aftershock_productivity_GCMT_ISC(defaultInpStr{:}, ...
                                           'PlateBoundaryClass','CTF', ...
                                           'PlateBoundaryDist',150);


plotprod = @(c) scatter(c.MSmag_appended_cat1, c.MSprod_appended_cat1,'filled','MarkerFaceAlpha',0.5);

plotYESNO = 1;
if plotYESNO
    plotmap  = @(c) scatterm(c.MSlat, c.MSlon, (c.MSmag_appended_cat1-4).^4,c.MSres_appended_cat1,'filled');
    worldmap_base
    plotmap(OTFCAT)
    
    worldmap_base
    plotmap(CTFCAT)
    
    plotmap  = @(c) scatterm(c.MSlat(I), c.MSlon(I), (c.MSmag_appended_cat1(I)-4).^4,c.MSres_appended_cat1(I),'filled');
    worldmap_base
    plotmap(CAT)
end

%% alt: CHANGE THE RESIDUAL BELOW!!
% tempCatFN           = 'temp_earthquake_cat_merge.mat';
% delete_temp         = @() system(sprintf('rm %s',tempCatFN));
% minMainshock = 6.8;
% defaultInpStr = {'eq_catalog.mat', ...
%     'SaveCatalog',tempCatFN,...
%     'MinMainshockMag',minMainshock, ...
%     'PlotYN', 'no', ...'ForeshockTimeWindow', 10      , ...
%     'DepthRange',          [0,55]  , ...
%     'ShowOverviewYN',      'no'}; % default input
% 
% % 1) 
% aftershock_productivity_kernel(defaultInpStr{:}); 
% CAT     = struct2table(load(tempCatFN)); delete_temp();
% 
% % 2)
% aftershock_productivity_kernel(defaultInpStr{:}, ...
%                                            'PlateBoundaryClass','OTF', ...
%                                            'PlateBoundaryDist',200); 
% OTFCAT  =  struct2table(load(tempCatFN)); delete_temp();
% 
% % 3)
% aftershock_productivity_kernel(defaultInpStr{:}, ...
%                                            'PlateBoundaryClass','CTF', ...
%                                            'PlateBoundaryDist',200);   
% CTFCAT  =  struct2table(load(tempCatFN)); delete_temp();

%% 

OTFCATSS= OTFCAT(OTFCAT.MSfms == 1,:);
CTFCATSS= CTFCAT(CTFCAT.MSfms == 1,:);

%% pair up earthquakes
dxCrit = 200;

ISS = CAT.MSfms == 1;
SSMS= CAT(ISS,:);
DSMS= CAT(~ISS,:);

XSS = get_coord(SSMS);
XDS = get_coord(DSMS);

n = 1;
clear SScolocCAT DScolocCAT
while true 
    [I,D] = knnsearch(XDS,XSS);
    
    [nextDistance, IminSS] = min(D); % get indices
    IminDS = I(IminSS);
    
    cond = nextDistance < dxCrit;
    if ~cond
        break
    end
    
    SScolocCAT(n,:) = SSMS(IminSS,:); %store eq into new array
    DScolocCAT(n,:) = DSMS(IminDS,:);
    
    [XSS, SSMS] = rem_row(IminSS, XSS, SSMS);
    [XDS, DSMS] = rem_row(IminSS, XDS, DSMS);
    n = n+1;
end

%% 
plot_output_map(SScolocCAT,DScolocCAT) 

%%
plotres({CAT.MSres_appended_cat1, ...
         OTFCATSS.MSres_appended_cat1, ...
         SSMS.MSres_appended_cat1, ...
         CTFCATSS.MSres_appended_cat1, ...
         SScolocCAT.MSres_appended_cat1, ...
         DScolocCAT.MSres_appended_cat1, ...
         DSMS.MSres_appended_cat1}, ...
         {'All','OTF-SS','Strike-slip','CTF-SS','Co-located strike-slip', ...
         'Co-located dip-slip', 'Dip-slip'})
%% 

function X = get_coord(CAT)
[x,y,z] = geodetic2ecef(wgs84Ellipsoid('kilometers'),CAT.MSlat,CAT.MSlon,-10*CAT.MSdepth);
X = [x,y,z];
end

function varargout = rem_row(I,varargin)
    for n = 1:length(varargin)
       arg = varargin{n};
       arg(I,:) = [];
       varargout{n} = arg;
    end
end

function plot_output_map(SScolocCAT,DScolocCAT) 

worldmap_base
ploteq = @(cat,c) scatterm(cat.MSlat,cat.MSlon,(cat.MSmag-4).^3,c,'filled');
ploteq(SScolocCAT,'b')
ploteq(DScolocCAT,'r')

end

function plotres(args,names)
assert(length(args) == length(names))

figure; hold on;
yTickLabel = cell(length(args),1);
for n = 1:length(args)
    plot_row(args{n},n)
    yTickLabel{n} = sprintf('%s, N=%g', names{n}, length(args{n}));
end
set(gca,'YTick',1:n,'YTickLabel',yTickLabel,'YLim',[0.1,n + 0.1],'Xlim',[-1,1])
xlabel('Residual productivity')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

    function plot_row(res,N)
        colors     = get(gca, 'ColorOrder');
        sz = 50;
        scatter(res,        N*ones(size(res)),  sz*2,'filled','MarkerFaceColor',colors(mod(N,7)+1,:), 'MarkerFaceAlpha',0.1)
        medRes = median(res);
        scatter(medRes,N,                 sz,  'filled','MarkerFaceColor',colors(mod(N,7)+1,:))
        plot(prctile(res,[40,60]),[N,N],'LineWidth',2,'Color',colors(mod(N,7)+1,:))
    end
end