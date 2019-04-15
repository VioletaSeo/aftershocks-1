% % robustly compare whether co located earthquakes have differences in
% % productivity based on focal mechanism
% 
% %% load catalog
% defaultInpStr = {'ForeshockTimeWindow', 10      , ...
%                  'DepthRange',          [0,55]  , ...
%                  'ShowOverviewYN',      'no'}; % default input
% 
% CAT      = aftershock_productivity_GCMT_ISC(defaultInpStr{:});
% OTFCAT   = aftershock_productivity_GCMT_ISC(defaultInpStr{:}, ...
%                                            'PlateBoundaryClass','OTF', ...
%                                            'PlateBoundaryDist',150);                                       
% CTFCAT  = aftershock_productivity_GCMT_ISC(defaultInpStr{:}, ...
%                                            'PlateBoundaryClass','CTF', ...
%                                            'PlateBoundaryDist',150);
% 
% 
% plotprod = @(c) scatter(c.MSmag_appended_cat1, c.MSprod_appended_cat1,'filled','MarkerFaceAlpha',0.5);
% 
% plotYESNO = 1;
% if plotYESNO
%     plotmap  = @(c) scatterm(c.MSlat, c.MSlon, (c.MSmag_appended_cat1-4).^4,c.MSres_appended_cat1,'filled');
%     worldmap_base
%     plotmap(OTFCAT)
%     
%     worldmap_base
%     plotmap(CTFCAT)
%     
%     plotmap  = @(c) scatterm(c.MSlat(I), c.MSlon(I), (c.MSmag_appended_cat1(I)-4).^4,c.MSres_appended_cat1(I),'filled');
%     worldmap_base
%     plotmap(CAT)
% end
% 
% %% alt: CHANGE THE RESIDUAL BELOW!!
% % tempCatFN           = 'temp_earthquake_cat_merge.mat';
% % delete_temp         = @() system(sprintf('rm %s',tempCatFN));
% % minMainshock = 6.8;
% % defaultInpStr = {'eq_catalog.mat', ...
% %     'SaveCatalog',tempCatFN,...
% %     'MinMainshockMag',minMainshock, ...
% %     'PlotYN', 'no', ...'ForeshockTimeWindow', 10      , ...
% %     'DepthRange',          [0,55]  , ...
% %     'ShowOverviewYN',      'no'}; % default input
% % 
% % % 1) 
% % aftershock_productivity_kernel(defaultInpStr{:}); 
% % CAT     = struct2table(load(tempCatFN)); delete_temp();
% % 
% % % 2)
% % aftershock_productivity_kernel(defaultInpStr{:}, ...
% %                                            'PlateBoundaryClass','OTF', ...
% %                                            'PlateBoundaryDist',200); 
% % OTFCAT  =  struct2table(load(tempCatFN)); delete_temp();
% % 
% % % 3)
% % aftershock_productivity_kernel(defaultInpStr{:}, ...
% %                                            'PlateBoundaryClass','CTF', ...
% %                                            'PlateBoundaryDist',200);   
% % CTFCAT  =  struct2table(load(tempCatFN)); delete_temp();
% %% 
% 
% 
%% entire catalog
function nearSSpairs()
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;
minMag  = 6.5;
% completensess = 4.3;
% magRange = [completeness,10];
maxDepth = 55;

defaultInpStr = {...
    CAT.time, ...
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
    'Completeness',4.3}; % default input

[ASinfo] = aftershock_productivity_kernel(defaultInpStr{:});
MSCat = CAT(ASinfo.ID,:); MSCat.MSres = ASinfo.MSres;

[ASinfo] = aftershock_productivity_kernel(defaultInpStr{:},'PlateBoundaryClass','OTF','PlateBoundaryDist',150);
OTFCAT = CAT(ASinfo.ID,:); OTFCAT.MSres = ASinfo.MSres; 


[ASinfo] = aftershock_productivity_kernel(defaultInpStr{:},'PlateBoundaryClass','CTF','PlateBoundaryDist',150);
CTFCAT = CAT(ASinfo.ID,:); CTFCAT.MSres = ASinfo.MSres;


OTFCATSS= OTFCAT(OTFCAT.fms == 1,:);
CTFCATSS= CTFCAT(CTFCAT.fms == 1,:);

%% pair up earthquakes
dxCrit = 200;

ISS = MSCat.fms == 1;
SSMS= MSCat(ISS,:);
DSMS= MSCat(~ISS,:);

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
% plot_output_map(SScolocCAT,DScolocCAT) 
XSS = get_coord(SScolocCAT);
XDS = get_coord(DScolocCAT);
r   = sqrt(sum((XSS-XDS).^2,2));
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);
figure; hold on
plot([-1.5,1.5],[-1.5,1.5],'--k','LineWidth',2);
scatter(SScolocCAT.MSres,DScolocCAT.MSres,[],r, 'Filled')
xlim([-1.5,1.5]); ylim([-1.5,1.5])
xlabel('Strike-slip mainshock N^*')
ylabel('Dip-slip mainshock N^*')
cH = colorbar;
ylabel(cH,'Distance (km)');
colormap('gray')
legend({'1:1','Co-located pairs'},'Location','southeast')
ftsz(gcf,12)
setsize(gcf,4.5,3.5)
%%
plotres({MSCat.MSres, ...
         OTFCATSS.MSres, ...
         SSMS.MSres, ...
         CTFCATSS.MSres, ...
         SScolocCAT.MSres, ...
         DScolocCAT.MSres, ...
         DSMS.MSres}, ...
         {'All','OTF-SS','Strike-slip','CTF-SS','Co-located strike-slip', ...
         'Co-located dip-slip', 'Dip-slip'})
%% 

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

function plot_output_map(SScolocCAT,DScolocCAT) 

worldmap_base
ploteq = @(cat,c) scatterm(cat.lat,cat.lon,(cat.M-4).^3,c,'filled');
ploteq(SScolocCAT,'b')
ploteq(DScolocCAT,'r')

end

function plotres(args,names)
assert(length(args) == length(names))

figure; hold on;
yTickLabel = cell(length(args),1);
for n = 1:length(args)
    plot_row(args{n},n)
    yTickLabel{n} = sprintf('%s', names{n});
end
set(gca,'YTick',1:n,'YTickLabel',yTickLabel,'YLim',[0.1,n + 0.1],'Xlim',[-1,1])
xlabel('Relative productivity')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

    function plot_row(res,N)
        colors     = get(gca, 'ColorOrder');
        sz = 50;
        scatter(res,        N*ones(size(res)),  sz*2,'filled','MarkerFaceColor',colors(mod(N,7)+1,:), 'MarkerFaceAlpha',0.01)
        medRes = median(res);
        scatter(medRes,N,                 sz,  'filled','MarkerFaceColor',colors(mod(N,7)+1,:)/1.4)
        plot(prctile(res,[25,75]),[N,N],'LineWidth',2,'Color',colors(mod(N,7)+1,:)/1.4)
    end
end