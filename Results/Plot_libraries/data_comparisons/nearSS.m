function nearSS(MSlat, MSlon, MSdepth, MSmag, MSfms, MSprod)

Mcutoff         = 5.5;  % magnitude cutoff
distanceCutoff  = 100; %km
depthCutoff     = 5;
fmsChoice       = 1;  % 1:SS 2:NS 3:RS

colors = {[     0    0.4470    0.7410],...
          [0.4660    0.6740    0.1880], ...
          [0.6350    0.0780    0.1840]};

if nargin == 0
    aftershock_productivity_kernel('eq_catalog.mat', ...
                                   'PlotYN', 'no', ...
                                   'DepthRange',[0,55], ...
                                   'SaveCatalog','temp_mainshocks.mat');
    load('temp_mainshocks.mat')
end

[I,I2] = is_in_fms_region(MSlat, MSlon, MSdepth, MSmag, MSfms, ...
    distanceCutoff, depthCutoff, Mcutoff, fmsChoice);
I      = I & MSfms' ~= fmsChoice;
ISS    = MSmag > Mcutoff & MSfms == fmsChoice;

% I2 is the location of strike slip earthquakes that are collocated with
% other dip slip earthquakes

[MSprodNearSS,  MSmagNearSS,MSresNearSS]= goodind(I,	MSprod,MSmag,MSres);
[MSprodSS,      MSmagSS,    MSresSS]	= goodind(I2 , 	MSprod,MSmag,MSres);
[MSprodSSall,   MSmagSSall, MSresSSall]	= goodind(ISS,	MSprod,MSmag,MSres);


figure

subplot(1,2,1)
hold on
worldmap({'pacific', 'USA','Indonesia'})
load coastlines
%plotm(coastlat, coastlon)
geoshow('landareas.shp', 'FaceColor', [0.9 0.9 0.9])
scatterm(MSlat,MSlon, (MSmag-4).^3.5,assign_color(MSfms,colors),'filled')
circlem(MSlat(ISS), MSlon(ISS), distanceCutoff,'edgecolor','none','facecolor','blue',  'facealpha',.3);
%circlem(MSlat(I2),  MSlon(I2),  distanceCutoff,'edgecolor','none','facecolor','green','facealpha',.3);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

subplot(1,2,2)
worldmap({'pacific', 'USA','Indonesia'})
load coastlines
%plotm(coastlat, coastlon)
geoshow('landareas.shp', 'FaceColor', [0.9 0.9 0.9])
scatterm(MSlat(I),MSlon(I), (MSmag(I)-4).^3.5,assign_color(MSfms(I),colors),'filled')
Sh = scatterm(MSlat,MSlon, (MSmag-4).^3.5,assign_color(MSfms,colors),'filled','MarkerFaceAlpha',0.1);
Sh.Children.MarkerFaceAlpha = 0.05;
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure; hold on;
plotprod(MSmag,         MSprod)
plotprod(MSmagNearSS,   MSprodNearSS)
plotprod(MSmagSS,       MSprodSS)
plotprod(MSmagSSall,    MSprodSSall)
set(gca,'YScale','log')
L = legend({'All earthquakes','',...
        'Dip slip earthuqakes near SS', '',...
        'SS with DS eqs', '',...
        'All SS eq', ''},'Location','southeast');
 

figure; hold on; 
Mth = 7;
rectangle('Position',[-1.5,1.5,3,2],'FaceColor',[colors{1},0.1],'EdgeColor','none')
fade_plot(MSres((MSfms~=1) & MSmag > Mth),  4,colors{3});
fade_plot(MSresNearSS(MSmagNearSS > Mth),   3,colors{3});
fade_plot(MSresSS(MSmagSS > Mth),           2,colors{1});
fade_plot(MSresSSall(MSmagSSall > Mth),     1,colors{1});
set(gca,'YTick',[0,1,2,3,4],'YTickLabel',{' ','Strike-slip','Co-located strike-slip','Co-located dip-slip','Dip-slip'},'YLim',[0.1,4.1])
xlabel('Residual productivity')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

    
!rm temp_mainshocks.mat

% figure;hold on
% scatter(MSmag,      MSprod,         'filled','MarkerFaceAlpha',0.3)
% scatter(MSmagNearSS,MSprodNearSS,   'filled','MarkerFaceAlpha',0.3)
% scatter(MSmagSS,    MSprodSS,       'filled','MarkerFaceAlpha',0.3)
% xlabel('M_w');
% ylabel('Number of Aftershocks');
% set(gca,'YScale','log')
% legend({'All Data','Near SS'})

end

function fade_plot(X,yPos,c)
sz     = 200;
scatter(X,repmat(yPos,1,length(X)),sz,'filled','MarkerFaceAlpha',0.03,'MarkerFaceColor',c)
    scatter(median(X),yPos,sz,'filled','MarkerFaceColor',c)
    plot(prctile(X,[35,65]),[yPos,yPos],'LineWidth',2,'Color',c)
end

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);
end

end
function plotprod(mag,prod)

[magArray,prodArray,K,A]  = productivity_law(mag,prod);
f = @(x) K*10.^(A*x);
mm = minmax(magArray);

plt = plot(mm,f(mm));
pltColor = get(plt, 'Color');
hold on
scatter(magArray,prodArray, 'filled', 'MarkerFaceColor', pltColor);
xlabel('M_w');
ylabel('Number of Aftershocks');

end
function [magArray,prodArray,K,A]            = productivity_law(MAGNITUDE,AFTERSHOCKS)

increment = 0.1; 
minmaxMag = minmax(MAGNITUDE');

magArray    = minmaxMag(1):increment:minmaxMag(2);
numMag      = length(magArray);
prodArray   = zeros(1,numMag);

for iM = 1:(numMag-1)
    I = MAGNITUDE>=magArray(iM) ...
        & ...
        MAGNITUDE<magArray(iM+1);
    prodArray(iM) = mean(AFTERSHOCKS(I));
end

[K,A] = getproductivity(magArray,prodArray);

end
function C = assign_color(X,c)
nX = length(X);
C = zeros(nX,3);
for n = 1:nX
    C(n,:) = c{X(n)};
end
end


% function bootstrapHist(MSmag,MSres,Mcutoff,nbtstrp)
% 
% I           = ~isinf(MSres) &  MSmag > Mcutoff;
% MSresNoInf  = MSres(I);
% nEqNoInf    = length(MSresNoInf);
% 
% hold on
% 
% for n = 1:3
%     Ipct = MSfmsNoInf == n;
%     numI = sum(Ipct);
%     [res] = bootstrp(floor(nbtstrp*numI/nEqNoInf),@mean, MSresNoInf(Ipct));
%     histogram(res,'BinWidth',0.01, 'FaceColor',colors{n},'EdgeColor','none')
% end
% title(sprintf('Bootstrap (n = %d) for M%d+ Mainshocks', nbtstrp, Mcutoff))
% xlabel('Residual Productivity')
% set(gca,'Ytick',[])
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% end