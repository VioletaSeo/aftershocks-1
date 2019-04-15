clear
close all

%%
SAVEFIG = 'yes';

%%
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;
minMag  = 6;
% completensess = 4.3;
% magRange = [completeness,10];
maxDepth = 55;
[ASinfo,k,alpha] = aftershock_productivity_kernel(...
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
    'Completeness',4.5);
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;
%%
M   = MSCat.M;
t   = MSCat.time;
res = MSCat.MSres;
prod= MSCat.MSprod;
fms = MSCat.fms;
E   = MSCat.BBtotalEnergy;
Enorm=E./mw2mo(M);
dur = MSCat.RuptureDuration;
durNorm = dur./(mw2mo(M)).^(1/3);

depth= MSCat.depth;
RT  = MSCat.RuptureDuration;
age = get_crust_age(MSCat.lat, MSCat.lon, MSCat.depth);

colors = {[0.7      0.7       0.7], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
c = assign_color(fms+1,colors);


%% 
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);
%% SENSITIVYT TESTS + OVERVIEWS
%% prod versus time (completeness)
figure 
scatter(t,res,'Filled','MarkerFaceAlpha',0.4);
datetick('x','yyyy','keeplimits')
xlabel('Time (yr)')
ylabel('Relative Productivity (N^*)')
savefigure(gcf,'completeness_with_time',SAVEFIG)

%% gutenberg richter
plotgr(M);
xlabel('M_w (declustered)')
ftsz(gcf,16);
setsize(gcf,2,2)
savefigure(gcf,'gr_plot',SAVEFIG)

irisCat = load('IRIS_DMC_with_FMS.mat');
plotgr(irisCat.CAT.M);
xlabel('M_w (ISC-NEIC)')
ftsz(gcf,16);
setsize(gcf,2,2)
savefigure(gcf,'gr_plot_declust',SAVEFIG)

%% sensitivity to time/space windows:
if 1
nT = 3;
timeSelAry = logspace(1,3,nT);
timeWinAry = 100/60 * timeSelAry;

nDX = 3;
distSelAry = logspace(0,2,nDX);
distWinAry = 5/4 * distSelAry;
count = 0;
figure
for n = 1:nT
    for m = 1:nDX
        count = count+1;
        subplot(nT,nDX,count)
        [ASinfoSENS,kSENS,alphaSENS] = aftershock_productivity_kernel(...
            CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',minMag, ...
            'DepthRange',[0, maxDepth], ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'TimeWindow', timeWinAry(n), ...
            'TimeSelectionWindow',timeSelAry(n), ...
            'SearchRadius',distWinAry(m),...
            'SelectionRadius',distSelAry(m), ...
            'Completeness', 4.3);

        SENScat = CAT(ASinfoSENS.ID,:);
        scatter(SENScat.M,ASinfoSENS.MSprod,'Filled','MarkerFaceAlpha',0.4); hold on
        p1 = plot(SENScat.M,kSENS*10.^(alphaSENS*SENScat.M));
       
        [ASinfoSENS,kSFL,alphaSFL] = aftershock_productivity_kernel(...
            CAT.time, ...
            CAT.lat, ...
            CAT.lon, ...
            CAT.depth, ...
            CAT.M, ...
            CAT.fms, ...
            'MinMainshockMag',minMag, ...
            'DepthRange',[0, maxDepth], ...
            'ReturnCatalog', 'yes', ...
            'SaveCatalog', 'no', ...
            'PlotYN','no', ...
            'TimeWindow', timeWinAry(n), ...
            'TimeSelectionWindow',timeSelAry(n), ...
            'SearchRadius',distWinAry(m),...
            'SelectionRadius',distSelAry(m), ...
            'Completeness', 4.3, ...
            'ShuffleYN','yes');
        SENScat = CAT(ASinfoSENS.ID,:);
        scatter(SENScat.M,ASinfoSENS.MSprod,'Filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[.4 .4 .4]);
        p2 = plot(SENScat.M,kSFL*10.^(alphaSFL*SENScat.M),'k' );
        
        legend([p1,p2],{sprintf('N=%0.1g10^{%0.1g M_w}',kSENS, alphaSENS) ...
                sprintf('N_{SFL}=%0.1g10^{%0.1g M_w}',kSFL, alphaSFL)}, ...
                'Location','southeast')
        
        title(['T = ',num2str(timeSelAry(n)), ' days, R = ',num2str(distSelAry(m)), ' R_{source}'])
        set(gca,'YScale','log')
        set(gca,'YLim',[0.1,10000])
        
    end
end

setsize(gcf,9,6)
ftsz(gcf,6)
savefigure(gcf,'space_time_sensitivity',SAVEFIG)

end
      
%% RESULTS
%% Productivity law

figure; 
subplot(2,3,[1 5]);hold on;
scatter(M,prod,30,c,'filled', 'MarkerFaceAlpha',0.3);
set(gca,'Yscale','log');
plot(minmax(M),k*10.^(alpha*minmax(M)),'--k','LineWidth',2);
legend({'Earthquake sequence',sprintf('N = %0.2g 10^{%0.2g M_w}',k,alpha)}, ...
    'Location','southeast');
ylabel('Number of Aftershocks (N)')
xlabel('M_w')

edges = min(res(~isinf(res))):0.2:max(res);
for n = [3,1,2]
    I = fms == n;
    subplot(2,3,3); hold on
    histogram(res(I),edges,'FaceColor',colors{n+1},'EdgeColor',[1 1 1])
    subplot(2,3,6); hold on
    fade_plot(res(I),n,colors{n+1});
end

subplot(2,3,3)
set(gca,'Xlim',[-2,2])
set(gca,'XTick',[])
ylabel('N')
subplot(2,3,6)
set(gca,'YTick',[0,1,2,3],'YTickLabel',{' ','Strike-slip','Normal','Reverse'},'YLim',[0.1,4.1])
xlabel('Relative productivity (N^*)')
set(gca,'Xlim',[-2,2]) 

ftsz(gcf,16)
setsize(gcf,6.5,5)
savefigure(gcf,'fms_prod',SAVEFIG)

%%
% p value for
res2 = res(~isinf(res));
fms2 = fms(~isinf(res));
[h,p] = kstest2(res2(fms2 == 1),res2(fms2 == 3)) % ... tiny!!
[h,p] = kstest2(res2(fms2 == 2),res2(fms2 == 3))

%% productivity by tectonic region
%% Productivity maps

worldmap_res('IRIS DMC');
h = colorbar;
ylabel(h,'Relative Productivity')
setsize(gcf,8,6)

savefigure(gcf,'worldmap_res',SAVEFIG)

%% prod vs depth

[ASWD,~,~] = aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'MinMainshockMag',minMag, ...
    'ReturnCatalog', 'yes', ...
    'SaveCatalog', 'no', ...
    'PlotYN','no');

figure; hold on
scatter(CAT.depth(ASWD.ID), ASWD.MSres,...
        [],assign_color(CAT.fms(ASWD.ID)+1,colors), ...
        'Filled','MarkerFaceAlpha',0.4);
ax = gca;
yLim = ax.YLim;
plot([maxDepth maxDepth],yLim,'LineStyle','--','Color',[colors{1},0.4],'LineWidth',4)
xlabel('Depth (km)')
ylabel('Relative Productivity (N^*)')
set(gca,'XScale','log')
xLim = ax.XLim;
ax.XLim = [3,xLim(2)];
legend({'Earthquake Sequences','Maximum depth Considered in this study'})

ftsz(gcf,12);
setsize(gcf,5,3);
savefigure(gcf,'prod_vs_depth',SAVEFIG)

%% NearSS
nearSSpairs;
ftsz(gcf,12);
setsize(gcf,3.8,1.7);
%%
savefigure(gcf,'nearSS',SAVEFIG);

%% plate boundary
plate_boundaries
ftsz(gcf,12);
setsize(gcf,2.5,1.8);
%%
savefigure(gcf,'prod_by_plate_boundary',SAVEFIG);

%% prod vs age
figure;
yyaxis right
hold on
numBin      = 30;
edges       = linspace(0,max(age),numBin);
ylim([0,200]); ylabel('No aftershocks (N)')
I = isinf(res);
histogram(age(I),edges,'EdgeColor','none');

yyaxis left; ylim([-2,2]);

scatter(age,res,25,c,'filled','MarkerFaceAlpha',0.2);
xlabel('Sea floor age')
ylabel('Relative Productivity')
[Mage,Mres,Merr] = mov_mean(age,res,150,30);
plot(Mage,Mres,'-','LineWidth',3,'Color',[1 0 0 0.5])
plot(Mage,Merr(:,1),'--','LineWidth',1,'Color',[1 0 0 0.5])
plot(Mage,Merr(:,2),'--','LineWidth',1,'Color',[1 0 0 0.5])
grid on
ftsz(gcf,12);
setsize(gcf,5,3);

%%
savefigure(gcf,'prod_vs_age',SAVEFIG)




%% SOURCE ATTRIBUTES
%% prod vs enery
figure
subplot(1,2,1)
scatter(E,prod,[],c,'Filled','MarkerFaceAlpha',0.4)
xlabel('Broad-band radiated energy (J)')
ylabel('Number of aftershocks (N)')
set(gca,'Yscale','log');
set(gca,'Xscale','log');

subplot(1,2,2)
scatter(Enorm, res,[],c,'Filled','MarkerFaceAlpha',0.4);
xlabel('Scaled Energy (E_{BB}/Mo)')
ylabel('Relative productivity (N^*)')
set(gca,'Xscale','log');
ftsz(gcf,12);
setsize(gcf,6,3);
savefigure(gcf,'prod_vs_energy',SAVEFIG)

%% prod vs rupture duration

figure
subplot(1,2,1)
scatter(dur,prod,[],c,'Filled','MarkerFaceAlpha',0.4)
xlabel('Rupture Duration [s]')
ylabel('Number of aftershocks (N)')
set(gca,'Yscale','log');
set(gca,'Xscale','log');

subplot(1,2,2)
scatter(durNorm, res,[],c,'Filled','MarkerFaceAlpha',0.4);
xlabel('Scaled duration, t_R/Mo^{1/3} [s/(Nm)^{1/3}]')
ylabel('Relative productivity (N^*)')
set(gca,'XScale','log');

%%
ftsz(gcf,12);
setsize(gcf,6,3);

savefigure(gcf,'prod_vs_dur',SAVEFIG)


%% FINITE FAULT INVERSIONS ANALYSIS
%% ML
hayesCat = load('Hayes_fully_merged_catalog.mat');
hcat = hayesCat.cat;



%% 

function C = assign_color(X,c)
nX = length(X);
C = zeros(nX,3);
for n = 1:nX
    C(n,:) = c{X(n)};
end
end

function fade_plot(X,yPos,c)
sz     = 200;
scatter(X,repmat(yPos,1,length(X)),sz,'filled','MarkerFaceAlpha',0.03,'MarkerFaceColor',c)
plot(prctile(X,[25,75]),[yPos,yPos],'LineWidth',3,'Color',c/1.4)  
scatter(median(X),yPos,sz,'filled','MarkerFaceColor',c/1.4)
end

function savefigure(fh,figureName,plotYN)
if strcmp(plotYN,'yes')
saveas(fh,figureName);
print(fh,[figureName,'.png'],'-dpng','-r300');
end
end