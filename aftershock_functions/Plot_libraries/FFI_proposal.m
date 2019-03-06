%% FINITE FAULT INVERSIONS ANALYSIS

hayesCat = load('Hayes_fully_merged_catalog.mat');
c = hayesCat.cat;

%%
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);

%% define input values:
fms   = c.fms_appended_cat1;
res   = c.MSres_appended_cat1;
prod  = c.MSprod_appended_cat1;

Mo    = c.Mo;
Wnorm = c.Width.*Mo.^(-1/3);
Lnorm = c.Length.*Mo.^(-1/3);
areaNorm=c.Length.*c.Width./(10.^(0.59*c.Mw));
logAnorm=log10(areaNorm);
H     = c.Heterogeneity;
A     = c.AspectRatio;
A(A<1) = 1./A(A<1);
logA  = log10(A);
V     = c.RuptureVelocity;
SD    = c.StressDrop;
logSD = log10(SD);
Y     = c.Young;
P     = c.Poisson;
Enorm = c.BBtotalEnergy_appended_cat1./Mo;
durNorm=c.RuptureDuration_appended_cat1;
age   = c.age_appended_cat1;
dip   = c.Dip;
rake  = c.Rake;
depth = c.Depth;

%%
colors = {[0.7      0.7       0.7], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
colorAry = assign_color(fms+1,colors);

%% r2 stem
X = [log10(areaNorm), Wnorm,Lnorm,H,A,V,log10(SD),Y,P,log10(Enorm),log10(durNorm)];
Y = c.MSres_appended_cat1;
SourceParameters = {'A*','W*','L*','H','A_r','V_r','\Delta\sigma','Y','P','E*', 't*'};
r2analysis(X,Y,SourceParameters,10000);
%%
setsize(gcf,6,2)
savefigure(gcf,'r2stems','yes')

%% Aspect ratio and stress drop
figure
subplot(1,2,1)
scatter(A,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
set(gca,'XScale','log')
xlabel('Aspect Ratio')
ylabel('Relative Productivity')

subplot(1,2,2)
scatter(SD,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
set(gca,'XScale','log')
xlabel('Stress Drop, \Delta \sigma_E [MPa]')
ylabel('Relative Productivity')

%%
ftsz(gcf,12)
setsize(gcf,2.8*2,2.5)
savefigure(gcf,'prod_vs_aspect','yes')
%%

goodCat = table(Enorm, depth,logA,logAnorm,logSD,dip,rake,age,res);

%%
load('SVMmodel.mat')
figure;
hold on; 
plot(minmax(res'),minmax(res'),'--','LineWidth',2,'Color','k');
scatter(res,trainedModel.predictFcn(goodCat),[],colorAry,'filled','MarkerFaceAlpha',0.5);
legend({'1:1','Mainshocks'},'Location','northwest')
xlabel('Relative Productivity')
ylabel('Prediction')
ftsz(gcf,11);
xlim([-1.2,1.2])
ylim([-1.2,1.2])
setsize(gcf,3,3)
%%
savefigure(gcf,'response','yes')

%%

function savefigure(fh,figureName,plotYN)
if strcmp(plotYN,'yes')
saveas(fh,figureName);
print(fh,[figureName,'.png'],'-dpng','-r300');
end
end

function C = assign_color(X,c)
nX = length(X);
C = zeros(nX,3);
for n = 1:nX
    C(n,:) = c{X(n)};
end
end