% This file (hopefully) contains manuscript-ready figures for the
% aftershock productivity project. The structure should be very similar to
% that of the proposal script
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preamble:
% clear
% close all
SAVEFIG = 'yes';
p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and format
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;

% MAGIC NUMBERS                                                 Justification
% -------------                                                 -----------
minMag      = 6.5; % Minimum magnituce of mainshocks considered   (Mc + 2.5)
maxDepth    = 55;  % Maximum mainshock depth considered         (No deep focus earthquakes)
completeness= 4.5; % Completeness of the catalog                (from KS test for Mid ocean rideges)
distance2pb = 400; % Maximum distance to plate                  ~ maxDepth/tand(10) 
startDate   = 1993;% start of the catalog                       big drop in MC around start of 1990s 

CAT     = CAT(CAT.time > datenum(startDate,01,01),:);

% aftershock statistics
[ASinfo,k,alpha] = aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'MinMainshockMag',minMag, ...    'DepthRange',[0,maxDepth], ...
    'ReturnCatalog', 'yes', ...
    'SaveCatalog', 'no', ...
    'PlotYN','no', ...
    'Completeness',completeness);
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;

wholeCat    = MSCat;
wholeCat.age= get_crust_age(wholeCat.lat, wholeCat.lon, wholeCat.depth);
wholeCat.pb = assign_PB_class(wholeCat.lat,wholeCat.lon, distance2pb,'yes',wholeCat.fms);

deepMSCat   = MSCat(MSCat.depth > maxDepth,:);
MSCat       = MSCat(MSCat.depth < maxDepth,:); % lazy way to also store deep eq

% shorter variable names
M       = MSCat.M;
lat     = MSCat.lat;
lon     = MSCat.lon;
depth   = MSCat.depth;
t       = MSCat.time;

% aftershock stats
res     = MSCat.MSres;
prod    = MSCat.MSprod;

% source properties
fms     = MSCat.fms;
E       = MSCat.BBtotalEnergy;
Enorm   = E./mw2mo(M);
dur     = MSCat.RuptureDuration;
durNorm = dur./(mw2mo(M)).^(1/3);
RT      = MSCat.RuptureDuration;

% assign ages
age = get_crust_age(MSCat.lat, MSCat.lon, MSCat.depth);


% plate boundaries
%                           see above   (enforce fms)                             
PB = assign_PB_class(lat,lon, distance2pb,'yes',fms);


% plotting colors (be focal mechanism)
colors = {[0.7      0.7       0.7], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
cAry = assign_color(fms+1,colors);

%%

hayesCat = load('Hayes_fully_merged_catalog.mat');
c = hayesCat.cat;
c = c(c.fms_appended_cat1 ~= 0,:);
c = c(abs(c.time_appended_cat1-c.Time+365)<1,:);
c = c(c.Time>datenum(1993,1,1),:);
c.DateTime = datetime(c.Time,'ConvertFrom','datenum');
c = c(~(c.MSres_appended_cat1 == max(c.MSres_appended_cat1)),:); % remove tuhoku foreshock

FFIfms   = c.fms_appended_cat1;
FFIres   = c.MSres_appended_cat1;
FFIprod  = c.MSprod_appended_cat1;
 
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
Yng     = c.Young;
P     = c.Poisson;
Enorm = c.BBtotalEnergy_appended_cat1./Mo;
durNorm=c.RuptureDuration_appended_cat1;
FFIage   = get_crust_age(c.Lat,c.Lon,0);
FFIdip   = c.Dip;
FFIrake  = c.Rake;
FFIdepth = c.Depth;
FFIlitho = get_litho(c.Lat,c.Lon);

colorAryFFI = assign_color(FFIfms+1,colors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main document figures:
% 1) Producivity law
% 2) World map (insets?)
% 3) Productivity by plate boundary
% 4) Prod vs age
% 5) Stem plots  (with AR and SD)
% 6) FMS dependence
% 7) KN vs SVM ML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% 1) productivity law 
plot_productivity_law(M,prod,res,cAry)
%%
ftsz(gcf,10)
setsize(gcf,4,3.5)
savefigure(gcf,'prod_law',SAVEFIG)

%% 
% 2)  world map:
plot_worldmap_prod(lat,lon,res)
%%
ftsz(gcf,10)
setsize(gcf,7.5,4)
savefigure(gcf,'worldmap_res',SAVEFIG)

%%
% 3) Productivity by plate boundary
plot_prod_by_pb(PB,res)
ftsz(gcf,12)
setsize(gcf,5,2)
savefigure(gcf,'prod_by_pb',SAVEFIG)

%%
% 4) Prod vs age
plot_prod_vs_age(age,res,M,cAry)
ftsz(gcf,12);
setsize(gcf,5,3);
%%
savefigure(gcf,'prod_vs_age',SAVEFIG)


%%
% 5) Stem
X = [log10(areaNorm), log10(Wnorm),log10(Lnorm),log10(H),(A),V,log10(SD),Yng,P,log10(Enorm),log10(durNorm),FFIdip]; Y = FFIres;
SourceParameters = {'A*','W*','L*','H','A_r','V_r','log(\Delta\sigma)','Yng','P','E*', 't*','dip'};
plot_r2_stems(X,Y,SourceParameters,colorAryFFI,SD,A,Wnorm)
ftsz(gcf,10);
setsize(gcf,6.5,4);

%%
savefigure(gcf,'stem_plot',SAVEFIG)

%%
% 6) FMS dependence


%%
% 7) Caltech plot
caltech_plot(wholeCat.depth,wholeCat.age,wholeCat.fms,wholeCat.pb,wholeCat.MSres, A, Wnorm, SD, FFIres) % FINITE FAULTS NEED TO HAVE THE SAME RES CALC!
ftsz(gcf,10);
setsize(gcf,6.5,4);
%%
savefigure(gcf,'cal_tech',SAVEFIG)

%%
% 8) kn vs SVM
% I = ~isnan(FFIage);
plot_ML([(FFIage),log10(Wnorm),(FFIdip),log10(areaNorm)],FFIres,c.Lat,c.Lon,c.Depth)
setsize(gcf,6.2,3)
ftsz(gcf,11);
savefigure(gcf,'response','yes')

%% Supplement:
% prod vs depth?
% space time window sensitivity
% sensitity to Mc
% co-variance matrix
% each parameter as a function of N*

%% functions

%%plots:
function plot_productivity_law(M,prod,res,c)
figure
subplot(3,3,[1:2,4:5])
hold on
[magArray,prodArray,K,A] = productivity_law(M,prod);
scatter(M,prod,40,c,'filled','MarkerFaceAlpha',0.2);
set(gca,'Yscale','log');
scatter(magArray,prodArray,'s','k')
plot(magArray,K*10.^(A*magArray),'--k','LineWidth',2);
ylabel('Number of Aftershocks (N)')
set(gca,'xlim', [min(M)-0.2,9.5],'XTickLabel',[])

lg = legend({'Mainshock','Median',sprintf('%0.2g 10^{%0.2g M_w}',K,A)}, ...
    'Location','southeast');

newPosition = [0.66 0.4 0.1 0.1];
newUnits = 'normalized';
set(lg,'Position', newPosition,'Units', newUnits);

subplot(3,3,7:8)
hold on
scatter(M,res,40,c,'filled','MarkerFaceAlpha',0.2);
set(gca,'xlim', [min(M)-0.2,9.5])
xlabel('M_w')
ylabel('N^*')


subplot(3,3,9)
histogram(res,25,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
view(90, -90)
set(gca,'XTickLabel',[])
end

function plot_worldmap_prod(lat,lon,res)

worldmap_base;
sh  = scatterm(lat,lon,[],res, 'filled');
sh.Children.MarkerFaceAlpha = .7;
h = colorbar;
h.Location = 'south'; 
set(h,'position',[.6,.125,.1,.05])
ylabel(h,'Relative Productivity')

end

function plot_prod_by_pb(PB,res)
figure; hold on;
colors     = get(gca, 'ColorOrder');
sz         = 50;
N          = 0;
PBopt      = unique(PB);
PBlongName = {'None', ...
    'Intraplate', ...
    'Subduction', ...
    'Ocean transform', ...
    'Speading ridge', ...
    'Ocean convenrgent', ...
    'Continental Transform',... 
    'Continental Rift',...
    'Contiental convergent'};

for n = 1:length(PBopt)
    I               = strcmp(PB,PBopt(n));
    ires            = res(I);
    numNegInf(n)    = sum(isinf(ires))/length(res);
    scatter(ires,        n*ones(size(ires)),  sz*2,'filled','MarkerFaceColor',  colors(mod(n,7)+1,:), 'MarkerFaceAlpha',0.1)
    scatter(median(ires),n,                   sz,  'filled','MarkerFaceColor',  colors(mod(n,7)+1,:)/1.4)
    plot([prctile(ires,25),prctile(ires,75)],[n,n],'LineWidth',2,'Color',       colors(mod(n,7)+1,:)/1.4)
    tickName{n} = sprintf('%s', PBlongName{end+1-n});  
end
    
set(gca,'YTick',1:n,'YTickLabel',tickName','YLim',[0.1,n + 0.1])
xlabel('Relative productivity')

end

function plot_prod_vs_age(age,res,M,c)

figure;
yyaxis right
hold on
numBin      = 20;
edges       = linspace(0,max((age)),numBin);
ylim([0,100]); ylabel('No aftershocks (N)')
I = isinf(res);
histogram(age(I),edges,'EdgeColor','none','FaceAlpha',0.5);

yyaxis left; ylim([-2,2]);

scatter((age),res, (M-4).^4,c,'filled','MarkerFaceAlpha',0.2);
xlabel('Sea floor age')
ylabel('Relative Productivity')
[Mage,Mres,Merr] = mov_mean((age),res,ceil(range(age)/20),10);
plot(Mage,Mres,'-','LineWidth',3,'Color',[1 0 0 0.5])
plot(Mage,Merr(:,1),'--','LineWidth',1,'Color',[1 0 0 0.7])
plot(Mage,Merr(:,2),'--','LineWidth',1,'Color',[1 0 0 0.7])
grid on

end

function plot_r2_stems(X,res,SourceParameters,colorAry, SD,A,Wnorm)

figure
subplot(2,3,1:3);
r2analysis(X,res,SourceParameters,10000);
set(gca,'YTickLabel','','YTick',[])

subplot(2,3,4)
scatter(SD,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
set(gca,'XScale','log')
xlabel('Stress Drop, \Delta \sigma_E [MPa]')
ylabel('Relative Productivity')

subplot(2,3,5)
scatter(Wnorm,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
set(gca,'XScale','log')
set(gca,'Yticklabel','')
xlabel('Normalized width, W*')

subplot(2,3,6);
scatter(A,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
%set(gca,'XScale','log')
xlabel('Aspect Ratio')
set(gca,'Yticklabel','')


end

function caltech_plot(depth,age,fms,PB,res, A, Wnorm, SD, FFIres)

% yougn age
% intermediate age
% plate boundaries
% fms
% FFI stress drop
% FFI ASPECT top third botton third
% FFI Width top third botton third
% earthquake depth (shallow, intermediate, and deep)

PBnames= flip({    'None', ...
    'Intraplate', ...
    'Subduction', ...
    'Ocean transform', ...
    'Speading ridge', ...
    'Ocean convenrgent boundary', ...
    'Continental Transform',... 
    'Continental Rift',...
    'Contiental convergent boundary',...
    'High aspect ratio',...
    'Low aspect ratio',...
    'High stress drop', ...
    'Low stress drop',...
    'High normalized width',...
    'Low normalized width'});

rowNames = {'Shallow (<55km)', ...
    'Intermediate (55-300km)', ...
    'Deep (>300km)', ...
    'Strike-slip', ...
    'Normal', ...
    'Reverse',...
    'Yougn crust (<40Ma)',...
    'Old crust (>100Ma)',...
    PBnames{:}};

I55 = depth < 55;
depthShallow        = res(depth < 55);
depthIntermediate   = res(depth < 300 & depth > 55);    
depthDeep           = res(depth >300);
SS                  = res(fms == 1 & I55);
RF                  = res(fms == 3 & I55);
NF                  = res(fms == 2 & I55);
ageYougn            = res(age < 40 & I55);
ageOld              = res(age > 100& I55);

highAspect          = FFIres(A > prctile(A,80));
lowAspect           = FFIres(A < prctile(A,20));
highSD              = FFIres(SD> prctile(SD,80));
lowSD               = FFIres(SD< prctile(SD,20));
highWnorm           = FFIres(Wnorm > prctile(Wnorm,80));
lowWnorm            = FFIres(Wnorm < prctile(Wnorm,20));



PBopt      = unique(PB);

PBres = cell(length(PBopt),1);
for n = 1:length(PBopt)
    I               = strcmp(PB,PBopt(n));
    PBres{n}            = res(I  & I55);
end

featureArray = {depthShallow,depthIntermediate,depthDeep,SS,RF,NF,ageYougn, ...
    ageOld, highAspect, lowAspect, highSD, lowSD, highWnorm, lowWnorm, PBres{:}};
Nfeat = length(featureArray);
[errUp,errDown,meds] = deal(zeros(Nfeat,1));

for n = 1:Nfeat
    err = prctile(featureArray{n},[25,75]);
    errUp(n)    = err(1);
    errDown(n)  = err(2);
    meds(n) = median(featureArray{n});
end


T = table(rowNames',featureArray',meds,errUp,errDown,'VariableNames',{'Subset','Data','Median','ErrDown', 'ErrUp'});
T = sortrows(T,'Median');

figure; hold on;
sz = 50;
for n = 1:height(T) 
    ires            = T.Data{n};
    scatter(ires,        n*ones(size(ires)),  sz*2,'filled','MarkerFaceColor',[0.5,0.5,0.5], 'MarkerFaceAlpha',0.1)
    scatter(median(ires),n,                   sz,  'filled','MarkerFaceColor', 'k')
    
    errBar = [prctile(ires,25),prctile(ires,75)];
    if isinf(errBar(1))
        newErrBar = linspace(min(ires(~isinf(ires))),errBar(2),10)';
        alpha_values = linspace(0,1,10)';
        patch(newErrBar,repmat(n,10,1),'red','EdgeColor','k',...
            'FaceVertexAlphaData',alpha_values,'AlphaDataMapping','none',...
            'EdgeAlpha','interp', 'LineWidth',2);
    else
        plot(errBar,[n,n],'LineWidth',2,'Color', 'k')
    end
     
end
set(gca,'YTick',1:n,'YTickLabel',T.Subset,'YLim',[0.1,n + 0.1])
ytickangle(60)

xlabel('Relative productivity')
view(90, -90)
grid on


end

function plot_ML(predictors,res,lat,lon,depth)

[resP,R,mdl] = train_SVM(predictors,res);
R;

mm = minmax(res');

figure
subplot(1,2,1)
xval_knn_prediction(lat,lon,depth,res)
xlim(mm)
ylim(mm)
title({'k-nearest neighbor';'site effects'}, 'FontWeight','normal')

subplot(1,2,2)
% worldmap_base;
% scatterm(lat,lon,100,abs(resP-res),'filled')

hold on; 
plot(mm,mm,'--','LineWidth',2,'Color','k');
scatter(res,resP,[],'filled','MarkerFaceColor',[0.5 0.5 0.5], 'MarkerFaceAlpha',0.5);
legend({'1:1','Mainshocks'},'Location','southeast')
xlim(mm)
ylim(mm)

xlabel('Relative Productivity')
title({'SVM'; 'site and source effects'},'FontWeight','normal')

end 

%% other functions

function C = assign_color(X,c)
nX = length(X);
C = zeros(nX,3);
for n = 1:nX
    C(n,:) = c{X(n)};
end
end

function    [magArray,prodArray,K,A]            = productivity_law(MAGNITUDE,AFTERSHOCKS)

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