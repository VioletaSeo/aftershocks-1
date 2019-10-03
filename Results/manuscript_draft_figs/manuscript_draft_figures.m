 % This file (hopefully) contains manuscript-ready figures for the
% aftershock productivity project. The structure should be very similar to
% that of the proposal script
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preamble:
clear
close all
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
minMag      = 6.3; % Minimum magnituce of mainshocks considered (Mc + 2)
maxDepth    = 55;  % Maximum mainshock depth considered         (No deep focus earthquakes)
completeness= 4.5; % Completeness of the catalog                  (from KS test for Mid ocean rideges)
distance2pb = 400; % Maximum distance to plate                  ~ maxDepth/tand(10) 
startDate   = 1990;% start of the catalog                       big drop in MC around start of 1990s 

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
            'Completeness', completeness};

% aftershock statistics
[ASinfo,k,alpha] = aftershock_productivity_kernel(inputCell{:});

%%
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;

MSCat = MSCat(MSCat.fms ~= 0,:); % should I do this??????????????????????????

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
prod    = MSCat.MSprod;
res     = get_res(M,prod);

% source properties
fms     = MSCat.fms;
E       = MSCat.BBtotalEnergy;
Ehf     = MSCat.HFtotalEnergy;
Enormall= E./mw2mo(M);
dur     = MSCat.RuptureDuration;
durNormall = dur./(mw2mo(M)).^(1/3);
RT      = MSCat.RuptureDuration;

% assign ages
age = get_crust_age(MSCat.lat, MSCat.lon, 0);
rate= get_crust_rate(MSCat.lat, MSCat.lon, 0);

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

load('FSPcat_manuscript_draft.mat')

mscatTimeShifted = table(t,lat,lon,M,fms,res,prod,E,Ehf,RT);
FSPformat = FSPcat.Properties.VariableNames([2:3,5]); 
MSformat = mscatTimeShifted.Properties.VariableNames(2:4);
[mergedCat,I] = merge_eq_catalog(FSPcat, FSPformat, [1,1,1.5], ...
                                 mscatTimeShifted,MSformat);
c = mergedCat(I,:);
c = c(c.fms_appended_cat1 ~= 0,:);
c = c(abs(c.t_appended_cat1+365-c.Time)<1,:);
c = c(c.Time>datenum(startDate,1,1),:);
c.DateTime = datetime(c.Time,'ConvertFrom','datenum');
c = c(~(c.res_appended_cat1 == max(c.res_appended_cat1)),:); % remove tuhoku foreshock
 
Mo       = c.Mo;
Wnorm    = c.Width.*Mo.^(-1/3);
Lnorm    = c.Length.*Mo.^(-1/3);
areaNorm = c.Length.*c.Width./(10.^(0.59*c.Mw));
logAnorm = log10(areaNorm);
A        = c.AspectRatio;
A(A<1)   = 1./A(A<1);
logA     = log10(A);
V        = c.RuptureVelocity;
SD       = c.StressDrop;
logSD    = log10(SD);
Yng      = c.Young;
P        = c.Poisson;
H        = c.Heterogeneity./Mo.*(Yng./2.*(1+P));
Enorm    = c.E_appended_cat1./Mo;
durNorm  = c.RT_appended_cat1;
FFIage   = get_crust_age(c.Lat,c.Lon,0);
FFIdip   = c.Dip;
FFIrake  = c.Rake;
FFIdepth = c.Depth;
FFIlitho = get_litho(c.Lat,c.Lon);
FFIfms   = c.fms_appended_cat1;
FFIprod  = c.prod_appended_cat1;
FFIres   = c.res_appended_cat1;

colorAryFFI = assign_color(FFIfms+1,colors);

load('fsp_multi_rupt_manuscript_draft.mat')
cM = array2table(zeros(length(fspMultRupt),4),'VariableNames',{'t','lat','lon','M'});
for n = 1:length(fspMultRupt)
    eq = fspMultRupt(n);
    cM.t(n)     = eq.DateNum;
    cM.lat(n)   = eq.Location.LAT;
    cM.lon(n)   = eq.Location.LON;
    cM.M(n)     = eq.Size.Mw;
end


cMformat = {'lat','lon','M'}; 
MSformat = mscatTimeShifted.Properties.VariableNames(2:4);
[cM,I] = merge_eq_catalog(cM,   cMformat, [1,1,1.5], ...
                          mscatTimeShifted,MSformat);
cM = cM(I,:);
cM = cM(cM.fms_appended_cat1 ~= 0,:);
cM = cM(abs(cM.t_appended_cat1+365-cM.t)<1,:);
cM = cM(cM.t>datenum(startDate,1,1),:);
cM.DateTime = datetime(cM.t,'ConvertFrom','datenum');

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

%% 1) productivity law 
plot_productivity_law(M,prod,res,cAry)
ftsz(gcf,8)
setsize(gcf,4.5,3.5)
savefigure(gcf,'prod_law_mw5',SAVEFIG)

%% 2) Space-time sensitivity
plot_ST_sensitiviety([inputCell,{'DepthRange',[0,maxDepth]}], CAT, 100)
ftsz(gcf,8)
setsize(gcf,6,6)
%%

savefigure(gcf,'sensitivity',SAVEFIG)

%% 3) World map
plot_worldmap_prod(wholeCat.lat, wholeCat.lon,wholeCat.MSres)
ftsz(gcf,8)
setsize(gcf,7.5,4.2)
savefigure(gcf,'worldmap_res',SAVEFIG)

%% 4) Prod vs depth
% to do
plot_prod_vs_depth(wholeCat.depth,wholeCat.MSres, assign_color(wholeCat.fms+1,colors),maxDepth)
ftsz(gcf,8)
setsize(gcf,3.5,2)
savefigure(gcf,'prod_vs_depth',SAVEFIG)

%% 5) regional maps
plot_regions(lat,lon,res)
ftsz(gcf,8)
setsize(gcf,6.6,2.5)
savefigure(gcf,'regions',SAVEFIG)

%% 6) Productivity by plate boundary
plot_prod_by_pb(PB,res)
ftsz(gcf,8)
setsize(gcf,2.5,2.5)

savefigure(gcf,'prod_by_pb',SAVEFIG)

%% 7) Prod vs age
plot_prod_vs_age(age,res,M,cAry)
ftsz(gcf,9);
setsize(gcf,3.5,2.3);
savefigure(gcf,'prod_vs_age',SAVEFIG)

%% 8-9) Stem - covariance matrix
X = [log10(areaNorm), log10(Wnorm),log10(Lnorm),log10(H),(A),V,log10(SD),Yng,P,log10(Enorm),log10(durNorm),FFIdip]; Y = FFIres;
SourceParameters = {'Area^*','Width^*','Length^*','log(H)','Aspect','Rupt. Vel.','log(\Delta\sigma)','Young','Poisson','Energy^*', 'Duration^*','Dip'};
%% 
plot_r2_stems(X,Y,SourceParameters,colorAryFFI,SD,A,Wnorm)
ftsz(gcf,8);
setsize(gcf,6.5,4);
%%
savefigure(gcf,'stem_plot',SAVEFIG)

%% 
plot_cov_matrix([X,FFIdepth,log10(FFIage),Y],[SourceParameters,{'Depth','Crust age','\Delta log(N)'}])
ftsz(gcf,10);
setsize(gcf,3.8,3.4);

%%
savefigure(gcf,'covariance_plot',SAVEFIG)

%% 10) FMS dependence

fms_pairs(MSCat,res, 200,colors)
ftsz(gcf,8);
setsize(gcf,4.5,2.2);
savefigure(gcf,'fmspairs',SAVEFIG)


%% 11) Caltech plot
caltech_plot(wholeCat.depth,wholeCat.age,wholeCat.fms,wholeCat.pb,wholeCat.MSres, A, Wnorm, SD, FFIres, cM.res_appended_cat1)
ftsz(gcf,10);
setsize(gcf,5,4.5);
savefigure(gcf,'cal_tech',SAVEFIG)

%% 12) kn vs SVM
I = ~isnan(FFIage);
pred = [(FFIage),(FFIdip),log10(areaNorm)];
plot_ML(pred(I,:),FFIres(I),c.Lat(I),c.Lon(I),c.Depth(I), ...
    lat,lon,depth,res)
setsize(gcf,5.2,2.4)
ftsz(gcf,8);
savefigure(gcf,'response','yes')

%% SUPPLEMENT:
% scaling
% prod vs depth?
% sensitity to Mc
% co-variance matrix
% each parameter as a function of N*
%   rate
%   energy
%   time
%   ...

SAVEFIG = 'no';


%% scaling 
plot_scaling_comparison(CAT.M,CAT.fms,ASinfo.ID,ASinfo.MSprod,inputCell)
ftsz(gcf,10);
setsize(gcf,5.8,3)
savefigure(gcf,'scaling_comparison','yes')

%% plate age by faulting mech.
faultingStyle = {'Strike-slip','Normal','Reverse'};
for n = 1:3
I = fms == n;
plot_prod_vs_age(age(I),res(I),M(I),cAry(I,:))
ftsz(gcf,12);
setsize(gcf,5,3);
savefigure(gcf,['prod_vs_age_',faultingStyle{n}],SAVEFIG)
end

%% completeness (run seperately with Mc = 5)
% will be done manually

%% covariance matrix

%% N* as a function each parameter
it(-1); % it is an iterator that adavances everytime it is called ;)
pf = @(x,nm) pr(x,res,cAry,nm,[5,3],it);
pf(age              ,'age (Ma)');
pf(rate             ,'rate (mm/yr)');
pf(Enormall         ,'Energy*'); set(gca,'xscale','log');
pf(durNormall       ,'Duration*');set(gca,'xscale','log');

pf = @(x,nm) pr(x,FFIres,colorAryFFI,nm,[5,3],it);
pf(P                ,'Poisson''s ratio');
pf(Yng              ,'Young''s');
pf(Lnorm            ,'Length*');            set(gca,'xscale','log');
pf(Wnorm            ,'Width*');             set(gca,'xscale','log');
pf(A                ,'Aspect');
pf(Enorm            ,'Energy^*');           set(gca,'xscale','log');
pf(durNorm          ,'Duration^*');         set(gca,'xscale','log');
pf(FFIdip           ,'Dip ^o');
pf(FFIrake          ,'Rake ^o');
pf(H                ,'Heterogeneity');      set(gca,'xscale','log');
pf(V                ,'Rupture Velocity (km/s)'); 
ftsz(gcf,8)
setsize(gcf,6,7)
savefigure(gcf,'misc',SAVEFIG)


%% covariance
plot_cov(FFIage,log10(areaNorm),FFIdip,FFIres,FFIfms,colorAryFFI,ftsz,setsize,@savefigure, SAVEFIG)

%% PLOTTING FUNCTIONS:

%%plots:
function plot_productivity_law(M,prod,res,c)

tpos = [-0.15,1];
ptsz = 25;
alpha= 0.15;

figure

% plot a)
subplot(3,3,[1:2,4:5])
t = title('a)'); set(t,'Position',tpos,'Units','normalized')
hold on

[magArray,prodArray,K,A] = productivity_law(M,prod);
scatter(M,prod,ptsz,c,'filled','MarkerFaceAlpha',alpha);
set(gca,'Yscale','log');
squareH = scatter(magArray,prodArray,'s','k');
blankH  = scatter([],[],'k');
lineH   = plot(magArray,K*10.^(A*magArray),'--k','LineWidth',2);
ylabel('Number of Aftershocks, N')
set(gca,'xlim', [min(M)-0.2,9.5],'XTickLabel',[])

kexp = log10(K);
lg = legend([blankH,squareH,lineH],{'Mainshock','Median',sprintf('10^{%0.1f M %0.1f}',A,kexp)}, ...
    'Location','southeast');

newPosition = [0.48 0.42 0.1 0.1];
newUnits = 'normalized';
set(lg,'Position', newPosition,'Units', newUnits);

% plot b)
subplot(3,3,7:8)
t = title('b)'); set(t,'Position',tpos,'Units','normalized')

hold on
scatter(M,res,ptsz,c,'filled','MarkerFaceAlpha',alpha);
set(gca,'xlim', [min(M)-0.2,9.5])
xlabel('M_w')
ylabel({'Relative', 'productivity','\Delta log(N)'})

% plot c)
subplot(3,3,9); hold on
t = title('c)'); set(t,'Position',tpos,'Units','normalized')
histogram(res,30,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
view(90, -90)
set(gca,'XTickLabel',[])
ylabel('Number of Mainshocks')
end

function plot_ST_sensitiviety(inputCell, CAT, Nit)


nT = 3;
timeSelAry = [10,60,300];
timeWinAry = 100/60 * timeSelAry;

nDX = 3;
distSelAry = [1,3,10];
distWinAry = 4/3 * distSelAry;
count = 0;
figure
for n = 1:nT
    for m = 1:nDX
        count = count+1;
        %subplot(nT,nDX,count)
        subplot('Position',[(m-1)*0.9/nDX+0.07+0.02 (nT-n)*0.9/nT+0.07 1/nDX*0.85 1/nT*0.85])
        timeWin = {...
            'TimeWindow', timeWinAry(n), ...
            'TimeSelectionWindow',timeSelAry(n), ...
            'SearchRadius',distWinAry(m),...
            'SelectionRadius',distSelAry(m)};
        
        % get stats & plot
        [ASinfoSENS,~,~] = aftershock_productivity_kernel(inputCell{:},timeWin{:});
        SENScat = CAT(ASinfoSENS.ID,:);
        badDATA = SENScat.fms == 0;
        SENScat = SENScat(~badDATA,:);
        ASinfoSENS = ASinfoSENS(~badDATA,:); 
        scatter(SENScat.M,ASinfoSENS.MSprod,10,'Filled','MarkerFaceAlpha',0.1); hold on
        [magArray,prodArray,kSENS,alphaSENS] = productivity_law(SENScat.M,ASinfoSENS.MSprod);
        scatter(magArray,prodArray,'s','Filled','MarkerFaceColor',[0    0.4470    0.7410], 'MarkerFaceAlpha',0.5)
        p1 = plot(SENScat.M,kSENS*10.^(alphaSENS*SENScat.M),'Color',  [0    0.4470    0.7410], 'LineWidth',2);
        
        text(8.6,1.5,sprintf('N_{x}:%g',sum(ASinfoSENS.MSprod==0)),'Color','r')
        
        % shuffle test
        nMag =  length(magArray);
        SFLprod = zeros(Nit,nMag);
        [SFLk,SFLalpha] = deal(zeros(Nit,1));
        
        %if Nit > 1; parforArg = 0; else; parforArg = Inf; end
        for nSFL = 1: Nit %for (nSFL = 1: Nit, parforArg)
            [ASinfoSENS,~,~] = aftershock_productivity_kernel( inputCell{:},timeWin{:},'ShuffleYN','yes');
            SENScat = CAT(ASinfoSENS.ID,:);
            % scatter(SENScat.M,ASinfoSENS.MSprod,10,'Filled','MarkerFaceAlpha',0.05,'MarkerFaceColor',[.4 .4 .4]);
            [~,SFLprod(nSFL,:),SFLk(nSFL),SFLalpha(nSFL)] = productivity_law(SENScat.M,ASinfoSENS.MSprod);
        end
        prodArray = median(SFLprod);
        k = median(SFLk);
        a = median(SFLalpha);
        
        scatter(magArray,prodArray,'sk','MarkerEdgeAlpha',0.5)
        p2 = plot(SENScat.M,k*10.^(a*SENScat.M),'k','LineWidth',2);
        
        % legends
        leg = legend([p1,p2],{sprintf('N=10^{%0.1f M %0.1f}',alphaSENS, log10(kSENS)), ...
                sprintf('N_{SFL}=10^{%0.1f M %0.1f}',a, log10(k))}, ...
                'Location','northwest','Box','off');
        
        leg.ItemTokenSize = [15,18];
        
        if n == 1
            title([num2str(distSelAry(m)), ' R_{source}'])
        end
        
        if m == 1
            if n == 2
                ylabel({'Number of Aftershocks,', [num2str(timeSelAry(n)), ' days']})
            else
                ylabel([num2str(timeSelAry(n)), ' days'])
            end
        else
            set(gca,'YTickLabel',[])
        end
            
        if n == 3 && m == 2 
            xlabel('Mainshock magnitude, M_w')
        elseif n ~= 3
            set(gca,'XTickLabel',[])
        end
        
        set(gca,'YScale','log')
        if n==2 && m==2
           set(gca,'Color',[0.92 1 0.9]) 
        end 
        
        xlim([min(SENScat.M),max(CAT.M)])
        ylim([0.8,10000])
    end
end

end

function plot_worldmap_prod(lat,lon,res)

worldmap_base;
sh  = scatterm(lat,lon,20,res, 'filled');

% colormap(rgb)
set(gca,'CLim',[-1, 1])
sh.Children.MarkerFaceAlpha = .6;
h = colorbar;
h.Location = 'south'; 
set(h,'position',[.87,.122,.08,.03])

ylabel(h,{'Relative productivity','\Delta log(N)'})

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

end

function plot_prod_vs_depth(depth,res,c,maxDepth)
figure; hold on
scatter(depth, res,...
        10,c, ...
        'Filled','MarkerFaceAlpha',0.4);
ax = gca;
yLim = ax.YLim;
blankH  = scatter([],[],'k');
lineH = plot([maxDepth maxDepth],yLim,'LineStyle','--','Color',[0 0 0 0.4],'LineWidth',2);
ylim(yLim)
xlabel('Depth (km)')
ylabel('Relative Productivity (\Delta log(N))')
set(gca,'XScale','log')
xLim = ax.XLim;
ax.XLim = [3,xLim(2)];
legend([blankH,lineH],{'Mainshock','Depth cutoff'},'Location','northeast')

end

function plot_regions(lat,lon,res)

inp = {{'USA','Mexico','Ecuador'},'no'};

figure
ax1 = subplot('Position',[0.005,0,0.496,1]);
worldmap_base(inp{:});
sh  = scatterm(lat,lon,[],res, 'filled');
sh.Children.MarkerFaceAlpha = .7;
colormap(ax1,'parula')
h = colorbar;
h.Location = 'south'; 
h.Limits = [-1.5,1.5];
set(h,'position',[.37,.8,.06,.03])

ylabel(h,'\Delta log(N)')
t = title('a)'); 
set(t,'Units','normalized','Position',[0.05,0.9])

textm(37,-118.8, 'SAF')
textm(54,-127, 'QCF')
textm(11,-100,'Cocos')
textm(39,-138,'MTJ') 
textm(45,-179,'Aleutian Arc')
textm(27,-125,'GOC')

ax2 = subplot('Position',[0.506,0,0.496,1]);
worldmap_base(inp{:})
load('age.mat'); 
ageAry = age_3_6(~isnan(age_3_6(:,3)),:);
scatterm(ageAry(:,2),ageAry(:,1),2,ageAry(:,3))
colormap(ax2,'hsv')
h = colorbar;
h.Limits = [0,150];
h.Location = 'south'; 
set(h,'position',[.9,.8,.06,.03])
ylabel(h,'age (Ma)')
t = title('b)'); 
set(t,'Units','normalized','Position',[0.05,0.9])

end

function plot_prod_by_pb(PB,res)
% figure; hold on;
% colors     = get(gca, 'ColorOrder');
% sz         = 50;
% N          = 0;
% PBopt      = unique(PB);
% PBlongName = {'None', ...
%     'Intraplate', ...
%     'Subduction', ...
%     'Ocean transform', ...
%     'Speading ridge', ...
%     'Ocean convenrgent', ...
%     'Continental Transform',... 
%     'Continental Rift',...
%     'Contiental convergent'};
% 
% for n = 1:length(PBopt)
%     I               = strcmp(PB,PBopt(n));
%     ires            = res(I);
%     numNegInf(n)    = sum(isinf(ires))/length(res);
%     scatter(ires,        n*ones(size(ires)),  sz*2,'filled','MarkerFaceColor',  colors(mod(n,7)+1,:), 'MarkerFaceAlpha',0.1)
%     scatter(median(ires),n,                   sz,  'filled','MarkerFaceColor',  colors(mod(n,7)+1,:)/1.4)
%     plot([prctile(ires,25),prctile(ires,75)],[n,n],'LineWidth',2,'Color',       colors(mod(n,7)+1,:)/1.4)
%     tickName{n} = sprintf('%s', PBlongName{end+1-n});  
% end
%     
% set(gca,'YTick',1:n,'YTickLabel',tickName','YLim',[0.1,n + 0.1])
% xlabel('Relative productivity')

figure; hold on;

% subplot(1,3,1:2)
% worldmap_base

rowNames= flip({    'Intraplate', ...
    'Intraplate^*', ...
    'Subduction', ...
    'Ocean transform', ...
    'Speading ridge', ...
    'Ocean convergent', ...
    'Continental transform',... 
    'Continental rift',...
    'Continental convergent'});

PBopt      = unique(PB);

PBres = cell(length(PBopt),1);
for n = 1:length(PBopt)
    I               = strcmp(PB,PBopt(n));
    PBres{n}            = res(I);
end

featureArray = PBres;
Nfeat = length(featureArray);
[errUp,errDown,meds,neq,pVal] = deal(zeros(Nfeat,1));
ni = @(x) x(~isinf(x));

for n = 1:Nfeat
    iRes = featureArray{n};
    err = prctile(featureArray{n},[25,75]);
    errUp(n)    = err(1);
    errDown(n)  = err(2);
    meds(n) = median(featureArray{n});
    neq(n)  = length(featureArray{n});
    [~,pVal(n)]    = kstest2(ni(iRes),ni(res));
end

T = table(rowNames',featureArray,meds,errUp,errDown,neq,pVal,'VariableNames',{'Subset','Data','Median','ErrDown', 'ErrUp','N_eq','p'});
T = sortrows(T,'Median');
T = T(T.N_eq>3,:);

% subplot(1,3,3); hold on
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
ytickangle(50)

xlabel({'Relative productivity','\Delta log(N)'})
view(90, -90)
grid on

end

function plot_prod_vs_age(age,res,M,c)

fh = figure;
left_color = [0.8 0.1 0.1];
right_color = [0 0 0];
set(fh,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis right
ylim([0,250]);
yticks([0 50 100])
hold on
edges       = 0:10:max(age);
ylabel({'% no aftershocks'})
I = isinf(res);
[Nno,~] = histcounts(age(I) ,edges);
[Ntot,~] = histcounts(age,edges);
frac = Nno./(Ntot);
Inan = isnan(frac);
frac(Inan) = 0;
frac(Ntot<5) = 0;
frac = frac*100;

histogram('BinEdges',edges,'BinCounts',frac,'EdgeColor','none','FaceAlpha',0.5);

yyaxis left; ylim([-3,2]);
yticks([-2 -1 0 1 2])

scatter((age),res, (M-4.5).^4,c,'filled','MarkerFaceAlpha',0.2);
xlabel('Sea floor age (Ma)')
ylabel({'Relative Productivity','\Delta log(N)'})
r = 10;
[Mage,Mres,Merr] = mov_mean((age),res,ceil(range(age)/(2*r)),r);
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
t = title('a)'); 
set(t,'Units','normalized','Position',[-0.04,0.9])

subplot(2,3,4)
scatter(SD,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
set(gca,'XScale','log')
xlabel('Stress Drop, \Delta \sigma_E [MPa]')
ylabel({'Relative Productivity','\Delta log(N)'})
set(gca,'Xtick',[10^-1,1,100])
t = title('b)'); 
set(t,'Units','normalized','Position',[-0.12,0.98])

subplot(2,3,5)
scatter(Wnorm,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
set(gca,'XScale','log')
set(gca,'Yticklabel','')
xlabel('Normalized width, W*')
t = title('c)'); 
set(t,'Units','normalized','Position',[-0.12,0.98])

subplot(2,3,6);

scatter(A,res,[],colorAry,'filled','MarkerFaceAlpha',0.4)
xlabel('Aspect Ratio')
set(gca,'Yticklabel','')
t = title('d)'); 
set(t,'Units','normalized','Position',[-0.12,0.98])


end

function plot_cov_matrix(X,SourceParameters)
% make a plot of the covariance among the predictors for relative
% productivity. This is done by calculating the covariance matrix of X and
% plotting the matrix as an image

X = X(~isnan(mean(X,2)),:);
C = abs(corrcoef(X));
C = tril(C);
C(C==0) = nan;
C = [C;C(end,:)];
SourceParameters = [SourceParameters,'(Target)'];
figure;
imagesc(C,'AlphaData',~isnan(C));
set(gca,...
        'Xtick', 1:length(SourceParameters),...
        'Ytick', 1:length(SourceParameters),...
        'XtickLabel',SourceParameters, ...
        'YtickLabel',SourceParameters)
xtickangle(50)
ch = colorbar;
ch.Location = 'north';
xlabel(ch,{'Pearson Correlation Coef.','$|PCC| = |\frac{cov(p_i,p_j)}{\sigma_i\sigma_j}|$'},'Interpreter','latex')
ch.Position = [0.52,0.82,0.3,0.06];
colormap(cubehelix([],0.7,-0.7,1,1.8,[0.1,0.9],[0.1,0.95]))
set(gca,'TickLength',[0 0])
hold on
rectangle('Position',[0.6,14.5,length(SourceParameters)-0.1-1,1.9],...
    'EdgeColor','r',...
     'LineWidth',1.5)

end

function fms_pairs(MSCat,res,dxCrit,c)
MSCat.MSres = res;
tpos = [-0.20,0.9];

figure

subplot('Position',[0.15,0.25,0.2,0.65]); hold on
for n = [3,1,2]
    I = MSCat.fms == n;
    fade_plot(MSCat.MSres(I),n,c{n+1});
end
ylim([0.1,3.9])
xlim([-1.5,1.5])
yticks(1:3)
xticks([-1, 0, 1])
yticklabels({'Strike-slip','Normal','Reverse'})
xlabel({'Relative productivity','\Delta log(N)'})
view(90, -90)
grid on
ytickangle(60)
t = title('a)'); 
set(t,'Units','normalized','Position',tpos)


subplot('Position',[0.5,0.25,0.4,0.65]); hold on
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
    
    SScolocCAT(n,:) = SSMS(IminSS,:); %#ok<AGROW> %store eq into new array
    DScolocCAT(n,:) = DSMS(IminDS,:); %#ok<AGROW>
    
    [XSS, SSMS] = rem_row(IminSS, XSS, SSMS);
    [XDS, DSMS] = rem_row(IminSS, XDS, DSMS);
    n = n+1;
end

XSS = get_coord(SScolocCAT);
XDS = get_coord(DScolocCAT);
r   = sqrt(sum((XSS-XDS).^2,2));

plot([-1.5,1.5],[-1.5,1.5],'--k','LineWidth',2);
scatter(SScolocCAT.MSres,DScolocCAT.MSres,[],r, 'Filled','MarkerFaceAlpha',0.8)
xlim([-1.5,1.5]); ylim([-1.5,1.5])
xticks([-1, 0, 1])
yticks([-1, 0, 1])
xlabel({'Strike-slip relative productivity', '\Delta log(N)'})
ylabel({'Dip-slip relative productivity','\Delta log(N)'})
cH = colorbar;
ylabel(cH,'Distance (km)');
colormap('gray')
t = title('b)'); 
set(t,'Units','normalized','Position',tpos-[0.2 0])

%legend({'1:1','Pairs'},'Location','southeast')


% subplot(3,3,1:6)
% worldmap_base
% ploteq = @(cat,c) scatterm(cat.lat,cat.lon,(cat.M-4).^3,c,'filled');
% ploteq(SScolocCAT,'b')
% ploteq(DScolocCAT,'r')



end

function caltech_plot(depth,age,fms,PB,res, A, Wnorm, SD, FFIres,multRes)

% yougn age
% intermediate age
% plate boundaries
% fms
% FFI stress drop
% FFI ASPECT top third botton third
% FFI Width top third botton third
% earthquake depth (shallow, intermediate, and deep)

PBnames= flip({'Intraplate', ...
    'Intraplate^*', ...
    'Subduction', ...
    'Ocean transform', ...
    'Speading ridge', ...
    'Convergent (oceanic)', ...
    'Continental transform',... 
    'Continental rift',...
    'Convergent (continental)'});

rowNames = {'Shallow (<55km)', ...
    'Intermediate (55-300km)', ...
    'Deep (>300km)', ...
    'Strike-slip', ...
    'Normal', ...
    'Reverse',...
    'Young crust (<40Ma)',...
    'Old crust (>100Ma)',...
    'High aspect ratio',...
    'Low aspect ratio',...
    'High stress drop', ...
    'Low stress drop',...
    'High normalized width',...
    'Low normalized width',...
    'Multi-segment rupture',...
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

featureArray = {depthShallow,depthIntermediate,depthDeep,SS,NF,RF,ageYougn, ...
    ageOld, highAspect, lowAspect, highSD, lowSD, highWnorm, lowWnorm, multRes, PBres{:}};
Nfeat = length(featureArray);
[errUp,errDown,meds,neq] = deal(zeros(Nfeat,1));

for n = 1:Nfeat
    err = prctile(featureArray{n},[25,75]);
    errUp(n)    = err(1);
    errDown(n)  = err(2);
    meds(n) = median(featureArray{n});
    neq(n)  = length(featureArray{n});
end

T = table(rowNames',featureArray',meds,errUp,errDown,neq,'VariableNames',{'Subset','Data','Median','ErrDown', 'ErrUp','N_eq'});
T = sortrows(T,'Median');
T = T(T.N_eq>3,:);

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

xlabel({'Relative productivity','\Delta log(N)'})
view(90, -90)
grid on


end

function plot_ML(predictors,res,lat,lon,depth,latAll,lonAll,depthAll,resAll)

I = ~isinf(res);
mm = minmax(res(I)')*1.5;

% KNN prediction
figure
subplot(1,2,1)
wgs84   = wgs84Ellipsoid('kilometer');
loc     = geodetic2ecef(wgs84,latAll, lonAll, -depthAll);
qloc    = geodetic2ecef(wgs84,lat,lon,-depth); 
resPknn = knn_regress(loc,resAll,qloc,10);
hold on; 
pH = plot(mm,mm,'--','LineWidth',2,'Color','k');
scatter(res,resPknn,[],'filled','MarkerFaceColor',[0.5 0.5 0.5], 'MarkerFaceAlpha',0.5);
I2 = ~isinf(resPknn);
I3 = I & I2;
R       = sqrt(sum((res(I3)-resPknn(I3)).^2)/length(res(I3)));
legend(pH,{sprintf('RMSE: %0.2g',R)},'Location','southeast','Box','off')
xlim(mm)
ylim(mm)
ylabel({'Predicted relative productivity','\Delta log(N)'})
xlabel('Relative productivity, \Delta log(N)')
title({'a) k-nearest neighbor (site effects)'}, 'FontWeight','normal');


% SVM prediction
subplot(1,2,2)
predictors = predictors(I,:);
res = res(I);
[resP,R,~] = train_SVM(predictors,res);
hold on; 
pH = plot(mm,mm,'--','LineWidth',2,'Color','k');
scatter(res,resP,[],'filled','MarkerFaceColor',[0.5 0.5 0.5], 'MarkerFaceAlpha',0.5);
legend(pH,{sprintf('1:1, RMSE: %0.2g',R)},'Location','southeast','Box','off')
xlim(mm)
ylim(mm)

xlabel('Relative productivity, \Delta log(N)')
title({'b) SVM (site and source effects)'},'FontWeight','normal');



%  worldmap_base;
%  I = ~isnan(resP);
%  scatterm(lat(I),lon(I),100,abs(resP(I)-res(I)),'filled')
end 

function pr(x,res,c,xName,DIM,N)%, ftsz, setsize, savefigure, SAVEfig)

subplot(DIM(1),DIM(2),N)
scatter(x,res,5,c,'Filled','MarkerFaceAlpha',0.5)
xlabel(xName)

t = title([char(96+N),')']); 
set(t,'Units','normalized','Position',[-0.15,0.96])
ylim([-1.5,1.5])

yyaxis right
ylim([0,250]);
yticks([0 50 100])
hold on
edges       = linspace(nanmin(x),nanmax(x),10);
%ylabel({'Fraction with', 'no aftershocks'})
I = isinf(res);
[Nno,~] = histcounts(x(I) ,edges);
[Ntot,~] = histcounts(x,edges);
frac = Nno./(Ntot)*100;
Inan = isnan(frac);
frac(Inan) = 0;
frac(Ntot<5) = 0;

histogram('BinEdges',edges,'BinCounts',frac,'EdgeColor','none','FaceAlpha',0.5);

switch N
    case 7
        yyaxis left
        ylabel({'Relative Productivity','\Delta log(N)'})
    case 9 
        yyaxis right
        ylabel('% no aftershocks')
end

switch mod(N,DIM(2))
    case 1
        yyaxis right
        set(gca,'Ytick',[])
    case 2 
        yyaxis right
        set(gca,'Ytick',[])
        yyaxis left
        set(gca,'Ytick',[])
    case 3
        yyaxis left
        set(gca,'Ytick',[])
end
end

function plot_cov(age,logAnorm,dip,res,fms,colorMatrix)

figure
[h,ax,~] = gplotmatrix([age,logAnorm,dip,res],[],fms,colorMatrix,[], [],'off');
set(ax(4,:),'Color',[0.95 0.95 0.95])
set(ax(:,4),'Color',[0.95 0.95 0.95])
for n = 1:4
    for m = 1:3
        set(h(n,n,m),'DisplayStyle','bar','FaceColor',colorMatrix(m,:))
    end
end
ax(1,1).YLabel.String = 'Age (Ma)';
ax(2,1).YLabel.String = 'A_{norm}';
ax(3,1).YLabel.String = 'Dip';
ax(4,1).YLabel.String = '\Delta log(N)';

ax(4,1).XLabel.String = 'Age (Ma)';
ax(4,2).XLabel.String = 'A_{norm}';
ax(4,3).XLabel.String = 'Dip';
ax(4,4).XLabel.String = '\Delta log(N)';

ftsz(gcf,8);
setsize(gcf,6,5)
savefigure(gcf,'plotmatrix',SAVEFIG)

end

function plot_scaling_comparison(M,fms,MSID,MSprod,INPUTCELL)
% compare productivity measurements for difference scaling relationships
% inputs:
% MSID: id for productivity measurements


alpha = 0.09;
markerSize = 10;

[ASinfo,~,~] = aftershock_productivity_kernel(INPUTCELL{:},'Scaling','B19');
[I,i1,i2] = intersect(MSID,ASinfo.ID);

prod1 = MSprod(i1);
prod2 = ASinfo.MSprod(i2);
M = M(I);

[prod1,prod2,M] = goodind(fms(I) ~= 0, prod1,prod2,M);

jit1 = 1; %+ (rand(length(prod1),1)-0.5)/2;
jit2 = 1; %+ (rand(length(prod1),1)-0.5)/2;
jit3 = 1; % + (rand(length(prod1),1)-0.5)/10;


figure;
subplot(1,2,1)
scatter(prod1.*jit1,prod2.*jit2,markerSize,'filled','MarkerFaceAlpha',alpha)
set(gca,'XScale','log','YScale','log')
xlabel({'N_{WC94}', 'R_{source} \sim 10^{0.59M_W}'})
ylabel({'N_{B19}', 'R_{source} \sim 10^{0.4951M_W}'})
hold on
Nrange = [1,max([prod1;prod2])];
plot(Nrange,Nrange,'--k','linewidth',1.5)
legend({'Mainshocks','1:1'},'Location','northwest','Box','off')

subplot(1,2,2)
prodR = log10(prod1./prod2);
scatter(M,(prodR.*jit3),markerSize,'filled','MarkerFaceAlpha',alpha)

% mm = (min(M):0.1:max(M))-0.05;
% for n = 1:length(mm)
%     p025(n) = prctile(prodR(M>mm(n) & M < mm(n)+0.1),25);
%     p975(n) = prctile(prodR(M>mm(n) & M < mm(n)+0.1),75);
% end
%     
% hold on
% plot(mm,p025,'--k')
% plot(mm,p975,'--k')
xlabel('Mainshock Magnitude (M_W)')
ylabel('log(N_{WC94}/N_{B19})')
legend(sprintf('Mean abs. error: %0.1g',nanmean(abs(prodR(~isinf(prodR))))), ...
    'Box','off', ...
    'Location','southeast');
% set(gca,'Yscale','log')
disp(sprintf('%0.2g\5 agreement', sum(prod1==prod2)/length(prod1)))

end

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

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);  %#ok<AGROW>
end

end
