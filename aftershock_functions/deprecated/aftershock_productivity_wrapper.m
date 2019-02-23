function aftershock_productivity_wrapper()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTENT:
% entire catalog
% iterate over the distance from the plate boundary
% plot different plate boundaries
% iterate over dept zones
% fms vs PB
% intraplate vs plateboundary
% compare crust thicknesses
% compare focal mechanisms within a single plate boundary

%% SENSITIVITY TESTS
% time window

%% DATA ANALYSIS
% aggregate fit and residuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% entire catalog 
figure
aftershock_productivity_kernel('DepthRange', [0,55]);

%% plot different focal mechanisms
% 
% figure
% 
% maxDepth = 55;
% aftershock_productivity_kernel('FocalMechanism','strike slip', ...
%     'DepthRange', [0,maxDepth]);
% aftershock_productivity_kernel('FocalMechanism','reverse', ...
%     'DepthRange', [0,maxDepth]);
% aftershock_productivity_kernel('FocalMechanism','normal', ...
%     'DepthRange', [0,maxDepth ]);
% 
% title('Mean number of aftershocks as a function of Magnitude for each focal mechanism');

%% plot different plate boundaries
% 
% maxDepth = 55;
% plateBoundaryDist = 100;
% 
% aftershock_productivity_kernel(...
%     'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        'transform', ...
%     'PlateBoundaryDist',    plateBoundaryDist);
% aftershock_productivity_kernel(...
%     'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        'convergent', ...
%     'PlateBoundaryDist',    plateBoundaryDist);
% crustThickness = [0 25 300];
% maxDepth            = 80; % km
% 
% figure
% for n = 1:(length(crustThickness)-1)
%     aftershock_productivity_kernel(...
%                                    'DepthRange',            [0, maxDepth] ,...
%                                    'CrustThicknessRange',   [crustThickness(n),crustThickness(n+1)]);
% end  
% title('Mean number of aftershocks as a function of available crust');

% aftershock_productivity_kernel(...
%     'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        'divergent', ...
%     'PlateBoundaryDist',    plateBoundaryDist);
% 
% title('Mean number of aftershocks as a function of Magnitude for each plate boundary');

 
%% iterate over the distance from the plate boundary
% This is plot of looking wether the productivity of earthquakes depends on
% the distance to the plate boundaries

% N = 4;
% distance2PB = linspace(10,1000,N);
% 
% for n = 1:N
%     aftershock_productivity_kernel('PlateBoundary','transform','PlateBoundaryDist',distance2PB(n))
% end  

% title('Mean number of aftershocks as a function of Magnitude for increasing distance from plate boundary');

%% iterate over dept zones
% This is a plot comparing the productivity of earthquakes in different
% depth zones at convergent boundaries

% depth = [0 15 35 55];
% 
% figure
% for n = 1:(length(depth)-1)
%     aftershock_productivity_kernel('PlateBoundary','concat1.lon(Idx)vergent',        ...
%                                    'FocalMechanism','reverse',          ...
%                                    'PlateBoundaryDist',1000,             ...
%                        % resAllData  = allAS(MAG,NUMEQ);
% resRange    = minmax(resAllData(~isnan(resAllData) & ~isinf(resAllData)));
% 
% binEdges    = linspace(resRange(1),resRange(2),10);
% % binEdges    = logspace(log10(resRange(1)),log10(resRange(2)),10);
% 
% figure
% histogram(resAllData,binEdges,'FaceColor', [0.8 0.8 0.8],'EdgeColor',[1 1 1])
% hold on
% 
% faceColor = get(h,'CData');
% for n = 1:length(xData)
%     mag     = xData{n};
%     numEq   = yData{n};
%     res     = allAS(mag,numEq);  
%     
%     %offset = n*diff(binEdges(1:2))*0.2;
%     offset = 0;
%     histogram(res,binEdges+offset,'FaceColor',faceColor{n},'EdgeColor',[1 1 1])
%     
% end
% 
% xlabel('\Delta log N')
% ylabel('Count')
% 
% title('Residual Analysis')            'DepthRange',[depth(n), depth(n+1)]  );
% end  
% title('Mean number of aftershocks as a function of Magnitude in zones A B and C');
    
%% fms vs PB
% This is a plot comparing the productivity of earthquakes that have
% different focal mechanisms versus those in transform plate boundaries

% figure
% aftershock_productivity_kernel('PlateBoundary','transform',        ...
%                                'PlateBoundaryDist',30, ...
%                                'DepthRange', [0,25]);
% aftershock_productivity_kernel('FocalMechanism','strike slip', ...
%                                'DepthRange', [0,25]);
% title('Mean number of aftershocks as a function of Magnitude for transform vs strike slip eq');
 
%% intraplate vs plateboundary for each focal mechanism
% this compares intra plate earthquakes to plate boundary eartquakes (in
% the sense that they are close or far away from plate boundaries
% 
maxDepth            = 55; % km

aftershock_productivity_kernel(  quakes.otime, ...
                                quakes.lat, ...
                                quakes.lon, ...
                                quakes.depth, ...
                                quakes.mag, ...
                                'PlateBoundary', 'plate boundaries'     , ...
                                'PlateBoundaryDist',100                 , ...
                                'DepthRange', [0,maxDepth]              , ...
                                'MinMainshockMag', 5.5);
aftershock_productivity_kernel(  quakes.otime, ...
                                quakes.lat, ...
                                quakes.lon, ...
                                quakes.depth, ...
                                quakes.mag, ...
                                'PlateBoundary', 'none'                 , ...
                                'PlateBoundaryDist', 500                , ...
                                'DepthRange', [0,maxDepth]              , ...
                                'MinMainshockMag', 5.5);
%                             
%  title('Mean number of aftershocks as a function of Magnitude (intraplate vs PB)')

%% compare crust thicknesses:

% This is a plot comparing the productivity of earthquakes in different
% ranges of crust thicknessaftershock_productivity_kernel('ShowOverviewYN','yes',...
%                                'NewPlotYN','yes');

% 
% crustThickness = [0 7 55];
% maxDepth            = 80; % km
% 
% figure
% for n = 1:(length(crustThickness)-1)
%     aftershock_productivity_kernel(...
%                                    'DepthRange',            [0, maxDepth] ,...
%                                    'CrustThicknessRange',   [crustThickness(n),crustThickness(n+1)]);
% end  
% title('Mean number of aftershocks as a function of available crust');MSres = allAS(MSmag,MSprod);

%% compare focal mechanisms within a single plate boundary

% maxDepth = 55;
% plateBoundaryType = 'convnumberOfIterations = 10;

% plateBoundaryDist = 500;
% 
% aftershock_productivity_kernel('FocalMechanism','strike slip', ...
%     'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        plateBoundaryType, ...
%     'PlateBoundaryDist',    plateBoundaryDist);
% aftershock_productivity_kernel('FocalMechanism','reverse', ...
%     'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        plateBoundaryType, ...
%     'PlateBoundaryDist',    plateBoundaryDist);
% aftershock_productivity_kernel('FocalMechanism','normal', ...
%     'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        plateBoundaryType, ...
%     'PlateBoundaryDist',    plateBoundaryDist);
% 
% title('Mean number of aftershocks as a function of Magnitude');

%% strike slip vs strike slip transform

% maxDepth = 55;
% plateBoundaryDist = 400;
% 
% aftershock_productivity_kernel(...
%  aftershock_productivity_kernel(...
%                                'DepthRange', [0,55]);
%    'DepthRange',           [0,maxDepth], ...
%     'PLateBoundary',        'transform', ...
%     'PlateBoundaryDist',    plateBoundaryDist, ...
%     'FocalMechanism',       'strike slip', ...
%     'DepthRange',           [0,maxDepth]);
% 
% 
% aftershock_productivity_kernel('FocalMechanism','strike slip', 'PLateBoundary','convergent',...
%     'DepthRange', [0,maxDepth]);
% 
% aftershock_productivity_kernel('FocalMechanism','strike slip', 'PLateBoundary','divergent',...
%     'DepthRange', [0,maxDepth]);
% aftershock_productivity_kernel('FocalMechanism','strike slip', 'PLateBoundary','none',...
%     'DepthRange', [0,maxDepth]);

% -> interpretation: the plate boundary environment has effectively no
% influence on the productivity behavior of the 

 
%% SENSITIVTY TESTS
%% time selection window

% maxDepth            = 55; % km
% 
% N = 5;
% timeWindow = logspace(-2,2,N);
% 
% for n = 1:N
%     aftershock_productivity_kernel('TimeSelectionWindow',timeWindow(n), ...
%         'DepthRange', [0,maxDepth]);
% end  
% 
% title('Sensitivity test: time selection window')

%% time search window
% 
% maxDepth            = 55; % km

% 
% N = 5;
% timeWindow = logspace(-2,2,N);
% 
% for n = 1:N
%     aftershock_productivity_kernel('TimeWindow',timeWindow(n), ...
% end  
% 
% title('Sensitivity test: time search window')function    [p,notUsingDefault] = parse_input(input)


%% shuffle times

% numberOfIterations = 10;
% for n = 1:numberOfIterations
%     aftershock_productivity_kernel('ShuffleYN','yes')
% end


%% randomize timeload('mainshock_info.mat')


% numberOfIterations = 5;
% for n = 1:numberOfIterations
%      aftershock_productivity_kernel('SyntheticTimeYN','yes');
% end114.11042944785275

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AGGREGATE DATA ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('mainshock_info.mat')
MSres= getMSres(MSmag,MSprod);

%% Residuals
% 
% resAllData  = allAS(MAG,NUMEQ);
% resRange    = minmax(resAllData(~isnan(resAllData) & ~isinf(resAllData)));
% 
% binEdges    = linspace(resRange(1),resRange(2),10);
% % binEdges    = logspace(log10(resRange(1)),log10(resRange(2)),10);
% 
% figure
% histogram(resAllData,binEdges,'FaceColor', [0.8 0.8 0.8],'EdgeColor',[1 1 1])
% hold on
% 
% faceColor = get(h,'CData');
% for n = 1:length(xData)
%     mag     = xData{n};
%     numEq   = yData{n};
%     res     = allAS(mag,numEq);  
%     
%     %offset = n*diff(binEdges(1:2))*0.2;
%     offset = 0;
%     histogram(res,binEdges+offset,'FaceColor',faceColor{n},'EdgeColor',[1 1 1])
%     
% end
% 
% xlabel('\Delta log N')
% ylabel('Count')
% 
% title('Residual Analysis')

%% Main shock decription
% produce a map of the main shock productivity (in residual form)

% Mcrit = 6.5;
% make_productivity_anomaly_map(MSlat,MSlon,MSmag, MSt, MSfms, Mcrit, MSres);


%% compare productivity and mainshock area (Hayes data)
cat = merge_cat('mainshock_info.mat','eq_catalog_moderate_to_large.mat');
I       = ~isnan(cat.area) & MSmag > 7.5;

% make a maps of the data being analysed below
make_eq_map(I,MSlat,MSlon,MSmag,MSt, MSfms)

% in residual form
areaRes = compare_residuals(cat, MSmag, MSfms, MSres, MSt, I);

% geometrical model fit - area
geometric_model_area(cat,MSmag,MSfms,MSprod,I)

% geometrical model fit - perimeter
geometric_model_perimeter(cat,MSmag,MSfms,MSprod,I)

end

%% blocks:

function MSRES = getMSres(MSMAG,MSPROD)



h = findobj(gca,'Type','scatter');
xData=get(h,'Xdata');
yData=get(h,'Ydata');

if iscell(xData)
    MAG     = [xData{:}];
    NUMEQ 	= [yData{:}];
else
    MAG     = xData;
    NUMEQ   = yData;
end

[prefactor,exponent] = getproductivity(MAG,NUMEQ);
plot(minmax(MAG),prefactor*10.^(exponent*minmax(MAG)),'--','LineWidth',2);


allAS       = @(x,y) log10(y) - log10(prefactor*10.^(exponent*x));
MSRES       = allAS(MSMAG,MSPROD);

end

function make_productivity_anomaly_map(MSlat,MSlon,MSmag,MSt, MSfms, Mcrit, MSres)

critMag = Mcrit;
Imag    = MSmag>critMag;

% go though each focal mechannism
shapes = {'o','s','d'};

% select a certain percentile range
pctThreshold    = 95;
resTheshold     = prctile(MSres(Imag), pctThreshold);
Ipercentile    = MSres > resTheshold;

% create map:
worldmap_base

for n = 1:3
    
    I = (MSfms==n) & Imag & Ipercentile;
    [tempLat,tempLon,tempMag,tempRes,tempT] = goodind(I,MSlat,MSlon,MSmag,MSres,MSt);
    scatterm(tempLat,tempLon,(tempMag-4).^3.5,tempRes,'filled',shapes{n});
    show_eq_info_map(tempLat,tempLon,tempMag,tempT);
    
end


end

function make_eq_map(I,MSlat,MSlon,MSmag,MSt, MSfms)
    worldmap_base
    [MSlat,MSlon,MSmag,MSt,MSfms] = goodind(I,MSlat,MSlon,MSmag,MSt,MSfms);
    scatterm(MSlat,MSlon,(MSmag-5).^4,MSfms,'filled','MarkerFaceAlpha',0.75);
    show_eq_info_map(MSlat,MSlon,MSmag,MSt);
end

function AREARES = compare_residuals(cat, MSmag,MSfms, MSres, MSt, I)

d       = polyfit(MSmag(I),log10(cat.area(I)),1);
func    = @(M,A) log10(A) - log10(10.^(d(1)*M+d(2)));
areaRes = func(MSmag(I), cat.area(I));

figure
scatter(MSres(I), areaRes,(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75)
show_eq_info(MSres(I), areaRes,MSmag(I),MSt(I));

xlabel('residuals in productivity scaling: log(N_i/k10^{\alpha M_i})')
ylabel('residuals in area scaling: log(A_i/b10^{bM_i})')
colormap jet

grid on

hold on 

msres = MSres(I);
msfms = MSfms(I);

I2 = ~isinf(msres);
dres = polyfit(msres(I2), areaRes(I2),1);
plot(msres(I2),dres(1)*msres(I2)+dres(2),'-')

%...directly
figure
scatter(MSmag(I),cat.area(I),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75)
set(gca,'xscale','log','yscale','log')
xlabel('M_w');
ylabel('area m^2');
colormap jet

AREARES = areaRes;

end

function geometric_model_area(cat,MSmag,MSfms,MSprod,I)

Xmat        = [ones(length(cat.stress_drop(I)),1) ,log10(cat.stress_drop(I)')];
Ymat        = log10(MSprod(I)./(cat.area(I)));
infInd      = isinf(Ymat);
Xmat        = Xmat(~infInd,:);
Ymat        = Ymat(~infInd)'; % Caution!! This will affect the fit!

BETAmat     = (Xmat'*Xmat)\(Xmat'*Ymat);

prodModelFit= 10^BETAmat(1).*cat.area(I).*cat.stress_drop(I).^BETAmat(2);

figure
subplot(1,3,1)
scatter(prodModelFit,MSprod(I),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75)
set(gca,'xscale','log','yscale','log')
colormap jet

xLabelStr = sprintf('10^{%g} A SD ^{%g}', [BETAmat(1),BETAmat(2)]);
xlabel(xLabelStr);
ylabel('Number of aftershocks');

grid on


subplot(1,3,2)
scatter(MSmag(I),log10(prodModelFit./MSprod(I)),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75);

xlabel('M_w');
ylabel('Model Residual');

grid on


subplot(1,3,3)
scatter(MSprod(I),log10(prodModelFit./MSprod(I)),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75);

set(gca,'xscale','log')
xlabel('Number of Aftershocks');
ylabel('Model Residual');

grid on

end

function geometric_model_perimeter(cat,MSmag,MSfms,MSprod,I)

Xmat        = [ones(length(cat.stress_drop(I)),1) ,log10(cat.stress_drop(I)')];
Ymat        = log10(MSprod(I)./(cat.rlen(I)+cat.rwid(I)).^2);
infInd      = isinf(Ymat);
Xmat        = Xmat(~infInd,:);
Ymat        = Ymat(~infInd)'; % Caution!! This will affect the fit!

BETAmat     = (Xmat'*Xmat)\(Xmat'*Ymat);

prodModelFitPerim= 10^BETAmat(1).*(cat.rlen(I)+cat.rwid(I)).^2.*cat.stress_drop(I).^BETAmat(2);

figure
subplot(1,3,1)
scatter(prodModelFitPerim,MSprod(I),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75);
set(gca,'xscale','log','yscale','log')
colormap jet

% xLabelStr = sprintf('%g P L SD ^{2/3}', perimModelFitC);
% xlabel(xLabelStr);
xLabelStr = sprintf('10^{%g} (a1+a2)^2 SD ^{%g}', [BETAmat(1),BETAmat(2)]);
xlabel(xLabelStr);
ylabel('Number of aftershocks');

grid on

subplot(1,3,2)
scatter(MSmag(I),log10(prodModelFitPerim./MSprod(I)),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75);

xlabel('M_w');
ylabel('Model Residual');

grid on


subplot(1,3,3)
scatter(MSprod(I),log10(prodModelFitPerim./MSprod(I)),(MSmag(I)-5).^4,MSfms(I),'filled','MarkerFaceAlpha',0.75);

set(gca,'xscale','log')
xlabel('Number of Aftershocks');
ylabel('Model Residual');

grid on

end

%% accessory functions

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);
end

end

function show_ind(x,y,I)
% plot the index number on the pllot as a legend of sorts
indLoc = find(I);
numInd = length(indLoc);

for n = 1:numInd
    text(x(n),y(n),num2str(n));
end

end

function show_ind_map(x,y,I)
% plot the index number on the pllot as a legend of sorts (for maps)
indLoc = find(I);
numInd = length(indLoc);

for n = 1:numInd
    textm(x(n),y(n),num2str(n));
end

end

function show_eq_info_map(lat, lon, Mw, t)
    dim = size(Mw);
    if dim(2) > 1
       Mw = Mw'; 
    end
    
    dateStr = datestr(t);
    yearStr = dateStr(:,8:12);
    eqStr   = [yearStr,num2str(Mw)];
    textm(lat, lon, eqStr, 'FontSize',6,'Rotation',45,'HorizontalAlignment','left');

end

function show_eq_info(x, y, Mw, t)
    dim = size(Mw);
    if dim(2) > 1
       Mw = Mw'; 
    end
    
    dateStr = datestr(t);
    yearStr = dateStr(:,8:12);
    eqStr   = [yearStr,num2str(Mw)];
    text(x, y, eqStr, 'FontSize',6,'Rotation',45,'HorizontalAlignment','left');

end