load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;
minMag  = 6;
[ASinfo,k,alpha] = aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'MinMainshockMag',minMag, ...
    'DepthRange',[0,55], ...
    'ReturnCatalog', 'yes', ...
    'SaveCatalog', 'no', ...
    'PlotYN','no');
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
depth= MSCat.depth;
RT  = MSCat.RuptureDuration;
age = get_crust_age(MSCat.lat, MSCat.lon, MSCat.depth);

colors = {[0.7      0.7       0.7], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
c = assign_color(fms+1,colors);

ftsz = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
      
%%
figure; hold on;
scatter(M,prod,30,c,'filled', 'MarkerFaceAlpha',0.3);
set(gca,'Yscale','log');
plot(minmax(M),k*10.^(alpha*minmax(M)),'--k','LineWidth',2);
ftsz(gcf,12)


%%
figure;hold on

edges = min(res(~isinf(res))):0.2:max(res);

for n = [3,1,2]
    I = fms == n;
    subplot(2,1,1); hold on
    histogram(res(I),edges,'FaceColor',colors{n+1},'EdgeColor',[1 1 1])
    subplot(2,1,2); hold on
    fade_plot(res(I),n,colors{n+1});
end
subplot(2,1,1)
set(gca,'Xlim',[-2,2])
set(gca,'XTick',[])
ylabel('N')

subplot(2,1,2)
set(gca,'YTick',[0,1,2,3],'YTickLabel',{' ','Strike-slip','Normal','Reverse'},'YLim',[0.1,4.1])
xlabel('Residual productivity')
set(gca,'Xlim',[-2,2]) 

ftsz(gcf,16)

%%
% p value for
res2 = res(~isinf(res));
fms2 = fms(~isinf(res));
[h,p] = kstest2(res2(fms2 == 1),res2(fms2 == 3)) % ... tiny!!
[h,p] = kstest2(res2(fms2 == 2),res2(fms2 == 3))

%% prod vs depth
figure
scatter(depth,res,'Filled','MarkerFaceAlpha',0.4);

%% prod vs age
figure
scatter(age,res,[],c,'Filled','MarkerFaceAlpha',0.4);

%% prod vs 

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
