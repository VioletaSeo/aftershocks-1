load('IRIS_DMC_with_FMS.mat')
minMag = 6;
[MSCat,k,alpha] = aftershock_productivity_kernel(...
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

%%
M   = MSCat.MSmag;
t   = MSCat.MSt;
res = MSCat.MSres;
prod= MSCat.MSprod;
fms = MSCat.MSfms;

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


%%
figure;hold on
for n = 1:3
    I = fms == n;
    fade_plot(res(I),n,colors{n+1});
end
set(gca,'YTick',[0,1,2,3],'YTickLabel',{' ','Strike-slip','Normal','Reverse'},'YLim',[0.1,4.1])
xlabel('Residual productivity')
ftsz(gcf,12)

%%
% p value for 
[h,p] = kstest2(res(fms == 1),res(fms == 3))
[h,p] = kstest2(res(fms == 2),res(fms == 3))

%%

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
