inpStr = {'ForeshockTimeWindow', 10      , ...
                 'DepthRange',          [0,55]  , ...
                 'ShowOverviewYN',      'no'}; % default input
glo_CAT = aftershock_productivity_GCMT_ISC(inpStr{:});

otf_CAT = aftershock_productivity_GCMT_ISC(inpStr{:}, ...
                        'PlateBoundaryClass',   'OTF');
worldmap_base
plotres = @(c) scatterm(c.MSlat,c.MSlon,(c.MSmag-4).^3,c.MSres_appended_cat1);
plotres(otf_CAT)
                    
ctf_CAT = aftershock_productivity_GCMT_ISC(inpStr{:}, ...
    'PlateBoundaryClass',   'CTF');
worldmap_base
plotres(otf_CAT)

%% 
Mc      = 6;

I       = glo_CAT.MSmag > Mc & glo_CAT.MSfms ~= 1;
gloMag  = glo_CAT.MSmag(I);
gloProd = glo_CAT.MSprod_appended_cat1(I);
gloRes  = glo_CAT.MSres_appended_cat1(I);

ISS     = otf_CAT.MSfms == 1 & otf_CAT.MSmag > Mc;
otfMag  = otf_CAT.MSmag(ISS);
otfProd = otf_CAT.MSprod_appended_cat1(ISS);
otfRes  = otf_CAT.MSres_appended_cat1(ISS);

ISS     = ctf_CAT.MSfms == 1 & ctf_CAT.MSmag > Mc;
ctfMag  = ctf_CAT.MSmag(ISS);
ctfProd = ctf_CAT.MSprod_appended_cat1(ISS);
ctfRes  = ctf_CAT.MSres_appended_cat1(ISS);

%%

figure; hold on;
plotProd = @(mag,prod) scatter(mag,prod,'filled','MarkerFaceAlpha',0.5);

plotProd(gloMag,gloProd);
plotProd(otfMag,otfProd);
plotProd(ctfMag,ctfProd);
set(gca,'YScale','log')
legend({'Global Catalog','Oceanic Transform Fault', 'Continental Transform Fault'});

%%
figure; hold on
colors = {[     0    0.4470    0.7410],...
          [0.6350    0.0780    0.1840],...
          [     0    1         0     ]};
sz     = 100;
N = 1;
scatter(otfRes,N*ones(size(otfRes)),sz,'filled','MarkerFaceAlpha',0.05,'MarkerFaceColor',colors{N})
scatter(median(otfRes),N,sz,'filled','MarkerFaceColor',colors{N})
plot(prctile(otfRes,[35,65]),[N,N],'LineWidth',2,'Color',colors{N})
N = 2;
scatter(ctfRes,N*ones(size(ctfRes)),sz,'filled','MarkerFaceAlpha',0.05,'MarkerFaceColor',colors{N})
scatter(median(ctfRes),N,sz,'filled','MarkerFaceColor',colors{N})
plot(prctile(ctfRes,[35,65]),[N,N],'LineWidth',2,'Color',colors{N})
N = 3;
scatter(gloRes,N*ones(size(gloRes)),sz,'filled','MarkerFaceAlpha',0.05,'MarkerFaceColor',colors{N})
scatter(median(gloRes),N,sz,'filled','MarkerFaceColor','k')
plot(prctile(gloRes,[35,65]),[N,N],'LineWidth',2,'Color','k')

set(gca,'YTick',[0,1,2,3],'YTickLabel',{' ','OTF','CTF','Dip Slip'},'YLim',[0.1,3.1])
xlabel('Residual productivity')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

