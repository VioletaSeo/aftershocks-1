% reproduce Nadav's figure

fm = mergedCat.FocalMechanism;
figure


SD = mergedCat.StressDrop;
res= mergedCat.MSres_appended_cat1;

pct= prctile(log10(SD),[33,66]);
g1 = log10(SD) < pct(1);
g3 = log10(SD) > pct(2);

figure; hold on
scatter(SD,res)
set(gca,'XScale','log')
xlabel('Stress Drop')
ylabel('Productivity')

figure; hold on
scatter((SD(g1)),       (res(g1)),      'b')
scatter((SD(g3)),       (res(g3)),      'r')
scatter((SD(~(g3|g1))),  (res(~(g3|g1))),'k')
hold on
scatter(mean(SD(g1)),nanmean(res(g1)),'b','filled')
scatter(mean(SD(g3)),nanmean(res(g3)),'r','filled')
set(gca,'XScale','log')
xlabel('Stress Drop')
ylabel('Productivity')

I  = strcmp(fm,'Reverse');
figure
scatter(SD(I),res(I))
hold on
scatter(mean(SD(g1 & I)),nanmean(res(g1 & I)),'b','filled')
scatter(mean(SD(g3 & I)),nanmean(res(g3 & I)),'r','filled')
set(gca,'XScale','log')

% using ye data
[mergedMSYe,I] = merge_eq_catalog(ye12,{'t','Lat','Lon','Hc','Mw'},MScat,MSformat);
SD      = mergedMSYe.StressDrop;
res     = mergedMSYe.MSres_appended_cat1;

pct= prctile(log10(SD),[33,66]);
g1 = log10(SD) < pct(1);
g3 = log10(SD) > pct(2);

figure
scatter(SD,res)
hold on
scatter(mean(SD(g1)),nanmean(res(g1)),'b','filled')
scatter(mean(SD(g3)),nanmean(res(g3)),'r','filled')
set(gca,'XScale','log')