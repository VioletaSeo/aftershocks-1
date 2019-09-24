CAT = aftershock_productivity_GCMT_ISC('LatRange',[20,45],'LonRange',[55,90]);


%%
hAry   = cell(3,1);
figure; hold on
for n = 1:3
    I = CAT.MSfms == n & CAT.MSmag_appended_cat1 > 6.5;
    gloRes = CAT.MSres_appended_cat1(I);
    hAry{n} = scatter(gloRes,repmat(n,1,sum(I)),50, ...
        'filled','MarkerFaceAlpha',0.3);
    scatter(median(gloRes),n,50,'filled','MarkerFaceColor','k')
    plot(prctile(gloRes,[18,82]),[n,n],'LineWidth',2,'Color','k')

end
legend([hAry{:}],{'Strike-Slip','Normal','Reverse'});

