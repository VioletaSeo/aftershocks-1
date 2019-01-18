% compare background productivity to the productivity of the mainshock
% better to rerun the kernel to avoid mistakenly using the wrong catalog

[k,alpha] = aftershock_productivity_kernel('eq_catalog.mat','DepthRange',[0,55]);
load('mainshock_info.mat');
MSback = get_background_rate();

figure; hold on;
L = {};
Mrange = diff(minmax(MSmag'));
IzeroAS= MSprod == 0;
for iM = 7.5
    I = MSmag == iM;
    if sum(I) > 5
        L = [L,{num2str(iM)}];
        [bootStat] = bootstrp(10000,@(x) mean(x), [MSback(I),MSprod(I)]);
        backArray = bootStat(:,1);
        prodArray = bootStat(:,2);
        resArray = log10(prodArray./mean(MSprod(I)));
        scatter(backArray,resArray,6,'r','filled', 'MarkerFaceAlpha', 0.1)
        scatter(MSback(I),log10(MSprod(I)/mean(MSprod(I))))
        
        stem(MSback(I&IzeroAS),0.5*ones(1,sum(I&IzeroAS)))
    end
end
set(gca,'XScale','log')
legend(L)
