function plotgr(M)
%%% Plot the Gutenberg-Richter relationship based on a 1 by N vector of
%%% earthquake magnitudes (M). This function certainly can be improved
%%% (b-val, completeness).
inc     = 0.1; 
mags    = min(M):inc:(max(M)-0.1);
nMags   = length(mags);
gr      = zeros(nMags,1);

for n = 1:nMags
    gr(n) = sum(M > mags(n));    
end

mc = KS_completeness(M,'yes');

figure
hold on
plot(mags,gr,'.-','LineWidth',2,'MarkerSize',10)
xlabel('Magnitude')
ylabel('N (M>M_i)')
set(gca,'Yscale','log')

utsu_correction = inc/2;
b = 1/log(10) * 1/mean(M(M>=mc)-mc+utsu_correction);
I = mags >= mc;
a = log10(sum(M>mc)); %mean(log10(gr(I).*10.^(b*mags(I)')));
plot(mags(I),10.^(a-b*(mags(I)-mc)),'--')

legend({'',sprintf('b = %0.2g, M_c = %0.2g',b,mc)})

end