function plotgr(M)
%%% Plot the Gutenberg-Richter relationship based on a 1 by N vector of
%%% earthquake magnitudes (M). This function certainly can be improved
%%% (b-val, completeness).
mags    = min(M):0.1:max(M);
nMags   = length(mags);
gr      = zeros(nMags,1);

for n = 1:nMags
    gr(n) = sum(M > mags(n));    
end
figure
plot(mags,gr,'.-','LineWidth',2,'MarkerSize',10)
xlabel('Magnitude')
ylabel('N (M>M_i)')
set(gca,'Yscale','log')
end