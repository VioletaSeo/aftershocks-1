function b_val_timeseries(T,M,Twin,Dt,Nmin,Mc)

inc = 0.1;
mags    = min(M):inc:(max(M)-0.1);
nMags   = length(mags);

T0 = min(T);
Tmax = max(T);
Tarray = T0:Dt:Tmax;
WinN = ceil(Twin/Dt);
Nt = length(Tarray);

gr      = zeros(nMags,Nt);
[a,b,err]   = deal(zeros(1,Nt)); 
for iT = WinN:Nt
    
    Mwin = M(T>(Tarray(iT)-Twin) & T<Tarray(iT));
    
    for n = 1:nMags
        gr(n,iT) = sum(Mwin > mags(n));
    end
    
    gr(:,iT) = gr(:,iT)/max(gr(:,iT)); 
    
    if Mc == true
       mc = KS_completeness(Mwin,'no');       
    else
       mc = Mc;
    end
    
    utsu_correction = inc/2;
    b(iT) = 1/log(10) * 1/mean(Mwin(Mwin>=mc)-mc+utsu_correction);
    a(iT) = (sum(Mwin>mc));
    
    if a(iT) > Nmin
        cumDiff = sum((Mwin(Mwin>=mc)-mean(Mwin(Mwin>=mc))).^2);
        err(iT) = 2.3*b(iT)^2 * ...
            sqrt(cumDiff/(a(iT)*(a(iT)-1)));
    end

end

grDiff = gr./(mean(gr,2)*ones(1,Nt));

figure;

subplot(4,1,1)
scatter(T,M,(M+1).^2,'filled','MarkerFaceAlpha',0.5)
xlabel('Time')
ylabel('Magnitude')

subplot(4,1,2)
imagesc([0,1],[min(M),max(M)],log10(grDiff))
ylabel('Magnitude')
ch = colorbar;
ylabel(ch,'log(N(M>M_i))')
xticks([])

subplot(4,1,3)
plot(Tarray,b,'k');
ylabel('b-values')
xlabel('Time')

subplot(4,1,4)
plot(Tarray,a);
set(gca,'yscale','log')
ylabel('a-values')
xlabel('Time')

end