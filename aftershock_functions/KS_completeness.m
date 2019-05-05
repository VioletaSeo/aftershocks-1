function MCOUT = KS_completeness(M,plotYN)

Mmax = max(M);
Mmin = min(M);
inc  = 0.1;
Mvec = Mmin:inc:(Mmax-4*inc);
nM   = length(Mvec);
[k] = deal(zeros(nM,1));

for n = 1:nM
    
    Mc = Mvec(n);
    iM = M(M>Mc);
    %iM = iM + randn(length(iM),1)*0.05;
    
    xmin=10^Mc;
    vX_tmp=10.^sort(iM);
   
    alpha=length(iM)/sum(log(vX_tmp./xmin));
    obsCumul=(0:(length(iM)-1))'/length(iM);
    modCumul=1-(xmin./vX_tmp).^alpha;
    k(n)=max(abs(obsCumul-modCumul));
    
%     % fit b value
%     b(n) = get_b(iM,Mc);
%     
%     % generate catalog
%     iMsynth = geneq(b(n),Mc,Mmax,length(iM));
    
%     % perform ks test
%     try
%         [~,p(n),k(n)] = kstest2(iMsynth,iM);
%     catch
%         [k(n),p(n)] = deal(nan);
%     end
end

% plot output if so desired
if nargin == 2
    if strcmp(plotYN,'yes')
        figure
        plot(Mvec,k)

    end
end

% mc_candidates = find(islocalmin(k));
MCOUT = Mvec(k == min(k));
    
end
% 
% function b = get_b(M,Mc)
% 
% mags    = min(M):0.1:max(M);
% nMags   = length(mags)-1;
% gr      = zeros(nMags,1);
% 
% for n = 1:nMags
%     gr(n) = sum(M > mags(n));
% end
% d = polyfit(mags(1:end-1)',log10(gr),1);
% b = -d(1);
% %b = 1/log(10)*1./(mean(M)-Mc)
% end
% 
% function MAG = geneq(B,MMIN,MMAX,N)
% % inverse CDF of the normalized GR relationship from MMIN to MMAX
% 
% u   = rand(1,N);
% MAG = -1/B * log10(-u * (10^(-B*MMIN)-10^(-B*MMAX)) + 10^(-B*MMIN));
% 
% end