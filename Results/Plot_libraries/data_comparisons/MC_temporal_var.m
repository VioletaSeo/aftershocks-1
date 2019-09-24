function MC_temporal_var(M,t)
Nyr = 4;
tVec = linspace(min(t),max(t),Nyr);
for n = 1:Nyr-1
    mc(n) = KS_completeness(M(t>tVec(n) & t<tVec(n+1)) ,'yes');
end
% figure
% plot(tVec(1:n),mc)

end