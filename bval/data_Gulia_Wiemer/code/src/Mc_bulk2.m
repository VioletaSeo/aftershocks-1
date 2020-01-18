
function [Mc,cntnumb]=Mc_bulk2(catS, Mrange,corr)


cntnumb=hist(catS(:,6),Mrange);
[~,ind]=max(cntnumb);
Mc=Mrange(ind)+corr;


    




