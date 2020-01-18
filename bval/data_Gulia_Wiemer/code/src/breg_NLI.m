% function breg_NLI will compute the regional pre-mainchosk b-value in the
% case that not enough are detected in the box. 
% Author: Laura Gulia, 2019
% 

function [Mcreg, breg, areg,bprevalue,cat_Reg]=breg_NLI(cat_Reg,mbin,Tms,Mms,corr)


b_orig=cat_Reg;
N=length(b_orig(:,1));
Tmax=Tms;
Tmin=min(b_orig(:,3));
T=Tmax-Tmin;

Nmin=50;
Mrange=min(b_orig(:,6)):mbin:max(b_orig(:,6));
cntnumb_orig=hist(b_orig(:,6),Mrange);
cntnumbA_orig=cntnumb_orig./T;
Mtarg=Mms;
bv=nan; av_ann=nan; av=nan; pr=nan; result_flag=nan;

b=b_orig;
cntnumb=hist(b(:,6),Mrange);
[maxn,ind]=max(cntnumb);
magco=Mrange(ind)+corr;


ll=b(:,6)>=magco;
b=b(ll,:);
N_b=sum(ll);


if N_b>=Nmin

% Call the function to compute value and uncertainties, checking for non linertity also     
[bestmc,bv,result_flag]=NLIndex_version1_sigma(b,magco,'PreDefinedMc');
    
   if or(result_flag==3,result_flag==2)
        
        magco=bestmc;
        ll=b(:,6)>=magco;
        N_b=sum(ll);
        
        av_ann=log10(N_b/T)+ bv*magco;
        av=log10(N_b)+ bv*magco;
        tr=1/10.^(av_ann-bv*Mtarg);
        pr=1-exp(-1/tr);
    else
        bv=nan;
        magco=bestmc;
        av_ann=nan;
        pr=nan;
    end
    
end



if isfinite(magco)
    Numbh=cntnumb_orig(end:-1:1);
    Ncumh=cumsum(Numbh);
    Ncum=Ncumh(end:-1:1);
    N1=10^(av_ann-bv*magco);
    N2=10^(av_ann-bv*max(Mrange));
    
end
breg=bv;
areg=av_ann;
Mcreg=magco;
prreg=pr;


bprevalue=isfinite(breg);

ll=b_orig(:,6)>=Mcreg;
cat_Reg=b_orig(ll,:);

