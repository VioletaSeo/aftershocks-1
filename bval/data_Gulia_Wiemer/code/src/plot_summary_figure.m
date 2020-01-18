% this script will plots a summary figure of all resulst of the time series
% analysis of pre and post-mainshco analysis. 
% Author: Laura Gulia, 2019 
figure(summary)


subplot(2,5,[1,2,6,7])
Mrange=min(catS(:,6)):mbin:max(catS(:,6));

cntnumb=hist(catSpre_orig(:,6),Mrange);
Numbh=cntnumb(end:-1:1);
Ncumh=cumsum(Numbh);
Ncumpre=Ncumh(end:-1:1);
Tcatpre=catSpre_orig(end,3)-Tmin;
semilogy(Mrange,Ncumpre./Tcatpre,'ob')
hold on
ll=catSpre_orig(:,6)>=Mcpre_orig;
catSpre_orig=catSpre_orig(ll,:);
if length(catSpre_orig(:,6))>=50
    bv=(1/(mean(catSpre_orig(:,6))-(Mcpre_orig-mbin/2)))*log10(exp(1));
    av_ann=log10(length(catSpre_orig(:,1))/Tcatpre)+ bv*Mcpre_orig;
    N1=10^(av_ann-bv*Mcpre_orig);
    N2=10^(av_ann-bv*max(Mrange));
    semilogy([Mcpre_orig,max(Mrange)],[N1,N2],'b');


    
    %Shi and Bold, 1982
    sig1 = (sum((catSpre_orig(:,6)-mean(catSpre_orig(:,6))).^2))/(length(catSpre_orig(:,6))*(length(catSpre_orig(:,6))-1));
        sig1 = sqrt(sig1);
        sig1_pre = 2.30*sig1*bv^2
end

title('Annualized FMDs - BK-AFT')
xlim([1, Mms+0.2])
ylabel('Cumulative Number')
xlabel('M')
bvplot=round(bv,2);
text(Mcpre_orig+0.3,N2+0.1*N2,[num2str(floor(catSpre_orig(1,3))),'/',num2str(catSpre_orig(1,4)),...
        '/',num2str(catSpre_orig(1,5)),' - ', num2str(floor(catSpre_orig(end,3))), '/',...
        num2str(catSpre_orig(end,4)), '/',num2str(catSpre_orig(end,5)),'  b=',num2str(bvplot,2)],'color','blue','FontSize',14)
   


Mcpost_orig2=Mcpost_orig+corr_post;
ll=catS_post(:,6)>=Mcpost_orig2;
catS_post=catS_post(ll,:);

cntnumb=hist(catS_post(:,6),Mrange);
Numbh=cntnumb(end:-1:1);
Ncumh=cumsum(Numbh);
Ncumpost=Ncumh(end:-1:1);
Tcatpost=catS_post(end,3)-catS_post(1,3);
semilogy(Mrange,Ncumpost./Tcatpost,'o','color',[0.5 0.5 0.5])
hold on
if length(catS_post(:,6))>=50
    bv=(1/(mean(catS_post(:,6))-(Mcpost_orig2-mbin/2)))*log10(exp(1));
    av_ann=log10(length(catS_post(:,1))/Tcatpost)+ bv*Mcpost_orig2;
    N3=10^(av_ann-bv*Mcpost_orig2);
    N4=10^(av_ann-bv*max(Mrange));
    semilogy([Mcpost_orig2,max(Mrange)],[N3,N4],'k');
    
        %Shi and Bold, 1982
    sig1 = (sum((catS_post(:,6)-mean(catS_post(:,6))).^2))/(length(catS_post(:,6))*(length(catS_post(:,6))-1));
        sig1 = sqrt(sig1);
        sig1_post = 2.30*sig1*bv^2

end

bvplot=round(bv,2);
text(Mcpost_orig+0.3,N4+0.1*N4,[num2str(floor(catS_post(1,3))),'/',num2str(catS_post(1,4)),...
        '/',num2str(catS_post(1,5)),' - ', num2str(floor(catS_post(end,3))), '/',...
        num2str(catS_post(end,4)), '/',num2str(catS_post(end,5)),'  b=',num2str(bvplot,2)],'color',[0.5 0.5 0.5],'FontSize',14)
   

TLS_post=(bvplot*100)/breference;

    
    
 ll=catS_post2(:,3)<=tpostmax2;
    catmax=catS_post2(ll,:);
    
    if length(catmax(:,1))>Npost
        catmax=catmax(end:-1:end-Npost+1,:);
        catmax=sortrows(catmax,3);
    end
    size(catmax)
    
    Tdiff3=tpostmax2-min(catmax(:,3));
    cntnumb=hist(catmax(:,6),Mrange);
    
    Numbh=cntnumb(end:-1:1);
    Ncumh=cumsum(Numbh);
    Ncum3=Ncumh(end:-1:1);
    
    N5=10^(apostmax2-bpostmax2*magcopostmax2);
    N6=10^(apostmax2-bpostmax2*max(Mrange));
    semilogy([magcopostmax2,max(Mrange)],[N5,N6],'m');
    semilogy(Mrange,Ncum3/Tdiff3,'om')
    bpostmaxplot2=round(bpostmax2,2);
    text(Mcpost_orig+0.3,N5+0.1*N6,[num2str(floor(catS_post(1,3))),'/',num2str(catS_post(1,4)),...
        '/',num2str(catS_post(1,5)),' - ', num2str(floor(catS_post(end,3))), '/',...
        num2str(catS_post(end,4)), '/',num2str(catS_post(end,5)),'  b=',num2str(bpostmaxplot2,2)],'color','magenta','FontSize',14)

    
            %Shi and Bold, 1982
    sig1 = (sum((catmax(:,6)-mean(catmax(:,6))).^2))/(length(catmax(:,6))*(length(catmax(:,6))-1));
        sig1 = sqrt(sig1);
        sig1_post2 = 2.30*sig1*bpostmaxplot2^2

    
    TLS_post1=(bpostmaxplot2*100)/breference;

    
    
subplot(2,5,3:5)
plot([Tms Tms], [0.1 2.4],'r-') % mainshock occurrence
hold on
plot([Tmin max(catS_post2(:,3))], [breference breference],'b-')
ylabel('b-value')
xlim([2016, Tmax])
ylim([0.4, (2)])
plot(tpost,bpost,'.','color',[0.5 0.5 0.5])
plot(tpost2,bpost2,'m.')
plot([Tms2 Tms2], [0.1 2.4],'m-') % mainshock occurrence
title('b time-series')
text(2000,(max(bpost)+0.15),'b time-series')


% Define the TLS level
subplot(2,5,[9,10])
title('TLS level')
xlim([1, 2])
ylim([1, 2])
axis off;
tt1 = text(0.5,1.6,['TLS post 1st event = ',num2str(TLS_post,4), 4],'Fontsize',[16])
tt2 = text(0.5,1.4,['TLS post 2nd event = ',num2str(TLS_post1,4),4],'Fontsize',[16])


if TLS_post<90
  tt1.Color = 'red'
elseif TLS_post>110
  tt1.Color = 'green'
else
  tt1.Color = 'yellow'
end



if TLS_post1<90
  tt2.Color = 'red'
elseif TLS_post1>110
  tt2.Color = 'green'
else
  tt2.Color = 'yellow'
end

