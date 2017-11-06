% density plot;
% called inside MeasDist_all;

% set up density and plot


MM=[Rec.Mag];

MM=round(MM*100)/100;
iMag=find((MM>=3)&(MM<4));
DD=[Rec(iMag).DD];
DDfor=[Rec(iMag).Dfor];

LW=2;
% get aftershock density
[Dsort Isort]=sort(DD);
D_dx=Dsort(2:end)-Dsort(1:end-1);
D_mid=(Dsort(2:end)+Dsort(1:end-1))/2;
dens=1./D_dx;
dens=dens(dens<1e12);
D_mid=D_mid(dens<1e12);
ha=loglog(D_mid,dens,'x','LineWidth',LW);
hold on 
%ha=scatter(D_mid,dens,30,tsort(Isort(1:end-1)),'x');


%%
%foreshock density
[Dsortfor Isortfor]=sort(DDfor);
Dfor_dx=Dsortfor(2:end)-Dsortfor(1:end-1);
Dfor_mid=(Dsortfor(2:end)+Dsortfor(1:end-1))/2;

densfor=1./Dfor_dx;

% cut off outliers
densfor=densfor(densfor<1e12);
Dfor_mid=Dfor_mid(densfor<1e12);

hf=loglog(Dfor_mid,densfor,'ro','LineWidth',LW); 
%hf=scatter(Dfor_mid,densfor,30,tsort(Isortfor(1:end-1)),'o');
%hold off
set(gca,'FontSize',18);
legend([ha,hf],'Aftershocks','Foreshocks');
xlabel('Distance from Target Earthquake (km)');
ylabel('Linear Density (Earthquakes/km)');
%pause;

Iafter=(D_mid>0.5)&(D_mid<40)&(isfinite(dens));
Ifor=((Dfor_mid>0.5)&(Dfor_mid<40));

D_mid_fit=D_mid(Iafter);
dens_fit=dens(Iafter);
densfor_fit=densfor(Ifor);
Dfor_mid_fit=Dfor_mid(Ifor);

Ndensa=length(dens_fit);
Ndensf=length(densfor_fit);
for ibstrp=1:1000,
    pba=rand(1,Ndensa);
    Iba=round(pba.*Ndensa+0.5);
    pbf=rand(1,Ndensf);
    Ibf=round(pbf.*Ndensf+0.5);
    %Voutafter=lsqcurvefit(@(V,x) V(1)*x.^V(2),[400;-1.4],D_mid_fit(Iba),dens_fit(Iba));
    %Voutfore=lsqcurvefit(@(V,x) V(1)*x.^V(2),[25;-1.4],Dfor_mid_fit(Ibf),densfor_fit(Ibf));
   % Voutafter=lsqcurvefit(@(V,x) V(1)+x.*V(2),[400;-1.4],log10(D_mid_fit(Iba)),log10(dens_fit(Iba)));
   % Voutfore=lsqcurvefit(@(V,x) V(1)+x.*V(2),[25;-1.4],log10(Dfor_mid_fit(Ibf)),log10(densfor_fit(Ibf)));
    Voutafter=polyfit(log10(D_mid_fit(Iba)),log10(dens_fit(Iba)),1);
    Voutfore=polyfit(log10(Dfor_mid_fit(Ibf)),log10(densfor_fit(Ibf)),1);
    run(ibstrp).Voutafter2=Voutafter(1);
    run(ibstrp).Voutafter1=Voutafter(2);
    if (~isnan(Voutfore))
        run(ibstrp).Voutfore2=Voutfore(1);
        run(ibstrp).Voutfore1=Voutfore(2);  
    end
end

Voutafter(1)=10.^mean([run.Voutafter1]);
Voutafter(2)=mean([run.Voutafter2]);
Voutafter_err=std([run.Voutafter2]);

Voutfore(1)=10.^mean([run.Voutfore1]);
Voutfore(2)=mean([run.Voutfore2]);
Voutfore_err=std([run.Voutfore2]);

Dsynth=D_mid_fit;
hold on;
loglog(Dsynth,Voutfore(1).*Dsynth.^Voutfore(2),'r')
loglog(Dsynth,Voutafter(1).*Dsynth.^Voutafter(2),'b')
hold off;

['After :' num2str(Voutafter(2)) '+/-' num2str(Voutafter_err)]
['Fore :' num2str(Voutfore(2)) '+/-' num2str(Voutfore_err)]

%ylim([1e-3 1e7])
xlim([0.1 1e3])
sum(Iafter)/sum(Ifor)
[sum(Iafter) sum(Ifor)]
