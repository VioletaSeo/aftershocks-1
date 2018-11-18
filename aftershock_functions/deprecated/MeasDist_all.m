% calculate earthquake density
clear
clf;

A=load('SCSN_1984-2008.cleaned');
t=datenum(A(:,1:6))';
lat=A(:,8)';
lon=A(:,9)';
depth=A(:,10)';
M=A(:,7)';

%%%%%%%generic code below here
x=lat;
y=lon;
z=depth;

I=find((y<-114)&(x<40)&(x>32.5)&(z<20)&(z>0));
t=t(I);
M=M(I);
x=x(I);
y=y(I);
z=z(I);

%mag=round(M*10)/10;
mag=M;tsort=t;

mth=2; % complete to M=2.1 (this is accurate for ANSS)


Imag=find(mag>mth);
mag=mag(Imag);
x=x(Imag);
y=y(Imag);
z=z(Imag);
tsort=tsort(Imag);

%[tsort,I]=sort(tsort);
%mag=mag(I);
[mag,I]=sort(mag,'descend');
tsort=tsort(I);
x=x(I);
y=y(I);
z=z(I);

Dmaxbig=100;
Dtbig=3; % or 0.5 day before

Dt=(15/60)/24;
Dtfor=-Dt;

minmag=3;

flag=tsort*0;
ijk=0;

% define structure
Rec(1).DD=[];
Rec(1).Mag=[];
Rec(1).Dfor=[];
Rec(1).C=[];
Rec(1).NumMain=[];
Nmag=round((max(mag)-minmag)*10);
Rec(Nmag)=Rec(1);
maxmag=7.3;

for mm=maxmag:-0.1:minmag
    mm
    ijk=ijk+1; 
    Rec(ijk).Mag=mm;
    NumMain=0;
    I1=find(abs(mag-mm)<0.05);
    Nparfor=length(I1);   
    if (Nparfor>0),
        clear Recm;
        Recm(Nparfor).DD=[];
        Recm(Nparfor).Dfor=[];
        Recm(Nparfor).C=[];
        for i=1:Nparfor,
          t0=tsort(I1(i));
          lat0=x(I1(i));
          lon0=y(I1(i));
          z0=z(I1(i));
          
          trel=tsort-t0;
          % kf code
%          y = (lat2 - lat1) * 111.2;
 %         meanlat = mean([lat2 lat1],2);
 %         x = ((lon2 - lon1).*111.2)./cos((meanlat./360).*(2*pi));
 %         dist = x.^2 + y.^2;
  %       dist = sqrt(dist);

          %dx=(x-lat0)*111.2;
          %meanlat=(x+lat0)./2;
          %dy=((y-lon0).*111.2).*cos((meanlat./180).*(pi));
          dz=z-z0;
 %         D=sqrt(dx.^2+dy.^2+dz.^2);
             D=distance(x,y,lat0,lon0);
%          D=D*111.1;
          
          
          D=deg2km(D);
          %D=sqrt(D.^2+(z-z0).^2);
          
          Itbig=((trel>-0.5)&(trel<Dtbig));
          Idbig=(abs(D)<Dmaxbig);
          Ibig=and(Itbig,Idbig);
        
          if (flag(I1(i))==0)
              I=((trel>0)&(trel<Dt));%&(abs(D)<100));
              Ifor=((trel<0)&(trel>Dtfor));%&(abs(D)<100));
              I(I1(i))=0; % guard against including mainshock
              Ifor(I1(i))=0;
          %    if ((sum([flag(I) flag(Ifor)]))<1)   % none of the aftershocks or foreshocks are flagged                                            
                 NumMain=NumMain+1;
                   Recm(i).DD=D(I);   
                   Recm(i).Dfor=D(Ifor);
          %	end 
	   end
    
    flag(Ibig)=1;
    end   
    % Stop when everything is blocked out
   % if (sum(flag)>(length(flag)-5)), break, end
        
    Rec(ijk).DD=[Recm.DD];
    Rec(ijk).Dfor=[Recm.Dfor];
    end    
    Rec(ijk).C=length(Rec(ijk).DD)/Nparfor;
    Rec(ijk).NumMain=NumMain;
    I1old=I1;
end

ijklast=ijk;
Rec=Rec(1:ijklast);
PlotDens;
save sorted;

