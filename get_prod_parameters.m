function [prefactor,alpha] = get_prod_parameters(t, lat,lon,depth,M,mth)
x=lat;
y=lon;
z=depth;

mag=M;tsort=t;



% remove earthquakes below detection
Imag=find(mag>mth);
mag=mag(Imag);
x=x(Imag);
y=y(Imag);
z=z(Imag);
tsort=tsort(Imag);

% sort earthquakes with respect to size
%[tsort,I]=sort(tsort);
%mag=mag(I);
[mag,I]=sort(mag,'descend');
tsort=tsort(I);
x=x(I);
y=y(I);
z=z(I);

Dmaxbig=100; % max search radius
Dtbig=3; % or 0.5 day before % this is the time after the show;

Dt=(30/60)/24;
Dtfor=-Dt;

minmag=mth; %just search for all magnitudes

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
maxmag=max(mag);

for mm=maxmag:-0.2:minmag
    
    ijk             = ijk+1;
    Rec(ijk).Mag    = mm;
    NumMain         = 0;
    I1              = find(abs(mag-mm)<0.05); %0.05?
    Nparfor         = length(I1);
    
    
    
    if (Nparfor>0)
        clear Recm;
        Recm(Nparfor).DD    =[];
        Recm(Nparfor).Dfor  =[];
        Recm(Nparfor).C     =[];
        
        for i =1:Nparfor
            
            t0    = tsort(I1(i));
            lat0  = x(I1(i));
            lon0  = y(I1(i));
            z0    = z(I1(i));
            
            trel  = tsort-t0;
            
            D     =distance(x,y,lat0,lon0);
            D     =deg2km(D);
            
            Itbig =((trel>-0.5)&(trel<Dtbig));
            Idbig =(abs(D)<Dmaxbig);
            Ibig  =and(Itbig,Idbig);
            
    
    
            if (flag(I1(i))==0)
                I             =(trel>0)&(trel<Dt)&(abs(D)<10)';
                Ifor          =(trel<0)&(trel>Dtfor)&(abs(D)<10)';
                I(I1(i))      =0; % guard against including mainshock
                Ifor(I1(i))   =0;
                %    if ((sum([flag(I) flag(Ifor)]))<1)   % none of the aftershocks or foreshocks are flagged
                NumMain=NumMain+1;
                Recm(i).DD=D(I);
                Recm(i).Dfor=D(Ifor);
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

figure
plot([Rec.Mag],[Rec.C], 'o')
set(gca,'YScale','log')
title('Aftershock Productivity')
xlabel('Magnitude')
ylabel('Average Number of Aftershocks')


% get the productivity: 
% to do this I model the aftershock production (plotted aboves) as two
% logarythmic segments: aftershock production and background productivity.
% The fit is accomplished in peicewise fashion is a free knot.

[prefactor,alpha] = getproductivity([Rec.Mag], [Rec.C]);
disp('alpha:')
disp(alpha)
disp('prefactor:')
disp(prefactor)

hold on
minmaxMag = minmax([Rec.Mag]);
plot(minmaxMag,prefactor*10.^(alpha*minmaxMag),'--');

end
