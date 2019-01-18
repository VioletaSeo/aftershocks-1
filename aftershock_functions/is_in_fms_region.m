function [YN,varargout] = is_in_fms_region...
    (LAT,LON,DEPTH,MAG,FMS,distCutoff,depthCutoff,Mth, fmsQuerry)

I  = MAG > Mth & FMS == fmsQuerry;
loc= find(I);

wgs84 = wgs84Ellipsoid('kilometer');
[X,Y,Z] = geodetic2ecef(wgs84, LAT, LON, -DEPTH');

[Xf,Yf, Zf, DEPTHf] = goodind( I, ...
    X, Y, Z, DEPTH);

YN      = zeros(1,length(LAT));
Icumm   = zeros(length(Xf),1);

parfor n = 1:length(LAT)
    YN(n) = 0;
    pt = [X(n),Y(n),Z(n)];
    Iin  = distCutoff > sqrt(sum((pt-[Xf,Yf,Zf]).^2,2)) ...
        & ...
        depthCutoff > abs(DEPTH(n)-DEPTHf);
    if sum(Iin) > 1
        YN(n) = true;
        Icumm = Icumm | Iin;        
    end
end

YN = logical(YN);
varargout{1} = loc(Icumm);

end

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);
end

end