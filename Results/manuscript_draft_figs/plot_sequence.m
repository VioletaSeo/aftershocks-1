function plot_sequence(t0,CAT,varargin)

I = CAT.time > t0 & CAT.time < t0 +40;
IFS = CAT.time > t0-40 & CAT.time < t0;
IMS = CAT.time == t0;

figure
worldmap_base;
scatterm(CAT.lat(I),  CAT.lon(I),  (CAT.M(I)-2).^3,   CAT.time(I)-t0)
scatterm(CAT.lat(IFS),CAT.lon(IFS),(CAT.M(IFS)-2).^3, 'r')
scatterm(CAT.lat(IMS),CAT.lon(IMS),(CAT.M(IMS)-2).^4, 'r','Filled')

title("Foreshocks and aftershocks 10 before the mainshock to 60 after")

end