function SM = fsp2slipmap(FSP)

SM = slipmaps;

SM.EID      = FSP.Tag.
SM.FileName = FSP.FileName;
SM.XGrid    = xGrid; % to do
SM.YGrid    = yGrid; % to do
SM.SlipGrid = slipMap; % to do
SM.PointSpacing = [xPtSp, yPtSp];
SM.SlipArray= [xF',yF',slip];
SM.El       = El;
SM.Kode     = Kode;
SM.Poisson 	= Poisson;
SM.Young    = Yougn * 10^5; % bar-> N/M^2
SM.Cdepth   = Cdepth;
SM.Fric     = Fric;
SM.RStress  = RStress;
SM.Mo       = Mo/10^7; % dyne cm -> Nm
SM.Mw       = Mw;
SM.EllipseError     = elpseError;
SM.EllipseErrorNorm = elpseErrorNorm;
SM.Dip      = mean(El(:,7));
SM.Rake     = rakeMap;
SM.StressDrop = [];
SM.Perimeter= [];


end