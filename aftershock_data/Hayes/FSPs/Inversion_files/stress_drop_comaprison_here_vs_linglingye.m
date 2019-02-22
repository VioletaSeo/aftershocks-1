load('ye12.mat');
ye12.t = datenum(   [char(ye12.Date), ...
                    repmat(' ',height(ye12),1), ...
                    char(ye12.Time)]);
load('FSPcat.mat');
FSPye = FSPcat;
numeq = height(FSPye);
FSPye.StressDropYe = nan(numeq,1);
for n = 1:numeq
    [dt,I] = min(abs(FSPye.Time(n) - datenum(ye12.t)));
    if dt < 2
        FSPye.StressDropYe(n) = ye12.StressDrop(I);
    end
end

figure
scatter(FSPye.StressDrop, FSPye.StressDropYe/10^6,FSPye.Mo/(2*10^19))
set(gca,'XScale','log','YScale','log','Xlim',[0.3,30],'Ylim',[0.3,30])
title(sprintf('n = %g',sum(~isnan( FSPye.StressDropYe))))

