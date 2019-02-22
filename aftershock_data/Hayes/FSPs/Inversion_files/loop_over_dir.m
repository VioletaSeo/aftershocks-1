function fspInvObjArray = loop_over_dir(dirQ)

dirInfo = dir(dirQ);
nF  = length(dir);
eqs = repmat(fsp_inversion,nF,1);
loc = zeros(1,nF);

for n = 1:nF
    fn = dirInfo(n).name;
                              % remove duplicate inversions
    if contains(fn,'.fsp') && ~contains(fn,'_2.fsp')
        eqs(n) = fsp_inversion(fn);
        loc(n) = true;
    end
end

% remove blank entries
loc = logical(loc);
eqs = eqs(loc);

% try removing multisegment inversions
MS  = zeros(length(eqs), 1);
for n = 1:length(eqs)
    MS(n) = eqs(n).Parameters.Nsg;
end

eqSingleRupt = eqs(MS == 1);
fspInvObjArray = eqSingleRupt;

save('fsp.mat','fspInvObjArray')

end