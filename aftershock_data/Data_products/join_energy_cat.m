% script used to merge ISC and energy catalogs into one

load eq_energy.mat
load IRIS_DMC_with_FMS.mat
nC                  = CAT;
nEq                 = height(nC);
[nC.BBtotalEnergy, nC.RuptureDuration, nC.HFtotalEnergy] = deal(nan(nEq,1));
nMerged = 0;

for n = 1:height(nC)
    I = nC.EVENT_ID(n) == eqenergy.EVENT_ID;
    if sum(I) == 1
        nMerged = nMerged + 1;
        nC.HFtotalEnergy(n)     = eqenergy.HfTotalEnergy(I);
        nC.BBtotalEnergy(n)     = eqenergy.BbTotalEnergy(I);
        nC.RuptureDuration(n)   = eqenergy.BbRuptureDuration(I);
        
    end
end

disp([num2str(nMerged), ' out of ', num2str(height(eqenergy)), ' earthquakes were successfully matched']);

iris_dmc_cat_with_fms_and_energy = nC;
save('IRIS_DMC_with_FMS_and_energy.mat','iris_dmc_cat_with_fms_and_energy');
