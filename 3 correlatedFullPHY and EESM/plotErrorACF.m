load('snrPer_CBW20_Model-D_4-by-2_MCS4Tp2503rd.mat')
snrIdx = 2;
perInsAbsVec = results{snrIdx}.perInsAbsVec;
packetErrorAbs = rand(1,length(perInsAbsVec))<=perInsAbsVec;
packetErrorAbs = double(packetErrorAbs);
perStore = results{snrIdx}.perStore;
figure
subplot(2,1,1)
autocorr(perStore,'NumLags',40)
hold on
autocorr(packetErrorAbs,'NumLags',40)
subplot(2,1,2)
parcorr(perStore,'NumLags',40)
hold on
parcorr(packetErrorAbs,'NumLags',40)