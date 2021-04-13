% Plot 3 curves:
% 1) Full PHY packet error state (0 or 1)
% 2) EESM predicted instantaneous PER
% 3) Random realization of EESM predicted packet error state (0 or 1)
function plotValidation(snrs,results,snrIdx,maxNumPackets,attri)
txPeriod = results{snrIdx}.txPeriod; % Packet Tx period in sec
t = txPeriod*[1:maxNumPackets];
if attri == 'PER'
    perStore = results{snrIdx}.perStore;
    perInsVec = results{snrIdx}.perInsAbsVec;
    effSnrVec = results{snrIdx}.effSnrVec;
    packetErrorAbs = rand(1,maxNumPackets)<=perInsVec;
    elementwiseError = (packetErrorAbs ~= perStore);
    relativeError = sum(elementwiseError)/maxNumPackets;
    plot(t, perStore,'*','DisplayName','Full PHY packet error state (0 or 1)')
    hold on
    plot(t, perInsVec,'o','DisplayName','EESM predicted instantaneous PER');
    hold on
    plot(t, packetErrorAbs,'.','DisplayName','Random realization of EESM predicted packet error state (0 or 1)');
    legend
    grid on
    xlabel('Time (sec)')
    ylabel('Packet error state or PER')
    title(['Coherence time: 0.978s, RX SNR = ' ,num2str(snrs(snrIdx)), 'dB']);
elseif attri == 'effSnr' 
    effSnrVec = results{snrIdx}.effSnrVec;
    plot(t, effSnrVec)
    grid on
    xlabel('Time (sec)')
    ylabel('Effective SNR (dB)')
    title(['Coherence time: 0.978s, RX SNR = ' ,num2str(snrs(snrIdx)), 'dB']);
end
end