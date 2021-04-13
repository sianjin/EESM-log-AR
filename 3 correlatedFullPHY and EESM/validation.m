%% mcs 4
clear all
tStart = tic;% Simulation Parameters
mcs = [4]; % Vector of MCS to simulate between 0 and 9
numTxRx = [4 2]; % Matrix of MIMO schemes, each row is [numTx numRx]
chan = "Model-D"; % String array of delay profiles to simulate
maxnumberrors = 200*1e3;  % The maximum number of packet errors at an SNR point 
maxNumPackets = 200*1e3; % The maximum number of packets at an SNR point
% maxnumberrors = 1e2;  % The maximum number of packet errors at an SNR point
% maxNumPackets = 1e2; % The maximum number of packets at an SNR point
beta = 7.2502; % EESM tuning parameter

% Fixed PHY configuration for all simulations
cfgHE = wlanHESUConfig;
cfgHE.ChannelBandwidth = 'CBW20'; % Channel bandwidth
bandwidth = cfgHE.ChannelBandwidth;
cfgHE.APEPLength = 1000;          % Payload length in bytes
cfgHE.ChannelCoding = 'LDPC';     % Channel coding

% Generate a structure array of simulation configurations. Each element is
% one SNR point to simulate.
simParams = getBox0SimParams(chan,numTxRx,mcs,cfgHE,maxnumberrors,maxNumPackets,beta);
snrs = [simParams.SNR];

% Simulate each configuration
results = cell(1,numel(simParams));
parfor isim = 1:numel(simParams)  % Use 'parfor' to speed up the simulation
    results{isim} = box0Simulation(simParams(isim));
end
tEnd = toc(tStart);
fname_I = sprintf('snrPer_%s_%s_%s-by-%s_MCS%s.mat',bandwidth,char(chan),num2str(numTxRx(1)),num2str(numTxRx(2)),num2str(mcs));
save(fname_I,'results','mcs','numTxRx','chan','cfgHE','maxNumPackets','snrs','tEnd')

snrIdx = 1;
plotValidation(snrs,results,snrIdx,maxNumPackets,"effSnr");
plotPERvsSNR(simParams,results)
