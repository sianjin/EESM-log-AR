function out = box0Simulation(simParams)
% box0Simulation Example helper function

%   Copyright 2019 The MathWorks, Inc.

% Extract configuration
cfgHE = simParams.Config;
substreamidx = simParams.RandomSubstream;
maxNumPackets = simParams.MaxNumPackets;
maxNumErrors = simParams.MaxNumErrors;
snr = simParams.SNR;
% Create an NDP packet with the correct number of space-time streams to
% generate enough LTF symbols
cfgNDP = wlanHESUConfig('APEPLength',0,'GuardInterval',0.8); % No data in an NDP
cfgNDP.ChannelBandwidth = cfgHE.ChannelBandwidth;
cfgNDP.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
cfgNDP.NumSpaceTimeStreams = cfgHE.NumTransmitAntennas;

% Set random substream index per iteration to ensure that each
% iteration uses a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',99);
stream.Substream = substreamidx;
RandStream.setGlobalStream(stream);

% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgHE);

% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE);

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
% snr here refers to the SNR on active subcarriers
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
N0 = 10^(-awgnChannel.SNR/10);

% Get path filers for last channel (same for all channels)
tgaxChannel = simParams.Channel;
tgaxChannelInfo = info(tgaxChannel);
pathFilters = tgaxChannelInfo.ChannelFilterCoefficients; % [NP numChannelTaps]
chInfo = getChanInfoParams(tgaxChannel); % Get Tx and Rx antenna correlation matrices

% Create and configure comm.MIMOChannel same as wlanTGaxChannel
% Doppler spectrum: Jakes; This is a must if FadingTechnique is set to SOS
mimoChan = comm.MIMOChannel;
mimoChan.FadingTechnique = 'Sum of sinusoids'; % Set FadingTechnique to SOS
mimoChan.SampleRate = wlanSampleRate(cfgHE);
mimoChan.AveragePathGains = tgaxChannelInfo.AveragePathGains;
mimoChan.PathDelays = tgaxChannelInfo.PathDelays;
mimoChan.SpatialCorrelationSpecification = 'Separate Tx Rx';
mimoChan.TransmitCorrelationMatrix  = permute(chInfo.TxCorrelationMatrix,[3 2 1]); 
mimoChan.ReceiveCorrelationMatrix  = permute(chInfo.RxCorrelationMatrix,[3 2 1]);
wavelength = 3e8/tgaxChannel.CarrierFrequency;
tgaxDopplerShift = tgaxChannel.EnvironmentalSpeed*(5/18)/wavelength; % Change km/h to m/s
mimoChan.MaximumDopplerShift = tgaxDopplerShift; % For Jakes model
mimoChan.PathGainsOutputPort = true; 

% Create object to deal with abstraction
Abstraction = tgaxLinkPerformanceModel;

% Loop to simulate multiple packets
perStore = nan(maxNumPackets,1);
perAbsStore = nan(maxNumPackets,1);
perAbsRawStore = nan(maxNumPackets,1);
snreffStore = nan(maxNumPackets,1);
sinrStore = nan(ofdmInfo.NumTones,cfgHE.NumSpaceTimeStreams,maxNumPackets); % Nsc-by-Nsts-by-maxNumPackets
numPacketErrors = 0;
numPacketErrorsAbs = 0;
numPkt = 1; % Index of packet transmitted

while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    % reset(tgaxChannel); % Reset channel for different realization
    reset(mimoChan); % Reset channel for different realization
    
    % -------------------------------------------------------------------
    % Comment this to run simulation without beamforming
    % Generate NDP packet - with an empty PSDU as no data
    txNDP = wlanWaveformGenerator([],cfgNDP);
    % For each user STA, pass the NDP packet through the channel and calculate
    % the feedback channel state matrix by SVD.
    % Received waveform at user STA with 50 sample padding. No noise.
    rxNDP = tgaxChannel([txNDP; zeros(50,size(txNDP,2))]);
    
    % Get the full-band beamforming feedback for a user
    staFeedback = heUserBeamformingFeedback(rxNDP,cfgNDP);
    % For each RU, calculate the steering matrix to apply
    % Calculate the steering matrix to apply to the RU given the feedback
    steeringMatrix = heSUCalculateSteeringMatrix(staFeedback,cfgHE,cfgNDP);
    
    % Apply the steering matrix to each RU
    cfgHE.SpatialMapping = 'Custom';
    cfgHE.SpatialMappingMatrix = steeringMatrix;
    % -------------------------------------------------------------------
    
    % Generate a packet with random PSDU
    psduLength = getPSDULength(cfgHE); % PSDU length in bytes
    txPSDU = randi([0 1],psduLength*8,1,'int8');
    tx = wlanWaveformGenerator(txPSDU,cfgHE);
    % Add trailing zeros to allow for channel delay
    txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];

    % Pass through a fading indoor TGax channel
    
    [rx,pathGains] = mimoChan(txPad);

    % Get perfect timing offset and channel matrix for HE-LTF field
    heltfPathGains = pathGains(ind.HELTF(1):ind.HELTF(2),:,:,:,:);
    pktOffset = channelDelay(heltfPathGains,pathFilters);
    chan = helperPerfectChannelEstimate(heltfPathGains,pathFilters,ofdmInfo.FFTLength,ofdmInfo.CPLength,ofdmInfo.ActiveFFTIndices,pktOffset);

    % Calculate SINR using abstraction
    % As multiple symbols returned average over symbols and permute
    % for calculations
    % Get precoding matrix for abstraction
    Wtx = getPrecodingMatrix(cfgHE); % Include cyclic shift applied per STS
    Wtx = Wtx/sqrt(cfgHE.NumSpaceTimeStreams);
    Htxrx = permute(mean(chan,2),[1 3 4 2]); % Nst-by-Nt-by-Nr
    Ptxrx = 1; % Assume transmit power is 0dBW
    sinr = calculateSINR(Htxrx,Ptxrx,Wtx,N0);
    sinrStore(:,:,numPkt) = sinr;

    % Link performance model - estimate PER using abstraction
    [perAbs,effSINR] = estimateLinkPerformance(Abstraction,sinr,cfgHE);

    % Flip a coin for the abstracted PHY
    packetErrorAbs = rand(1)<=perAbs;
    numPacketErrorsAbs = numPacketErrorsAbs+packetErrorAbs;

    % Store outputs for analysis
    perAbsRawStore(numPkt) = perAbs;
    perAbsStore(numPkt) = packetErrorAbs;
    snreffStore(numPkt) = effSINR;

    % Pass the waveform through AWGN channel
    rx = awgnChannel(rx);

    % Demodulate data symbols
    rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
    demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE);

    % Extract data subcarriers from demodulated symbols and channel
    % estimate
    demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);

    % Get channel estimate from channel matrix (include spatial mapping
    % and cyclic shift)
    chanEst = heChannelToChannelEstimate(chan,cfgHE);
    chanEstAv = permute(mean(chanEst,2),[1 3 4 2]); % Average over symbols
    chanEstData = chanEstAv(ofdmInfo.DataIndices,:,:);

    % Calculate single stream pilot estimates per symbol and noise
    % estimate
    chanEstSSPilots = permute(sum(chanEst(ofdmInfo.PilotIndices,:,:,:),3),[1 2 4 5 3]);
    demodPilotSym = demodSym(ofdmInfo.PilotIndices,:,:);
    nVarEst = heNoiseEstimate(demodPilotSym,chanEstSSPilots,cfgHE);

    % Equalization and STBC combining
    [eqDataSym,csi] = heEqualizeCombine(demodDataSym,chanEstData,nVarEst,cfgHE);
    rxPSDU = wlanHEDataBitRecover(eqDataSym,nVarEst,csi,cfgHE);

    % Determine if any bits are in error, i.e. a packet error
    packetError = ~isequal(txPSDU,rxPSDU);
    perStore(numPkt) = packetError;
    numPacketErrors = numPacketErrors+packetError;

    numPkt = numPkt+1;
end

% Remove last increment
numPkt = numPkt-1;

% Calculate packet error rate (PER) at SNR point
packetErrorRate = numPacketErrors/numPkt;
packetErrorRateAbs = numPacketErrorsAbs/numPkt;

% Return results
out = struct;
out.packetErrorRate = packetErrorRate;
out.perStore = perStore;
out.numPkt = numPkt;
out.sinrStore = sinrStore;
out.packetErrorRateAbs = packetErrorRateAbs;
out.perAbsRawStore = perAbsRawStore;
out.perAbsStore = perAbsStore;

disp([char(simParams.DelayProfile) ' '...
      num2str(simParams.NumTransmitAntennas) '-by-' ...
      num2str(simParams.NumReceiveAntennas) ','...
      ' MCS ' num2str(simParams.MCS) ','...
      ' SNR ' num2str(simParams.SNR) ...
      ' completed after ' num2str(out.numPkt) ' packets,'...
      ' PER:' num2str(out.packetErrorRate)]);
  
end

%% Get spatial correlation prameters
% Copyright 2017-2019 The MathWorks, Inc.
% Faisal Darbari <fdarbari@mathworks.com>

function chInfo = getChanInfoParams(tgaxChannel)
%getChanInfoParams Get TGax spatial parameters

   modelConfig = struct( ...
              'NumTransmitAntennas',tgaxChannel.NumTransmitAntennas, ...
              'NumReceiveAntennas',tgaxChannel.NumReceiveAntennas, ...
              'TransmitAntennaSpacing',tgaxChannel.TransmitAntennaSpacing, ...
              'ReceiveAntennaSpacing',tgaxChannel.ReceiveAntennaSpacing, ...
              'DelayProfile',tgaxChannel.DelayProfile, ...
              'UserIndex',tgaxChannel.UserIndex, ...
              'ChannelBandwidth',tgaxChannel.ChannelBandwidth, ...
              'TransmitReceiveDistance',tgaxChannel.TransmitReceiveDistance, ...
              'CarrierFrequency',tgaxChannel.CarrierFrequency, ...
              'TransmissionDirection',tgaxChannel.TransmissionDirection, ...
              'NumPenetratedFloors',tgaxChannel.NumPenetratedFloors, ...
              'NumPenetratedWalls',tgaxChannel.NumPenetratedWalls, ...
              'WallPenetrationLoss',tgaxChannel.WallPenetrationLoss, ...
              'FormatType',class(tgaxChannel), ...
              'InputDataType','double');
          
    chInfo = spatialCorrelation(modelConfig);

end
