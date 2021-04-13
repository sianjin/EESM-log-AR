function [y,cpe] = heCommonPhaseErrorTracking(x,chanEst,cfg,varargin)
%heCommonPhaseErrorTracking HE common pilot phase tracking
%
%   [Y,CPE] = heCommonPhaseErrorTracking(X,CHANEST,CFGHE) performs common
%   pilot error phase tracking of the single user HE format input X.
%
%   Y is a complex Nst-by-Nsym-by-Nr array containing the common phase
%   corrected OFDM symbols. Nst is the number of occupied subcarriers, Nsym
%   is the number of symbols, and Nr is the number of receive antennas. Y
%   is the common phase error corrected symbols.
%
%   X is a complex Nst-by-Nsym-by-Nr array containing the received OFDM
%   symbols.
%
%   CHANEST is a complex Nst-by-Nsts-by-Nr array containing the channel
%   gains for all active subcarriers, or a Nsp-by-Nsts-by-Nr array
%   containing the channel gains for only pilot subcarriers. Nsts is the
%   number of space-time streams.
%
%   CFGHE is the format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>,
%   <a href="matlab:help('wlanHETBConfig')">wlanHETBConfig</a> or <a href="matlab:help('wlanHERecoveryConfig')">wlanHERecoveryConfig</a>.
%
%   [Y,CPE] = heCommonPhaseErrorTracking(X,CHANEST,CFGMU,RUNUMBER) performs
%   common pilot error phase tracking of the multi-user HE format input X.
%
%   CFGMU is the format configuration object of type <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a>.
%
%   RUNUMBER is the RU (resource unit) number.

%   Copyright 2017-2019 The MathWorks, Inc.

%#codegen

validateattributes(cfg,{'wlanHESUConfig','wlanHEMUConfig','wlanHETBConfig','wlanHERecoveryConfig'},{'scalar'},mfilename,'format configuration object');

if isa(cfg,'wlanHERecoveryConfig')
    pktFormat = cfg.PacketFormat;
    if strcmp(pktFormat,'HE-MU')
        numSpaceTimeStreamsPerRU = cfg.RUTotalSpaceTimeStreams;
        s = getSIGBLength(cfg);
        numHESIGB = s.NumSIGBSymbols;
    else % SU or EXT_SU
        numSpaceTimeStreamsPerRU = cfg.NumSpaceTimeStreams;
        numHESIGB = 0;
    end
    ruIdx = cfg.RUIndex;
    ruSize = cfg.RUSize;
else % wlanHESUConfig,wlanHEMUConfig,wlanHETBConfig
    pktFormat = packetFormat(cfg);
    allocInfo = ruInfo(cfg);
    if strcmp(pktFormat,'HE-MU')
        narginchk(4,4)
        ruNumber = varargin{1};
        sigbInfo = wlan.internal.heSIGBCodingInfo(cfg);
        numHESIGB = sigbInfo.NumSymbols;
        numSpaceTimeStreamsPerRU = allocInfo.NumSpaceTimeStreamsPerRU(ruNumber);
        ruSize = allocInfo.RUSizes(ruNumber);
        ruIdx = allocInfo.RUIndices(ruNumber);
    else
        % SU, EXT SU, TB
        numHESIGB = 0;
        ruIdx = allocInfo.RUIndices;
        numSpaceTimeStreamsPerRU = allocInfo.NumSpaceTimeStreamsPerRU;
        ruSize = allocInfo.RUSizes;
    end
end

numOFDMSym = size(x,2);
n = (0:numOFDMSym-1);
if strcmp(pktFormat,'HE-EXT-SU')
    numHESIGA = 4;
else % SU or MU
    numHESIGA = 2;
end

z = 2+numHESIGA+numHESIGB; % Pilot symbol offset
refPilots = wlan.internal.hePilots(ruSize,numSpaceTimeStreamsPerRU,n,z);

% Estimate CPE and phase correct symbols
info = wlanHEOFDMInfo('HE-Data',cfg.ChannelBandwidth,cfg.GuardInterval,[ruSize ruIdx]);

if numel(info.PilotIndices)==size(chanEst,1)
    % Assume channel estimate is only for pilots
    chanEstPilots = chanEst;
else
    % Otherwise extract pilots from channel estimate
    chanEstPilots = chanEst(info.PilotIndices,:,:);
end
cpe = wlan.internal.commonPhaseErrorEstimate(x(info.PilotIndices,:,:),chanEstPilots,refPilots);
y = wlan.internal.commonPhaseErrorCorrect(x,cpe);

end