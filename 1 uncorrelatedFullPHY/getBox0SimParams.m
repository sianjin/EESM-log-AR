function simParams = getBox0SimParams(chans,numTxRx,mcs,cfgHE,maxNumErrors,maxNumPackets)
% getBox0SimParams Example helper function

%   Copyright 2019 The MathWorks, Inc.

% These arrays define the value and order SNRs are defined
channelConfigs = ["Model-B","Model-D"];
anteannaSNRConfigs = [1 1; 2 1; 2 2; 4 2; 8 2];

snr = {
    % Model-B
    [ ...
    {... % 1x1
    [-3:4:9,11], ...   % MCS 0
    [1:4:13], ...     % MCS 1
    [2:4:18], ...     % MCS 2
    [7:4:19,21], ...     % MCS 3
    [9:4:25], ...     % MCS 4
    [13:4:29], ...    % MCS 5
    [14:4:30], ...    % MCS 6
    [16:4:32,34], ...    % MCS 7
    [18:4:34,36] ... % MCS 8
    [18:4:38] ...     % MCS 9
    }; ...
        {... % 2x1
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    }; ...
    	{... % 2x2
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    }; ...
        {... % 4x2
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    }; ...
        {... % 8x2
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    16:4:40, ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    }; ...
    ];

% Model-D
    [ ...
    {... % 1x1
    [-3 -2 -1 10], ...     % MCS 0
    0:4:12, ...     % MCS 1
    2:4:18, ...     % MCS 2
    8:4:20, ...     % MCS 3
    [11,15,19,21], ...    % MCS 4
    [14:4:26,28], ...    % MCS 5
    16:4:28, ...    % MCS 6
    18:4:30, ...    % MCS 7
    21.5:4:37.5 ... % MCS 8
    20:4:36 ...     % MCS 9
    }; ...
    	{... % 2x1
    [-3:1:-1,-0.5], ...  %[-3:1:-1,-0.5,0,0.5,0.75,1], ...   % MCS 0
    [0,1,1.5,2], ... %[0,1,1.5,2,2.5,2.75,3,3.25], ...   % MCS 1
    [4:1:6,6.5], ... %[4:1:6,6.5,7,7.5,8,8.25], ...   % MCS 2
    [5:1:7,7.5], ... %[5:1:7,7.5,8,8.5,9,9.25], ...  % MCS 3
    [12,15,18.5,20], ... %[11,11.5,12,12.5,13,13.5,13.75,14], ...  % MCS 4
    [12:1:15], ... %[12:1:15,15.5,16,16.5,16.75], ...  % MCS 5
    [15:1:17,17.5], ... %[15:1:17,17.5,18,18.5,19,19.25], ...  % MCS 6
    [16:1:18,18.5], ... %[16:1:18,18.5,19,19.5,20,20.25], ...  % MCS 7
    [18.5:1:21.5], ... %[18.5:1:21.5,22,22.5,23,23.5] ...   % MCS 8
    [22,23,23.5,24], ... %[22,23,23.5,24,24.5,25,25.25,25.5] ...   % MCS 9
    }; ...
    	{... % 2x2
    [-3:1:-1,-0.5], ...  %[-3:1:-1,-0.5,0,0.5,0.75,1], ...   % MCS 0
    [0,1,1.5,2], ... %[0,1,1.5,2,2.5,2.75,3,3.25], ...   % MCS 1
    [4:1:6,6.5], ... %[4:1:6,6.5,7,7.5,8,8.25], ...   % MCS 2
    [5:1:7,7.5], ... %[5:1:7,7.5,8,8.5,9,9.25], ...  % MCS 3
    [12,12.5,13,13.5], ... %[11,11.5,12,12.5,13,13.5,13.75,14], ...  % MCS 4
    [12:1:15], ... %[12:1:15,15.5,16,16.5,16.75], ...  % MCS 5
    [15:1:17,17.5], ... %[15:1:17,17.5,18,18.5,19,19.25], ...  % MCS 6
    [16:1:18,18.5], ... %[16:1:18,18.5,19,19.5,20,20.25], ...  % MCS 7
    [18.5:1:21.5], ... %[18.5:1:21.5,22,22.5,23,23.5] ...   % MCS 8
    [22,23,23.5,24], ... %[22,23,23.5,24,24.5,25,25.25,25.5] ...   % MCS 9
    }; ...
        {... % 2x2
    [-3:1:-1,-0.5], ...  %[-3:1:-1,-0.5,0,0.5,0.75,1], ...   % MCS 0
    [0,1,1.5,2], ... %[0,1,1.5,2,2.5,2.75,3,3.25], ...   % MCS 1
    [4:1:6,6.5], ... %[4:1:6,6.5,7,7.5,8,8.25], ...   % MCS 2
    [5:1:7,7.5], ... %[5:1:7,7.5,8,8.5,9,9.25], ...  % MCS 3
    [11,11.5,12,12.5], ... %[11,11.5,12,12.5,13,13.5,13.75,14], ...  % MCS 4
    [12:1:15], ... %[12:1:15,15.5,16,16.5,16.75], ...  % MCS 5
    [15:1:17,17.5], ... %[15:1:17,17.5,18,18.5,19,19.25], ...  % MCS 6
    [16:1:18,18.5], ... %[16:1:18,18.5,19,19.5,20,20.25], ...  % MCS 7
    [18.5:1:21.5], ... %[18.5:1:21.5,22,22.5,23,23.5] ...   % MCS 8
    [22,23,23.5,24], ... %[22,23,23.5,24,24.5,25,25.25,25.5] ...   % MCS 9
    }; ...
    	{... % 8x2
    1:3:17, ...  % MCS 0
    5:3:23, ...  % MCS 1
    9:4:30, ...  % MCS 2
    12:4:36, ... % MCS 3
    [11:1:14], ... % MCS 4
    20:4:43, ... % MCS 5
    22:4:46, ... % MCS 6
    24:4:48, ... % MCS 7
    26:4:50 ...  % MCS 8
    29:4:54 ...  % MCS 9
    }; ...
    ] ...
    };

% Create channel configuration
tgaxChannel = wlanTGaxChannel;
tgaxChannel.DelayProfile = 'Model-D';
tgaxChannel.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
tgaxChannel.NumReceiveAntennas = 1;
tgaxChannel.TransmitReceiveDistance = 15; % Distance in meters for NLOS
tgaxChannel.ChannelBandwidth = cfgHE.ChannelBandwidth;
tgaxChannel.LargeScaleFadingEffect = 'None';
fs = wlanSampleRate(cfgHE);
tgaxChannel.SampleRate = fs;
tgaxChannel.PathGainsOutputPort = true;
tgaxChannel.NormalizeChannelOutputs = false;

% Generate a structure array containing the simulation parameters,
% simParams. Each element contains the parameters for a simulation.
simParamsRef = struct('MCS',0,'SNR',0,'RandomSubstream',0,'Config',cfgHE, ...
    'MaxNumPackets',maxNumPackets,'MaxNumErrors',maxNumErrors, ...
    'NumTransmitAntennas',0,'NumReceiveAntennas',0,'DelayProfile',"Model-B",...
    'Channel',tgaxChannel);
simParams = repmat(simParamsRef,0,0);
% There must be a SNR cell for each channel
assert(all(numel(channelConfigs)==numel(snr)))
% There must be a SNR cell element for each MIMO configuration
assert(all(size(anteannaSNRConfigs,1)==cellfun(@(x)size(x,1),snr)))
for ichan = 1:numel(chans)
    channelIdx = chans(ichan)==channelConfigs;
    for itxrx = 1:size(numTxRx,1)
        numTxRxIdx = all(numTxRx(itxrx,:)==anteannaSNRConfigs,2);
        for imcs = 1:numel(mcs)
            snrIdx = mcs(imcs)+1;
            for isnr = 1:numel([snr{channelIdx}{numTxRxIdx,snrIdx}])
                % Set simulation specific parameters
                sp = simParamsRef;
                sp.MCS = mcs(imcs);
                sp.NumTransmitAntennas = numTxRx(itxrx,1);
                sp.NumReceiveAntennas = numTxRx(itxrx,2);
                sp.DelayProfile = chans(ichan);

                % Set random substream for reproducible results
                sp.RandomSubstream = isnr;
                
                % Setup PHY configuration
                sp.Config.MCS = mcs(imcs);
                sp.Config.NumTransmitAntennas = numTxRx(itxrx,1);
                sp.Config.NumSpaceTimeStreams = 2;
                sp.Config.SpatialMapping = 'Fourier';

                % Configure channel model
                sp.Channel = clone(tgaxChannel);
                sp.Channel.DelayProfile = chans(ichan);
                sp.Channel.NumTransmitAntennas = numTxRx(itxrx,1);
                sp.Channel.NumReceiveAntennas = numTxRx(itxrx,2);

                % Lookup SNR to simulate
                sp.SNR = snr{channelIdx}{numTxRxIdx,snrIdx}(isnr);

                % Append to other tests
                simParams = [simParams sp]; %#ok<AGROW>
            end
        end
    end
end

end