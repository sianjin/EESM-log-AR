function S = spatialCorrelation(simConfig)
%SPATIALCORRELATION between antenna pairs of TGn and TGac channel model
%   S = spatialCorrelation function returns the power delay profile and
%   spatial correlation properties of IEEE 802.11 delay profiles as
%   specified by the IEEE 802.11 Wireless LAN Working group [1,2,3]
%
%   %   References:
%   [1] Erceg, V. et al. "TGn Channel Models."  Doc. IEEE 802.11-03/940r4.
%   [2] Briet, G. et al. "TGac Channel Model Addendum." Doc. IEEE 802.11-09/0308r12.
%   [3] Porat. R, Sk. Yong, K. Doppler. "TGah Channel Model." IEEE 802.11-11/0968r4.

% Copyright 2015-2020 The MathWorks, Inc.

%#codegen

    % Setup input parameters
    modelCfg = channelModelParameters(simConfig);
    % Transmit and receive antenna correlation
    [txCorr,rxCorr] = antennaCorrelation(modelCfg,simConfig);

    S.AntennaCorrelation = coder.nullcopy(complex(zeros(size(txCorr,1)*size(txCorr,2)*size(rxCorr,2),size(txCorr,1)*size(txCorr,2)*size(rxCorr,2),simConfig.InputDataType)));

    for ii = 1:size(rxCorr,1)
        idx = (ii-1)*size(txCorr,2)*size(rxCorr,2)+1:1:ii*size(txCorr,2)*size(rxCorr,2);
        txCorrShaped = reshape(txCorr(ii,:,:),size(txCorr,2),size(txCorr,3));
        rxCorrShaped = reshape(rxCorr(ii,:,:),size(rxCorr,2),size(rxCorr,3));
        S.AntennaCorrelation(idx,idx) = chol(kron(txCorrShaped,rxCorrShaped))';
    end

    % Set outputs
    S.PathPower = modelCfg.PathPower;
    S.PathDelays = modelCfg.PathDelays;
    S.Pathloss = modelCfg.Pathloss;
    S.ShadowFading = modelCfg.ShadowFading;
    S.Kfactor = modelCfg.Kfactor;
    S.TxLOSRadians = modelCfg.TxLOSRadians;
    S.RxLOSRadians = modelCfg.RxLOSRadians;
    S.TxCorrelationMatrix = txCorr;
    S.RxCorrelationMatrix = rxCorr;
end

function [txCorrelation,rxCorrelation] = antennaCorrelation(modelCfg,simConfig)

    numTxAnt = simConfig.NumTransmitAntennas;
    numRxAnt = simConfig.NumReceiveAntennas;
    spacingTx = simConfig.TransmitAntennaSpacing;
    spacingRx = simConfig.ReceiveAntennaSpacing;
    powerPerAngle = 10.^(.1.*modelCfg.PowerPerAngle);

    txCorrelation = coder.nullcopy(complex(zeros(size(modelCfg.PathPower,2),numTxAnt,numTxAnt,simConfig.InputDataType)));
    rxCorrelation = coder.nullcopy(complex(zeros(size(modelCfg.PathPower,2),numRxAnt,numRxAnt,simConfig.InputDataType)));

    for ii = 1:size(modelCfg.PathPower,2)
        index = transpose(find(modelCfg.PowerPerAngle(:,ii) > -Inf));
        if (~isempty(index))
           txCorr = correlationGeneration(numTxAnt, ...
               linspace(0,(numTxAnt-1)*spacingTx,numTxAnt), ...
               size(index,2),powerPerAngle(index,ii), ...
               modelCfg.TxAoD(index,ii).',modelCfg.TxAS(index,ii).', ...
               180.*ones(1,size(index,2)));
           rxCorr = correlationGeneration(numRxAnt, ...
               linspace(0,(numRxAnt-1)*spacingRx, numRxAnt), ...
               size(index,2),powerPerAngle(index,ii), ...
               modelCfg.RxAoA(index,ii).',modelCfg.RxAS(index,ii).', ...
               180.*ones(1,size(index,2)));
            txCorrelation(ii,:,:) = txCorr;
            rxCorrelation(ii,:,:) = rxCorr;
        end
    end
end

function R = correlationGeneration(numAnt,antennaSpacing,numCluster,cluterAmplitude,powerPerAngle,azimuthSpread,angleOfDeparture)

   [Q,sigmaDegrees] = normalisationLaplacian(numCluster,cluterAmplitude,azimuthSpread,angleOfDeparture);
   if (numAnt>1)
        Rxx = besselj(0,2*pi*antennaSpacing)+Q*RxxLaplacian(antennaSpacing,powerPerAngle,sigmaDegrees,angleOfDeparture);
        Rxy = Q*RxyLaplacian(antennaSpacing, powerPerAngle,sigmaDegrees,angleOfDeparture);
        R = toeplitz(complex(Rxx,Rxy));
   else
        R = 1;
   end
end

function [Q,sigmaDegrees] = normalisationLaplacian(numclusters,cluterAmplitude,azimuthSpread,angleOfDeparture)

    sigmaDegrees = zeros(1,numclusters);
    for k=1:numclusters
        tmp = sigmaValues;
        pos = find(tmp(:,1)>=azimuthSpread(k));
        sigmaDegrees(k) = ((tmp(pos(1),1)-azimuthSpread(k))*tmp(pos(1)-1,2) + ...
            (azimuthSpread(k)-tmp(pos(1)-1,1))*tmp(pos(1),2))/(tmp(pos(1),1)-tmp(pos(1)-1,1));
    end

    sigmaRadians = sigmaDegrees.*(pi/180);
    angleOfDepartureDegreees = angleOfDeparture.*(pi/180);

    if numclusters==1
      Q = 1/(1-exp(((-1)*sqrt(2)*angleOfDepartureDegreees)/sigmaRadians));
    else
      A = zeros(numclusters);
      A(1:numclusters-1,1) = 1/(sigmaRadians(1)*cluterAmplitude(1));
      for k=2:numclusters
        A(k-1,k) = (-1)/(sigmaRadians(k)*cluterAmplitude(k));
      end
      A(numclusters,:) = ones(1,numclusters)-exp(((-1).*sqrt(2).*angleOfDepartureDegreees)./sigmaRadians);
      B = zeros(numclusters,1);
      B(numclusters,1) = 1;
      Q = (A\B).';
    end
end

function out = RxyLaplacian(antSpacing,phase,sigma,angleOfDeparture)

    D = 2*pi*antSpacing;
    phi_0_rad = deg2rad(phase).';
    sigma_rad = deg2rad(sigma).';
    aod_rad = deg2rad(angleOfDeparture).';
    sqrt2 = sqrt(2);

    mrxy = (0:100)';
    RxyVec = besselj(repmat(2*mrxy,1,size(D,2))+1,repmat(D,101,1));
    m = (0:100);
    out = ((sin((2*m+1).*phi_0_rad).*((sqrt2./sigma_rad)+(exp(-sqrt2.*aod_rad./sigma_rad) ...
        .*((2*m+1).*sin((2*m+1).*aod_rad)-sqrt2.*cos((2*m+1).*aod_rad)./sigma_rad)))) ...
        ./(sqrt2.*sigma_rad.*((sqrt2./sigma_rad).^2+(2*m+1).^2)))*4*RxyVec;
end

function out = RxxLaplacian(antSpacing,phase,sigma,angleOfDeparture)

    D = 2*pi*antSpacing;
    phi_0_rad = deg2rad(phase).';
    sigma_rad = deg2rad(sigma).';
    aod_rad = deg2rad(angleOfDeparture).';
    sqrt2 = sqrt(2);

    mrxx = (1:100)';
    RxxVec = besselj(repmat(2*mrxx,1,size(D,2)),repmat(D,100,1));
    m = (1:100);
    out = (cos(2*m.*phi_0_rad).*(sqrt2./sigma_rad+exp(-sqrt2*aod_rad./sigma_rad) ...
        .*(2*m.*sin(2*m.*aod_rad)-sqrt2*cos(2*m.*aod_rad)./sigma_rad)) ...
        ./(4*sqrt2.*sigma_rad.*m.^2+sqrt(8)./sigma_rad))*4*RxxVec;
end

function out = channelModelParameters(inp)

    modelType = inp.DelayProfile;
    distanceTxRx = inp.TransmitReceiveDistance;

    switch modelType
        case ('Model-A')
            % Flat fading with 0 ns RMS delay spread(one tap at 0 ns delay)
            % Average power in linear scale
            pathPower = 1;
            % Relative delay (ns)
            pathDelays = 0;
            % Power roll-off coefficients
            powerPerAngle = 0;
            % Tx
            txAoDDegrees = 45;
            txASDegrees  = 40;
            % Rx
            rxAoADegrees = 45;
            rxASDegrees  = 40;

            % Path loss break point
            breakPointDistance = 5;
            if (distanceTxRx < breakPointDistance)
                % LOS conditions
                kFactor = 0; % Rician K-factor
                shadowFading = 3; % Shadow fading standard deviation
            else
                % NLOS conditions
                kFactor = -100; % Rician K-factor
                shadowFading = 4; % Shadow fading standard deviation
            end
    case 'Model-B'
        % Typical residential environment, 15 ns RMS delay spread
        % Average power in linear scale
        pathPower = 10.^(.1*[0 -5.4287 -2.5162 -5.8905 -9.1603 -12.5105 -15.6126 -18.7147 -21.8168]);
        % Relative delay (ns)
        pathDelays = [0 10e-9 20e-9 30e-9 40e-9 50e-9 60e-9 70e-9 80e-9 ];
        % Power roll-off coefficients
        powerPerAngle = [0 -5.4287 -10.8574 -16.2860 -21.7147 -Inf -Inf -Inf -Inf;
                        -Inf -Inf -3.2042 -6.3063 -9.4084 -12.5105 -15.6126 -18.7147 -21.8168];
        % Tx
        txAoDDegrees = [225.1084.*ones(1,5) -Inf.*ones(1,4);
                       -Inf.*ones(1,2) 106.5545.*ones(1,7)];
        txASDegrees  = [14.4490.*ones(1,5) -Inf.*ones(1,4);
                       -Inf.*ones(1,2) 25.4311.*ones(1,7)];
        % Rx
        rxAoADegrees = [4.3943.*ones(1,5) -Inf.*ones(1,4);
                       -Inf.*ones(1,2) 118.4327.*ones(1,7)];
        rxASDegrees  = [14.4699.*ones(1,5) -Inf.*ones(1,4);
                       -Inf.*ones(1,2) 25.2566.*ones(1,7)];

        % Path loss break point
        breakPointDistance = 5;
        if (distanceTxRx < breakPointDistance)
            % LOS conditions
            kFactor = [0,(-100).*ones(1,8)]; % Rician K-factor
            shadowFading = 3; % Shadow fading standard deviation
        else
            % NLOS conditions
            kFactor = (-100).*ones(1,9); % Rician K-factor
            shadowFading = 4; % Shadow fading standard deviation
        end  
    case 'Model-C'
        % Typical residential or small office environment, 30 ns RMS delay spread
        % Average power in linear scale
        pathPower = 10.^(.1*[0 -2.1715 -4.3429 -6.5144 -8.6859 -10.8574 -4.3899 -6.5614 -8.7329 -10.9043 -13.7147 -15.8862 -18.0577 -20.2291]);
        % Relative delay (ns)
        pathDelays = [0 10e-9 20e-9  30e-9 40e-9 50e-9 60e-9 70e-9 80e-9 90e-9 110e-9 140e-9 170e-9 200e-9];

        % Power roll-off coefficients
        powerPerAngle = [0 -2.1715 -4.3429 -6.5144 -8.6859 -10.8574 -13.0288 -15.2003 -17.3718 -19.5433 -Inf -Inf -Inf -Inf;
                        -Inf -Inf -Inf -Inf -Inf -Inf -5.0288 -7.2003 -9.3718 -11.5433 -13.7147 -15.8862 -18.0577 -20.2291];
        % Tx
        txAoDDegrees = [13.5312.*ones(1,10) -Inf.*ones(1,4);
                       -Inf.*ones(1,6) 56.4329.*ones(1,8)];
        txASDegrees  = [24.7897.*ones(1,10) -Inf.*ones(1,4);
                       -Inf.*ones(1,6) 22.5729.*ones(1,8)];
        % Rx
        rxAoADegrees = [290.3715.*ones(1,10) -Inf.*ones(1,4);
                       -Inf.*ones(1,6) 332.3754.*ones(1,8)];
        rxASDegrees  = [24.6949.*ones(1,10) -Inf.*ones(1,4);
                       -Inf.*ones(1,6) 22.4530.*ones(1,8)];

        % Path loss break point
        breakPointDistance = 5;
        if (distanceTxRx < breakPointDistance)
            % LOS conditions
            kFactor = [0,(-100).*ones(1,13)]; % Rician K-factor
            shadowFading = 3; % Shadow fading standard deviation
        else
            % NLOS conditions
            kFactor = (-100).*ones(1,14); % Rician K-factor
            shadowFading = 5; % Shadow fading standard deviation
        end
    case 'Model-D'
        % Typical office environment, 50 ns RMS delay spread
        % Average power in linear scale
        pathPower = 10.^(.1*[0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18 -22.4 -26.7]);
        % Relative delay (ns)
        pathDelays = [0 10e-9 20e-9 30e-9 40e-9 50e-9 60e-9 70e-9 80e-9 90e-9 110e-9 140e-9 170e-9 200e-9 240e-9 290e-9 340e-9 390e-9];
        % Power roll-off coefficients
        powerPerAngle = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -9.0712046 -11.199064 -13.795428 -16.391791 -19.370991 -23.201722 -Inf -Inf;
                        -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -6.6756386 -9.5728825 -12.175385 -14.777891 -17.435786 -21.992788 -25.580689 -Inf;
                        -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -18.843300 -23.238125 -25.246344 -26.7];
        % Tx
        txAoDDegrees = [332.1027.*ones(1,16) -Inf.*ones(1,2);
                       -Inf.*ones(1,10) 49.3840.*ones(1,7) -Inf;
                       -Inf.*ones(1,14) 275.9769.*ones(1,4)];
        txASDegrees  = [27.4412.*ones(1,16) -Inf.*ones(1,2);
                       -Inf.*ones(1,10) 32.1430.*ones(1,7) -Inf;
                       -Inf.*ones(1,14) 36.8825.*ones(1,4)];
        % Rx
        rxAoADegrees = [158.9318.*ones(1,16) -Inf.*ones(1,2);
                       -Inf.*ones(1,10) 320.2865.*ones(1,7) -Inf;
                       -Inf.*ones(1,14) 276.1246.*ones(1,4)];
        rxASDegrees  = [27.7580.*ones(1,16) -Inf.*ones(1,2);
                       -Inf.*ones(1,10) 31.4672.*ones(1,7) -Inf;
                       -Inf.*ones(1,14) 37.4179.*ones(1,4)];

        % Path loss break point
        breakPointDistance = 10;
        if (distanceTxRx < breakPointDistance)
            % LOS conditions
            kFactor = [3,(-100).*ones(1,17)]; % Rician K-factor
            shadowFading = 3; % Shadow fading standard deviation
        else
            % NLOS conditions
            kFactor = (-100).*ones(1,18); % Rician K-factor
            shadowFading = 5; % Shadow fading standard deviation
        end
    case 'Model-E'
        % Large office environments, 100 ns RMS delay spread
        % Average power in linear scale
        pathPower = 10.^(.1*[-2.5 -3.0 -3.5 -3.9 0 -1.3 -2.6 -3.9 -3.4 -5.6 -7.7 -9.9 -12.1 -14.3 -15.4 -18.4 -20.7 -24.6]);
        % Relative delay (ns)
        pathDelays = [0 10e-9 20e-9 30e-9 50e-9 80e-9 110e-9 140e-9 180e-9 230e-9 280e-9 330e-9 380e-9 430e-9 490e-9 560e-9 640e-9 730e-9];
        % Power roll-off coefficients
        powerPerAngle = [-2.6 -3 -3.5 -3.9 -4.5644301 -5.6551533 -6.9751533 -8.2951533 -9.8221791 -11.785521 -13.985521 -16.185521 -18.385521 -20.585521 -22.985195 -Inf -Inf -Inf;
                         -Inf -Inf -Inf -Inf -1.8681171 -3.2849115 -4.5733656 -5.8619031 -7.1920408 -9.9304493 -10.343797 -14.353720 -14.767068 -18.776991 -19.982151 -22.446411 -Inf -Inf;
                         -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -7.9044978 -9.6851670 -14.260649 -13.812819 -18.603831 -18.192376 -22.834619 -Inf -Inf -Inf;
                         -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -20.673366 -20.574381 -20.7 -24.6];
        % Tx
        txAoDDegrees = [105.6434.*ones(1,15) -Inf.*ones(1,3);
                      -Inf.*ones(1,4) 293.1199.*ones(1,12) -Inf.*ones(1,2);
                      -Inf.*ones(1,8) 61.9720.*ones(1,7) -Inf.*ones(1,3);
                      -Inf.*ones(1,14) 275.7640.*ones(1,4)];
        txASDegrees  = [36.1176.*ones(1,15) -Inf.*ones(1,3);
                      -Inf.*ones(1,4) 42.5299.*ones(1,12) -Inf.*ones(1,2);
                      -Inf.*ones(1,8) 38.0096.*ones(1,7) -Inf.*ones(1,3);
                      -Inf.*ones(1,14) 38.7026.*ones(1,4)];
        % Rx
        rxAoADegrees = [163.7475.*ones(1,15) -Inf.*ones(1,3);
                      -Inf.*ones(1,4) 251.8792.*ones(1,12) -Inf.*ones(1,2);
                      -Inf.*ones(1,8) 80.0240.*ones(1,7) -Inf.*ones(1,3);
                      -Inf.*ones(1,14) 182.0000.*ones(1,4)];
        rxASDegrees  = [35.8768.*ones(1,15) -Inf.*ones(1,3);
                      -Inf.*ones(1,4) 41.6812.*ones(1,12) -Inf.*ones(1,2);
                      -Inf.*ones(1,8) 37.4221.*ones(1,7) -Inf.*ones(1,3);
                      -Inf.*ones(1,14) 40.3685.*ones(1,4)];

       % Path loss break point
       breakPointDistance = 20;
       if (distanceTxRx < breakPointDistance)
            % LOS conditions
            kFactor = [6,(-100).*ones(1,17)]; % Rician K-factor
            shadowFading = 3; % Shadow fading standard deviation
       else
            % NLOS conditions
            kFactor = (-100).*ones(1,18); % Rician K-factor
            shadowFading = 6; % Shadow fading standard deviation
       end
    otherwise % 'Model-F'
        % Large space (indoor and outdoor), 150 ns RMS delay spread
        % Average power in linear scale
        pathPower = 10.^(.1*[-3.3 -3.6 -3.9 -4.2 0 -0.9 -1.7 -2.6 -1.5 -3.0 -4.4 -5.9 -5.3 -7.9 -9.4 -13.2 -16.3 -21.2]);
        % Relative delay (ns)
        pathDelays = [0 10e-9 20e-9 30e-9 50e-9 80e-9 110e-9 140e-9 180e-9 230e-9 280e-9 330e-9 400e-9 490e-9 600e-9 730e-9 880e-9 1050e-9];          
        % Power roll-off coefficients
        powerPerAngle = [-3.3 -3.6 -3.9 -4.2 -4.6474101 -5.393095 -6.293095 -7.193095 -8.2370567 -9.5792981 -11.079298 -12.579298 -14.358636 -16.731146 -19.978784 -Inf -Inf -Inf;
                         -Inf -Inf -Inf -Inf -1.8241629 -2.8069486 -3.5527879 -4.4527879 -5.3075764 -7.4275717 -7.0894233 -10.30481 -10.44412 -13.837355 -15.762121 -19.940313 -Inf -Inf;
                         -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -5.7960011 -6.7737346 -10.475827 -9.6416705 -14.107182 -12.752335 -18.503266 -Inf -Inf -Inf;
                         -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -8.8824645 -13.319464 -18.733410 -Inf -Inf -Inf;
                         -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -12.947155 -14.233751 -Inf -Inf;
                         -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -16.3 -21.2];
        % Tx
        txAoDDegrees = [56.2139.*ones(1,15) -Inf.*ones(1,3);
                       -Inf.*ones(1,4) 183.7089.*ones(1,12) -Inf.*ones(1,2);
                       -Inf.*ones(1,8) 153.0836.*ones(1,7) -Inf.*ones(1,3);
                       -Inf.*ones(1,12) 112.5317.*ones(1,3) -Inf.*ones(1,3);
                       -Inf.*ones(1,14) 291.0921.*ones(1,2) -Inf.*ones(1,2);
                       -Inf.*ones(1,16) 62.3790.*ones(1,2)];
        txASDegrees  = [41.6936.*ones(1,15) -Inf.*ones(1,3);
                       -Inf.*ones(1,4) 55.2669.*ones(1,12) -Inf.*ones(1,2);
                       -Inf.*ones(1,8) 47.4867.*ones(1,7) -Inf.*ones(1,3);
                       -Inf.*ones(1,12) 27.2136.*ones(1,3) -Inf.*ones(1,3);
                       -Inf.*ones(1,14) 33.0126.*ones(1,2) -Inf.*ones(1,2);
                       -Inf.*ones(1,16) 38.0482.*ones(1,2)];
        % Rx
        rxAoADegrees = [315.1048.*ones(1,15) -Inf.*ones(1,3);
                       -Inf.*ones(1,4) 180.4090.*ones(1,12) -Inf.*ones(1,2);
                       -Inf.*ones(1,8) 74.7062.*ones(1,7) -Inf.*ones(1,3);
                       -Inf.*ones(1,12) 251.5763.*ones(1,3) -Inf.*ones(1,3);
                       -Inf.*ones(1,14) 68.5751.*ones(1,2) -Inf.*ones(1,2);
                       -Inf.*ones(1,16) 246.2344.*ones(1,2)];
        rxASDegrees  = [48.0084.*ones(1,15) -Inf.*ones(1,3);
                       -Inf.*ones(1,4) 55.0823.*ones(1,12) -Inf.*ones(1,2);
                       -Inf.*ones(1,8) 42.0885.*ones(1,7) -Inf.*ones(1,3);
                       -Inf.*ones(1,12) 28.6161.*ones(1,3) -Inf.*ones(1,3);
                       -Inf.*ones(1,14) 30.7745.*ones(1,2) -Inf.*ones(1,2);
                       -Inf.*ones(1,16) 38.2914.*ones(1,2)];

        % Path loss break point
        breakPointDistance = 30;
        if (distanceTxRx < breakPointDistance)
            % LOS conditions
            kFactor = [6,(-100).*ones(1,17)]; % Rician K-factor
            shadowFading = 3; % Shadow fading standard deviation
        else
            % NLOS conditions
            kFactor = (-100).*ones(1,18); % Rician K-factor
            shadowFading = 6; % Shadow fading standard deviation
        end
    end

    % The AoA and AoD of the LOS component are fixed to 45 degrees (pi/4 in
    % radians) in TGn channel model. Ref: Implementation note version
    % 3.2-May 2004, page-25, Table-4.
    angleLOS = 45; % Fixed AoA and AoD of the LOS component in degrees
    if inp.UserIndex > 0
        % MU-MIMO model parameters

        % The range of LOS/NLOS diversity for both Tx and Rx is 360
        % degrees. Ref:IEEE 802.11-09/0308r12, page 12
        range = 360; % Range for LOS/NLOS diversity for both Tx and Rx 

        % Ref:IEEE 802.11-09/0308r12, page 12
        seedAPLOS = 39161; % Seed for LOS D/L AoD or U/L AoA offset
        seedAPNLOS = 2803; % Seed for NLOS D/L AoD or U/L AoA offset
        seedSTALOS = 45191; % Seed for LOS D/L AoA or U/L AoD offset
        seedSTANLOS = 13367; % Seed for NLOS D/L AoA or U/L AoD offset

        % Assign model parameters based on transmission direction
        if strcmp(inp.TransmissionDirection,'Downlink')
            seedTxLOS = seedAPLOS;
            seedTxNLOS = seedAPNLOS;
            seedRxLOS = seedSTALOS;
            seedRxNLOS = seedSTANLOS;
        else
            seedRxLOS = seedAPLOS;
            seedRxNLOS = seedAPNLOS;
            seedTxLOS = seedSTALOS;
            seedTxNLOS = seedSTANLOS;
        end

        % Tx LOS
        s = RandStream('v4','Seed',seedTxLOS); % Initialize generator
        randSequence1 = rand(s,inp.UserIndex,1);
        txOffsetLOSDegrees = (randSequence1(end)-0.5)*range;

        % Tx NLOS
        s = RandStream('v4','Seed',seedTxNLOS); % Initialize generator
        randSequence2 = rand(s,inp.UserIndex,1);
        txOffsetNLOSDegrees = (randSequence2(end)-0.5)*range;

        % Rx LOS
        s = RandStream('v4','Seed',seedRxLOS); % Initialize generator
        randSequence3 = rand(s,inp.UserIndex,1);
        rxOffsetLOSDegrees = (randSequence3(end)-0.5)*range;

        % Rx NLOS
        s = RandStream('v4','Seed',seedRxNLOS); % Initialize generator
        randSequence4 = rand(s,inp.UserIndex,1);
        rxOffsetNLOSDegrees = (randSequence4(end)-0.5)*range;
    else % UserIndex = 0
        % User Index 0 is reserved and denotes a user with channel
        % parameters equivalent to the TGn channel model, e.g. zero angle
        % offsets. Ref:IEEE 802.11-09/0308r12, page 7.
        txOffsetLOSDegrees = 0;
        txOffsetNLOSDegrees = 0;
        rxOffsetLOSDegrees = 0;
        rxOffsetNLOSDegrees = 0;
    end
    % For LOS. Ref:IEEE 802.11-09/0308r12, page 13
    txLOSRadians = deg2rad(mod(angleLOS + txOffsetLOSDegrees,360)); % In radians
    rxLOSRadians = deg2rad(mod(angleLOS + rxOffsetLOSDegrees,360)); % In radians

    % For NLOS
    txAoDDegrees = mod(txAoDDegrees + txOffsetNLOSDegrees,360);
    rxAoADegrees = mod(rxAoADegrees + rxOffsetNLOSDegrees,360);

    txAoDDegrees(isnan(txAoDDegrees)) = -Inf;
    rxAoADegrees(isnan(rxAoADegrees)) = -Inf;

    % PDP Tap spacing as per doc: IEEE 802.11-09/0308r1212 (Table-1, page-3)
    switch inp.ChannelBandwidth
        case {'CBW1','CBW2','CBW4','CBW20','CBW40'}
            tapSpacing = 10; % In nsec
        case {'CBW8','CBW80'}
            tapSpacing = 5; % In nsec
        otherwise
            tapSpacing = 2.5; % In nsec
    end

    % Tap interpolation (TGac)
    bwFactor = 10/tapSpacing; % Bandwidth expansion factor

    if tapSpacing < 10 && ~strcmp(inp.DelayProfile,'Model-A')

        numTGnPaths = size(powerPerAngle);

        pathDelaysTGac = [];
        for cluster = 1:numTGnPaths(1)
            clusterTaps = find(powerPerAngle(cluster,:) > -Inf);
            for tap = clusterTaps(1): clusterTaps(end-1)
               tap_times_TGac = (0:bwFactor-1)*(tapSpacing/1e9) + pathDelays(1,tap);
               pathDelaysTGac = [pathDelaysTGac tap_times_TGac];%#ok<AGROW>
            end
            pathDelaysTGac = [pathDelaysTGac pathDelays(1,clusterTaps(end))]; %#ok<AGROW>
        end

        pathDelaysTGac = sort(unique(pathDelaysTGac));
        numTGacPaths = [numTGnPaths(1) size(pathDelaysTGac,2)];

        powerPerAngleTGac = -Inf*ones(numTGacPaths);

        txAoDTGac = -Inf*ones(numTGacPaths);
        txASTGac = -Inf*ones(numTGacPaths);
        rxAoATGac = -Inf*ones(numTGacPaths);
        rxASTGac = -Inf*ones(numTGacPaths);

        % Generate interpolated matrices
        for cluster = 1:numTGnPaths(1)
            validTapsTGn = find(powerPerAngle(cluster,:) > -Inf);
            pathDelaysTGn = pathDelays(1,:);
            pathDelaysTGn(validTapsTGn(1)) = pathDelaysTGn(validTapsTGn(1))*(1-(1e-9));
            pathDelaysTGn(validTapsTGn(end)) = pathDelaysTGn(validTapsTGn(end))*(1+(1e-9));
            powerPerAngleTGac(cluster,:) = interp1(pathDelaysTGn,powerPerAngle(cluster,:),pathDelaysTGac);

            validTapsTGac = ~isnan(powerPerAngleTGac(cluster,:));
            powerPerAngleTGac(cluster,~validTapsTGac) = -Inf;
            clusterAngle = max(txAoDDegrees(cluster,:));
            txAoDTGac(cluster,validTapsTGac) = clusterAngle;
            clusterAngle = max(rxAoADegrees(cluster,:));

            rxAoATGac(cluster,validTapsTGac) = clusterAngle;
            clusterAS = max(rxASDegrees(cluster,:));
            rxASTGac(cluster,validTapsTGac) = clusterAS;

            clusterAS = max(txASDegrees(cluster,:));
            txASTGac(cluster,validTapsTGac) = clusterAS;
        end
        kFactorTGac = (-100).*ones(1,numTGacPaths(2));
        kFactorTGac(1) = kFactor(1);

        powerPerAngleTGacLinear = 10.^(powerPerAngleTGac/10);
        pathPowerTGac = sum(powerPerAngleTGacLinear);
    end

    if tapSpacing < 10 && ~strcmp(inp.DelayProfile,'Model-A')
        out.PathPower = cast(pathPowerTGac,inp.InputDataType);
        out.PathDelays = cast(pathDelaysTGac,inp.InputDataType);
        out.PowerPerAngle = cast(powerPerAngleTGac,inp.InputDataType);
        out.TxAoD = cast(txAoDTGac,inp.InputDataType);
        out.RxAoA = cast(rxAoATGac,inp.InputDataType);
        out.TxAS = cast(txASTGac,inp.InputDataType);
        out.RxAS = cast(rxASTGac,inp.InputDataType);
        out.Kfactor = cast(kFactorTGac,inp.InputDataType);
        out.TxLOSRadians = cast(txLOSRadians,inp.InputDataType);
        out.RxLOSRadians = cast(rxLOSRadians,inp.InputDataType);
    else
        out.PathPower = cast(pathPower,inp.InputDataType);
        out.PathDelays = cast(pathDelays,inp.InputDataType);
        out.PowerPerAngle = cast(powerPerAngle,inp.InputDataType);
        out.TxAoD = cast(txAoDDegrees,inp.InputDataType);
        out.RxAoA = cast(rxAoADegrees,inp.InputDataType);
        out.TxAS = cast(txASDegrees,inp.InputDataType);
        out.RxAS = cast(rxASDegrees,inp.InputDataType);
        out.Kfactor = cast(kFactor,inp.InputDataType);
        out.TxLOSRadians = cast(txLOSRadians,inp.InputDataType);
        out.RxLOSRadians = cast(rxLOSRadians,inp.InputDataType);
    end

    % Floor attenuation loss defined in doc IEEE 802.11-11/0968r4 and IEEE
    % 802.11-14/0882r4 for 802.11ah and 802.11ax respectively.
    floorAttenuationLoss = 0; % Prefix to no attenuation
    wallAttenuationLoss = 0; % Prefix to no attenuation

    % For 802.11ah channel.
    % Always assume single floor for channel model-C, D, E and F. The floor
    % separation property is hidden in this case.
    if any(strcmp(inp.FormatType,{'wlanTGahChannel','hTGahChannel_legacy'})) && any(strcmp(modelType,{'Model-C','Model-D','Model-E','Model-F'}))
        inp.NumPenetratedFloors = 0;
    end

    if any(strcmp(inp.FormatType,{'wlanTGahChannel','hTGahChannel_legacy'})) % For TGah Channel
        if inp.NumPenetratedFloors == 0
            shadowFading = shadowFading-1;
        end
        % Multi-floor scenario is only valid for channel model A and B
        if (inp.NumPenetratedFloors > 0) && any(strcmp(modelType,{'Model-A','Model-B'}))
            switch inp.NumPenetratedFloors
                case 1
                    shadowFading = 7;
                    floorAttenuationLoss = 12.9;
                case 2
                    shadowFading = 2.8;
                    floorAttenuationLoss = 18.7;
                case 3
                    shadowFading = 1.7;
                    floorAttenuationLoss = 24.4;
                otherwise
                    shadowFading = 1.5;
                    floorAttenuationLoss = 27.7;
            end
        end

        if ~(distanceTxRx > breakPointDistance)
            % Floor attenuation loss is not applicable when Tx/Rx separation
            % is less than the breakpoint distance for 802.11ah standard.
            % Ref: doc IEEE 802.11-11/0968r4
           floorAttenuationLoss = 0;
        end
    end

    if any(strcmp(inp.FormatType,{'wlanTGaxChannel','hTGaxChannel_legacy'}))
        % Floor attenuation loss. Ref: IEEE 802.11-14/0882r4
        n = inp.NumPenetratedFloors;
        floorAttenuationLoss = 18.3*n.^((n+2)/(n+1)-0.46); % In dBs

        % Wall attenuation loss
        wallAttenuationLoss = inp.NumPenetratedWalls*inp.WallPenetrationLoss; % In dBs
    end

    % Path loss calculations
    if (distanceTxRx > breakPointDistance)
        pathloss = 20*log10(4*pi*inp.CarrierFrequency/3e8)+20*log10(breakPointDistance)+(35*log10(distanceTxRx/breakPointDistance))+floorAttenuationLoss+wallAttenuationLoss;
    else
        pathloss = 20*log10(4*pi*inp.CarrierFrequency/3e8)+20*log10(distanceTxRx)+floorAttenuationLoss+wallAttenuationLoss;
    end

    out.Pathloss = cast(pathloss,inp.InputDataType);
    out.ShadowFading = cast(shadowFading,inp.InputDataType);
end

function out = sigmaValues

   out = [ ...
   0.001000000000000   0.001000000000000;
   0.002000000000000   0.002000000000000;
   0.003000000000000   0.003000000000000;
   0.004000000000000   0.004000000000000;
   0.005000000000000   0.005000000000000;
   0.006000000000000   0.006000000000000;
   0.007000000000000   0.007000000000000;
   0.008000000000000   0.008000000000000;
   0.009000000000000   0.009000000000000;
   0.010000000000000   0.010000000000000;
   0.020000000000000   0.020000000000000;
   0.030000000000000   0.030000000000000;
   0.040000000000000   0.040000000000000;
   0.050000000000000   0.050000000000000;
   0.060000000000000   0.060000000000000;
   0.069999999999996   0.070000000000000;
   0.079999999999674   0.080000000000000;
   0.089999999989971   0.090000000000000;
   0.099999999846174   0.100000000000000;
   0.109999998574985   0.110000000000000;
   0.119999990949232   0.120000000000000;
   0.129999956976519   0.130000000000000;
   0.139999837054695   0.140000000000000;
   0.149999485261405   0.150000000000000;
   0.159998596411514   0.160000000000000;
   0.169996608592331   0.170000000000000;
   0.179992589520795   0.180000000000000;
   0.189985119823364   0.190000000000000;
   0.199972187875949   0.200000000000000;
   0.209951109180660   0.210000000000000;
   0.219918479364342   0.220000000000000;
   0.229870165049242   0.230000000000000;
   0.239801332199394   0.240000000000000;
   0.249706507829649   0.250000000000000;
   0.259579668505112   0.260000000000000;
   0.269414347858725   0.270000000000000;
   0.279203755214353   0.280000000000000;
   0.288940898034351   0.290000000000000;
   0.298618702019563   0.300000000000000;
   0.308230124021961   0.310000000000000;
   0.317768254292421   0.320000000000000;
   0.327226405848202   0.330000000000000;
   0.336598189831417   0.340000000000000;
   0.345877576608532   0.350000000000000;
   0.355058943029479   0.360000000000000;
   0.364137106739593   0.370000000000000;
   0.373107348745202   0.380000000000000;
   0.381965425604661   0.390000000000000;
   0.390707572681411   0.400000000000000;
   0.399330499881729   0.410000000000000;
   0.407831381230977   0.420000000000000;
   0.416207839537671   0.430000000000000;
   0.424457927269851   0.440000000000000;
   0.432580104634431   0.450000000000000;
   0.440573215715887   0.460000000000000;
   0.448436463401379   0.470000000000000;
   0.456169383699028   0.480000000000000;
   0.463771819946593   0.490000000000000;
   0.471243897310273   0.500000000000000;
   0.478585997887944   0.510000000000000;
   0.485798736657463   0.520000000000000;
   0.492882938447934   0.530000000000000;
   0.499839616059197   0.540000000000000;
   0.506669949611134   0.550000000000000;
   0.513375267168660   0.560000000000000;
   0.519957026659509   0.570000000000000;
   0.526416799079023   0.580000000000000;
   0.532756252958393   0.590000000000000;
   0.538977140059256   0.600000000000000;
   0.545081282247542   0.610000000000000;
   0.551070559492422   0.620000000000000;
   0.556946898931507   0.630000000000000;
   0.562712264940743   0.640000000000000;
   0.568368650146199   0.650000000000000;
   0.573918067314977   0.660000000000000;
   0.579362542063437   0.670000000000000;
   0.584704106322536   0.680000000000000;
   0.589944792502317   0.690000000000000;
   0.595086628300083   0.700000000000000;
   0.600131632099618   0.710000000000000;
   0.605081808911741   0.720000000000000;
   0.609939146809489   0.730000000000000;
   0.614705613814253   0.740000000000000;
   0.619383155192148   0.750000000000000;
   0.623973691122821   0.760000000000000;
   0.628479114705672   0.770000000000000;
   0.632901290271164   0.780000000000000;
   0.637242051967444   0.790000000000000;
   0.641503202594861   0.800000000000000;
   0.645686512663307   0.810000000000000;
   0.649793719649333   0.820000000000000;
   0.653826527432048   0.830000000000000;
   0.657786605888584   0.840000000000000;
   0.661675590631663   0.850000000000000;
   0.665495082873338   0.860000000000000;
   0.669246649400458   0.870000000000000;
   0.672931822648735   0.880000000000000;
   0.676552100863522   0.890000000000000;
   0.680108948336552   0.900000000000000;
   0.683603795708888   0.910000000000000;
   0.687038040331338   0.920000000000000;
   0.690413046674382   0.930000000000000;
   0.693730146780515   0.940000000000000;
   0.696990640752585   0.950000000000000;
   0.700195797272374   0.960000000000000;
   0.703346854144260   0.970000000000000;
   0.706445018859369   0.980000000000000;
   0.709491469176061   0.990000000000000;
   0.712487353713093   1.000000000000000;
   0.715433792552187   1.010000000000000;
   0.718331877847086   1.020000000000000;
   0.721182674436521   1.030000000000000;
   0.723987220458834   1.040000000000000;
   0.726746527966224   1.050000000000000;
   0.729461583536870   1.060000000000000;
   0.732133348883386   1.070000000000000;
   0.734762761456275   1.080000000000000;
   0.737350735041188   1.090000000000000;
   0.739898160349028   1.100000000000000;
   0.742405905597993   1.110000000000000;
   0.744874817086865   1.120000000000000;
   0.747305719758908   1.130000000000000;
   0.749699417755864   1.140000000000000;
   0.752056694961647   1.150000000000000;
   0.754378315535372   1.160000000000000;
   0.756665024433481   1.170000000000000;
   0.758917547920750   1.180000000000000;
   0.761136594070047   1.190000000000000;
   0.763322853250757   1.200000000000000] *1e2;
end