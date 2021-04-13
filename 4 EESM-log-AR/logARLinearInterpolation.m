% Does not deal with outliers
% Does not remove trend before processing
% Gaussian innovation AR
% Gaussian innovation AR can produce zero-mean residual
% Standarized residual has standard devidation 1
% Residual looks white
% Use KPSS and Leybourne-McCabe stationarity test 
% in econometricModeler to check stationarity;
% choose no trend when checking stationarity.
% Both PHY layer RPs: effSnrVecLog_d and perStore are stationary
%% Load parameters at gamma_1 and gamma_2
arLags = 10;
innovdist = struct('Name',"Gaussian");
load('snrPer_CBW20_Model-D_4-by-2_MCS4.mat')
% Load data at gamma_1
snrIdx = 1; 
effSnrVec1 = results{snrIdx}.effSnrVec;
effSnrVecdB1 = effSnrVec1';
effSnrVecLinear1 = 10.^(effSnrVecdB1/10); % Transfer dB into linear
effSnrVecLog1 = log(effSnrVecLinear1); % Transfer linear into log domain
Mdl1 = arima('ARLags',1:arLags,'Distribution',innovdist);
EstMdl1 = estimate(Mdl1,effSnrVecLog1)
constLI1 = EstMdl1.Constant;
varLI1 = EstMdl1.Variance;
arLI1 = EstMdl1.AR;
% Load data at gamma_2
snrIdx = 3; 
effSnrVec2 = results{snrIdx}.effSnrVec;
effSnrVecdB2 = effSnrVec2';
effSnrVecLinear2 = 10.^(effSnrVecdB2/10); % Transfer dB into linear
effSnrVecLog2 = log(effSnrVecLinear2); % Transfer linear into log domain
Mdl2 = arima('ARLags',1:arLags,'Distribution',innovdist);
EstMdl2 = estimate(Mdl2,effSnrVecLog2)
constLI2 = EstMdl2.Constant;
varLI2 = EstMdl2.Variance;
arLI2 = EstMdl2.AR;
%% Parameter linear interpolation (LI)
constLI = (constLI1 + constLI2)/2;
varLI = (varLI1 + varLI2)/2;
for m = 1:arLags
    arLI{m} = (arLI1{m} + arLI2{m})/2;
end
%% Loading effective SNR in dB
load('snrPer_CBW20_Model-D_4-by-2_MCS4.mat')
snrIdx = 2;
effSnrVec = results{snrIdx}.effSnrVec;
effSnrVecdB = effSnrVec';
effSnrVecLinear = 10.^(effSnrVecdB/10); % Transfer dB into linear
effSnrVecLog = log(effSnrVecLinear); % Transfer linear into log domain
%% AR model with white Gaussian innovation
innovdist = struct('Name',"Gaussian");
Mdl = arima('ARLags',1:arLags,'Distribution',innovdist);
EstMdl = estimate(Mdl,effSnrVecLog)
EstMdl.Constant = constLI;
EstMdl.Variance = varLI;
EstMdl.AR = arLI;
%% Analysis of residuals: whiteness test
% Source: https://www.mathworks.com/help/econ/residual-diagnostics.html#btb6li_
% In time series models, the innovation process is assumed to be uncorrelated
% As an informal check,you can plot the sample ACF and PACF.
% If either plot shows significant autocorrelation in the residuals, 
% you can consider modifying your model to include additional AR or MA terms.
% Sample ACF and PACF plot
[res,condVar] = infer(EstMdl,effSnrVecLog);
% stdr = res/sqrt(EstMdl.Variance);
figure
subplot(2,1,1)
autocorr(res,'NumLags',40)
subplot(2,1,2)
parcorr(res,'NumLags',40)
sgtitle('Whiteness test')
% More formally, you can conduct a Ljung-Box Q-test on the residual series.
% This tests the null hypothesis of jointly zero autocorrelations up to lag m, 
% against the alternative of at least one nonzero autocorrelation.
% If you obtain res by fitting a model to data, then you should reduce the degrees of freedom (the argument DoF) by the number of estimated coefficients, excluding constants. 
% For example, if you obtain res by fitting an ARMA(p,q) model, 
% set DoF to L−p−q, where L is Lags.
% Ljung-Box Q-test for residual autocorrelation
lags = 30;         
dof  = lags - arLags; % One autoregressive parameter
[h_LB,p_LB] = lbqtest(res,'Lags',lags,'DOF',dof)
%% Analysis of residuals: heteroscedasticity (nonconstant variance) test
% Source: https://www.mathworks.com/help/econ/residual-diagnostics.html#btb6li_
% A white noise innovation process has constant variance
% After fitting a model, you can infer residuals and check them for heteroscedasticity (nonconstant variance)
% As an informal check, you can plot the sample ACF and PACF of the squared residual series
% If either plot shows significant autocorrelation, 
% you can consider modifying your model to include a conditional variance process.
% ACF and PACF plot
figure
subplot(2,1,1)
autocorr(res.^2,'NumLags',40)
subplot(2,1,2)
parcorr(res.^2,'NumLags',40)
sgtitle('Nonconstant variance test')
% More formally, you can conduct an Engle’s ARCH test on the residual series.
% This tests the null hypothesis of no ARCH effects against the alternative ARCH model with k lags.
%% Residual distribution fitting: PDF fitting
figure
sampleRes = linspace(min(res),...
    max(res),10000).';
histogram(res,'Normalization','pdf')
xlabel('Residual')
ylabel('PDF')
% s_stdr = skewness(stdr);
% k_stdr = kurtosis(stdr);
% title(sprintf('\\bf Residual Distribution (Skewness $\\gamma$ = %4f, Kurtosis $\\tau$ = %4f)',s_stdr,k_stdr),'Interpreter','latex')
grid on
hold on
normalFit = fitdist(res,'Normal');
normalPDF = pdf(normalFit,sampleRes);
plot(sampleRes,normalPDF,'LineWidth',2)
legend({'Standardized residual','Normal Fitting'})
hold off
%% Residual distribution fitting: CDF fitting
figure
[stdResCDF,queryRes] = ecdf(res); % Empirical CDF
normalFitCDF = cdf(normalFit,sampleRes); % Fitted normal Distribution
stairs(queryRes,stdResCDF,'LineWidth', 1.5)
hold on
plot(sampleRes,normalFitCDF,'LineWidth', 1.5)
xlabel('Residual')
ylabel('CDF')
title('Empirical and Fitted CDFs of Residuals')
grid on
legend({'Residual','Normal Fitting'},...
        'Location','southeast')
normalFitCDF = cdf(normalFit,queryRes);  % Fitted cdf with respect to empirical cdf query points         
% Obtain maximum discrepancy of normal fitting
[maxDiscrepNormal,maxPosNormal] = max(abs(stdResCDF - normalFitCDF));
resAtMaxDiffNormal = queryRes(maxPosNormal);
% Plot maximum discrepancy of normal fitting
plot(resAtMaxDiffNormal*ones(1, 2),...
    [stdResCDF(maxPosNormal) normalFitCDF(maxPosNormal)],'k','LineWidth',2,...
    'DisplayName',['Max Discrepancy of Normal Fitting (%): ',num2str(100*maxDiscrepNormal,'%.2f')])        
%% Residual distribution fitting: Test
% % One-sample Kolmogorov-Smirnov test
% test_normal_cdf = [stdr,cdf('Normal',stdr,normalFit.mu,normalFit.sigma)];
% [h_KS,p_KS] = kstest(stdr,'CDF',test_normal_cdf)
% % Chi-square goodness-of-fit test
% [h_Chi,p_Chi] = chi2gof(stdr,'CDF',normalFit)
%% Plot effective SNR (dB) VS samples, and obtain average PER
% Setup for obtaining PER
channelCoding = cfgHE.ChannelCoding; 
dataLength = 1000;  
format = 'HE_SU'; % hard code for SUConfig
BW = cfgHE.ChannelBandwidth;
abstraction = tgaxEESMLinkPerformanceModel;
eesmPer = results{snrIdx}.packetErrorRateAbs
% Sample interval
sampleInterval = 250/10^3;
numSteps = length(effSnrVecdB);
sampleTime = [1:numSteps]*sampleInterval;
% 1) Original effective SNR
figure
plot(sampleTime,effSnrVecdB);
hold on
% 2) AR model with white Gaussian innovation
tStart = tic;
% arEffSnrVecLog = simulate(EstMdl,size(effSnrVecLog,1),'E0',res,'V0',condVar,'Y0',effSnrVecLog);
arEffSnrVecLog = simulate(EstMdl,size(effSnrVecLog,1));
arEffSnrVecLinear = exp(arEffSnrVecLog);
arEffSnrVecdB = 10*log10(arEffSnrVecLinear);
[arPerVec,arPerPL0Vec,L0,lut] = estimatePER(abstraction,arEffSnrVecdB,format,mcs,channelCoding,dataLength);
arPer = mean(arPerVec)
tEnd = toc(tStart);
plot(sampleTime,arEffSnrVecdB)
% Label and legend
xlabel('Time (sec)')
ylabel('Effective SNR (dB)')
legend({'EESM','EESM-log-AR'})
%% Plot ACF and PACF of AR-estimated PER
arPacketError = rand(length(arPerVec),1)<=arPerVec;
arPacketError = double(arPacketError);
figure
subplot(2,1,1)
autocorr(arPacketError,'NumLags',40)
subplot(2,1,2)
parcorr(arPacketError,'NumLags',40) 
sgtitle('ACF and PACF of AR-estimated PER')