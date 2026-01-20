const codes = {

exp7_ofdm: `
clc;
clear;
close all;

%% Parameters
numBits = 1e5;                        % Number of bits
modSchemes = {'QPSK', '16QAM', '64QAM'}; % Modulation schemes
M_vals = [4, 16, 64];                 % Modulation orders for QPSK, 16QAM, 64QAM
numSNR = 40;                          % Number of SNR points
SNR_dB = linspace(0, 40, numSNR);     % SNR range in dB
numTrials = 100;                       % Monte Carlo trials

% OFDM parameters
N = 64;                               % Number of subcarriers (FFT size)
CP_len = 16;                          % Cyclic prefix length
numOFDMsymbols = floor(numBits ./ log2(M_vals(1)) / N); % Number of OFDM symbols

BER = zeros(length(modSchemes), numSNR);

% Store time-domain signals for comparison plot
txSignalsIFFT = cell(1,length(modSchemes));
txSignalsIFFT_CP = cell(1,length(modSchemes));

for mIndex = 1:length(modSchemes)
    M = M_vals(mIndex);             % Current modulation order
    k = log2(M);                   % Bits per symbol
    
    % Generate data bits for current modulation
    curNumBits = numOFDMsymbols * N * k;
    
    % Map constellation once for TX plots
    data = randi([0 1], curNumBits, 1);
    if M == 4
        modSymbols = pskmod(bi2de(reshape(data, k, []).','left-msb'), M, pi/4);
    else
        modSymbols = qammod(bi2de(reshape(data, k, []).','left-msb'), M, 'gray');
    end
    
    % Transmit Constellation Diagram (before IFFT)
    figure;
    scatterplot(modSymbols(1:5000));
    title(['Transmit Constellation - ' modSchemes{mIndex}]);
    grid on;
    
    % Reshape symbols into OFDM symbols
    modSymbolsMatrix = reshape(modSymbols, N, []);
    
    % Perform IFFT (OFDM modulation)
    txOFDM = ifft(modSymbolsMatrix, N, 1);
    
    % Save time-domain OFDM signal without CP (for plotting)
    txSigIFFT = txOFDM(:);
    txSignalsIFFT{mIndex} = txSigIFFT(1:5*N);
    
    % Add cyclic prefix
    cyclicPrefix = txOFDM(end-CP_len+1:end, :);
    txOFDM_CP = [cyclicPrefix; txOFDM];
    
    % Save time-domain OFDM signal with CP (for plotting)
    txSigIFFT_CP = txOFDM_CP(:);
    txSignalsIFFT_CP{mIndex} = txSigIFFT_CP(1:5*(N+CP_len));
    
    % ---------------- MONTE CARLO BER-SNR ---------------- %
    for snrIdx = 1:numSNR
        errCount = 0;
        bitCount = 0;
        
        for trial = 1:numTrials
            % New random bits for each trial
            data = randi([0 1], curNumBits, 1);

            % Modulate
            if M == 4
                modSymbols = pskmod(bi2de(reshape(data, k, []).','left-msb'), M, pi/4);
            else
                modSymbols = qammod(bi2de(reshape(data, k, []).','left-msb'), M, 'gray');
            end
            
            % Reshape into OFDM symbols
            modSymbolsMatrix = reshape(modSymbols, N, []);
            
            % IFFT
            txOFDM = ifft(modSymbolsMatrix, N, 1);
            
            % Add cyclic prefix
            cyclicPrefix = txOFDM(end-CP_len+1:end, :);
            txOFDM_CP = [cyclicPrefix; txOFDM];
            
            % Serialize for transmission
            txSignal = txOFDM_CP(:);
            
            % AWGN channel
            rxSignal = awgn(txSignal, SNR_dB(snrIdx), 'measured');
            
            % Reshape to OFDM symbols with CP
            rxOFDM_CP = reshape(rxSignal, N+CP_len, []);
            
            % Remove CP
            rxOFDM = rxOFDM_CP(CP_len+1:end, :);
            
            % FFT
            rxSymbolsMatrix = fft(rxOFDM, N, 1);
            
            % Serialize
            rxSymbols = rxSymbolsMatrix(:);
            
            % Demodulate
            if M == 4
                rxDec = pskdemod(rxSymbols, M, pi/4);
            else
                rxDec = qamdemod(rxSymbols, M, 'gray');
            end
            
            % Convert symbols to bits
            rxBitsMat = de2bi(rxDec, k, 'left-msb').';
            rxBits = rxBitsMat(:);
            
            % Count errors
            errCount = errCount + sum(rxBits ~= data);
            bitCount = bitCount + length(data);
        end
        
        % Monte Carlo BER
        BER(mIndex, snrIdx) = errCount / bitCount;
    end
    
    % Received Constellation plot
    figure;
    validSymbols = rxSymbols(~isnan(rxSymbols) & ~isinf(rxSymbols));
    scatterplot(validSymbols(1:min(5000,numel(validSymbols))));
    title(['Received Constellation - ' modSchemes{mIndex} ' OFDM']);
    grid on;
end

% BER vs SNR plot
figure;
semilogy(SNR_dB, BER(1,:), '-o', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, BER(2,:), '-s', 'LineWidth', 2);
semilogy(SNR_dB, BER(3,:), '-d', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for OFDM with Different Modulation Schemes (Monte Carlo)');
legend(modSchemes, 'Location', 'southwest');
axis([0 40 1e-5 1]);
hold off;

% Plot comparison of time-domain OFDM waveforms (with and without CP)
figure;
colors = ['b','r','k'];
for i = 1:length(modSchemes)
    subplot(3,1,i);
    plot(real(txSignalsIFFT{i}), 'Color', colors(i), 'LineWidth', 1.5);
    hold on;
    plot(real(txSignalsIFFT_CP{i}), '--', 'Color', colors(i), 'LineWidth', 1.5);
    hold off;
    title([modSchemes{i} ': OFDM Time-Domain Waveforms']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    legend('Without CP','With CP');
    grid on;
end

%% Command line output analysis:
fprintf('\nAnalysis:\n');
fprintf('Monte Carlo averaging over %d trials makes BER vs SNR smoother.\n', numTrials);
fprintf('Higher-order modulations (16-QAM, 64-QAM) give higher data rates but need more SNR.\n');
fprintf('QPSK is more robust at low SNR.\n');
fprintf('OFDM with CP combats ISI effectively.\n');

	
`,

exp8_cu_du: `
clc;
clear;
close all;

% Required Toolbox: 5G Toolbox, Communications Toolbox

%% Simulation Parameters
simTime = 1;             % Simulation duration in seconds
numFrames = 20;          % Number of 5G NR frames simulated
frameDuration = 10e-3;   % 10 ms frame duration in 5G NR
fronthaulLatencyRange = 0:0.5:5;  % Fronthaul latency values in ms

% Data rates (bps) assumptions for illustration
baseDataRateSA = 1e9;    % 1 Gbps baseline for Standalone
baseDataRateNSA = 1.2e9; % 1.2 Gbps baseline for Non-Standalone (due to LTE anchor)

% Latency assumptions (processing + transport) without fronthaul delay (in ms)
processingLatencyCU = 0.1;
processingLatencyDU = 0.1;

%% Initialize arrays to hold results
latencySA = zeros(size(fronthaulLatencyRange));
latencyNSA = zeros(size(fronthaulLatencyRange));
throughputSA = zeros(size(fronthaulLatencyRange));
throughputNSA = zeros(size(fronthaulLatencyRange));

%% Simulate varying fronthaul latency impact
for idx = 1:length(fronthaulLatencyRange)
    fhLatency = fronthaulLatencyRange(idx);
    
    % Total latency = fronthaul + CU/DU processing delays
    totalLatencySA = processingLatencyCU + processingLatencyDU + fhLatency;
    totalLatencyNSA = processingLatencyCU + processingLatencyDU + fhLatency / 2; % NSA assumed less latency
    
    latencySA(idx) = totalLatencySA;
    latencyNSA(idx) = totalLatencyNSA;
    
    % Throughput affected inversely by latency (simplified)
    throughputSA(idx) = baseDataRateSA * (1 - min(totalLatencySA/10, 0.9));
    throughputNSA(idx) = baseDataRateNSA * (1 - min(totalLatencyNSA/10, 0.9));
end

%% Plot Latency vs Fronthaul Latency
figure;
plot(fronthaulLatencyRange, latencySA, '-o', 'LineWidth', 2, 'DisplayName', 'Standalone (SA)');
hold on;
plot(fronthaulLatencyRange, latencyNSA, '-s', 'LineWidth', 2, 'DisplayName', 'Non-Standalone (NSA)');
grid on;
xlabel('Fronthaul Latency (ms)');
ylabel('Total Latency (ms)');
title('Latency vs Fronthaul Latency for SA and NSA 5G Modes');
legend('Location','northwest');
hold off;

%% Plot Throughput vs Fronthaul Latency
figure;
plot(fronthaulLatencyRange, throughputSA/1e9, '-o', 'LineWidth', 2, 'DisplayName', 'Standalone (SA)');
hold on;
plot(fronthaulLatencyRange, throughputNSA/1e9, '-s', 'LineWidth', 2, 'DisplayName', 'Non-Standalone (NSA)');
grid on;
xlabel('Fronthaul Latency (ms)');
ylabel('Throughput (Gbps)');
title('Throughput vs Fronthaul Latency for SA and NSA 5G Modes');
legend('Location','northeast');
hold off;

%% Protocol Stack Visualization

% Define layers for user plane and control plane stacks for SA and NSA
layers = {'RLC', 'PDCP', 'SDAP', 'RRC', 'NAS'}; % Common control/user plane layers

% Standalone mode stacks (CU + DU)
SA_UserPlane = {'PHY-DU', 'MAC-DU', 'RLC-CU', 'PDCP-CU', 'SDAP-CU'};
SA_ControlPlane = {'PHY-DU', 'MAC-DU', 'RLC-CU', 'PDCP-CU', 'RRC-CU', 'NAS'};

% Non-Standalone mode stacks (DU linked to LTE eNodeB)
NSA_UserPlane = {'PHY-DU', 'MAC-DU', 'RLC-LTE', 'PDCP-LTE', 'SDAP-LTE'};
NSA_ControlPlane = {'PHY-DU', 'MAC-DU', 'RLC-LTE', 'PDCP-LTE', 'RRC-LTE', 'NAS'};

% Create figure for protocol stack
figure('Name','Protocol Stack Visualization','Position',[100, 100, 900, 600]);

subplot(2,2,1);
barh(1:length(SA_UserPlane), ones(1,length(SA_UserPlane)), 'FaceColor', [0 0.4470 0.7410]);
set(gca,'yticklabel',SA_UserPlane,'ytick',1:length(SA_UserPlane));
title('SA User Plane Stack');
xlabel('Layer Activity');
xlim([0 1.5]);

subplot(2,2,2);
barh(1:length(SA_ControlPlane), ones(1,length(SA_ControlPlane)), 'FaceColor', [0.8500 0.3250 0.0980]);
set(gca,'yticklabel',SA_ControlPlane,'ytick',1:length(SA_ControlPlane));
title('SA Control Plane Stack');
xlabel('Layer Activity');
xlim([0 1.5]);

subplot(2,2,3);
barh(1:length(NSA_UserPlane), ones(1,length(NSA_UserPlane)), 'FaceColor', [0.4660 0.6740 0.1880]);
set(gca,'yticklabel',NSA_UserPlane,'ytick',1:length(NSA_UserPlane));
title('NSA User Plane Stack');
xlabel('Layer Activity');
xlim([0 1.5]);

subplot(2,2,4);
barh(1:length(NSA_ControlPlane), ones(1,length(NSA_ControlPlane)), 'FaceColor', [0.4940 0.1840 0.5560]);
set(gca,'yticklabel',NSA_ControlPlane,'ytick',1:length(NSA_ControlPlane));
title('NSA Control Plane Stack');
xlabel('Layer Activity');
xlim([0 1.5]);

sgtitle('5G Base Station Protocol Stack: CU-DU Split in SA vs NSA Modes');

%% Command window summary
fprintf('\n5G CU-DU Split Architecture Simulation Analysis:\n');
fprintf(' - Fronthaul latency linearly increases total latency in both SA and NSA modes.\n');
fprintf(' - NSA mode exhibits lower latency due to reliance on LTE core and optimized processing split.\n');
fprintf(' - Throughput decreases as fronthaul latency increases, with SA mode being slightly lower throughput\n');
fprintf('   due to more end-to-end 5G stack processing.\n');
fprintf(' - Protocol stack visualization shows additional layers involved in SA mode compared to NSA.\n');
fprintf(' - This simulation provides a simplified but insightful view of the impact of CU-DU split and deployment modes.\n\n');

	
`,

exp9_numerology: `
% Experiment 9: 5G Numerology and Frame Structure
% Prerequisites: 5G Toolbox
clear; close all; clc;

% Numerologies to test (μ = 0 to 4)
mu_values = 0:4;
carrierFreq = 3.5e9; % Carrier frequency 3.5 GHz
nRB = 52;            % Number of resource blocks (common for FR1)
nAntennas = 1;

for mu = mu_values
    fprintf('\n--- Numerology µ = %d ---\n', mu);

    % Subcarrier Spacing (SCS) in kHz
    scs = 15 * 2^mu;          % in kHz
    Tslot = 1e-3 / (2^mu);    % Slot duration in seconds
    fprintf('Subcarrier Spacing: %d kHz\n', scs);
    fprintf('Slot Duration: %.3f ms\n', Tslot*1e3);

    % Carrier and grid configuration
    carrier = nrCarrierConfig;
    carrier.SubcarrierSpacing = scs;
    carrier.NSizeGrid = nRB;
    carrier.NSlot = 0;
    carrier.NStartGrid = 0;

    % Get OFDM info (for sample rate)
    ofdmInfo = nrOFDMInfo(carrier);

    % Create empty resource grid
    grid = nrResourceGrid(carrier, nAntennas);
    [K, L, R] = size(grid);

    % Fill grid with dummy complex QPSK data (for visualization)
    dataSym = qammod(randi([0 3], K*L,1), 4, 'UnitAveragePower', true);
    grid(:,:,1) = reshape(dataSym, K, L);

    % Time-domain signal generation
    txWaveform = nrOFDMModulate(carrier, grid);

    % --- Plot time-domain waveform ---
    figure;
    t = (0:length(txWaveform)-1)/ofdmInfo.SampleRate;
    plot(t*1e3, abs(txWaveform),'b');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    grid on;
    title(sprintf('Time-domain Signal | µ = %d, SCS = %d kHz', mu, scs));

    % --- Plot resource grid magnitude ---
    figure;
    imagesc(0:L-1, 0:K-1, abs(grid(:,:,1)));
    axis xy; colormap(jet); colorbar;
    xlabel('OFDM Symbol Index');
    ylabel('Subcarrier Index');
    title(sprintf('Resource Grid Magnitude | µ = %d, SCS = %d kHz', mu, scs));
end

	
`,

exp10_csi_rs: `

clc; clear; close all;

%% Parameters
N_tx = 4;              
N_rx = 1;              
N_sc = 256;            
SNR_demo = 10;         
M = 4;                 
SNR_dB_range = 0:5:30; 
N_realizations = 50;  

%% Generate CSI-RS pilot (QPSK)
data_symbols = randi([0 M-1],1,N_sc);
pilot = exp(1j*2*pi*data_symbols/M);  

%% Rayleigh channel and noise
H_demo = (randn(N_rx,N_tx,N_sc)+1j*randn(N_rx,N_tx,N_sc))/sqrt(2);
noise_demo = (randn(N_rx,N_sc)+1j*randn(N_rx,N_sc))/sqrt(2)*10^(-SNR_demo/20);

%% Received signal before beamforming
rx_signal_demo = zeros(N_rx,N_sc);
rx_signal_MRT_demo = zeros(N_rx,N_sc);
rx_signal_EGT_demo = zeros(N_rx,N_sc);
for sc=1:N_sc
    tx_signal_demo = repmat(pilot,N_tx,1);
    
    % Received without beamforming
    rx_signal_demo(:,sc) = H_demo(:,:,sc)*tx_signal_demo(:,sc) + noise_demo(:,sc);
    
    % MRT beamforming
    w_mrt = H_demo(:,:,sc)'; w_mrt = w_mrt/norm(w_mrt);
    rx_signal_MRT_demo(:,sc) = H_demo(:,:,sc)*w_mrt.*pilot(sc);
    
    % EGT beamforming
    w_egt = exp(1j*angle(H_demo(:,:,sc)')); w_egt = w_egt/norm(w_egt);
    rx_signal_EGT_demo(:,sc) = H_demo(:,:,sc)*w_egt.*pilot(sc);
end

%% Figures 1-6: Step-by-step signals

% Figure 1: Original Transmitted Signal
figure('Name','Figure 1: Transmitted Signal','NumberTitle','off');
plot(real(pilot),'b'); hold on; plot(imag(pilot),'r');
title('Original Transmitted CSI-RS Signal'); xlabel('Subcarrier'); ylabel('Amplitude');
legend('Real','Imag'); grid on;

% Figure 2: Noise
figure('Name','Figure 2: Noise','NumberTitle','off');
plot(real(noise_demo),'k'); hold on; plot(imag(noise_demo),'m');
title('Noise Signal'); xlabel('Subcarrier'); ylabel('Amplitude');
legend('Real','Imag'); grid on;

% Figure 3: Received Signal (Tx + Noise)
figure('Name','Figure 3: Received Signal','NumberTitle','off');
plot(real(rx_signal_demo),'b'); hold on; plot(imag(rx_signal_demo),'r');
title('Received Signal = Transmitted + Noise'); xlabel('Subcarrier'); ylabel('Amplitude');
legend('Real','Imag'); grid on;

% Figure 4: No Beamforming
figure('Name','Figure 4: No Beamforming','NumberTitle','off');
plot(real(rx_signal_demo),'b'); hold on; plot(imag(rx_signal_demo),'r');
title('Received Signal Without Beamforming'); xlabel('Subcarrier'); ylabel('Amplitude');
legend('Real','Imag'); grid on;

% Figure 5: MRT Beamforming
figure('Name','Figure 5: MRT Beamforming','NumberTitle','off');
plot(real(rx_signal_MRT_demo),'b'); hold on; plot(imag(rx_signal_MRT_demo),'r');
title('MRT Beamformed Signal'); xlabel('Subcarrier'); ylabel('Amplitude');
legend('Real','Imag'); grid on;

% Figure 6: EGT Beamforming
figure('Name','Figure 6: EGT Beamforming','NumberTitle','off');
plot(real(rx_signal_EGT_demo),'b'); hold on; plot(imag(rx_signal_EGT_demo),'r');
title('EGT (Phase-only) Beamformed Signal'); xlabel('Subcarrier'); ylabel('Amplitude');
legend('Real','Imag'); grid on;

%% 7. Constellation Plots
figure('Name','Constellations','NumberTitle','off');
subplot(1,3,1);
plot(pilot,'ko','MarkerFaceColor','k'); hold on;
plot(rx_signal_demo,'ro'); title('Before Beamforming'); grid on; axis equal;
xlabel('Re'); ylabel('Im'); legend('Ideal','Received');

subplot(1,3,2);
plot(pilot,'ko','MarkerFaceColor','k'); hold on;
plot(rx_signal_MRT_demo,'bo'); title('MRT Beamformed'); grid on; axis equal;
xlabel('Re'); ylabel('Im'); legend('Ideal','Received');

subplot(1,3,3);
plot(pilot,'ko','MarkerFaceColor','k'); hold on;
plot(rx_signal_EGT_demo,'mo'); title('EGT Beamformed'); grid on; axis equal;
xlabel('Re'); ylabel('Im'); legend('Ideal','Received');

%% 8. SER vs SNR
% Initialize SER arrays
SER_before = zeros(size(SNR_dB_range));
SER_MRT = zeros(size(SNR_dB_range));
SER_EGT = zeros(size(SNR_dB_range));

for idx_snr = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(idx_snr);
    errors_before = 0; errors_MRT = 0; errors_EGT = 0;
    
    for r = 1:N_realizations
        tx_signal = repmat(pilot,N_tx,1);
        H = (randn(N_rx,N_tx,N_sc)+1j*randn(N_rx,N_tx,N_sc))/sqrt(2);
        noise = (randn(N_rx,N_sc)+1j*randn(N_rx,N_sc))/sqrt(2)*10^(-SNR_dB/20);
        
        rx_signal = zeros(N_rx,N_sc); rx_signal_MRT = zeros(N_rx,N_sc); rx_signal_EGT = zeros(N_rx,N_sc);
        for sc = 1:N_sc
            rx_signal(:,sc) = H(:,:,sc)*tx_signal(:,sc) + noise(:,sc);
            w_mrt = H(:,:,sc)'; w_mrt = w_mrt/norm(w_mrt);
            rx_signal_MRT(:,sc) = H(:,:,sc)*w_mrt.*pilot(sc);
            w_egt = exp(1j*angle(H(:,:,sc)')); w_egt = w_egt/norm(w_egt);
            rx_signal_EGT(:,sc) = H(:,:,sc)*w_egt.*pilot(sc);
        end
        
        detected_before = pskdemod(rx_signal,M,0);
        detected_MRT = pskdemod(rx_signal_MRT,M,0);
        detected_EGT = pskdemod(rx_signal_EGT,M,0);
        
        errors_before = errors_before + sum(detected_before ~= data_symbols);
        errors_MRT = errors_MRT + sum(detected_MRT ~= data_symbols);
        errors_EGT = errors_EGT + sum(detected_EGT ~= data_symbols);
    end
    
    SER_before(idx_snr) = errors_before/(N_sc*N_realizations);
    SER_MRT(idx_snr) = errors_MRT/(N_sc*N_realizations);
    SER_EGT(idx_snr) = errors_EGT/(N_sc*N_realizations);
end

figure('Name','SER vs SNR','NumberTitle','off');
semilogy(SNR_dB_range,SER_before,'r-o','LineWidth',1.5); hold on;
semilogy(SNR_dB_range,SER_MRT,'b-s','LineWidth',1.5);
semilogy(SNR_dB_range,SER_EGT,'m-d','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Symbol Error Rate (SER)');
title('SER vs SNR with CSI-RS-based Beamforming');
legend('Before BF','MRT','EGT');


%% Monte Carlo averaging for SNR comparison
N_iter = 1000;
SNR_noBF_total = 0;
SNR_MRT_total = 0;
SNR_EGT_total = 0;

for iter = 1:N_iter
    % Random Rayleigh channel
    H = (randn(N_rx,N_tx,N_sc)+1j*randn(N_rx,N_tx,N_sc))/sqrt(2);
    noise = (randn(N_rx,N_sc)+1j*randn(N_rx,N_sc))/sqrt(2)*10^(-SNR_demo/20);
    % Total transmit power fixed
    tx_signal = repmat(pilot/sqrt(N_tx), N_tx, 1);  % divide by sqrt(N_tx)

    % Noise-free received signals
    rx_signal_nf = zeros(N_rx,N_sc); 
    rx_signal_MRT_nf = zeros(N_rx,N_sc); 
    rx_signal_EGT_nf = zeros(N_rx,N_sc);

    for sc = 1:N_sc
        % No Beamforming (sum of channels without noise)
        rx_signal_nf(:,sc) = H(:,:,sc)*tx_signal(:,sc);

        % MRT Beamforming
        w_mrt = H(:,:,sc)'; w_mrt = w_mrt/norm(w_mrt);
        rx_signal_MRT_nf(:,sc) = H(:,:,sc)*w_mrt.*pilot(sc);

        % EGT Beamforming
        w_egt = exp(1j*angle(H(:,:,sc)')); w_egt = w_egt/norm(w_egt);
        rx_signal_EGT_nf(:,sc) = H(:,:,sc)*w_egt.*pilot(sc);
    end

    % Compute SNR per iteration
    SNR_noBF_total = SNR_noBF_total + mean(abs(rx_signal_nf).^2)/mean(abs(noise).^2);
    SNR_MRT_total = SNR_MRT_total + mean(abs(rx_signal_MRT_nf).^2)/mean(abs(noise).^2);
    SNR_EGT_total = SNR_EGT_total + mean(abs(rx_signal_EGT_nf).^2)/mean(abs(noise).^2);
end

% Average over iterations and convert to dB
SNR_noBF_dB = 10*log10(SNR_noBF_total/N_iter);
SNR_MRT_dB = 10*log10(SNR_MRT_total/N_iter);
SNR_EGT_dB = 10*log10(SNR_EGT_total/N_iter);


% Display in command line
fprintf('\nAverage SNRs (Monte Carlo, %d iterations):\n', N_iter);
fprintf('MRT Beamforming: %.2f dB\n', SNR_MRT_dB);
fprintf('EGT Beamforming: %.2f dB\n', SNR_EGT_dB);
fprintf('No Beamforming: %.2f dB\n', SNR_noBF_dB);

% Bar graph in correct order MRT > EGT > No BF
figure('Name','Average SNR Comparison','NumberTitle','off');
bar([SNR_MRT_dB SNR_EGT_dB SNR_noBF_dB],'FaceColor',[0.2 0.6 0.8]);
set(gca,'XTickLabel',{'MRT','EGT','No BF'});
ylabel('SNR (dB)');
title('Average SNR Comparison (Monte Carlo)');
grid on;

	
`,

exp11_prach: `

clc; clear; close all;

%% PARAMETERS
N_subcarriers = 72;       % total subcarriers (6 PRBs x 12 SC each)
N_sym = 14;               % OFDM symbols per slot
SC_per_PRB = 12;
M = 4;                    % QPSK modulation
k = log2(M);              % bits per symbol
SNR_dB = 15;              % channel SNR

%% --- Initialize Resource Grid ---
grid = zeros(N_subcarriers, N_sym);        % 0=empty, 1=PRACH, 2=PUCCH, 3=PUSCH
tx_grid_symbols = complex(zeros(N_subcarriers, N_sym));

%% --- Map PRACH ---
prach_prbs = 2:3;  
prach_sc_start = (prach_prbs(1)-1)*SC_per_PRB+1;
prach_sc_end   = prach_prbs(end)*SC_per_PRB;
prach_syms = 1:4;
grid(prach_sc_start:prach_sc_end, prach_syms) = 1;
pr_len = length(prach_sc_start:prach_sc_end)*length(prach_syms);
tx_grid_symbols(grid==1) = exp(1j*2*pi*rand(pr_len,1)); % PRACH pilots

%% --- Map PUCCH ---
pucch_syms = N_sym-1:N_sym; 
pucch_prbs = [1,6]; 
for prb = pucch_prbs
    sc = (prb-1)*SC_per_PRB+1 : prb*SC_per_PRB;
    grid(sc,pucch_syms) = 2;
end
pucch_len = sum(grid(:)==2);
tx_grid_symbols(grid==2) = pskmod(randi([0 1], pucch_len,1),2); % BPSK ACK/NACK

%% --- Map PUSCH (remaining symbols) ---
grid(grid==0) = 3;  
pusch_locs = find(grid==3);

N_pusch = length(pusch_locs);                  
tx_bits_pusch = randi([0 1], N_pusch*k,1); 

bits_reshaped = reshape(tx_bits_pusch,k,[]).'; 
symbols_idx = bits_reshaped(:,1)*2 + bits_reshaped(:,2); 
tx_symbols_pusch = pskmod(symbols_idx, M, pi/4); 

% Map to grid
tx_grid_symbols(pusch_locs) = tx_symbols_pusch;

%% --- Plot 1: Resource Allocation Grid ---
figure('Name','Resource Allocation Grid','NumberTitle','off');
imagesc(1:N_sym,1:N_subcarriers,grid);
colormap([1 0 0; 0 0 1; 0 1 0]); caxis([1 3]);
set(gca,'YDir','normal'); 
xlabel('OFDM Symbols (time)'); ylabel('Subcarrier Index (frequency)');
title('5G NR Uplink Resource Allocation (PRACH, PUCCH, PUSCH)');
grid on; box on;

% Legend
hold on;
h1 = plot(NaN,NaN,'s','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','r'); % PRACH
h2 = plot(NaN,NaN,'s','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','b'); % PUCCH
h3 = plot(NaN,NaN,'s','MarkerFaceColor',[0 1 0],'MarkerEdgeColor','g'); % PUSCH
legend([h1 h2 h3], {'PRACH','PUCCH','PUSCH'}, 'Location','bestoutside');
hold off;

%% --- Plot 2: Transmitted Symbols ---
figure('Name','Transmitted Symbols','NumberTitle','off');
plot(tx_symbols_pusch,'bo'); grid on; axis square;
title('UE Transmitted QPSK Symbols'); xlabel('In-Phase'); ylabel('Quadrature');

%% --- Channel: AWGN ---
rx_grid_symbols = awgn(tx_grid_symbols,SNR_dB,'measured');

%% --- Plot 3: Received Signal Magnitude Grid ---
figure('Name','Received Signal Magnitude Grid','NumberTitle','off');
imagesc(1:N_sym,1:N_subcarriers,abs(rx_grid_symbols));
set(gca,'YDir','normal'); colorbar;
xlabel('OFDM Symbols (time)'); ylabel('Subcarrier Index (frequency)');
title(['Received Signal Magnitude at gNB (SNR = ' num2str(SNR_dB) ' dB)']);
grid on; box on;

%% --- Receiver: Demodulate PUSCH ---
rx_pusch_symbols = rx_grid_symbols(pusch_locs);
rx_idx = pskdemod(rx_pusch_symbols, M, pi/4); 
rx_bits_mat = de2bi(rx_idx, k, 'left-msb'); 
rx_bits_pusch_rx = rx_bits_mat.'; 
rx_bits_pusch_rx = rx_bits_pusch_rx(:);

%% --- BER calculation ---
numErr = sum(tx_bits_pusch ~= rx_bits_pusch_rx);
BER = numErr/length(rx_bits_pusch_rx);
fprintf('Bit Error Rate (BER) at %d dB SNR = %.4f\n',SNR_dB, BER);

%% --- Plot 4: Constellations ---
figure('Name','Constellations','NumberTitle','off');
subplot(1,2,1);
plot(tx_symbols_pusch(1:200),'bo'); grid on; axis square;
title('Transmitted Constellation (QPSK)');
xlabel('In-Phase'); ylabel('Quadrature');

subplot(1,2,2);
plot(rx_pusch_symbols(1:200),'ro'); grid on; axis square;
title(['Received Constellation (SNR = ' num2str(SNR_dB) ' dB)']);
xlabel('In-Phase'); ylabel('Quadrature');

	
`,

exp12_mimo: `
% Initialize system constants
rng(2014);
gc = helperGetDesignSpecsParameters();

% Tunable parameters
tp.txPower = 9;           % watt
tp.txGain = -8;           % dB
tp.mobileRange = 2750;    % m
tp.mobileAngle = 3;       % degrees
tp.interfPower = 1;       % watt
tp.interfGain = -20;      % dB
tp.interfRange = 9000;    % m
tp.interfAngle =   20;    % degrees
tp.numTXElements = 8;
tp.steeringAngle = 0;     % degrees
tp.rxGain = 108.8320 - tp.txGain; % dB

numTx= tp.numTXElements;

[encoder,scrambler,modulatorOFDM,steeringvec,transmitter,...
    radiator,pilots,numDataSymbols,frmSz] = helperMIMOTxSetup(gc,tp);
txBits = randi([0, 1], frmSz,1);
coded = encoder(txBits);
bitsS = scrambler(coded);
tx = qammod(bitsS,gc.modMode,'InputType','bit','UnitAveragePower',true);

ofdm1 = reshape(tx, gc.numCarriers,numDataSymbols);

ofdmData = repmat(ofdm1,[1, 1, numTx]);
txOFDM = modulatorOFDM(ofdmData, pilots);
%scale
txOFDM = txOFDM * ...
    (gc.FFTLength/sqrt(gc.FFTLength-sum(gc.NumGuardBandCarriers)-1));

% Amplify to achieve peak TX power for each channel
for n = 1:numTx
    txOFDM(:,n) = transmitter(txOFDM(:,n));
end

radiator.CombineRadiatedSignals = false;

wR = steeringvec(gc.fc,[-tp.mobileAngle;0]);

wT = steeringvec(gc.fc,[tp.steeringAngle;0]);
weight = wT.* wR;

txOFDM = radiator(txOFDM,repmat([tp.mobileAngle;0],1,numTx),conj(weight));

[channel,interferenceTransmitter,toRxAng,spLoss] = ...
    helperMIMOEnvSetup(gc,tp);
[sigFade, chPathG] =  channel(txOFDM);
sigLoss = sigFade/sqrt(db2pow(spLoss(1)));

% Generate interference and apply gain and propagation loss
numBits = size(sigFade,1);
interfSymbols = wgn(numBits,1,1,'linear','complex');
interfSymbols = interferenceTransmitter(interfSymbols);
interfLoss = interfSymbols/sqrt(db2pow(spLoss(2)));


[collector,receiver,demodulatorOFDM,descrambler,decoder] = ...
    helperMIMORxSetup(gc,tp,numDataSymbols);

rxSig = collector([sigLoss interfLoss],toRxAng);

% Front-end amplifier gain and thermal noise
rxSig = receiver(rxSig);

rxOFDM = rxSig * ...
    (sqrt(gc.FFTLength-sum(gc.NumGuardBandCarriers)-1)) / (gc.FFTLength);

% OFDM Demodulation
rxOFDM = demodulatorOFDM(rxOFDM);

% Channel estimation
hD = helperIdealChannelEstimation(gc,  numDataSymbols, chPathG);

% Equalization
rxEq = helperEqualizer(rxOFDM, hD, numTx);

% Collapse OFDM matrix
rxSymbs = rxEq(:);

rxBitsS = qamdemod(rxSymbs,gc.modMode,'UnitAveragePower',true,...
    'OutputType','bit');
rxCoded = descrambler(rxBitsS);
rxDeCoded = decoder(rxCoded);
rxBits = rxDeCoded(1:frmSz);

ber = comm.ErrorRate;
measures = ber(txBits, rxBits);
fprintf('BER = %.2f%%; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1)*100,measures(3), measures(2));

constdiag = comm.ConstellationDiagram('SamplesPerSymbol', 1,...
    'ReferenceConstellation', [], 'ColorFading',true,...
    'Position', gc.constPlotPosition);
% Display received constellation
constdiag(rxSymbs);


tp.steeringAngle = tp.mobileAngle;

% Steer the transmitter main lobe
wT = steeringvec(gc.fc,[tp.steeringAngle;0]);

[txBits, rxBits,rxSymbs] = helperRerunMIMOBeamformingExample(gc,tp,wT);

reset(ber);
measures = ber(txBits, rxBits);
fprintf('BER = %.2f%%; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1)*100,measures(3), measures(2));

constdiag(rxSymbs);


% analog phase shifter with quantization effect
release(steeringvec);
steeringvec.NumPhaseShifterBits = 4;
wTq = steeringvec(gc.fc,[tp.steeringAngle;0]);

[txBits, rxBits,rxSymbs] = helperRerunMIMOBeamformingExample(gc,tp,wTq);

reset(ber);
measures = ber(txBits, rxBits);
fprintf('BER = %.2f%%; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1)*100,measures(3), measures(2));

constdiag = comm.ConstellationDiagram('SamplesPerSymbol', 1,...
    'ReferenceConstellation', [], 'ColorFading',true,...
    'Position', gc.constPlotPosition);
constdiag(rxSymbs);
`,

evolution_modulation: `
% BER vs SNR Comparison for 1G to 5G using Monte Carlo Simulation
% Representative modulation schemes:
% 1G - Analog (approximated with BPSK)
% 2G - GSM (approximated with BPSK/GMSK)
% 3G - QPSK
% 4G - 16QAM
% 5G - 64QAM

clc; clear; close all;

% Simulation parameters
Nsym = 1e5;                   % Symbols per trial
Ntrials = 100;                 % Monte Carlo trials
SNR_dB = 0:2:30;              % SNR range
BER = zeros(5,length(SNR_dB));

% Modulation schemes
modSchemes = {'BPSK','BPSK','QPSK','16QAM','64QAM'};
titles = {'1G (Analog~BPSK)','2G (GSM~BPSK)','3G (QPSK)','4G (16QAM)','5G (64QAM)'};

for gen = 1:5
    % Select constellation size
    switch modSchemes{gen}
        case 'BPSK'
            M = 2;
        case 'QPSK'
            M = 4;
        case '16QAM'
            M = 16;
        case '64QAM'
            M = 64;
    end
    
    for i = 1:length(SNR_dB)
        errCount = 0;
        bitCount = 0;
        
        for t = 1:Ntrials
            % Random data
            data = randi([0 M-1],Nsym,1);
            
            % Modulation
            txSig = qammod(data,M,'UnitAveragePower',true);
            
            % AWGN channel
            rxSig = awgn(txSig,SNR_dB(i),'measured');
            
            % Demodulation
            rxData = qamdemod(rxSig,M,'UnitAveragePower',true);
            
            % Error counting
            errCount = errCount + sum(rxData~=data);
            bitCount = bitCount + length(data)*log2(M);
        end
        
        % Monte Carlo BER
        BER(gen,i) = errCount/bitCount;
    end
end

% Plot results
figure;
semilogy(SNR_dB,BER(1,:),'-o','LineWidth',1.5); hold on;
semilogy(SNR_dB,BER(2,:),'-s','LineWidth',1.5);
semilogy(SNR_dB,BER(3,:),'-d','LineWidth',1.5);
semilogy(SNR_dB,BER(4,:),'-^','LineWidth',1.5);
semilogy(SNR_dB,BER(5,:),'-v','LineWidth',1.5);
grid on; axis([0 30 1e-5 1]);
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title('BER vs SNR for 1G → 5G');
legend(titles,'Location','southwest');

	
`

};

function copyCode(key) {
  navigator.clipboard.writeText(codes[key]);
  alert("Experiment code copied!");
}
