const codes = {

exp7_ofdm: `
PASTE CODE FROM:
ofdm_waveform with diff modulation_exp7
`,

exp8_cu_du: `
PASTE CODE FROM:
cu_du_split_5g_exp8
`,

exp9_numerology: `
PASTE CODE FROM:
exp9_diff5g_numerology
`,

exp10_csi_rs: `
PASTE CODE FROM:
exp10_csi_rs
`,

exp11_prach: `
PASTE CODE FROM:
exp11_prach
`,

exp12_mimo: `
PASTE CODE FROM:
exp12_mimo
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
