clear;clc;
% ========================================================================
%               Question 1 : comaparing between FFT and DFT 
% ========================================================================

fprintf('Running Part 2.1: FFT/DFT Simulation...\n');
L= 4096;
X_N = rand(1,L);
% count time taken for DFT 
tic;
X_K_DFT= DFT_function(X_N);
taken_time_DFT=toc;
% count timw taken for fft 
tic;
X_K_fft= fft(X_N);
taken_time_fft=toc;
% print time taken for DFT and fft 
fprintf("DFT excution time : %.6f seconds\n",taken_time_DFT);
fprintf("fft excution time : %.6f seconds\n",taken_time_fft);

% ---------------- Bar Graph ----------------
times = [taken_time_DFT, taken_time_fft];
x=["DFT" "FFT"];
figure;
b = bar(x, times);   % <-- store bar handle
grid on;
ylabel('Execution Time (seconds)');
xlabel('type of transformer');
title('Execution Time Comparison: DFT vs FFT');
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);

text(xtips, ytips, labels, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom');

function X_k = DFT_function(X_n)
N=length(X_n);
X_k=zeros(1,N);
for k = 1:N-1
    for n = 1:N-1
        X_k(k+1) = X_k(k+1) + X_n(n+1) * exp((-1i*2*pi*n)/N); 
    end
end
end

%% =================================================================================================
%   Question 2 : Bit-error rate performance for BPSK and 16-QAM over Rayleigh flat fading channel
% =================================================================================================
% ==============================
% Single-File Rayleigh Fading Simulator
% Supports: BPSK, QPSK, 16-QAM
% Repetition factor fixed at 5
% Computes BER: normal vs repetition
% ==============================

clear; clc; rng(0);

% ------------------------------
% PART 1: User Input / Parameters
% ------------------------------

modulation_type = input('Select modulation type (BPSK / QPSK / 16QAM): ','s');
repetition_factor = 5;          % Fixed repetition factor
SNR_db = -3:0.5:10;             % SNR range (dB)
Eb = 1;                         % Energy per bit

% Set number of bits based on modulation
switch upper(modulation_type)
    case 'BPSK'
        Num_bits = 5e6; 
        bits_per_symbol = 1;
    case 'QPSK'
        Num_bits = 5e5; 
        bits_per_symbol = 2;
    case '16QAM'
        Num_bits = 1.2e6; 
        bits_per_symbol = 4;
    otherwise
        error('Unsupported modulation type');
end

fprintf('Running %s simulation with repetition factor = %d...\n', modulation_type, repetition_factor);

% Generate original bits (column vector)
Tx_bits = randi([0 1], Num_bits, 1);

% Preallocate BER results
BER_Normal = zeros(size(SNR_db));
BER_Rep5   = zeros(size(SNR_db));

% ------------------------------
% PART 2: Main SNR Loop
% ------------------------------

for idx = 1:length(SNR_db)
    
    % ----- Noise power -----
    EbNo_linear = 10^(SNR_db(idx)/10);
    N0 = Eb / (EbNo_linear * bits_per_symbol); 
    
    % ----- BER without repetition (Normal) -----
    tx_symbols_normal = map_bits(Tx_bits, modulation_type);
    [rx_symbols_normal, h_normal] = rayleigh_channel(tx_symbols_normal, N0);
    z_normal = equalize(rx_symbols_normal, h_normal);
    rx_bits_normal = demap_symbols(z_normal, modulation_type);
    BER_Normal(idx) = compute_ber(Tx_bits, rx_bits_normal);

    % ----- BER with repetition (factor 5) -----
    Tx_bits_rep = encode_bits(Tx_bits, repetition_factor);
    tx_symbols_rep = map_bits(Tx_bits_rep, modulation_type);
    [rx_symbols_rep, h_rep] = rayleigh_channel(tx_symbols_rep, N0);
    z_rep = equalize(rx_symbols_rep, h_rep);
    rx_bits_rep = demap_symbols(z_rep, modulation_type);
    Rx_bits = decode_bits(rx_bits_rep, repetition_factor);
    BER_Rep5(idx) = compute_ber(Tx_bits, Rx_bits);

end

% ------------------------------
% PART 3: Plot Results
% ------------------------------
figure;
semilogy(SNR_db, BER_Normal, 'b-o','LineWidth',1.5); hold on;
semilogy(SNR_db, BER_Rep5, 'r-s','LineWidth',1.5); grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
legend('Normal','Repetition 5');
title(sprintf('%s Performance over Rayleigh Channel', modulation_type));

% PART 4: Local Functions
function bits_rep = encode_bits(bits, repetition_factor)
    % Repeat each bit 'repetition_factor' times
    bits_rep = repelem(bits, repetition_factor);
 end

function bits_dec = decode_bits(bits_rx, repetition_factor)
    % Majority vote decoding
    bits_matrix = reshape(bits_rx, repetition_factor, []);
    bits_dec = sum(bits_matrix,1) >= ceil(repetition_factor/2);
    bits_dec = bits_dec.'; % column vector
  end

function symbols = map_bits(bits, modulation)
    % Maps bits to symbols based on modulation
    switch upper(modulation)
        case 'BPSK'
            symbols = 2*bits - 1; % 0->-1, 1->1
        case 'QPSK'
            b1 = bits(1:2:end-1);
            b2 = bits(2:2:end);
            symbols = complex(2*b1-1, 2*b2-1);
        case '16QAM'
            % 16-QAM Gray Mapping
            bits = bits(:);
            b = reshape(bits,4,[]).';
            gray = zeros(size(b));
            gray(:,1) = b(:,1);
            gray(:,2) = xor(b(:,1),b(:,2));
            gray(:,3) = xor(b(:,2),b(:,3));
            gray(:,4) = xor(b(:,3),b(:,4));
            idx = bi2de(gray,'left-msb')+1;
            % 16-QAM constellation (normalized)
            I_levels = [-3 -1 3 1]; Q_levels = [-3 -1 3 1];
            const = zeros(16,1);
            for k=0:15
                bin_k = de2bi(k,4,'left-msb');
                bin_nat = zeros(size(bin_k));
                bin_nat(1) = bin_k(1);
                bin_nat(2) = xor(bin_nat(1), bin_k(2));
                bin_nat(3) = xor(bin_nat(2), bin_k(3));
                bin_nat(4) = xor(bin_nat(3), bin_k(4));
                i_idx = bi2de(bin_nat(1:2),'left-msb')+1;
                q_idx = bi2de(bin_nat(3:4),'left-msb')+1;
                const(k+1) = I_levels(i_idx)+1j*Q_levels(q_idx);
            end
            const = const/sqrt(mean(abs(const).^2));
            symbols = const(idx).';
        otherwise
            error('Unsupported modulation type');
    end
end

function bits = demap_symbols(symbols, modulation)
    % Demap symbols back to bits
    switch upper(modulation)
        case 'BPSK'
            bits = real(symbols) >= 0;
        case 'QPSK'
            bits = zeros(2*length(symbols),1);
            bits(1:2:end-1) = real(symbols)>=0;
            bits(2:2:end) = imag(symbols)>=0;
       case '16QAM'
    symbols = symbols(:); % ensure column vector
    % Define constellation
    I_levels = [-3 -1 3 1]; Q_levels = [-3 -1 3 1];
    const = zeros(16,1);
    for k=0:15
        bin_k = de2bi(k,4,'left-msb');
        bin_nat = zeros(size(bin_k));
        bin_nat(1) = bin_k(1);
        bin_nat(2) = xor(bin_nat(1), bin_k(2));
        bin_nat(3) = xor(bin_nat(2), bin_k(3));
        bin_nat(4) = xor(bin_nat(3), bin_k(4));
        i_idx = bi2de(bin_nat(1:2),'left-msb')+1;
        q_idx = bi2de(bin_nat(3:4),'left-msb')+1;
        const(k+1) = I_levels(i_idx)+1j*Q_levels(q_idx);
    end
    const = const/sqrt(mean(abs(const).^2));
    % Nearest neighbor detection
    [~, idx] = min(abs(symbols - const.'),[],2); % now safe
    gray = de2bi(idx-1,4,'left-msb');
    % Reverse Gray
    bits_mat = zeros(size(gray));
    bits_mat(:,1) = gray(:,1);
    bits_mat(:,2) = xor(bits_mat(:,1),gray(:,2));
    bits_mat(:,3) = xor(bits_mat(:,2),gray(:,3));
    bits_mat(:,4) = xor(bits_mat(:,3),gray(:,4));
    bits = reshape(bits_mat.',[],1);

        otherwise
            error('Unsupported modulation type');
    end
end

function [rx, h] = rayleigh_channel(tx, N0)
    tx = tx(:); % ensure column vector
    h = (randn(length(tx),1) + 1j*randn(length(tx),1))/sqrt(2); % Rayleigh
    n = sqrt(N0/2)*(randn(length(tx),1) + 1j*randn(length(tx),1)); % AWGN
    rx = h .* tx + n;
end
 
function z = equalize(rx,h)
    % Channel equalization
    z = rx ./ h;
end

function BER = compute_ber(tx,rx)
    % Compute Bit Error Rate
    BER = sum(tx~=rx)/length(tx);
end




%%
% ========================================================================
%               Question 3 : OFDM System Simulation
% ========================================================================
clear; clc;
% ===================== System Parameters =====================
cp = 16;                        % Cyclic prefix length
rate = 1/5;                     % Coding rate
N_fft = 256;                    % FFT size
Eb_N0_dB = 0:2:20;              % Eb/N0 range in dB
N_OFDM = 10000;                 % Number of OFDM symbols
modulation ='QPSK';             % Modulation type: (BPSK, QPSK, 16QAM)
energy_type = 'Same_Per_Trans'; % Energy normalization type : (Same_Per_Trans, Same_Per_Info)
channel_type = 'Flat';          % Channel type: Flat or Frequency Selective (Flat, FreqSel)

% ===================== Modulation Order =====================
if isequal(modulation, 'BPSK')
    m = 1;
elseif isequal(modulation, 'QPSK')  
    m = 2;
elseif isequal(modulation, '16QAM')
    m = 4;
end    

% Loop for uncoded (N_code=1) and coded (N_code=1/rate) systems
for N_code = [1 1/rate]
    for i = 1:length(Eb_N0_dB)
        symbol_error = [];
        total_bits = [];
        for j = 1:N_OFDM

        % ===================== Transmitter =====================
            
            % Generate random input bits
            data = generate_data(m, N_fft, N_code);

            % Encoder
            encoder_output = Encoder(data, m, N_fft, N_code);
        
            % Interleaver
            interleaver_output = Interleaver(modulation, encoder_output);
        
            
            % Mapper
            mapper_output = Mapper(interleaver_output, modulation);
            
            % ifft
            ifft_output = ifft(mapper_output, N_fft);

            % cyclic extension
            cyclic_output = [ifft_output(end-cp+1:end), ifft_output];
            
            % channel
            [channel_output, H] = Channel(cyclic_output, m, Eb_N0_dB(i), N_code, N_fft, cp, channel_type, energy_type);
        
        
        % ===================== Receiver =====================
        
            % Remove cyclic extension
            remove_cyclic_output = channel_output(cp+1:end);
            
            % fft
            fft_output = fft(remove_cyclic_output, N_fft);

            % channel equalization
            channel_equ = fft_output ./ H;
    
            % De-Mapper
            demapper_output = Demapper(channel_equ, modulation);
            
            % De-Interleaver
            deinterleaver_output = DeInterleaver(modulation, demapper_output);
        
            % Decoder
            decoder_output = Decoder(deinterleaver_output, m, N_fft, N_code);
            
            % Compute BER for curren ofdm symbol
            symbol_error(j) = sum(decoder_output ~= data)/length(data);
        
        end
        % Average BER for current Eb/N0
        BER(i) = sum(symbol_error) / N_OFDM;
    end
    
    if N_code == 1
        BER_uncoded = BER;
    else
        BER_coded = BER;
    end
end
%
% ============================== Plot ==============================
figure('Name', [modulation ', ' channel_type ' Channel']);
semilogy(Eb_N0_dB , BER_uncoded ,'b-', 'LineWidth', 1.25) ;  
hold on
semilogy(Eb_N0_dB , BER_coded ,'b--', 'LineWidth', 1.25);  
hold off
grid on
xlabel('Eb/No (dB)');  
ylabel('BER');
ylim([1e-5 1])
legend([modulation,'(Uncoded)'], [modulation,'(coded 1/5)']) ;

title([modulation ', ' channel_type ' Channel'])

% ===================== Functions =====================
function R = generate_data(m, N_fft, N_code)
    R = randi([0 1], 1, floor(N_fft/N_code)*m);
end


function R = Encoder(data, m, N_fft, N_code)
    encoded_data = repelem(data, N_code);
    R = [encoded_data, zeros(1, (m*N_fft)-length(encoded_data))];
end

function R = Decoder(data, m, N_fft, N_code)
    padding_size = length(data) - floor(N_fft/N_code)*m*N_code;
    data_noPad = data(1:end-padding_size);
    
    % Majority vote decoder
    R = sum(reshape(data_noPad, N_code, []), 1) >= ceil(N_code/2); 
end

function [Channel_Output, H] = Channel(data, m, Eb_N0_dB, N_code, N_fft, cp, channeltype, energy_type)
    % Noise
    segma = sqrt((1/2)./10^(Eb_N0_dB/10));
    Eb_actual = (sum(abs(data).^2)/length(data))/m;
    noise = sqrt(Eb_actual) * segma  * (randn(size(data)) + 1j*randn(size(data)));

    if isequal(energy_type, 'Same_Per_Info')
        noise = sqrt(N_code) * noise;
    end

    if isequal(channeltype, 'Flat')
        % Flat fading channel
        H = sqrt(1/2)*(randn(1) + 1j*randn(1));
        Channel_Output = H*data + noise;
    else
        % Frequency selective channel
        H = sqrt(1/2)*(randn(1, N_fft) + 1j*randn(1, N_fft));

        % Remove cyclic extension
        data_cp_removed = data(cp+1:end);

        % Frequency domain
        data_freq = fft(data_cp_removed, N_fft);
        
        % Applay Freqency Selctive Channel
        ch_out_freq = data_freq .* H; 
        
        % Time domain
        ch_out = ifft(ch_out_freq, N_fft);
        
        % Add cyclic extension
        ch_out = [ch_out(end-cp+1:end), ch_out];

        %Add Noise
        Channel_Output = ch_out + noise;
    end
end



function R = Interleaver(mod, data)
    if isequal(mod, 'BPSK')
       reshape_tx_data = reshape(data, 8, 32)';
    elseif isequal(mod, 'QPSK')
       reshape_tx_data = reshape(data, 16, 32)';
    elseif isequal(mod, '16QAM')
       reshape_tx_data = reshape(data, 32, 32)';
    end
    R = reshape_tx_data(:)';
end

function R = DeInterleaver(mod, data)
    if isequal(mod, 'BPSK')
       reshape_rx_data = reshape(data, 32, 8);       
    elseif isequal(mod, 'QPSK')
       reshape_rx_data = reshape(data, 32, 16);
    elseif isequal(mod, '16QAM')
       reshape_rx_data = reshape(data, 32, 32);
    end
    reshape_rx_data = reshape_rx_data';
    R = reshape_rx_data(:)';
end



function R = Mapper(data, mod)
    if isequal(mod, 'BPSK')
            R = data*2-1;
    elseif isequal(mod, 'QPSK')
        data_reshaped = reshape(data, 2, []);
        % 00 -> -1 -1 => -1 - j
        % 01 -> -1 +1 => -1 + j
        % 10 -> +1 -1 => +1 - j
        % 11 -> +1 +1 => +1 + j
        I = 2*data_reshaped(1,:) - 1;      % MSB
        Q = 2*data_reshaped(2,:) - 1;      % LSB
        data_mapped = I + 1j * Q;
        R = data_mapped;
    elseif isequal(mod, '16QAM')  
        R = mod16(data);
    end
end

function R = Demapper(data, mod)
    if isequal(mod, 'BPSK')
        % if data > 0 it represents 1 else represents 0 
        R = data > 0;
    elseif isequal(mod, 'QPSK')  
        demapped_I = real(data) > 0;
        demapped_Q = imag(data) > 0;
        R= reshape([demapped_I; demapped_Q], 1, []);
    elseif isequal(mod, '16QAM')  
        R = demod16(data);
    end
end



function [rxsig] = mod16(txbits) 
    psk16mod=[1+j*1 3+j*1 1+j*3 3+j*3 1-j*1 3-j*1 1-j*3 3-j*3 -1+j*1 -3+j*1 -1+j*3 -3+j*3 -1-j*1 -3-j*1 -1-j*3 -3-j*3]; 
    sigham=txbits; 
    m=4; 
    sigqam16=reshape(sigham,m,length(sigham)/m); 
    rxsig=(psk16mod(bi2de(sigqam16')+1)); 
end


function [rxbits]=demod16(rxsig) 
    m=4; 
    psk16demod=[15 14 6 7 13 12 4 5 9 8 0 1 11 10 2 3]; 
    rxsig(find(real(rxsig)>3))=3+j*imag(rxsig(find(real(rxsig)>3))); 
    rxsig(find(imag(rxsig)>3))=real(rxsig(find(imag(rxsig)>3)))+j*3; 
    rxsig(find(real(rxsig)<-3))=-3+j*imag(rxsig(find(real(rxsig)<-3))); 
    rxsig(find(imag(rxsig)<-3))=real(rxsig(find(imag(rxsig)<-3)))-j*3; 
    rxdemod=round(real((rxsig+3+j*3)/2))+j*round(imag((rxsig+3+j*3)/2)); 
    rxdebi=real(rxdemod)+4*(imag(rxdemod)); 
    sigbits=de2bi(psk16demod(rxdebi+1));
    rxbits= reshape(sigbits',1,size(sigbits, 1)*m);
end


