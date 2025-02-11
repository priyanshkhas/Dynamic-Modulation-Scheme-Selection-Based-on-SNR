% Main script - save this as 'adaptive_modulation.m'

% Parameters
num_symbols = 1000;  % Number of symbols to transmit
snr_range = -5:2:30;  % SNR range in dB
target_ber = 1e-3;    % Target Bit Error Rate

% Generate random binary data
data = randi([0 1], 1, num_symbols * 4);  % Generate enough bits for 16-QAM

% Initialize arrays to store results
ber_bpsk = zeros(1, length(snr_range));
ber_qpsk = zeros(1, length(snr_range));
ber_16qam = zeros(1, length(snr_range));
throughput = zeros(1, length(snr_range));
selected_scheme = cell(1, length(snr_range));

% Process for each SNR value
for idx = 1:length(snr_range)
    snr = snr_range(idx);
    
    % BPSK
    bpsk_symbols = bpsk_modulate(data(1:num_symbols));
    noise_bpsk = (randn(size(bpsk_symbols)) + 1j*randn(size(bpsk_symbols))) / sqrt(2);
    received_bpsk = bpsk_symbols + noise_bpsk * 10^(-snr/20);
    demod_bpsk = bpsk_demodulate(received_bpsk);
    ber_bpsk(idx) = sum(data(1:num_symbols) ~= demod_bpsk) / num_symbols;
    
    % QPSK
    qpsk_symbols = qpsk_modulate(data(1:num_symbols*2));
    noise_qpsk = (randn(size(qpsk_symbols)) + 1j*randn(size(qpsk_symbols))) / sqrt(2);
    received_qpsk = qpsk_symbols + noise_qpsk * 10^(-snr/20);
    demod_qpsk = qpsk_demodulate(received_qpsk);
    ber_qpsk(idx) = sum(data(1:num_symbols*2) ~= demod_qpsk) / (num_symbols * 2);
    
    % 16-QAM
    qam_symbols = qam16_modulate(data(1:num_symbols*4));
    noise_qam = (randn(size(qam_symbols)) + 1j*randn(size(qam_symbols))) / sqrt(2);
    received_qam = qam_symbols + noise_qam * 10^(-snr/20);
    demod_qam = qam16_demodulate(received_qam);
    ber_16qam(idx) = sum(data(1:num_symbols*4) ~= demod_qam) / (num_symbols * 4);
    
    % Select modulation scheme based on BER
    if ber_bpsk(idx) <= target_ber
        selected_scheme{idx} = 'BPSK';
        throughput(idx) = 1;
    elseif ber_qpsk(idx) <= target_ber
        selected_scheme{idx} = 'QPSK';
        throughput(idx) = 2;
    elseif ber_16qam(idx) <= target_ber
        selected_scheme{idx} = '16-QAM';
        throughput(idx) = 4;
    else
        selected_scheme{idx} = 'No transmission';
        throughput(idx) = 0;
    end
end

% Plotting results
figure('Position', [100, 100, 1200, 800]);

% Plot 1: BER vs SNR
subplot(2,2,1);
semilogy(snr_range, ber_bpsk, 'b-*', 'LineWidth', 2);
hold on;
semilogy(snr_range, ber_qpsk, 'r-o', 'LineWidth', 2);
semilogy(snr_range, ber_16qam, 'g-s', 'LineWidth', 2);
semilogy(snr_range, target_ber*ones(size(snr_range)), 'k--', 'LineWidth', 1);
grid on;
legend('BPSK', 'QPSK', '16-QAM', 'Target BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER vs SNR for Different Modulation Schemes');

% Plot 2: Selected Modulation Scheme
subplot(2,2,2);
scheme_numbers = zeros(size(snr_range));
for i = 1:length(snr_range)
    switch selected_scheme{i}
        case 'BPSK'
            scheme_numbers(i) = 1;
        case 'QPSK'
            scheme_numbers(i) = 2;
        case '16-QAM'
            scheme_numbers(i) = 3;
        otherwise
            scheme_numbers(i) = 0;
    end
end
plot(snr_range, scheme_numbers, 'b-o', 'LineWidth', 2);
grid on;
yticks([0 1 2 3]);
yticklabels({'No Tx', 'BPSK', 'QPSK', '16-QAM'});
xlabel('SNR (dB)');
ylabel('Selected Scheme');
title('Adaptive Modulation Scheme Selection');

% Plot 3: Throughput
subplot(2,2,3);
plot(snr_range, throughput, 'r-*', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bits per Symbol');
title('System Throughput');

% Plot 4: Constellation Diagram Example (for SNR = 20dB)
subplot(2,2,4);
snr_idx = find(snr_range == 20);
if ~isempty(snr_idx)
    switch selected_scheme{snr_idx}
        case 'BPSK'
            scatter(real(bpsk_symbols), imag(bpsk_symbols), 'b*');
        case 'QPSK'
            scatter(real(qpsk_symbols), imag(qpsk_symbols), 'r*');
        case '16-QAM'
            scatter(real(qam_symbols), imag(qam_symbols), 'g*');
    end
    grid on;
    axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(['Constellation Diagram at SNR = 20dB (', selected_scheme{snr_idx}, ')']);
end

% Adjust layout
sgtitle('Adaptive Modulation System Analysis');

% Function definitions
function mod_symbols = bpsk_modulate(bits)
    mod_symbols = 2 * bits - 1;  % 0 -> -1, 1 -> 1
end

function demod_bits = bpsk_demodulate(symbols)
    demod_bits = real(symbols) > 0;
end

function mod_symbols = qpsk_modulate(bits)
    % Group bits in pairs and map to QPSK symbols
    I = 2 * bits(1:2:end) - 1;
    Q = 2 * bits(2:2:end) - 1;
    mod_symbols = (I + 1j*Q) / sqrt(2);
end

function demod_bits = qpsk_demodulate(symbols)
    % Convert QPSK symbols back to bits
    I = real(symbols) > 0;
    Q = imag(symbols) > 0;
    demod_bits = reshape([I; Q], 1, []);
end

function mod_symbols = qam16_modulate(bits)
    % Custom implementation of binary to decimal conversion
    bits_reshaped = reshape(bits, 4, [])';
    decimal = zeros(size(bits_reshaped, 1), 1);
    
    % Manual binary to decimal conversion
    for i = 1:size(bits_reshaped, 1)
        decimal(i) = bits_reshaped(i,1)*8 + bits_reshaped(i,2)*4 + ...
                     bits_reshaped(i,3)*2 + bits_reshaped(i,4);
    end
    
    % 16-QAM constellation mapping
    constellation = [-3-3j, -3-1j, -3+3j, -3+1j, ...
                    -1-3j, -1-1j, -1+3j, -1+1j, ...
                     3-3j,  3-1j,  3+3j,  3+1j, ...
                     1-3j,  1-1j,  1+3j,  1+1j] / sqrt(10);
    
    mod_symbols = constellation(decimal + 1);
end

function demod_bits = qam16_demodulate(symbols)
    % 16-QAM constellation
    constellation = [-3-3j, -3-1j, -3+3j, -3+1j, ...
                    -1-3j, -1-1j, -1+3j, -1+1j, ...
                     3-3j,  3-1j,  3+3j,  3+1j, ...
                     1-3j,  1-1j,  1+3j,  1+1j] / sqrt(10);
    
    % Find nearest constellation point
    [~, idx] = min(abs(repmat(symbols(:), 1, 16) - repmat(constellation, length(symbols), 1)), [], 2);
    
    % Manual decimal to binary conversion
    idx = idx - 1;  % Zero-based index
    demod_bits = zeros(1, length(idx) * 4);
    
    for i = 1:length(idx)
        val = idx(i);
        bit_pos = (i-1)*4 + 1;
        % Convert to 4 bits
        demod_bits(bit_pos:bit_pos+3) = [floor(val/8), floor(mod(val,8)/4), ...
                                        floor(mod(val,4)/2), mod(val,2)];
    end
end