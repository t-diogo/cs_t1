close all;
clear all

%% ADC parameters
Vdd  = 1;
Vcm  = Vdd/2;
Vref = 2*Vdd;
num_stages = 3;
gain  = [4.0,4.0,4.0];
nbits = [4,4,6];

nbits_adc = sum(nbits) - length(nbits) + 1;

%% Component errors - standard deviation
sigma_caps = 0.002;              % capacitor mismatch error - 0.2% of the unit value
sigma_gain = 0.001;              % gain error - 0.1% of the unit value
sigma_comp = 8e-3 / sqrt(2) ;    % comparator - offset voltage error - 8mV RMS

% ideal case
% sigma_caps = 0.0;                  
% sigma_gain = 0.0;                   
% sigma_comp = 0.0;                   

% Cells to store components parameters
Gain = cell(1,num_stages);
Voff = cell(1,num_stages);
Cup  = cell(1,num_stages);
Cun  = cell(1,num_stages);

Cup_total = cell(1,num_stages);
Cun_total = cell(1,num_stages);
Cbinp = cell(1, num_stages);
Cbinm = cell(1, num_stages);


%% Monte Carlo Simulation
mc_simulations = 10000;

% Accuracy and linearity variables
codes   = zeros(1, mc_simulations);
lin_inl = zeros(1, mc_simulations);
lin_dnl = zeros(1, mc_simulations);


%% Input signal - sweep transfer function
deltaVin = linspace(-Vdd, Vdd, 2^(nbits_adc + 3) + 1);  % 3 => 2^3 = 8 points per Vlsb
np = max(size(deltaVin));

tic
for mc = 1:mc_simulations
    disp(mc);

    % Generate components for each stage
    for i = 1:num_stages

        Voff{i} = randn*sigma_comp;
        Gain{i} = gain(i) + randn*sigma_gain;
        Cup {i} = 1 + randn(1,2^nbits(i))*sigma_caps;
        Cun {i} = 1 + randn(1,2^nbits(i))*sigma_caps;

        % Total capacitance of the capacitor array
        Cup_total{i} = sum(Cup{i});
        Cun_total{i} = sum(Cun{i});

        % Generate equivalent capacitors
        Cbinp{i} = zeros(1,nbits(i));
        Cbinm{i} = zeros(1,nbits(i));
        for j = 0:nbits(i)-1
            Cbinp{i}(j+1) = sum (Cup{i}((2^j):(2^(j+1)-1)));
            Cbinm{i}(j+1) = sum (Cun{i}((2^j):(2^(j+1)-1)));
        end
        % Flip to get the right order to simplify further access
        Cbinp{i} = fliplr(Cbinp{i});
        Cbinm{i} = fliplr(Cbinm{i});
    end

    % Simulate pipeline
    [dout] = simulate_pipeline(num_stages, nbits, Vref, Vcm, deltaVin, np, Cbinp, Cbinm, Cup_total, Cun_total, Gain, Voff);

    % Find Voltage transitions
    Vt = deltaVin(find(dout(2:end)-dout(1:end-1)) + 1);
    ncodes = max(size(Vt));

    % INL e DNL
    Vlsb_real = (Vt(end)-Vt(1)) / (ncodes-1);
    inl = (Vt-(0:ncodes-1) * Vlsb_real - Vt(1)) / Vlsb_real;
    dnl = ((Vt(2:end) - Vt(1:end-1)) / Vlsb_real) - 1;

    codes(mc) = ncodes;
    lin_inl(mc) = nbits_adc + log(max(inl)-min(inl)) / log(2);
    lin_dnl(mc) = nbits_adc + log(max(dnl)-min(dnl)) / log(2);
end
toc

%% Plots

if mc_simulations == 1
    % Transfer Function
    figure(1)
    plot(deltaVin, dout, 'LineWidth', 1.5)
    grid on
    title(['Transfer Function with ' num2str(ncodes) ' codes'], 'FontSize', 14)
    xlabel('\DeltaVin (Input Voltage)', 'FontSize', 12)
    ylabel('Dout (Digital Output Code)', 'FontSize', 12)
    set(gca, 'FontSize', 12)

    figure(2)
    subplot(2,1,1)
    plot(inl, 'LineWidth', 1.5)
    grid on
    title('INL (Integral Non-Linearity)', 'FontSize', 14)
    xlabel('Code', 'FontSize', 12)
    ylabel('INL (LSB)', 'FontSize', 12)
    set(gca, 'FontSize', 12)

    subplot(2,1,2)
    plot(dnl, 'LineWidth', 1.5)
    grid on
    title('DNL (Differential Non-Linearity)', 'FontSize', 14)
    xlabel('Code', 'FontSize', 12)
    ylabel('DNL (LSB)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
end
    figure(3)
    subplot(2,1,1)
    hist(lin_inl, 100)
    grid on
    title('Linearity INL', 'FontSize', 14)
    xlabel('INL (LSB)', 'FontSize', 12)
    ylabel('Frequency', 'FontSize', 12)
    set(gca, 'FontSize', 12)

    subplot(2,1,2)
    hist(lin_dnl, 100)
    grid on
    title('Linearity DNL', 'FontSize', 14)
    xlabel('DNL (LSB)', 'FontSize', 12)
    ylabel('Frequency', 'FontSize', 12)
    set(gca, 'FontSize', 12)

    figure(4)
    hist(2^nbits_adc - codes, 100)
    grid on
    title('Number of Missing Codes', 'FontSize', 14)
    xlabel('Missing Codes', 'FontSize', 12)
    ylabel('Frequency', 'FontSize', 12)
    set(gca, 'FontSize', 12)
