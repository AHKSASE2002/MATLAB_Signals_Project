% Define constants
fs = 1000; % Sampling frequency (Hz)
t = 0:1/fs:0.5; % Time vector for half a second
f0 = 440; % Base frequency for DO
alpha = 2; % Alpha value for frequency calculation

% Generate the four signals for DO, RE, MI, and FA
n_values = [-9, -7, -5, -4]; % Corresponding values for n
notes = {'DO', 'RE', 'MI', 'FA'};

% Create a matrix to store signals with the length of the time vector as rows and columns as the number of notes
x = zeros(length(t), length(notes));

% Create a variable to store individual energies
energies = zeros(1, length(notes));

% Create the combined signal
combined_signal = [];

for i = 1:4
    fn = f0 * (alpha ^ (n_values(i) / 12)); % Calculate the frequency
    x(:, i) = cos(2 * pi * fn * t); % Generate and store the signal in the matrix
    energies(i) = integral(@(t) (cos(2 * pi * fn * t)).^2, 0,0.5 ); % Calculate the energy for each note using the integral function
    
    % Plot each signal
    subplot(2, 2, i);
    plot(t, x(:, i));
    title(['Signal for ' notes{i}]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
end

% Play the signals sequentially
combined_signal = reshape(x, 1, []); % Reshape the matrix to a row vector
audiowrite('Combined_Wave.wav', combined_signal, fs);

% Calculate the total energy of the combined signal
energy_x = sum(combined_signal.^2) / fs;
disp(['Energy of the combined signal: ', num2str(energy_x)]);

% FFT and Spectrum Plot
% Compute the FFT of the combined signal
N = length(combined_signal); % Number of sample points
fft_signal = fft(combined_signal);

% Compute the two-sided spectrum and then convert to a one-sided spectrum
two_sided_spectrum = abs(fft_signal / N);
one_sided_spectrum = two_sided_spectrum(1:N/2+1);
f = fs * (0:(N/2)) / N;

% Plot the frequency spectrum
figure;
plot(f, one_sided_spectrum);
title('Frequency Spectrum of the Combined Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');% Calculate the energy of the FFT-transformed signal
energy_fft = sum(abs(fft_signal).^2) / (2 * fs) / fs;
disp(['Energy of the FFT-transformed signal: ', num2str(energy_fft)]);% Step 8
filter_order = 20;
% since Mi frequency starts at 329 it would be suitable to cutoff frequency
% at 300
cutoff = (300/(fs/2));
[b, a] = butter(filter_order, cutoff,"low");

% Step 9
freqz(b,a,[],fs)

%step 10
y1 = filter(b, a, combined_signal);

%step 11
audiowrite('Low_Pass_Filter.wav',y1,fs);

%step 12
t_y1 = linspace(0, length(y1)/fs, length(y1));
figure;
plot(t_y1, y1);
title('Filtered Signal y1(t)');
xlabel('Time (s)');
ylabel('Amplitude');

%step 13
y1_energy = sum(y1.^2)/fs;
disp(['Energy of y1 signal: ',num2str(y1_energy)]);

%step 14
N_y1 = length(y1);
fft_y1 = fft(y1);

%step 15
two_sided_spectrum_y1 = abs(fft_y1 / N_y1);
one_sided_spectrum_y1 = two_sided_spectrum_y1(1:N_y1/2+1);
f_y1 = fs * (0:(N_y1/2)) / N_y1;

% Step 15: Plot the magnitude of the frequency spectrum of the filtered signal
figure;
plot(f_y1, abs(one_sided_spectrum_y1));
title('Magnitude Spectrum of Filtered Signal y1(t)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-fs/2, fs/2]);

%step 16
energy_y1_frequency_domain = sum(abs(fft_y1).^2) / (2*fs) / fs;
disp(['Energy of y1 signal in frequency domain: ',num2str(energy_y1_frequency_domain)]);

%step 17
cutoff_2 = (329/(fs/2));
[c,d] = butter(filter_order,cutoff_2, "high");

%step 18
freqz(c,d,[],fs);

%step 19
y2 = filter(c,d,combined_signal);

%step 20
audiowrite('High_Pass_Filter.wav',y2,fs);

%step 21
t_y2 = linspace(0, length(y2)/fs, length(y2));
figure;
plot(t_y2, y2);
title('Filtered Signal y2(t)');
xlabel('Time (s)');
ylabel('Amplitude');

% step 22
y2_energy = sum(y2.^2)/fs;
disp(['Energy of y2 signal: ',num2str(y2_energy)]);

% step 23
N_y2 = length(y2);
fft_y2 = fft(y2);

% step 24
two_sided_spectrum_y2 = abs(fft_y2 / N_y2);
one_sided_spectrum_y2 = two_sided_spectrum_y1(1:N_y2/2+1);
f_y2 = fs * (0:(N_y2/2)) / N_y2;

% step 24

% Combine the frequency spectra for y1 and y2
combined_one_sided_spectrum = [one_sided_spectrum_y1, one_sided_spectrum_y2];
combined_frequencies = [-fs/2 + f_y1, f_y2];

% Plot the combined frequency spectrum
figure;
plot(combined_frequencies, abs(combined_one_sided_spectrum));
title('Combined Magnitude Spectrum of Filtered Signals');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-fs/2, fs/2]);

% step 25
energy_y2_frequency_domain = sum(abs(fft_y2).^2) / (2*fs) / fs;
disp(['Energy of y2 signal in frequency domain: ',num2str(energy_y2_frequency_domain)]);
