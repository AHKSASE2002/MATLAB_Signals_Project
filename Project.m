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
audiowrite('X.wav', combined_signal, fs);

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
ylabel('Magnitude');

% Calculate the energy of the FFT-transformed signal
energy_fft = sum(abs(fft_signal).^2) / (2 * fs) / fs;
disp(['Energy of the FFT-transformed signal: ', num2str(energy_fft)]);

% Step 8
filter_order = 20;
% since Mi frequency starts at 329 it would be suitable to cutoff frequency
% at 300
cutoff = (300/(fs/2));
[b, a] = butter(filter_order, cutoff);


