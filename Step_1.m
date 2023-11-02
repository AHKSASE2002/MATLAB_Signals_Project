% Define constants
fs = 44100; % Sampling frequency (Hz)
t = 0:1/fs:0.5; % Time vector for half a second
f0 = 440; % Base frequency for DO
alpha = 2; % Alpha value for frequency calculation

% Generate the four signals for DO, RE, MI, and FA
n_values = [-9, -7, -5, -4]; % Corresponding values for n
notes = {'DO', 'RE', 'MI', 'FA'};

% Create a matrix to store signals with length of time vector as rows & couluns as no. of notes
x = zeros(length(t), length(notes));  

for i = 1:4
    fn = f0 * (alpha ^(n_values(i)/12)); % Calculate the frequency
    x(:, i) = cos(2 * pi * fn * t); % Generate and store the signal in the matrix
end

% Plot the signals
for i = 1:4
    subplot(2, 2, i);
    plot(t, x(:, i));
    title(['Signal for ' notes{i}]);
    xlabel('Time (s)');
    ylabel('Amplitude');
end

% Play the signals sequentially
combined_signal = reshape(x, 1, []); % Reshape the matrix to a row vector
sound(combined_signal, fs);

% Save the combined signal as an audio file
audiowrite('combined_notes.wav', combined_signal, fs);

energy_x = sum(combined_signal.^2); % Energy is the sum of squared samples

fprintf('Energy of the combined signal: %f\n', energy_x);
