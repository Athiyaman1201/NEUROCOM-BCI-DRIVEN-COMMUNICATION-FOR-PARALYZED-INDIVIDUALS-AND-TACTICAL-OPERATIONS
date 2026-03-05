% MATLAB: EEG 5 bands + FFT bins with detailed preprocessing

filename = 'synthetic_eeg and test1.mp3';	% Replace with MP3 file

% Read audio
[y, Fs] = audioread(filename);

% Convert stereo to mono if needed if size(y,2) == 2
y = mean(y,2);
end

% --- EEG Preprocessing Parameters ---
segLength = 2 * Fs;	% 2-second segments
numBins = 23;	% first 23 FFT bins

% Band-pass filter (0.5 - 100 Hz)
[b_bp, a_bp] = butter(4, [0.5 100]/(Fs/2), 'bandpass');

% Notch filter (50 Hz)
[b_notch, a_notch] = iirnotch(50/(Fs/2), 50/(Fs/2)/35);

% Preallocate matrix: 5 band powers + 23 FFT bins = 28 columns numSegments = floor(length(y)/segLength);
data_matrix = zeros(numSegments, 5 + numBins);

% EEG bands in Hz delta_band = [0.5 4];
theta_band = [4 8];
alpha_band = [8 13];
beta_band = [13 30];
gamma_band = [30 100];

for i = 1:numSegments
idx_start = (i-1)*segLength + 1;

idx_end = idx_start + segLength - 1; segment = y(idx_start:idx_end);

% --- Step 1: Band-pass filtering --- segment = filter(b_bp, a_bp, segment);

% --- Step 2: Notch filtering ---
segment = filter(b_notch, a_notch, segment);

% --- Step 3: Baseline correction / detrending --- segment = detrend(segment);

% --- Step 4: Artifact removal (simple threshold method example) ---
% Remove extreme spikes beyond 3 std deviations threshold = 3 * std(segment);
segment(segment > threshold) = 0; segment(segment < -threshold) = 0;

% --- FFT computation --- Y = fft(segment);
Y_mag = abs(Y(1:floor(length(Y)/2))) / length(Y);

% Frequency vector
freqs = (0:length(Y_mag)-1) * (Fs/segLength);

% --- Compute band powers ---

delta_power = sum(Y_mag(freqs >= delta_band(1) & freqs < delta_band(2)).^2);

theta_power = sum(Y_mag(freqs >= theta_band(1) & freqs < theta_band(2)).^2);

alpha_power = sum(Y_mag(freqs >= alpha_band(1) & freqs < alpha_band(2)).^2);

beta_power  = sum(Y_mag(freqs >= beta_band(1)  & freqs < beta_band(2)).^2);

gamma_power = sum(Y_mag(freqs >= gamma_band(1) & freqs <= gamma_band(2)).^2);



% First 5 columns = band powers (Alpha, Beta, Gamma, Delta, Theta)
data_matrix(i, 1:5) = [alpha_power, beta_power, gamma_power, delta_power, theta_power];
% Column headers
header = [ {'Alpha Power','Beta Power','Gamma Power','Delta Power','Theta Power'}, ...
strcat('FFT_', string(1:numBins)) ];

% Write CSV with headers fid =
fopen('EEG_FFT_2sec_segments.csv','w'); fprintf(fid,'%s,', header{1:end-1}); fprintf(fid,'%s\n', header{end}); fclose(fid);

dlmwrite('EEG_FFT_2sec_segments.csv', data_matrix, '-append');

disp('CSV saved: EEG_FFT_2sec_segments.csv with 5 band powers first, then 23 FFT bins.');

