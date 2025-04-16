%% *Method 1: SWT-SVD with Arnold Transform*

% Clear workspace and command window
clc; clear; close all;
%%
audioNames = {};
berArray = [];
snrArray = [];
odgArray = [];
elapsedTimeEmbedArray = [];
elapsedTimeExtractArray = [];
memoryEmbedArray = [];
memoryExtractArray = [];

% Define Excel file name
excelFileName = 'audio_watermarking1_results.xlsx';
%%
% Define the different audio options in a struct
audioOptions = struct(...
    'Male',        struct('original', 'Male Voice.wav',        'watermarkaudio', 'Reconstructed_Male1.wav',       'segmentDuration', 2), ...
    'Female',      struct('original', 'Female Voice.wav',      'watermarkaudio', 'Reconstructed_Female1.wav',     'segmentDuration', 2), ...
    'Classical',   struct('original', 'Classical.wav',         'watermarkaudio', 'Reconstructed_Class1.wav',      'segmentDuration', 1), ...
    'Electronic',  struct('original', 'Electronic.wav',        'watermarkaudio', 'Reconstructed_Elec1.wav',       'segmentDuration', 1), ...
    'Rock',        struct('original', 'Rock.wav',              'watermarkaudio', 'Reconstructed_Rock1.wav',       'segmentDuration', 1), ...
    'Jazz',        struct('original', 'Jazz.wav',              'watermarkaudio', 'Reconstructed_Jazz1.wav',       'segmentDuration', 1), ...
    'Pop',         struct('original', 'Pop.wav',               'watermarkaudio', 'Reconstructed_Pop1.wav',        'segmentDuration', 1) ...
);

% Get the list of all audio option names (fieldnames of the struct)
audioFields = fieldnames(audioOptions);

% Loop over each audio option
for idx = 1:length(audioFields)
    
    selectedOption = audioFields{idx};  % Get the current audio option name
    audioData = audioOptions.(selectedOption);  % Get the audio data for the current option
    baseFolderPath='C:\Users\Jem\Desktop\Audio Watermarking\Method1\';
    folderPath= fullfile(baseFolderPath,selectedOption);

    % Now, 'audioData' contains all the necessary values
    original = audioData.original;
    watermarkaudio = fullfile('C:\Users\Jem\Desktop\Audio Watermarking\Method1', audioData.watermarkaudio);
    segmentDuration = audioData.segmentDuration;



    % Display the current audio being processed
    disp(['Processing audio: ', original]);
    disp(['Watermarked audio: ', watermarkaudio]);
    disp(['Segment Duration: ', num2str(segmentDuration)]);
%%
    tic;
    profile on;  % Start profiling

     % Step 1: Load the audio file
    [audio, fs] = audioread(original);  % Load audio file
    segmentSamples = round(segmentDuration * fs);  % Number of samples per segment (ensure it's an integer)
    numSegments = floor(length(audio) / segmentSamples);  % Number of complete segments
    
    img  = imread('60x60.bmp'); %Get the input image
    I = rgb2gray(img);
    po=imbinarize(I);
    w1=imbinarize(I);%Convert to grayscale image
    iter=2;
    minus=1;
    
    w1=reshape(w1,1800,2);
    n = 3; k = 2; % A (3,2) cyclic code
    code1 = encode(double(w1),n,k,'cyclic/binary');   % 2 level of encoding
    code2 = encode(code1(:,1:2),n,k,'cyclic/binary');
    code3 = encode(code1(:,2:3),n,k,'cyclic/binary');
    part1=reshape(code2(:,1:2),1800,2);
    part2=reshape(code2(:,2:3),1800,2);
    part3=reshape(code3(:,1:2),1800,2);
    part4=reshape(code3(:,2:3),1800,2);
    a1=horzcat(part1,part2);
    b1=horzcat(part3,part4);
    z=vertcat(a1,b1);
    w1=reshape(z,120,120);
    
    w1=arnold(w1,iter);
    
    [cA_wm, cH_wm, cV_wm, cD_wm] = swt2(w1, 1, 'db6');
    
    [Uwm, Swm, Vwm] = svd(cA_wm);
    Swatermark=Swm(Swm>0);
    Swatermark=Swatermark(1:numSegments-minus);
    Swatermark=Swatermark*0.01;
    nonzero_indices = find(Swm>0);
    nonzero_indices=nonzero_indices(1:numSegments-minus);

    % Pad the audio with zeros if necessary to ensure the length is divisible by numSegments
    paddedLength = segmentSamples * numSegments;  % Calculate the padded length
    audio_padded = [audio; zeros(paddedLength - length(audio), 1)];  % Pad audio with zeros
    
    % Step 2: Preallocate a cell array to store low-frequency approximation (cA)
    cA_array = cell(numSegments, 1);  % Preallocate a cell array for storing cA
    cL_array = cell(numSegments, 1);
    % First for-loop: Perform SWT and store cA for each segment
    for i = 1:numSegments
        % Extract each segment
        startIdx = (i-1)*segmentSamples + 1;
        endIdx = i*segmentSamples;
        segment = audio_padded(startIdx:endIdx);  % Extract segment
        
        % Perform SWT (Stationary Wavelet Transform)
        [cA, cL] = swt(segment, 1, 'haar');  % Use 1-level SWT with Daubechies wavelet ('db1')
        
        % Store the low-frequency approximation (cA)
        cA_array{i} = cA;
        cL_array{i} = cL;
    end
    
    % Step 3: Preallocate cell arrays to store SVD results
    Uaudio_array = cell(numSegments, 1);  % Preallocate a cell array for U
    Saudio_array = cell(numSegments, 1);  % Preallocate a cell array for S
    Vaudio_array = cell(numSegments, 1);  % Preallocate a cell array for V
    
    % Second for-loop: Perform SVD on each stored cA and extract singular values
    for i = 1:numSegments
        % Retrieve the stored cA
        cA = cA_array{i};  % Use cell array indexing
        
        % Perform SVD on the low-frequency approximation (cA)
        [Uaudio, Saudio, Vaudio] = svd(cA, 'econ');  % Perform SVD with economy mode
        
        % Store SVD results
        Uaudio_array{i} = Uaudio;
        Saudio_array{i} = Saudio;
        Vaudio_array{i} = Vaudio;
    end
    Saudio_array=cell2mat(Saudio_array);
    
    % Step 4: Preallocate an array to store the reconstructed segments
    reconstructed_audio = zeros(paddedLength, 1);  % Initialize reconstructed audio array
    
    % Ensure that S_positive is the same length as allLowFreqCoeffs
    S_padded = [Swatermark; zeros(length(Saudio_array) - length(Swatermark), 1)];
    
    Snew_array = cell(numSegments, 1);
    % For-loop: Reconstruct each segment with fuzzy logic adjustments
    for i = 1:numSegments
        % Retrieve U, S, V for the current segment
        Uaudio = Uaudio_array{i};
        Saudio = Saudio_array(i);
        Vaudio = Vaudio_array{i};
        cL = cL_array{i};
        Swater = S_padded(i);

        % Update the singular value matrix
        Snew = Saudio + Swater;
        Snew_array{i} = Snew;
        % Reconstruct the low-frequency approximation (cA) using SVD components
        cA_reconstructed = Uaudio * Snew * Vaudio';

        % Perform inverse SWT (Stationary Wavelet Transform) to reconstruct the segment
        segment_reconstructed = iswt(cA_reconstructed, cL, 'haar');

        % Place the reconstructed segment into the full audio signal
        startIdx = (i-1) * segmentSamples + 1;
        endIdx = min(i * segmentSamples, length(audio));
        reconstructed_audio(startIdx:endIdx) = segment_reconstructed(1:(endIdx - startIdx + 1));
    end
    Snew_array=cell2mat(Snew_array);
    % Trim any padding added to match the original length of the audio
    reconstructed_audio = reconstructed_audio(1:length(reconstructed_audio));
    
    % Step 5: Save or play the reconstructed audio
    audiowrite(watermarkaudio, reconstructed_audio, fs);  % Save the reconstructed audio file
    
    % Sample Data
    x1 = 1:30; % X-axis data for the first dataset
    x2 = 1:30; % X-axis data for the second dataset
    y1 = Saudio_array; % Example Y-axis data for the first dataset (random numbers scaled)
    y2 = Snew_array; % Example Y-axis data for the second dataset (random numbers scaled)

    % Create the Scatter Plot
    figure; % Create a new figure
    scatter(x1, y1, 'red', 'filled'); % First dataset in red
    hold on; % Hold the plot for the second dataset
    scatter(x2, y2, 'blue', 'filled'); % Second dataset in blue

    % Customize Plot Appearance
    title('Raw vs Watermarked Audio Singular Values'); % Add a title
    xlabel('INDEX'); % Label for the x-axis
    ylabel('VALUES'); % Label for the y-axis

    % Position Legend Below the Title
    legend({'Raw', 'Watermarked'}, 'Location', 'north', 'Orientation', 'horizontal'); % Adjust legend position
    grid on; % Turn on the grid
    hold off; % Release the plot

    exportgraphics(gcf, fullfile(folderPath, [num2str(selectedOption), '*.png']));

    elapsedTimeEmbedding = toc;
    profile off; % Stop profiling
    %profile viewer;

    % Memory usage during embedding
    [usr, sys] = memory;
    memoryEmbed = usr.MemUsedMATLAB / 1e6;
    %elapsedTime = toc;
    %[usr, sys] = memory;
    %fprintf('Embedding Elapsed time: %.4f seconds\n', elapsedTime);
    %fprintf('Memory used by MATLAB Embedding: %.2f MB\n', usr.MemUsedMATLAB / 1e6);

%% 
% *Watermark Attack Test*
% 
% _a. Low-pass filtering_


% Example for a lowpass filter design in MATLAB
Fs = 22050; % Sampling frequency (Hz)
Fc = 9000;  % Cutoff frequency (Hz)
order = 4;  % Filter order

% Design the low-pass filter
[b, a] = butter(order, Fc/(Fs/2), 'low');

% Load the watermarked audio signal
[x, Fs] = audioread(watermarkaudio);

% Apply the low-pass filter
y_lpf = filter(b, a, x);

% Save or play the filtered audio
audiowrite('audio_lpf.wav', y_lpf, Fs);
%% 
% _b. Band-pass filtering_

% Band-pass Filtering (BPF)
Fs = 44100; % Sampling frequency (Hz)
lowCutoff = 100; % Low cutoff frequency (Hz)
highCutoff = 9000; % High cutoff frequency (Hz)
order = 4; % Filter order

% Design the band-pass filter
[b, a] = butter(order, [lowCutoff highCutoff]/(Fs/2), 'bandpass');

% Load the watermarked audio signal
[x, Fs] = audioread(watermarkaudio);

% Apply the band-pass filter
y_bpf = filter(b, a, x);

% Save or play the filtered audio
audiowrite('audio_bpf.wav', y_bpf, Fs);
%% 
% _c. MP3 Lossy Compression_

% MATLAB code to simulate MP3 lossy compression attack on an existing WAV file

% Specify the path to your existing WAV file
wav_filename = watermarkaudio;

% Perform MP3 compression at 96 kbps
system(['ffmpeg -y -i "' wav_filename '" -b:a 96k audio_96kbps.mp3']);

% Perform MP3 compression at 64 kbps
system(['ffmpeg -y -i "' wav_filename '"  -b:a 64k audio_64kbps.mp3']);

%% 
% _d. Noise Addition_

% Noise Addition
Fs = 44100; % Sampling frequency (Hz)
SNR_30dB = 30; % Signal-to-Noise Ratio in dB
SNR_10dB = 10; % Signal-to-Noise Ratio in dB

% Load the watermarked audio signal
[x, Fs] = audioread(watermarkaudio);

% Add noise with SNR of 30 dB
y_30dB = awgn(x, SNR_30dB, 'measured');

% Add noise with SNR of 10 dB
y_10dB = awgn(x, SNR_10dB, 'measured');

% Save or play the noisy audio
audiowrite('audio_30dB_noise.wav', y_30dB, Fs);
audiowrite('audio_10dB_noise.wav', y_10dB, Fs);
%% 
% _e. Advance: Cropping Attack_

% Cropping Attack
% Load the watermarked audio signal
[x, Fs] = audioread(watermarkaudio);

% Calculate the number of samples to remove
num_samples = length(x);
num_remove = round(0.1 * num_samples);

% Randomly select indices to remove
indices = randperm(num_samples, num_remove);

% Remove the selected samples
x_cropped = x;
x_cropped(indices, :) = []; 

% Save or play the cropped audio
audiowrite('audio_cropped.wav', x_cropped, Fs);
%% 
% _f. Advance: Time Shifting_

% Time Shifting Attack (Expand by 5%)
% Load the watermarked audio signal
[x, Fs] = audioread(watermarkaudio);

% Expansion factor (5%)
factor = 1 / (1 + 0.05);

% Resample the audio
x_expanded_5 = resample(x, round(Fs * factor), Fs);

% Save or play the expanded audio
audiowrite('audio_expanded_5.wav', x_expanded_5, round(Fs * factor));
%%

% Load the original signal (Classical.wav)
[orig, fs] = audioread(original);

% Load the reconstructed signal (Reconstructed_Classical.wav)
[watermarked, fs] = audioread(watermarkaudio);

% Ensure both audio arrays have the same length
min_length = min(length(orig), length(watermarked));
orig = orig(1:min_length);  % Trim orig to match the shorter length
watermarked = watermarked(1:min_length);  % Trim watermarked to match the shorter length
    
%%
   % Preallocate an array to store Swmextract_array for each mode
allSwmextract = cell(1, 9); % Create a cell array with 8 elements (one for each mode)

% Loop through modes 1 to 9
for choice = 1:9
    % Display the current mode being processed
    fprintf('Processing mode %d\n', choice);

    % Assign the appropriate filename based on the current mode
    switch choice
        case 1
            mode = watermarkaudio; % Original Watermarked Audio
        case 2
            mode = 'audio_lpf.wav';      % Low-Pass Filtering
        case 3
            mode = 'audio_bpf.wav';      % Band-Pass Filtering
        case 4
            mode = 'audio_96kbps.mp3';   % MP3 Compression (96 kbps)
        case 5
            mode = 'audio_64kbps.mp3';   % MP3 Compression (64 kbps)
        case 6
            mode = 'audio_30dB_noise.wav'; % Noise Addition (30 dB)
        case 7
            mode = 'audio_10dB_noise.wav'; % Noise Addition (10 dB)
        case 8
            mode = 'audio_cropped.wav';  % Cropping Attack
        case 9
            mode = 'audio_expanded_5.wav'; % Time Shifting Attack (Expand by 5%)
        otherwise
            error('Invalid mode selection.');
    end

    % WATERMARK EXTRACTION
    % Step 1: Load the watermarked audio file
    [audio, fs] = audioread(mode);  % Load watermarked audio file
    t = (0:length(audio)-1) / fs; % Time vector

    % Define the segment of interest
    startTime = 6; % Start time in seconds
    if audioData.segmentDuration==2
        endTime = 6.02; % End time in seconds
    else
        endTime=6.02;
    end


    % Find indices corresponding to the time range
    startIndex = round(startTime * fs) + 1; % Convert to index (1-based)
    endIndex = round(endTime * fs); % Convert to index (1-based)

    % Ensure the indices are within the bounds of the audio signal
    if endIndex > length(audio)
        endIndex = length(audio); % Adjust to prevent indexing out of bounds
    end

    % Extract the audio segment for analysis
    audioSegment = audio(startIndex:endIndex);
    tSegment = t(startIndex:endIndex);

    % Frequency-Domain Analysis
    N = length(audioSegment); % Number of samples in the segment
    fAxis = linspace(-fs/2, fs/2, N); % Frequency axis for plotting

    % Apply a window function to reduce spectral leakage
    windowedSegment = audioSegment .* hamming(N);

    % Compute Fourier Transform and shift it for centered frequency
    spectrum = fftshift(fft(windowedSegment) / N); % Normalize the FFT and shift it
    spectrum_dB = 20 * log10(abs(spectrum) + eps); % Convert to dB scale, add eps to avoid log(0)

    % Plot only the positive half of the spectrum (if audio is real-valued)
    halfIdx = floor(N / 2) + 1;
    fAxisPos = fAxis(halfIdx:end); % Positive frequency axis
    spectrum_dBPos = spectrum_dB(halfIdx:end); % Positive spectrum in dB

    % Define frequency ranges
    lowFreqLimit = 500;   % Upper bound of low frequencies (Hz)
    midFreqLimit = 2000;  % Upper bound of middle frequencies (Hz)

    % Index ranges for low, middle, and high frequencies
    lowFreqIdx = (fAxisPos < lowFreqLimit);
    midFreqIdx = (fAxisPos >= lowFreqLimit) & (fAxisPos <= midFreqLimit);
    highFreqIdx = (fAxisPos > midFreqLimit);

    % Plotting
    figure;

    % Time-Domain Plot
    subplot(4,1,1); % Create a subplot for time-domain signal
    plot(tSegment, audioSegment);
    title('Time-Domain Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    % Low-Frequency Spectrum Plot (below 500 Hz)
    subplot(4,1,2);
    plot(fAxisPos(lowFreqIdx), spectrum_dBPos(lowFreqIdx));
    title('Low-Frequency Spectrum (Below 500 Hz)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([0 lowFreqLimit]); % Limit x-axis to low-frequency range
    grid on;

    % Middle-Frequency Spectrum Plot (500 Hz to 2000 Hz)
    subplot(4,1,3);
    plot(fAxisPos(midFreqIdx), spectrum_dBPos(midFreqIdx));
    title('Middle-Frequency Spectrum (500 Hz to 2000 Hz)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([lowFreqLimit midFreqLimit]); % Limit x-axis to middle-frequency range
    grid on;

    % High-Frequency Spectrum Plot (above 2000 Hz)
    subplot(4,1,4);
    plot(fAxisPos(highFreqIdx), spectrum_dBPos(highFreqIdx));
    title('High-Frequency Spectrum (Above 2000 Hz)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([midFreqLimit fs/2]); % Limit x-axis to high-frequency range
    grid on;

    % Save the figure as an image
    saveas(gcf, fullfile(folderPath, [num2str(choice), '.png'])); % Saves the image as a PNG file

    % Pad the audio with zeros if necessary to ensure the length is divisible by numSegments
    paddedLength = segmentSamples * numSegments;  % Calculate the padded length
    audio_padded = [audio; zeros(paddedLength - length(audio), 1)];  % Pad audio with zeros

    % Step 2: Preallocate a cell array to store low-frequency approximation (cA)
    cAwma_array = cell(numSegments, 1);  % Preallocate a cell array for storing cA

    % First for-loop: Perform SWT and store cA for each segment
    for i = 1:numSegments
        % Extract each segment
        startIdx = (i-1)*segmentSamples + 1;
        endIdx = i*segmentSamples;
        segment = audio_padded(startIdx:endIdx);  % Extract segment

        % Perform SWT (Stationary Wavelet Transform)
        [cAwma, ~] = swt(segment, 1, 'haar');  % Use 1-level SWT with Daubechies wavelet ('db1')

        % Store the low-frequency approximation (cA)
        cAwma_array{i} = cAwma;
    end

    % Step 3: Preallocate cell arrays to store SVD results
    Swmaudio_array = cell(numSegments, 1);  % Preallocate a cell array for S

    % Second for-loop: Perform SVD on each stored cA and extract singular values
    for i = 1:numSegments
        % Retrieve the stored cA
        cAwm = cAwma_array{i};  % Use cell array indexing

        % Perform SVD on the low-frequency approximation (cA)
        [~, Swmaudio, ~] = svd(cAwm, 'econ');  % Perform SVD with economy mode

        % Store SVD results
        Swmaudio_array{i} = Swmaudio;
    end

    Swmextract_array = cell(numSegments, 1);
    % For-loop: Reconstruct each segment with fuzzy logic adjustments
    for i = 1:length(Swatermark)
        % Retrieve U, S, V for the current segment
        Swmaudio = Swmaudio_array{i};
        Saudio = Saudio_array(i);

        SmodVal = abs(Swmaudio - Saudio);

        if SmodVal == 0
            Swmextract = 0;  % Or any other indication that no watermark is present
        else
            Swmextract = Swmaudio-Saudio;
        end

        % Store the result in the cell array
        Swmextract_array{i} = Swmextract;
    end
    % Convert the cell array to a matrix and filter out non-positive values
    Swmextract_array = cell2mat(Swmextract_array);
    Swmextract_array=Swmextract_array/0.01;
    % Store Swmextract_array for the current mode
    allSwmextract{choice} = Swmextract_array



    if all(Swmextract_array == 0)
        disp('The loaded audio does not contain the watermark.');
    else
        % Create a new matrix with the same size as S, filled with zeros
        S_reconstructed= zeros(size(Swm));
        % Assign the non-zero singular values back to their original positions
        S_reconstructed(nonzero_indices) = Swmextract_array;

        % Correct SVD reconstruction of w2
        ecA_wm= Uwm*S_reconstructed*Vwm';
        w2=iswt2(ecA_wm,cH_wm, cV_wm, cD_wm,'db6');
        w2=iarnold(w2,iter);
        w2=abs(w2);
        w2=imbinarize(w2);

        z2=reshape(w2,3600,4);
        sub1=z2(1:1800,:);
        sub2=z2(1801:3600,:);
        sub11=sub1(1:1800,1:2);
        sub12=sub1(1:1800,3:4);
        sub21=sub2(1:1800,1:2);
        sub22=sub2(1:1800,3:4);
        new1=horzcat(sub11,sub12);
        new1(:,3)=[];
        new2=horzcat(sub21,sub22);
        new2(:,3)=[];


        % Convert to binary if not already
        new1_binary = double(new1 > 0); % Or use other thresholding method if needed
        new2_binary = double(new2 > 0);

        % Decode with binary inputs
        newmsg1 = decode(new1_binary, n, k, 'cyclic');
        newmsg2 = decode(new2_binary, n, k, 'cyclic');

        newmsg=horzcat(newmsg1,newmsg2);
        newmsg(:,3)=[];
        sbr = decode(newmsg,n,k,'cyclic');
        sbr=reshape(sbr,60,60);
        sbr=logical(sbr);

        figure;
        imshow(sbr);
        title('Extracted Watermark');
        % Save the figure as an image
        imwrite(sbr, fullfile(folderPath, [num2str(choice), 'wm.png'])); % Saves the image as a PNG file
%%
        % Calculate the noise and SNR
        noise = watermarked - orig;
        snr_value = snr(orig, noise);
        BER=biterr(po,sbr);
        
        % Set the reference and test audio file names
        Fref = original;
        Ftest = watermarkaudio;
        
        % Set the starting and ending sample positions (optional)
        StartS = [0, 0];  % Process both files from the start
        EndS = [];        % Process until the end of both files
        
        % Call the function to compute the Objective Difference Grade (ODG)
        ODG = PQevalAudio(Fref, Ftest, StartS, EndS);




        disp(['BER for ', selectedOption, ' audio: ', num2str(BER)]);
        disp(['SNR for ', selectedOption, ' audio: ', num2str(snr_value), ' dB']);
        % The ODG is now stored in the variable 'ODG' and can be printed
        disp(['The Objective Difference Grade (ODG) is: ', num2str(ODG)]);
    end
end
%%

% Sample Data
x = 1:29; % X-axis values (common for all datasets)
y_no_attack = allSwmextract{1}; % Data for 'No attack'
y_lpf = allSwmextract{2};       % Data for 'LPF'
y_bpf = allSwmextract{3};       % Data for 'BPF'
y_96kbps = allSwmextract{4};    % Data for '96 kbps'
y_64kbps = allSwmextract{5};    % Data for '64 kbps'
y_30db_noise = allSwmextract{6}; % Data for '30db noise'
y_10db_noise = allSwmextract{7}; % Data for '10db noise'
y_cropping = allSwmextract{8};  % Data for 'Cropping'
y_time_shifting = allSwmextract{9}; % Data for 'Time shifting'

% Create the Scatter Plot
figure; % Create a new figure

% Plot each dataset with its specified color and symbol
scatter(x, y_no_attack, 50, 'b', 'o', 'filled', 'DisplayName', 'No Attack'); % Blue circle
hold on; % Keep the plot open for more data
scatter(x, y_lpf, 50, 'r', '^', 'filled', 'DisplayName', 'LPF'); % Red triangle
scatter(x, y_bpf, 50, [1 0.5 0], 's', 'filled', 'DisplayName', 'BPF'); % Orange square
scatter(x, y_96kbps, 50, [0 0.5 0], 'd', 'filled', 'DisplayName', '96 kbps'); % Bluegreen diamond
scatter(x, y_64kbps, 50, 'b', '*', 'DisplayName', '64 kbps'); % Blue star
scatter(x, y_30db_noise, 50, 'g', 'x', 'DisplayName', '30db Noise'); % Green cross
scatter(x, y_10db_noise, 50, [0.678, 0.847, 0.902], 'p', 'filled', 'DisplayName', '10db Noise'); % Light blue pentagon
scatter(x, y_cropping, 50, 'b', 'h', 'DisplayName', 'Cropping'); % Blue hexagon
scatter(x, y_time_shifting, 50, [1 0.647 0], 'p', 'DisplayName', 'Time Shifting'); % Orange pentagon

% Customize Plot Appearance
title('No attack vs Attacked Watermark Singular Values'); % Add a title
xlabel('INDEX'); % Label for the x-axis
ylabel('VALUES'); % Label for the y-axis
legend('Location', 'east', 'Orientation', 'vertical'); % Add a legend below the title
grid on; % Turn on the grid
hold off; % Release the plot

% Save the Plot
exportgraphics(gcf, fullfile(folderPath, [num2str(selectedOption), '.png']));
end
%%

% Create a table to hold all the results
resultsTable = table(audioNames', berArray', snrArray', odgArray', ...
    elapsedTimeEmbedArray', elapsedTimeExtractArray', memoryEmbedArray', memoryExtractArray', ...
    'VariableNames', {'Audio', 'BER', 'SNR', 'ODG', 'EmbeddingTime', 'ExtractionTime', 'MemoryEmbedding', 'MemoryExtraction'});

% Write the table to an Excel file
writetable(resultsTable, excelFileName);

disp(['Results saved to ', excelFileName]);
profile viewer;