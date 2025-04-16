%% *Method 2: SWT-SVD with Fuzzy Logic*

% Clear workspace and command window
clc; clear; close  all;
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
excelFileName = 'audio_watermarking2_results.xlsx';
%%
% Define the different audio options in a struct
audioOptions = struct(...
    'Male',        struct('original', 'Male Voice.wav',        'watermarkaudio', 'Reconstructed_Male2.wav',       'segmentDuration', 2), ...
    'Female',      struct('original', 'Female Voice.wav',      'watermarkaudio', 'Reconstructed_Female2.wav',     'segmentDuration', 2), ...
    'Classical',   struct('original', 'Classical.wav',         'watermarkaudio', 'Reconstructed_Class2.wav',      'segmentDuration', 1), ...
    'Electronic',  struct('original', 'Electronic.wav',        'watermarkaudio', 'Reconstructed_Elec2.wav',       'segmentDuration', 1), ...
    'Rock',        struct('original', 'Rock.wav',              'watermarkaudio', 'Reconstructed_Rock2.wav',       'segmentDuration', 1), ...
    'Jazz',        struct('original', 'Jazz.wav',              'watermarkaudio', 'Reconstructed_Jazz2.wav',       'segmentDuration', 1), ...
    'Pop',         struct('original', 'Pop.wav',               'watermarkaudio', 'Reconstructed_Pop2.wav',        'segmentDuration', 1) ...
);

% Get the list of all audio option names (fieldnames of the struct)
audioFields = fieldnames(audioOptions);

% Loop over each audio option
for idx = 1:length(audioFields)
    selectedOption = audioFields{idx};  % Get the current audio option name
    audioData = audioOptions.(selectedOption);  % Get the audio data for the current option
    baseFolderPath='C:\Users\Jem\Desktop\Audio Watermarking\Method2\';
    folderPath= fullfile(baseFolderPath,selectedOption);

    % Now, 'audioData' contains all the necessary values
    original = audioData.original;
    watermarkaudio = fullfile('C:\Users\Jem\Desktop\Audio Watermarking\Method2', audioData.watermarkaudio);
    segmentDuration = audioData.segmentDuration;

    % Display the current audio being processed
    disp(['Processing audio: ', original]);
    disp(['Watermarked audio: ', watermarkaudio]);
    disp(['Segment Duration: ', num2str(segmentDuration)]);
%%
    tic;
     % Step 1: Load the audio file
    [audio, fs] = audioread(original);  % Load audio file
    segmentSamples = round(segmentDuration * fs);  % Number of samples per segment (ensure it's an integer)
    numSegments = floor(length(audio) / segmentSamples);  % Number of complete segments
    
    img  = imread('60x60.bmp'); %Get the input image
    I = rgb2gray(img);
    po=imbinarize(I);
    w1=imbinarize(I);%Convert to grayscale image
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
    
    [cA_wm, cH_wm, cV_wm, cD_wm] = swt2(w1, 1, 'db6');
    
    [Uwm, Swm, Vwm] = svd(cA_wm);
    Swatermark=Swm(Swm>0);
    Swatermark=Swatermark(1:(numSegments-minus));
    nonzero_indices = find(Swm>0);
    nonzero_indices=nonzero_indices(1:(numSegments-minus));

    
    % Sort the array in descending order
    sorted_array = sort(Swatermark, 'descend');
        
    % Find the second highest value
    if length(sorted_array) >= 2
        MImportant=max(Swatermark);
        Important = sorted_array(2);
        Important2 = sorted_array(3);
    else
        error('Array does not have enough elements to determine the second highest value.');
    end

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
        min_value1=min(Saudio_array);
        max_value1=max(Saudio_array);
        min_value2=min(S_padded);
        max_value2=max(S_padded);
        
        Vaudio = Vaudio_array{i};
        cL = cL_array{i};
        Swater = S_padded(i);
    
        % Update the singular value matrix
        Snew=evaluateFuzzySystemAdjusted(Saudio,Swater,min_value1,max_value1,min_value2,max_value2);
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

    % Memory usage during embedding
    [usr, sys] = memory;
    memoryEmbed = usr.MemUsedMATLAB / 1e6;

    elapsedTime = toc;
    [usr, sys] = memory;
    fprintf('Embedding Elapsed time: %.4f seconds\n', elapsedTime);
    fprintf('Memory used by MATLAB Embedding: %.2f MB\n', usr.MemUsedMATLAB / 1e6);
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
    tic;
    % Step 1: Load the watermarked audio file
    [audio, fs] = audioread(watermarkaudio);  % Load watermarked audio file
    
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
    
        Swmextract = evaluateFuzzySystemReverse(i,SmodVal,Important,Important2,MImportant,length(Swmextract_array));
        
    
        % Store the result in the cell array
        Swmextract_array{i} = Swmextract;
    end
    % Convert the cell array to a matrix and filter out non-positive values
    Swmextract_array = cell2mat(Swmextract_array);
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

        elapsedTimeExtraction = toc;

        % Memory usage during extraction
        [usr, sys] = memory;
        memoryExtract = usr.MemUsedMATLAB / 1e6;


        %elapsedTime = toc;
        %[usr, sys] = memory;
        %fprintf('Extraction Elapsed time: %.4f seconds\n', elapsedTime);
        %fprintf('Memory used by MATLAB Extraction: %.2f MB\n', usr.MemUsedMATLAB / 1e6);
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

        % Store the results
        audioNames{end+1} = selectedOption;
        berArray(end+1) = BER;
        snrArray(end+1) = snr_value;
        odgArray(end+1) = ODG;
        elapsedTimeEmbedArray(end+1) = elapsedTimeEmbedding;
        elapsedTimeExtractArray(end+1) = elapsedTimeExtraction;
        memoryEmbedArray(end+1) = memoryEmbed;
        memoryExtractArray(end+1) = memoryExtract;


        disp(['BER for ', selectedOption, ' audio: ', num2str(BER)]);
        disp(['SNR for ', selectedOption, ' audio: ', num2str(snr_value), ' dB']);
        % The ODG is now stored in the variable 'ODG' and can be printed
        disp(['The Objective Difference Grade (ODG) is: ', num2str(ODG)]);
    end
end

% Create a table to hold all the results
resultsTable = table(audioNames', berArray', snrArray', odgArray', ...
    elapsedTimeEmbedArray', elapsedTimeExtractArray', memoryEmbedArray', memoryExtractArray', ...
    'VariableNames', {'Audio', 'BER', 'SNR', 'ODG', 'EmbeddingTime', 'ExtractionTime', 'MemoryEmbedding', 'MemoryExtraction'});

% Write the table to an Excel file
writetable(resultsTable, excelFileName);

disp(['Results saved to ', excelFileName]);