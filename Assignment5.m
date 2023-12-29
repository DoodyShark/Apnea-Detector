clear;
clc;

% Step 2: Reading File and Resampling
fprintf('We''ll be dealing with the a01 information\n');
FID = fopen('a01.dat');
if (FID == -1)
    error('Couldn''t open the file, please make sure the file is in your working directory\n');
else
    rawSignal = fread(FID, 'int16');
    fclose(FID);
    Fs1 = 100; % Sampling frequency of the original data
    AtoVCF = 200; % A/D to milliVolts conversion factor
    rawSignal = rawSignal(1:4*3600*Fs1); % Truncating the values for times after 4 hours
    rawSignal = rawSignal./AtoVCF; % Converting raw signal data into Volts
    
    fprintf('a01.dat successfully read. Graph of first 2 seconds generated...\n');

    % Plotting the first 2 seconds of the raw signal
    figure('name','Raw Signal');
    numSec = 2;
    t = 1/Fs1:1/Fs1:3600*4; % Time axis in seconds
    grid on
    plot(t(1:numSec*Fs1),rawSignal(1:numSec*Fs1));
    xlabel('Time (sec)');
    ylabel('Voltage (mV)');
    title('Raw Signal');
    hold off

    saveas(gcf,'RawSignal.fig');
    
    % Resampling to 500 Hz
    Fs2 = 500;
    t2 =  1/Fs2:1/Fs2:3600*4;
    resampledSignal1 = interp1(t, rawSignal, t2, 'spline');

    fprintf('Signal successfully resampled. Graph of first 2 seconds generated...\n');
    
    % Plotting the first 2 seconds of the resampled signal
    figure('name','Resampled Signal');
    numSec = 2;
    plot(t2(1:numSec*Fs2),resampledSignal1(1:numSec*Fs2), '.');
    grid on
    hold on
    plot(t(1:numSec*Fs1),rawSignal(1:numSec*Fs1), '.');
    xlabel('Time (sec)');
    ylabel('Voltage (mV)');
    title('Raw Signal VS Resampled Signal');
    legend({'Resampled Signal', 'Raw Signal'})
    hold off
    
    saveas(gcf,'ResampledSignal.fig');

    fprintf('Signal successfully filtered. Graph of first 2 seconds generated...\n');

    % Step 3: Applying the low pass filter
    lowPassFilter = designfilt('lowpassfir', 'Filterorder', 25, 'CutoffFrequency', 16, 'SampleRate', Fs2);
    filteredSignal = filter(lowPassFilter, resampledSignal1);
    
    % Plotting the first 60 seconds of the filtered signal
    figure('name','Filtered Signal');
    numSec = 60;
    grid on
    plot(t2(1:numSec*Fs2),filteredSignal(1:numSec*Fs2));
    hold on
    plot(t2(1:numSec*Fs2),resampledSignal1(1:numSec*Fs2));
    title('Resampled Signal VS Filtered Signal');
    legend({'Filtered Signal', 'Resampled Signal'})
    xlabel('Time (sec)');
    ylabel('Voltage (mV)');
    hold off
    
    saveas(gcf,'FilteredSignal.fig');

    % Step 4: Finding the Average. Minimum, and Maximum HR for the first
    % minute
    
    minute = filteredSignal(1:1*60*Fs2);
    RIndeces = findRIndices(minute);
    RTimes = RIndeces./Fs2;
    RR = diff(RTimes);
    RRav = mean(RR);
    MaxHR = 60/min(RR);
    MinHR = 60/max(RR);
    HRav = 60/RRav; 

    fprintf('First minute Heart-Rate information analyzed:\nMin Heart-Rate: %.0f\nMax Heart-Rate: %.0f\nAverage Heart-Rate: %.0f\n', MinHR, MaxHR, HRav);
    
    % Step 5: Finding the Average, Minimum, and Maximum HR for each minute
    HRav = [];
    MaxHR = [];
    MinHR = [];
    
    for i = 1:240
        minute = filteredSignal((i-1)*60*Fs2 + 1:i*60*Fs2);
        RIndeces = findRIndices(minute);
        RTimes = RIndeces./Fs2;
        RR = diff(RTimes);
        RRav = mean(RR);
        MaxHR = [ MaxHR, 60/min(RR)];
        MinHR = [MinHR, 60/max(RR)];
        HRav = [HRav, 60/RRav];
    end
    
    fprintf('First 4-hour Heart-Rate information analyzed. Graphs generated\n');

    % Plotting the average, minimum, & maximum heart rate for each minute over time
    figure('name','HR Variation');
    minutes = 1:240;
    grid on
    plot(minutes, HRav);
    title('Average HR Over the 4 Hr Period');
    xlabel('Time (min)');
    ylabel('HR (bpm)');
    
    saveas(gcf,'HRavVariation.fig');

    figure('name','Max and Min HR Variation');
    grid on
    plot(minutes, MaxHR);
    hold on
    plot(minutes, MinHR);
    legend({'Maximum HR', 'Minimum HR'});
    title('Min and Max HR Over the 4 Hr Period');
    xlabel('Time (min)');
    ylabel('HR (bpm)');
    hold off
    
    saveas(gcf,'HRVariation.fig');
    
    % Step 6: Performing a Fast-Fourier-Transform on the First Minute of Data
    
    Fs3 = 4;
    t3 = 1/Fs3: 1/Fs3:4*3600;
    resampledSignal2 = resample(filteredSignal, Fs3, Fs2);
    
    minute = resampledSignal2(1:1*60*Fs3);
    N = length(minute);
    fftOfMin = fft(minute);
    PSD = (1/(Fs3*N))* abs(fftOfMin).^2;
    PSD = PSD(1:N/2); % Removing the negative side
    
    fprintf('Spectral analysis of first minute completed. Graph generated\n');

    % Plotting the Power-Spectral-Density of the first minute
    fAxis = 1/N:Fs3/N:Fs3/2;
    figure('name','Power-Spectral-Density First Minute');
    plot(fAxis, PSD);
    grid on
    title ('Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power');
    hold off
    
    saveas(gcf,'PSDFirtsMinute.fig');


    % Step 7: Performing FFt on every minute and plotting the spectrogram
    
    % Performing FFT on every minute and saving PSD into PSDArray
    PSDArray = [];
    hold on
    for i = 1:240
        minute = resampledSignal2((i-1)*60*Fs3 + 1:i*60*Fs3);
        %t = i/Fs + 1:1/Fs:i*60;
        N = length(minute);
        fftOfMin = fft(minute);
        PSD = (1/(Fs3*N))* abs(fftOfMin).^2;
        PSD = PSD(1:N/2);
        PSDArray(:,i) = PSD;
    end
    
    fprintf('Spectral analysis of the 4 hours completed. Spectrogram generated\n');

    figure('name','Spectrogram');
    fAxis = 1/N:Fs3/N:Fs3/2;
    tAxis = 1:240;
    surf(tAxis, fAxis, PSDArray); % A spectrogram is a 3D plot of the PSD looked at from above
    hold on
    view(2); % To look at the data from above
    plot3(1:length(HRav), HRav./(HRav(1)/1.05), ones(length(HRav)), 'y', 'LineWidth',2)
    title('Spectrogram of ECG');
    xlabel('Time (min)');
    ylabel('Frequency (Hz)');
    colorbar % Shows the z values for every color, useful to determine the minimum
    % Peak threshold for the algorithm
    hold off
        
    saveas(gcf,'Spectrogram.fig');

    
    % Step 8: Apnea-Detecting Algorithm
    

    % Analyizng the data:
    % Uncomment this section if you'd like to try it but make sure to
    % include the text file a01Labled.txt that I attached in my submission
%     fprintf('Analyzing the data...\n');
%     figure('name','Spectrogram Analysis');
%     fAxis = 1/N:Fs3/N:Fs3/2;
%     tAxis = 1:240;
%     surf(tAxis, fAxis, PSDArray); % A spectrogram is a 3D plot of the PSD looked at from above
%     hold on
%     view(2); % To look at the data from above
%     plot3(1:length(HRav), HRav./(HRav(1)/1.05), ones(length(HRav)), 'y', 'LineWidth',2)
%     title('Analyzing the Spectrogram of ECG');
%     %legend({'', 'Average HR'});
%     xlabel('Time (min)');
%     ylabel('Frequency (Hz)');
%     colorbar
%     
%     FID = fopen("a01Labeled.txt");
%     raw = fread(FID);
%     labeled = [];
%     j = 1;
%     for i = 1:length(raw)
%         if (raw(i) == 65)
%             labeled(j) = 1;
%             j = j + 1;
%         end
%         if (raw(i) == 78)
%             labeled(j) = 0;
%             j = j + 1;
%         end
%     end
%     labeled = labeled(1:240);
%     plot3(1:length(labeled), labeled.*1.08, ones(length(labeled)), '_w')
%     clear raw;
%     hold off
%     
%     saveas(gcf,'SpectrogramAnalysis.fig');
% 
%     figure('name','Analyzing HR Variation');
%     plot(tAxis, HRav);
%     hold on
%     title('Analyzing Average HR')
%     xlabel('Time (min)');
%     ylabel('Heart Rate (bpm)');
%     plot(labeled.*68, '.');
%     hold off
% 
%     saveas(gcf,'HRAnalysis.fig');

    
    % After analyzing the data for a01.dat, a minimum frequency threshold is
    % about 1.1 Hz and minimum HR threshold is 68
    
    fprintf('Apnea detection in progress...\n');

    % Algorithm:
    
    predictions = [];
    for i = 1:length(HRav)
        % Seeing if there are any peaks in the PSD below the threshold
        % requency 
        belowLinePeaks = findpeaks(PSDArray(1:65, i), 'MinPeakHeight', 0.002); % Index 65 corresponds to the frequency 1.1 Hz which was decided on after analyzing the data
        [msg id] = lastwarn; 
        warning('off',id) % Suppressing the warning message when the findpeaks finds no peaks
        if (HRav(i) <= 68 || any(belowLinePeaks)) % This threshold for HRav was decided upon after visually analyzing the data
            predictions(i) = 1;
        else
            predictions(i) = 0;
        end
    end
    [msg id] = lastwarn;
    warning('off',id) % Suppressing the warning message when the findpeaks finds no peaks
    
    
    % Measuring Prediction Accuracy:
    % Uncomment this section and the other commented section above to try it out
%     figure
%     plot(predictions, '.');
%     hold on
%     plot(labeled.*2, '.');
%     hold off
%     numSame = 0;
%     for i = 1:240
%         if (labeled(i) == predictions(i))
%             numSame = numSame  + 1;
%         end
%     end
%     accuracy = (numSame / 240)  * 100;
%     fprintf('The accuracy of generated labels is %.2f%%\n', accuracy);

    % Printing into File:
    
    FID = fopen("ApneaPredictions.dat", 'w');
    fprintf(FID, 'Time(min): Prediction:\n');
    for i = 1:240
        if (predictions(i) == 1)
            prediction = 'A';
        else
            prediction = 'N';
        end
        fprintf(FID, '%d %c\n', i, prediction);
    end
    fclose(FID);
    fprintf('Successfully completed execution of program. All Graphs saved into your system...\nHave a great day!\n');
end
