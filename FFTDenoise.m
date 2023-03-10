function [dataBlock] = FFTDenoise(dataBlock,cutOut,keepResult,lowNotch,highNotch)
%Please read LICENSE.md before using, and attribute Robin J. Abel, PhD using appropriate citation.
%
%Current Version 1.3, 2020-02-25
%
%A generic function to perform Fourier fast transform filtering to remove
%periodic (and other) chromatographic noise using notch filters via FFT,
%then reconstitute the denoised data matrix. Optional visualisation is
%included to assess the output result.
%
%Written to accept indexed data structure from the LoadPEG.m function, or
%an equivalent. Easily adaptable to other data types.
%
%Originally designed to remove periodic noise from mains & other
%electronic noise at 50 Hz and 60 Hz from LECO ChromaTOF data to resolve
%inconsistent peak detection observed in the ChromaTOF software due to
%peaks occurring over noise occurring at difference phases between
%replicate samples.
%
%Usage: [dataout] = FFTDenoise(dataBlock,cutOut,keepResult,lowNotch,highNotch)
%
%FFTDenoise - Remove noise from a LoadPEG dataframe.
%
%      dataBlock    LoadPEG/LoadCSV formatted data structure or equivalent
%      cutOut       (optional) boolean, remove desired frequencies?
%      keepResult   (optional) boolean, force keep resulting denoised matrix
%      lowNotch     (optional) numerical, specify low notch frequency (Hz)
%      highNotch    (optional) numerical, specify high notch frequency (Hz)
%
%   returns:
%      dataBlock    Input data structure, returned either unaltered, or with an added 'FFT' flag and the denoised data if the result is kept

%% Version History
% Version 1.0, around 2018-09-13, first version
% Version 1.1, 2019-02-04
%   updated displays and function to provide notching between user-supplied
%   frequencies. Can be used iteratively to notch several zones, and can be
%   used to bandpass as well.
% Version 1.2, 2019-02-05
%   added ability to specify low & high notch for batch processing
% Version 1.3, 2020-02-25
%   corrected residual plotting when keepResult is enabled. This allows
%   silent calling from within functions remain silent when desired.



%need to separate the data
if isa(dataBlock, 'struct')

    %Perform FFT on TIC & spectral data
    %just perform the FFT as usual
    tic_ffted = fft(dataBlock.tic);
    for page = 1:size(dataBlock.specdata,2)
        specdata_ffted(:,page) = fft(dataBlock.specdata(:,page));
    end

elseif isa(dataBlock, 'double')
    
    %generate a TIC
    ticData = sum(dataBlock, 2);
    
    %I am lazy so just make the exact same structure.
    dataBlock.tic = ticData;
    dataBlock.specData = dataBlock;

    %%%%NEEED TO FIX%%%%
    %The scan rate of the instrument needs to be known for other data
    %perhaps require the user to input data as a structure...
    %for now assume the datarate is 50Hz for a quad (pretty quick)
    dataBlock.dataRate = 50;
    
    %do the same as before
    tic_ffted = fft(dataBlock.tic);
    for page = 1:size(dataBlock.specdata,2)
        specdata_ffted(:,page) = fft(dataBlock.specdata(:,page));
    end

end

N = length(dataBlock.tic);
fs = dataBlock.dataRate;
X_mags = abs(tic_ffted);
fax_bins = [0 : N-1];
fax_hz = fax_bins*fs/N; %frequency axis in Hz
N_2 = floor(N/2);

if ~exist('keepResult','var')
    f = figure('Name','userFig');
    plot(fax_hz(1:N_2), X_mags(1:N_2));
    xlim([0 fax_hz(N_2)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Single-sided Magnitude spectrum');
    axis tight
    mess = msgbox('Please review the magnitude spectrum plot, then press enter when done.');
    pause();
    if isvalid(mess)
        close(mess);
    end
end
    
if ~exist('cutOut','var') || cutOut == 1
    if ~exist('lowNotch','var') || ~exist('highNotch','var')
        defaultHz = {'0','55'};
        if exist('lowNotch','var') defaultHz{1} = char(string(lowNotch)); end
        if exist('highNotch','var') defaultHz{2} = char(string(highNotch)); end
        notchType = questdlg('Do you wish to surround noise or signal?','Notch noise or signal?','Surround noise','Surround signal','Surround signal');
        notchPts = inputdlg({'Low notch freq (Hz):','High notch freq (Hz):'},'Define notch margins',[1 35],defaultHz);
        try
            lowNotch = str2num(notchPts{1});
            highNotch = str2num(notchPts{2});
        catch
            error('Notch limits must be numeric values');
        end
    end
    
    if isempty(lowNotch) || lowNotch < min(fax_hz)
        lowNotch = min(fax_hz);
        lowidx = min(find(round(fax_hz,3) == lowNotch));
    else
        lowidx = min(find(round(fax_hz,3) == lowNotch));
    end
    
    if isempty(highNotch) || highNotch <= lowNotch || highNotch > max(fax_hz)
        highNotch = max(fax_hz);
        highidx = max(find(round(fax_hz,3) == highNotch));
    else
        highidx= max(find(round(fax_hz,3) == highNotch));
    end
else
    lowidx = 1;
    highidx = length(fax_bins);
    notchType = 'Surround Signal';
end

close(findobj('type','figure','name','userFig'));

if exist('notchType','var') && strcmp(notchType, 'Surround signal')
    tic_filtered = zeros(size(tic_ffted));
    tic_filtered(lowidx:highidx) = tic_ffted(lowidx:highidx);
    specdata_filtered = zeros([size(dataBlock.specdata)]);
    specdata_filtered(lowidx:highidx,:) = specdata_ffted(lowidx:highidx,:);
else
    tic_filtered = zeros(size(tic_ffted));
    tic_filtered(1:lowidx) = tic_ffted(1:lowidx);
    tic_filtered(highidx:end) = tic_ffted(highidx:end);
    specdata_filtered = zeros([size(dataBlock.specdata)]);
    specdata_filtered(1:lowidx,:) = specdata_ffted(1:lowidx,:);
    specdata_filtered(highidx:end,:) = specdata_ffted(highidx:end,:);
end

if ~exist ('keepResult','var')
    f = figure('Name','userFig');
    ax1 = subplot(3,1,1);
    plot(dataBlock.tic);
    title('TIC, Pre-FFT Denoise')
    xlabel('Retention (s)')
    ylabel('Abundance')

    ax2 = subplot(3,1,2);
    plot(abs(ifft(tic_filtered)));
    title(strcat('TIC, Post-FFT Denoise (notched from ',num2str(fax_hz(lowidx)),' to ',num2str(fax_hz(highidx)),' Hz)'))
    xlabel('Retention (s)')
    ylabel('Abundance')
    
    linkaxes([ax1,ax2],'xy');
    
    subplot(3,1,3);
    X_mags = abs(tic_filtered);
    plot(fax_hz(1:N_2), X_mags(1:N_2));
    xlim([0 fax_hz(N_2)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Single-sided Magnitude spectrum');
    axis tight
    
    keepResult = input('Do you wish to keep the result? Y/N [N] - ','s');

    if ~isempty(findobj(f)) close(f); end
end

if isempty(keepResult) == 0 && (keepResult == 'Y' || keepResult == 'y' || keepResult == 1)
    dataBlock.tic = abs(ifft(tic_filtered));
    clearvars dataBlock.specdata specdata_ffted %attempt to reduce memory usage
    dataBlock.specdata = abs(ifft(specdata_filtered));
    dataBlock.denoised = 'FFT';
end

end