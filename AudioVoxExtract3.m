% attempt 3 to extract centre channel from stereo mix:
% go frame by frame and use euclidian distance between left and right
% coeffecients as a measure of their similarity, apply higher gain to more
% similar coeffecients (i.e. as per attempt except using euclidean distance
% between complex coeffecients rather than just phase difference)

clear

% import audio file and split to L and R channels
[in, Fs] = audioread('Yesterday.flac');
length = size(in);
length = length(1)
left  = 0.5*in(:,1);    % reduce volume to give extra headroom
right = 0.5*in(:,2);

% calculate number of frames needed for given size
L = 2048                            % number of samples per frame
T_L = L/44.1                        % frame length (mS)
numFrames = ceil(length/(L/2))-1    % number of frames (with 50% overlap)

% pad audio tracks with zeros to make length multiple of L/2
paddedLength = (numFrames+1)*(L/2);
extraSamples = paddedLength-length;
left  = [left;  zeros(extraSamples,1)];
right = [right; zeros(extraSamples,1)];

% define hann window function
hann = transpose((sin(pi.*(0:L-1)/(L-1))).^2);

% empty L x numFrames arrays to store STFT vectors
leftFT  = zeros(L,numFrames);
rightFT = zeros(L,numFrames);

% calculate hann windowed STFT and store in vectors
for index = 1:numFrames
    % take slices for window
    startSample = (index-1)*(L/2) + 1;
    endSample   = startSample + L-1;
    
    leftSlice  = left(startSample:endSample).*hann;
    rightSlice = right(startSample:endSample).*hann;
    
    % calculate ffts and store in vectors
    currentFFT = fft(leftSlice);
    leftFT(:, index) = currentFFT;
    
    currentFFT = fft(rightSlice);
    rightFT(:, index) = currentFFT;
end

% -------------- audio processing goes here --------------

% matrix of Euclidean distance between left and right coeffecients
differenceMatrix = leftFT - rightFT;
sumMatrix        = leftFT + rightFT;
normDiffs = abs(differenceMatrix);
normSums  = abs(sumMatrix);
normMatrix = normDiffs./normSums;

% get gain from gaussian function (mean = 0) applied to norms
%sigmaGain = 0.25
%gainMatrix = exp(-normMatrix.^2/(2*sigmaGain^2))/(sigmaGain*sqrt(2*pi));

% frequency dependent sigma - we want sigma to be higher in the range of
% the vocals and close to zero elsewhere (i.e. at LF)

% sigma as a function of frequency (1st order high pass)
sigmaGain = 0.3             % sigma value for f > cutoffFreq
cutoffFreq = 300           % corner freq of HPF


gainMatrix = zeros(size(normMatrix));
fftFreqs = fftshift(ceil(-L/2:L/2-1)/(1/Fs)/L); % frequency of nth fft bin

% go row by row (i.e. by freq) of normMatrix and calculate gain
% n.b. running though 1:L gives NaN for f=0
for index = 2:L
    % calculate gain vs norm gaussian for the current freq
    currentFreq = abs(fftFreqs(index));
    
    % old function to calculate sigma for current freq
    while false
    if currentFreq >= cutoffFreq
        currentSigma = sigmaGain;
    else
        currentSigma = sigmaGain*currentFreq/cutoffFreq;
    end
    end
    
    % new function to calculate sigma for current freq
    if currentFreq >= cutoffFreq
        currentSigma = sigmaGain;
    elseif currentFreq >= cutoffFreq/2
        currentSigma = 2*sigmaGain*currentFreq/cutoffFreq;
    else 
        currentSigma = 0.00000000000001;  % using 0.0 gives NaNs
    end
    
    sigmaFun = @(norm) exp(-norm^2/(2*currentSigma^2))/(currentSigma*sqrt(2*pi));
    %sigmaFun = @(norm) exp(-norm^2/(2*sigmaGain^2))/(sigmaGain*sqrt(2*pi));             % test
    
    % apply gain vs norm function to norms
    currentNormsRow = normMatrix(index, :);
    currentGainsRow = arrayfun(sigmaFun, currentNormsRow);
    %currentGainsRow = exp(-currentNormsRow.^2/(2*sigmaGain^2))/(sigmaGain*sqrt(2*pi));    % test
    
    % save results in gain matrix
    gainMatrix(index, :) = currentGainsRow;
end


% gain matrix smoothing parameters
sigmaMs = 75               % st. dev. of smoothing kernel in mS
sigmaFrames = sigmaMs/T_L;
kernelLen = 50;

% smooth gain matrix
kernel = exp(-(-kernelLen:kernelLen).^2/(2*sigmaFrames^2))/(sigmaFrames*sqrt(2*pi));
gainMatrix = conv2(gainMatrix,kernel,'same');

% copy fourier coeffecients with gain multiplier
leftProcFT  = gainMatrix.*leftFT;
rightProcFT = gainMatrix.*rightFT;

% sanity check short circuit
%leftProcFT  = leftFT;
%rightProcFT = rightFT;

% ----------------- end audio processing -----------------

% inverse FT

% empty arrays for output audio
leftProc  = zeros(paddedLength, 1);
rightProc = zeros(paddedLength, 1);

for index = 1:numFrames
    % inverse ffts for indexed frame
    leftFrame  = ifft(leftFT(:, index));
    rightFrame = ifft(rightFT(:, index));
    
    % processesed audio outputs
    leftProcFrame = ifft(leftProcFT(:, index));
    rightProcFrame = ifft(rightProcFT(:, index));
    
    % add the frame audio to main stream
    startSample = (index-1)*(L/2) + 1;
    endSample   = startSample + L - 1;
    
    leftProc(startSample:endSample) = leftProc(startSample:endSample) + leftProcFrame;
    rightProc(startSample:endSample) = rightProc(startSample:endSample) + rightProcFrame;
end

% sum left and right to mono
mono = left + right;
monoProc = leftProc + rightProc;

% original players
lp = audioplayer(left, Fs);
rp = audioplayer(right, Fs);
mp = audioplayer(mono, Fs);

% processed players
lpp = audioplayer(leftProc, Fs);
rpp = audioplayer(rightProc, Fs);
mpp = audioplayer(monoProc, Fs);

