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
L = 4096                            % number of samples per frame
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
sigma = 0.25
gainMatrix = exp(-normMatrix.^2/(2*sigma^2))/(sigma*sqrt(2*pi));

% gain matrix smoothing
sigmaMs = 150               % st. dev. of guassian kernel in mS
sigmaFrames = sigmaMs/T_L;
kernelLen = 50;

kernel = exp(-(-kernelLen:kernelLen).^2/(2*sigmaFrames^2))/(sigmaFrames*sqrt(2*pi));
gainMatrix = conv2(gainMatrix,kernel,'same');

% copy fourier coeffecients with gain multiplier
leftProcFT  = gainMatrix.*leftFT;
rightProcFT = gainMatrix.*rightFT;

% ----------------- end audio processing -----------------

% inverse FT

% empty arrays for output audio
leftNew = zeros(paddedLength, 1);
rightNew = zeros(paddedLength, 1);

leftProcNew = zeros(paddedLength, 1);
rightProcNew = zeros(paddedLength, 1);

for index = 1:numFrames
    % inverse ffts for indexed frame
    leftFrame  = ifft(leftFT(:, index));
    rightFrame = ifft(rightFT(:, index));
    
    % processesed audio outputs
    leftProcFrame = ifft(leftProcFT(:, index));
    rightProcFrame = ifft(rightProcFT(:, index));
    
    % add the frame audio to main part
    startSample = (index-1)*(L/2) + 1;
    endSample   = startSample + L - 1;
    
    leftNew(startSample:endSample) = leftNew(startSample:endSample) + leftFrame;
    rightNew(startSample:endSample) = rightNew(startSample:endSample) + rightFrame;

    leftProcNew(startSample:endSample) = leftProcNew(startSample:endSample) + leftProcFrame;
    rightProcNew(startSample:endSample) = rightProcNew(startSample:endSample) + rightProcFrame;
end

% sum left and right to mono
mono = left + right;
monoPhaseNew = leftProcNew + rightProcNew;

% original players
lp = audioplayer(left, Fs);
rp = audioplayer(right, Fs);
mp = audioplayer(mono, Fs);

% processed players
lpp = audioplayer(leftProcNew, Fs);
rpp = audioplayer(rightProcNew, Fs);
mpp = audioplayer(monoPhaseNew, Fs);

