% attempt to extract centre channel from stereo mix

clear all

% import audio file and split to L and R channels
[in, Fs] = audioread('Yesterday.flac');
Fs
length = size(in);
length = length(1)
left  = 0.9*in(:,1);    % reduce volume to give extra headroom
right = 0.9*in(:,2);

% calculate number of frames needed for given size
L = 2048                         % number of samples per frame
numFrames = ceil(length/(L/2))-1 % number of frames (with 50% overlap)

% pad audio tracks with zeros to make length multiple of L/2
paddedLength = (numFrames+1)*(L/2)
extraSamples = paddedLength-length
left  = [left;  zeros(extraSamples,1)];
right = [right; zeros(extraSamples,1)];

% define hann window function
hann = transpose((sin(pi.*(0:L-1)/(L-1))).^2);

% empty L x numFrames arrays to store STFT vectors
leftSTFT  = zeros(L,numFrames);
rightSTFT = zeros(L,numFrames);

% calculate hann windowed STFT and store in vectors
for index = 1:numFrames
    % take slices for window
    startSample = (index-1)*(L/2) + 1;
    endSample   = startSample + L-1;
    
    leftSlice  = left(startSample:endSample).*hann;
    rightSlice = right(startSample:endSample).*hann;
    
    % calculate ffts and store in vectors
    currentFFT = fft(leftSlice);
    leftSTFT(:, index) = currentFFT;
    
    currentFFT = fft(rightSlice);
    rightSTFT(:, index) = currentFFT;
end

% -------------- audio processing goes here --------------

% multiplying together FTs per frame -- gives distorted mess
convFT = (0.5*leftSTFT).*(0.5*rightSTFT);

% algorithm idea - for each coeffecient, keep in output if phase is roughly
% the same between left and right channels

k = 0.8  % phase difference cutoff in radians

% matrix of phase difference between each STFT coeffecient
phaseDiffs = angle(leftSTFT) - angle(rightSTFT);

% binary array of which components to keep
indicies = abs(phaseDiffs)<k;

% only copy across components with left to right phase diff less than k
leftPhaseFT  = indicies.*leftSTFT;
rightPhaseFT = indicies.*rightSTFT;

% ----------------- end audio processing -----------------

% inverse FT

% empty arrays for output audio
leftNew = zeros(paddedLength, 1);
rightNew = zeros(paddedLength, 1);
convNew = zeros(paddedLength, 1);
leftPhaseNew = zeros(paddedLength, 1);
rightPhaseNew = zeros(paddedLength, 1);

for index = 1:numFrames
    % inverse ffts for indexed frame
    leftFrame  = ifft(leftSTFT(:, index));
    rightFrame = ifft(rightSTFT(:, index));
    
    % processesed audio outputs
    convFrame = ifft(convFT(:, index));
    leftPhaseFrame = ifft(leftPhaseFT(:, index));
    rightPhaseFrame = ifft(rightPhaseFT(:, index));
    
    % add the frame audio to main part
    startSample = (index-1)*(L/2) + 1;
    endSample   = startSample + L - 1;
    
    leftNew(startSample:endSample) = leftNew(startSample:endSample) + leftFrame;
    rightNew(startSample:endSample) = rightNew(startSample:endSample) + rightFrame;
    
    convNew(startSample:endSample) = convNew(startSample:endSample) + convFrame;
    leftPhaseNew(startSample:endSample) = leftPhaseNew(startSample:endSample) + leftPhaseFrame;
    rightPhaseNew(startSample:endSample) = rightPhaseNew(startSample:endSample) + rightPhaseFrame;
end

lp = audioplayer(leftNew, Fs);
rp = audioplayer(rightNew, Fs);
cp = audioplayer(convNew, Fs);
lpp = audioplayer(leftPhaseNew, Fs);
rpp = audioplayer(rightPhaseNew, Fs);

