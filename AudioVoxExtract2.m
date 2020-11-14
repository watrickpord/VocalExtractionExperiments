% attempt 2 to extract centre channel from stereo mix:
% using correlations in phase over a number of STFT frames to determine
% is a source is left, centre, or right panned
while true
clear

% import audio file and split to L and R channels
[in, Fs] = audioread('Yesterday.flac');
length = size(in);
length = length(1)
left  = 0.6*in(:,1);    % reduce volume to give extra headroom
right = 0.6*in(:,2);

% calculate number of frames needed for given size
L = 512                            % number of samples per frame
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

leftPhaseMatrix = angle(leftFT);
rightPhaseMatrix = angle(rightFT);

% we want to take slices of a given length around each entry in the 
% phase matrix and make a r value matrix based on L-R correlation
rMatrix = zeros(size(leftPhaseMatrix));

sliceLen = 6
currentFrame = sliceLen/2 + 1;

while currentFrame+(sliceLen/2) <= numFrames
    % frame no. of start and end of current slice
    sliceStart = currentFrame - sliceLen/2;
    sliceEnd   = currentFrame + sliceLen/2;
    
    % calculate r values for each row in current column (i.e. frame)
    for rowIndex = 1:L
        currentR = corrcoef(leftPhaseMatrix(rowIndex, sliceStart:sliceEnd), rightPhaseMatrix(rowIndex, sliceStart:sliceEnd));
        rMatrix(rowIndex, currentFrame) = currentR(2);
    end
    
    % increment frame counter
    currentFrame = currentFrame + 1;
end

end

% hard cutoff for r value - very digital and warbly
%leftProcFT = leftFT.*(rMatrix>0.5);
%rightProcFT = rightFT.*(rMatrix>0.5);

% guassian cutoff as r falls away from 1
sigma = 0.25
gainFunction = @(r) exp(-(r-1)^2/(2*sigma^2))/(sqrt(2*pi)*sigma);
gainMatrix = arrayfun(gainFunction, rMatrix);

% smooth gain matrix
sigmaSmooth = 4
kernel = exp(-(-30:30).^2/(2*sigmaSmooth^2))/(sigmaSmooth*sqrt(2*pi))
smoothGainMatrix = conv2(gainMatrix,kernel,'same');

leftProcFT  = leftFT.*smoothGainMatrix;
rightProcFT = rightFT.*smoothGainMatrix;


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

