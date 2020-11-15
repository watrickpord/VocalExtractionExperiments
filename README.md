# Vocal Extraction Experiments
Experiments with extracting centre panned vocals from a stereo mix.

AudioVoxExtract.m - short-term fourier transform is taken of left and right channel, phase difference between each left and right component used to determine gain (i.e. closer to 0 degree the higher the gain) with gaussian gain vs. phase difference centered around zero.

AudioVoxExtract2.m - using correlation in phase between slices to determine gain. Theory is that if a particular frequency is dominated by a centre panned sound then there will be correlation in time between phases of left and right components at that frequency. Thus, slices are taken of several frames and these create an r value for each pair of left and right channels, which then determines gain (closer to 1 = centre panned = higher gain).

AudioVoxExtract3.m - go frame by frame and use euclidian distance between left and right coeffecients as a measure of their similarity, apply higher gain to more similar coeffecients (i.e. as per attempt except using pythagorean distance between complex coeffecients rather than just phase difference). Work very well compared to other two approaches. Frequency dependent sigma was added and also improves performance, once tuned in.

SumAmplitudeRandomPhaseDiff.m - plot of expected value for amplitude of two sines (of equal amplitude) added together with a random phase shift determined by a Gaussian probability distrobution - i.e. plot of sum amplitude vs. std. deviation of the phase offset. This is probably related to amplitudes of coherent vs. incoherent summing (6 dB vs 3 dB respectively), but can't figure out how.
