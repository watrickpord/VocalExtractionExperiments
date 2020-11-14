# Vocal Extraction Experiments
Experiments with extracting centre panned vocals from a stereo mix.

AudioVoxExtract - short-term fourier transform is taken of left and right channel, phase difference between each left and right component used to determine gain (i.e. closer to 0 degree the higher the gain) with gaussian gain vs. phase difference centered around zero.

AudioVoxExtract2 - using correlation in phase between slices to determine gain. Theory is that if a particular frequency is dominated by a centre panned sound then there will be correlation in time between phases of left and right components at that frequency. Thus, slices are taken of several frames and these create an r value for each pair of left and right channels, which then determines gain (closer to 1 = centre panned = higher gain).
