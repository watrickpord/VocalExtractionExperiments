# Vocal Extraction Experiments
Experiments with extracting centre panned vocals from a stereo mix.

AudioVoxExtract - short-term fourier transform is taken of left and right channel, phase difference between each left and right component used to determine gain (i.e. closer to 0 degree the higher the gain) with gaussian gain vs. phase difference centered around zero.
