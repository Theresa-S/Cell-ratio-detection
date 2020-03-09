# Cell-ratio-detection

This is a skript for evaluation of fluorescence microscopy images. It counts cells in two channels and calculates the ratio of both. 
Both channels are processed independently on a grid with a defined tile size. A Gaussian blur is applied, then cell numbers are counted with a Maxima Finder.
As a result you receive two images, one with said ratio and one stack containing the fluorescence intensities and cell numbers of each channel plus the ratio.

It is highly recommended to optimize the parameters (as indicated in the code) and validate the results in randomly chosen regions.
