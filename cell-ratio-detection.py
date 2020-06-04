# This script takes a folder of CZI image files containing two channels. It detects
# spots (ideally nuclei) in both channels and generates a relative map of number of 
# spots s1 and s2, for channel 1 and 2 respectively. The pixel value in the 
# resulting map corresponds to the ratio r = s1/s2. Finally, it generates a new image 
# containing 5 channels: the original channels 1 and 2 downsampled, a channel 
# representing s1, s2 and r. The resulting image is stored in the same folder with
# a different file ending.
# 
# 
# Copyright 2018, Robert Haase, MPI-CBG Dresden, rhaase@mpi-cbg.de
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this 
#    list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors 
#    may be used to endorse or promote products derived from this software without 
#    specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
#
####################################################################################
# 
# Configuration

# where are the files located?
foldername = "//insert_your_folderpath_here/take_care_to_use_slashes/no_backslashes/no_space_characters"

# smoothing parameter, needs to be optimized
spotDetection_GaussianBlur_sigma_GFP= 1;
spotDetection_GaussianBlur_sigma_DAPI= 1;

# local grey value threshold, needs to be optimized
spotDetection_findMaxima_noise_GFP = 600;
spotDetection_findMaxima_noise_DAPI = 200;
# this value might be changing from image to image, it is the local height a maximum 
# must have to be considered as maximum

# size of the tiles to be analysed in pixels
# choose according to your research question and image size
tileSize = 128;

#minimal number of detected cells to be included into evaluation
cutOff_GFP = 2;
cutOff_DAPI = 2;

#do you want to see it run? True/False
entertain = False;

####################################################################################


from ij import IJ
from ij.gui import Roi
from ij.gui import NewImage
from ij.plugin import Duplicator
from ij.plugin import RGBStackMerge

from java.io import File

def main():

	# go through directory and open all CZI files
	directory = File(foldername)
	for file in directory.listFiles():
		if (file.toString().endswith("czi")):
			processImage(file.toString());
	
def processImage(filename):
	# many zeiss images have a pyramid format. In order to speed up the evaluation, we chose only the second highest resolution (=series_2). 
	# depending on your computer/workstation this parameter can be optimized (i.e. chose series_1 if you have a fast computer or series_3 for slow ones) 
	IJ.run("Bio-Formats", "open=" + filename + " color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2");
	imp = IJ.getImage()
	imp.show();
	if (not entertain):
		imp.show();

	# preparing memory for saving the results
	original1 = NewImage.createFloatImage("Original 1", imp.getWidth() / tileSize + 1, imp.getHeight() / tileSize + 1, 1, NewImage.FILL_BLACK); 
	original1Map = original1.getProcessor();

	original2 = NewImage.createFloatImage("Original 2", imp.getWidth() / tileSize + 1, imp.getHeight() / tileSize + 1, 1, NewImage.FILL_BLACK); 
	original2Map = original2.getProcessor();
	
	spotCount1 = NewImage.createFloatImage("Spot count 1", imp.getWidth() / tileSize + 1, imp.getHeight() / tileSize + 1, 1, NewImage.FILL_BLACK); 
	spotCount1Map = spotCount1.getProcessor();

	spotCount2 = NewImage.createFloatImage("Spot count 2", imp.getWidth() / tileSize + 1, imp.getHeight() / tileSize + 1, 1, NewImage.FILL_BLACK); 
	spotCount2Map = spotCount2.getProcessor();
	
	ratio = NewImage.createFloatImage("Ratio", imp.getWidth() / tileSize + 1, imp.getHeight() / tileSize + 1, 1, NewImage.FILL_BLACK); 
	ratioMap = ratio.getProcessor();
	
	if (entertain):
		ratio.show();

	# go through all tiles
	for x in range(0, ratioMap.getWidth()-1):
		for y in range(0, ratioMap.getHeight()-1):
			# crop out the tile from the original images
			imp.setRoi(Roi(x * tileSize, y * tileSize, tileSize, tileSize));
			channel1 = Duplicator().run(imp, 1, 1, 1, 1, 1, 1);
			channel2 = Duplicator().run(imp, 2, 2, 1, 1, 1, 1);
		
			# spot detection
			spots1 = detectSpots(channel1, spotDetection_GaussianBlur_sigma_GFP,  spotDetection_findMaxima_noise_GFP);
			spots2 = detectSpots(channel2, spotDetection_GaussianBlur_sigma_DAPI, spotDetection_findMaxima_noise_DAPI);
			
			# pixel statistics
			statistics1 = channel1.getStatistics()
			statistics2 = channel2.getStatistics()
		
			# calculate ratio if spots were found
			s1 = 0;
			s2 = 0;
			r = 0;
			if (spots1 is not None and spots2 is not None):
				fp1 = spots1.getFloatPolygon();
				fp2 = spots2.getFloatPolygon();	
				s1 = fp1.npoints;
				s2 = fp2.npoints;
				if (s2 > cutOff_DAPI and s1 > cutOff_GFP):
					r = 1.0 * s1 / s2;

			# fill result memory
			original1Map.setf(x, y, statistics1.mean);
			original2Map.setf(x, y, statistics2.mean);
			spotCount1Map.setf(x, y, s1);
			spotCount2Map.setf(x, y, s2);
			ratioMap.setf(x, y, r);
			
		if (entertain):
			# show current result image
			ratio.updateAndDraw();
			IJ.run(ratio, "Enhance Contrast", "saturated=0.35");
	
	# put all results image channels together to one image
	images = [];
	images.append(original1);
	images.append(original2);
	images.append(spotCount1);
	images.append(spotCount2);
	images.append(ratio);
	resultMap = RGBStackMerge.mergeChannels(images, False);

	# fix pixel size
	# factor is multiplied by 2 because ImageJ seems to have a problem when using .czi file series of lower resolution (i.e. series_2); please check for individual cases!
	factor = (imp.getWidth() / resultMap.getWidth()) * 2;
	IJ.run(resultMap, "Properties...", "channels=5 slices=1 frames=1 unit=" + imp.getCalibration().getUnit() + " pixel_width=" + str(imp.getCalibration().pixelWidth * factor) + " pixel_height=" + str(imp.getCalibration().pixelHeight * factor) + " voxel_depth=1.0000");
	IJ.run(ratio, "Properties...", "channels=1 slices=1 frames=1 unit=" + imp.getCalibration().getUnit() + " pixel_width=" + str(imp.getCalibration().pixelWidth * factor) + " pixel_height=" + str(imp.getCalibration().pixelHeight * factor) + " voxel_depth=1.0000");
	
	# visualisation
	resultMap.setC(1);
	IJ.run(resultMap, "Green", "");
	IJ.run(resultMap, "Enhance Contrast", "saturated=0.35");
	resultMap.setC(2);
	IJ.run(resultMap, "Blue", "");
	IJ.run(resultMap, "Enhance Contrast", "saturated=0.35");
	resultMap.setC(3);
	IJ.run(resultMap, "mpl-inferno", "");
	IJ.run(resultMap, "Enhance Contrast", "saturated=0.35");
	resultMap.setC(4);
	IJ.run(resultMap, "mpl-inferno", "");
	IJ.run(resultMap, "Enhance Contrast", "saturated=0.35");
	resultMap.setC(5);
	IJ.run(resultMap, "mpl-inferno", "");
	resultMap.show();
	IJ.resetMinAndMax(resultMap);
	resultMap.setDisplayMode(IJ.COLOR);

	IJ.run(ratio, "mpl-inferno", "");
	IJ.setMinAndMax(ratio, 0, 1);
	
	# save result
	IJ.saveAsTiff(resultMap, filename + "_suitable-name_map.tif");
	IJ.saveAsTiff(ratio, filename + "_suitable-name_ratio.tif");
	
	IJ.run("Close All", "");


# local maxima detection after Gaussian blurring an image
def detectSpots(imp, sigma, noise):
	IJ.run(imp, "Gaussian Blur...", "sigma=" + str(sigma));
	IJ.run(imp, "Find Maxima...", "noise=" + str(noise) + " output=[Point Selection]");
	return imp.getRoi();

main();
print("Bye");
