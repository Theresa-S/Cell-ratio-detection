//
// Installation instructions:
// * Activate the CLIJ and CLIJ2 update sites in Fiji
//
//
// Configuration

// where are the files located?
foldername = "C:/structure/data/Theresa/";

// smoothing parameter, needs to be optimized
spotDetection_GaussianBlur_sigma_GFP= 1;
spotDetection_GaussianBlur_sigma_DAPI= 1;

// local grey value threshold, needs to be optimized
spotDetection_findMaxima_noise_GFP = 600;
spotDetection_findMaxima_noise_DAPI = 200;
// this value might be changing from image to image, it is the local height a maximum 
// must have to be considered as maximum

// size of the tiles to be analysed in pixels
// choose according to your research question and image size
tileSize = 256;

// minimal number of detected cells to be included into evaluation
cutOff_GFP = 2;
cutOff_DAPI = 2;

// do you want to see it run? True/False
entertain = true;

run("Close All");
run("CLIJ2 Macro Extensions", "cl_device=");

filelist = getFileList(foldername);

for (i = 0; i < lengthOf(filelist); i++) {
    if (endsWith(filelist[i], ".czi")) { 
        filename = foldername + File.separator + filelist[i];
		run("Bio-Formats", "open=" + filename + " color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2");
		
		input_stack = getTitle();
		Ext.CLIJ2_clear();
		Ext.CLIJ2_push(input_stack);
		Ext.CLIJ2_getDimensions(input_stack, input_width, input_height, input_depth);

		Ext.CLIJ2_create2D(original1, input_width / tileSize + 1, input_height / tileSize + 1, 32); 
		Ext.CLIJ2_create2D(original2, input_width / tileSize + 1, input_height / tileSize + 1, 32); 
		Ext.CLIJ2_create2D(spotCount1, input_width / tileSize + 1, input_height / tileSize + 1, 32); 
		Ext.CLIJ2_create2D(spotCount2, input_width / tileSize + 1, input_height / tileSize + 1, 32); 
		Ext.CLIJ2_create2D(ratio, input_width / tileSize + 1, input_height / tileSize + 1, 32); 

		Ext.CLIJ2_getDimensions(ratio, width, height, depth);

		original1_array = newArray(width);
		original2_array = newArray(width);
		spotCount1_array = newArray(width);
		spotCount2_array = newArray(width);	
		ratio_array = newArray(width);
			
		for (y = 0; y < height; y++) {	
			
				
			for (x = 0; x < width; x++) {
				if (entertain) {
					makeRectangle(x * tileSize, y * tileSize, tileSize, tileSize);
				}
				
				Ext.CLIJ2_crop3D(input_stack, cropped_stack, x * tileSize, y * tileSize, 0, tileSize, tileSize, 2);
				Ext.CLIJ2_copySlice(cropped_stack, channel1, 0);
				Ext.CLIJ2_copySlice(cropped_stack, channel2, 1);

				// pixel statistics
				Ext.CLIJ2_getMeanOfAllPixels(channel1, mean_channel1);
				Ext.CLIJ2_getMeanOfAllPixels(channel2, mean_channel2);

				// count spots
				s1 = detectSpots(channel1, spotDetection_GaussianBlur_sigma_GFP,  spotDetection_findMaxima_noise_GFP);
				s2 = detectSpots(channel2, spotDetection_GaussianBlur_sigma_DAPI, spotDetection_findMaxima_noise_DAPI);


				// calculate ratio if spots were found
				s1 = 0;
				s2 = 0;
				r = 0;
				
				if (s2 > cutOff_DAPI && s1 > cutOff_GFP) {
					r = 1.0 * s1 / s2;
				}

				// 
				original1_array[x] = mean_channel1;
				original2_array[x] = mean_channel2;
				spotCount1_array[x] = s1;
				spotCount2_array[x] = s2;
				ratio_array[x] = r;

			}

			// 
			Ext.CLIJ2_pushArray(original1_stripe, original1_array, width, 1, 1);
			Ext.CLIJ2_pushArray(original2_stripe, original2_array, width, 1, 1);
			Ext.CLIJ2_pushArray(spotCount1_stripe, spotCount1_array, width, 1, 1);
			Ext.CLIJ2_pushArray(spotCount2_stripe, spotCount2_array, width, 1, 1);
			Ext.CLIJ2_pushArray(ratio_stripe, ratio_array, width, 1, 1);

			// 
			Ext.CLIJ2_paste2D(original1_stripe, original1, 0, y);
			Ext.CLIJ2_paste2D(original2_stripe, original2, 0, y);
			Ext.CLIJ2_paste2D(spotCount1_stripe, spotCount1, 0, y);
			Ext.CLIJ2_paste2D(spotCount2_stripe, spotCount2, 0, y);
			Ext.CLIJ2_paste2D(ratio_stripe, ratio, 0, y);	

			// break;
			showProgress(y, height);
		}

    	// 
		Ext.CLIJ2_create3D(result_stack, width, height, 5, 16);
		Ext.CLIJ2_copySlice(original1_stripe, result_stack, 0);
		Ext.CLIJ2_copySlice(original2_stripe, result_stack, 1);
		Ext.CLIJ2_copySlice(spotCount1_stripe, result_stack, 2);
		Ext.CLIJ2_copySlice(spotCount2_stripe, result_stack, 3);
		Ext.CLIJ2_copySlice(ratio_stripe, result_stack, 4);

		// 
		Ext.CLIJ2_pull(result_stack);
    }
        // save result
	dot = indexOf(input_stack, "."); 
	title = substring(input_stack, 0, dot); 
	saveAs(".tiff", foldername + "/" + title + ".tif");
}
	run("Close All");
	
	
function detectSpots(image, sigma, noise) {
	setBatchMode(true);
	Ext.CLIJ2_gaussianBlur2D(image, blurred, sigma, sigma);
	Ext.CLIJ2_pull(blurred);
	run("Find Maxima...", "noise=" + noise + " output=[Single Points]");
	spots = getTitle();
	Ext.CLIJ2_push(spots);
	close();
	close();
	Ext.CLIJ2_getSumOfAllPixels(spots, spot_count);
	setBatchMode(false);
	return spot_count / 255;	
}

	
