# List of included files

## Data
### getHor1data.sh
Run this file to download the MAST data set for Horologium I into the data directory. This will be in a nested structure, but the orgFLC.py script will clean this up.

To get other targets, go to MAST, select what you need, and copy the API script with the "curl" commands. Paste those in the terminal while in the data/ directory.

### mastQuery.ipynb
Alternatively, if you're not a fan of using curl/would like to do things without having to go to MAST directly, this Jupyter Notebook shows how to use astroquery.mast to get the relevant data. If you want other targets for this program, simply change the "tname."

### From within the run.py...
If you use the run.py file in the scripts directory, it will run getObs.py, which will use astroquery.mast to get the observations.

## Scripts
Use run.py and the config.json file to download the images, drizzle them, and perform aperture photometry using photutils.  

The scripts currently in the run.py file are:  
#### 1. getObs.py  
Uses astroquery.mast to fetch HST FLC observations of the desired target.
#### 2. moveFLCs.py  
Moves the downloaded FLCs from their MAST folders to our functional directory structure.
#### 3. alignGaia.py  
Aligns the FLCs to Gaia. Possibly an unnecessary step, but I added it at some point because my WCS values were a bit off.
#### 4. drizIt.py  
Drizzles the images with the customized (optimized) configuration file. Read the Drizzlepac documentation in case you need to make your own tweaks. Notably, the final_pixfrac is 0.8 and the final_scale is 0.035 (arcsec/pixel). This creates the image we will now refer to as the DRC.
#### 5. maskACS.py  
Finds large sources (e.g., galaxies, saturated stars, chip gap) in the DRC and replaces the pixel values with a very large number. Moves the relevant file from the data directory to the results directory.

From within aperPhotCodes directory:  

#### 6. getZPTsky.py  
Fetches the zeropoint for the date that the image was taken. Additionally finds the sky sigma value (using code from dnidever) and outputs the zeropoint and sky sigma out to a file for later reference.  
#### 7. runPhotUtils.py  
Runs the aperture photometry routine from photutils.aperture on the images from both filers.  
#### 8. getDRCfiltRef.py  
Finds bright stars to use as references for the linear transformation from filter 2 to filter 1.  
#### 9. drcFiltLinTrans.py  
Uses the reference stars to find the 6D transformation from filter 2 to filter 1.  
#### 10. matchDRCfilt.py  
Matches sources between the two filters and creates a combined catalog (in observational space, not de-reddened).  
#### 11. makeCMD.py  
Creates a CMD from the matched DRC photometry. Also shows which sources were flagged as "likely stars" in at least one filter.
