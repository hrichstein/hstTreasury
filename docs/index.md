# List of included files

## Data
### getHor1data.sh
Run this file to download the MAST data set for Horologium I into the data directory. This will be in a nested structure, but the orgFLC.py script will clean this up.

To get other targets, go to MAST, select what you need, and copy the API script with the "curl" commands. Paste those in the terminal while in the data/ directory.

## Scripts
Order to run:
1. orgFLC.py
2. alignGaia.py
3. drizIt.py

### alignGaia.py
Aligns the FLCs to Gaia. Probably an unnecessary step, but I added it at some point because my WCS values were a bit off.

### drizIt.py
Drizzles the images with the customized configuration file. Read the Drizzlepac documentation in case you need to make your own tweaks. Notably, the final_pixfrac is 0.8 and the final_scale is 0.035 (arcsec/pixel).

### orgFLC.py
Moves the FLC files from the MAST downloaded folder in the data directory to new subfolders for the specific target.