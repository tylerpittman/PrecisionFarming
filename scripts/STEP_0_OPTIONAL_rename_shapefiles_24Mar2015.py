# This renames every shapefile in a specified directory
# 24 March 2015
# python -v /Users/tylerpittman/Farm/STEP_0_OPTIONAL_rename_shapefiles_24Mar2015.py
import re, glob, os

def renamer(files, pattern, replacement):
	for pathname in glob.glob(files):
		basename= os.path.basename(pathname)
		new_filename= re.sub(pattern, replacement, basename)
	if new_filename != basename:
		os.rename(pathname, os.path.join(os.path.dirname(pathname), new_filename))
			  
##########################################################################################

ws = "/Users/tylerpittman/Farm/input"
os.chdir(ws)
# batch rename all files in directory to truncate after _
renamer("*", '(?<=_2015_boundary)\w+',  r"2014_boundary")

ws = "/Users/tylerpittman/Farm/input"
os.chdir(ws)
# batch rename all files in directory to truncate before .
renamer("*", '\.',  r"_2014_yield.")

'''
Place the field boundaries and combined yield data (after converting from Greenstar to .shp with Apex software)
to the input directory for STEP 2.  Should have filenames similar to: 900001_2015_yield.shp and 900001_2015_boundary.shp
'''

