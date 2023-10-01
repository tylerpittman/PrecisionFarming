# This renames every shapefile in a specified directory
# 5 March 2016
# python /Users/tylerpittman/Farm/STEP_5_OPTIONAL_rename_shapefile_directory_5March2016.py
import glob, os

ws = "/Users/tylerpittman/Farm/contour_yields_ll"
os.chdir(ws)
for fn in glob.glob('*.shp'):
	k = fn.split('_')
	#print k[0]
	newName = k[0] + ".shp"
	print newName
	os.rename(fn, newName)
for fn in glob.glob('*.shx'):
	k = fn.split('_')
	#print k[0]
	newName = k[0] + ".shx"
	print newName
	os.rename(fn, newName)
for fn in glob.glob('*.dbf'):
	k = fn.split('_')
	#print k[0]
	newName = k[0] + ".dbf"
	print newName
	os.rename(fn, newName)		