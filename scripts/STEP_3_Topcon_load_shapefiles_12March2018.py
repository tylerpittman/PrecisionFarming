# This puts boundary and VRC shapefile data in format loadable to Topcon X30 monitor
# 12 March 2018
# python /Users/tylerpittman/Farm/STEP_3_Topcon_load_shapefiles_12March2018.py

import re, glob, os, errno, shutil

def renamer(files, pattern, replacement):
	for pathname in glob.glob(files):
		basename= os.path.basename(pathname)
		new_filename= re.sub(pattern, replacement, basename)
	if new_filename != basename:
		os.rename(pathname, os.path.join(os.path.dirname(pathname), new_filename))

ws = "/Users/tylerpittman/Farm/contour_yields_ll"
os.chdir(ws)
for fn in os.listdir(ws):
	if fn.endswith("2018.shp") and os.path.isfile(fn):
		#print fn
		k = fn.split('_')
		#print k[0]
        	newDir = os.path.join('/Users/tylerpittman/Farm/Clients/Agventure Farms Ltd/Agventure Farms Ltd', k[0], 'VRC', k[0])
		filename = os.path.join(newDir, fn)
        	directory = os.path.dirname(filename)
       		print directory
        	if not os.path.exists(directory):
                	os.makedirs(directory)
		for file in glob.glob(os.path.join(ws, k[0]+'*.*')):
			print file		
			shutil.copy2(file, os.path.join(newDir))
			renamer(directory + "/" + "*.shp", '_.*',  r".shp")
			renamer(directory + "/" + "*.shx", '_.*',  r".shx")
			renamer(directory + "/" + "*.dbf", '_.*',  r".dbf")
			myFile = os.path.basename(file)
			fileNew = os.path.join(newDir, myFile)
			print fileNew
			try:
    				os.remove(fileNew)
			except OSError:
    				pass
 			
		#m = re.search('\w+(?<=_)\w+', fn)
		#if m:
    		#	found = m.group(0)
		#	print found

ws = "/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon2017"
os.chdir(ws)
for fn in os.listdir(ws):
        if fn.endswith(".shp") and os.path.isfile(fn):
                #print fn
                k = fn.split('.')
                #print k[0]
                #newDir = os.path.join('/Users/tylerpittman/Farm/Clients/Agventure Farms Ltd/Agventure Farms Ltd', k[0], 'Boundaries', k[0])
                newDir = os.path.join('/Users/tylerpittman/Farm/Clients/Agventure Farms Ltd/Agventure Farms Ltd', k[0], 'Boundaries', 'ShapefilesForBoundary', k[0])
                #newDir = os.path.join('/Users/tylerpittman/Farm/Clients/Agventure Farms Ltd/Agventure Farms Ltd', k[0], 'Boundaries', k[0], 'ShapefilesForBoundary')
                #newDir = os.path.join('/Users/tylerpittman/Farm/Clients/Agventure Farms Ltd/Agventure Farms Ltd', k[0], 'BoundaryShapefiles')
                filename = os.path.join(newDir, fn)
                directory = os.path.dirname(filename)
                print directory
                if not os.path.exists(directory):
                        os.makedirs(directory)
                for file in glob.glob(os.path.join(ws, k[0]+'*.*')):
                        print file
                        shutil.copy2(file, os.path.join(newDir))
			#renamer(directory + "/" + "*.idx", '_.*',  r".idx")
			#renamer(directory + "/" + "*.ids", '_.*',  r".ids")
                        #renamer(directory + "/" + "*.shp", '_.*',  r".shp")
                        #renamer(directory + "/" + "*.shx", '_.*',  r".shx")
                        #renamer(directory + "/" + "*.dbf", '_.*',  r".dbf")
			#myFile = os.path.basename(file)
                        #fileNew = os.path.join(newDir, myFile)
                        #print fileNew
                        #try:
                                #os.remove(fileNew)
                        #except OSError:
                                #pass

		
'''
DO NOT DELETE EXISTING FOLDER /Users/tylerpittman/Farm/Clients/Agventure Farms Ltd
Place the field boundaries and combined yield data (after converting from Greenstar to .shp with Apex software)
to the input directory for STEP 2.  Should have VRC filenames similar to: 900001_2016_yield.shp and boundary filenames similar to 900001.shp
'''

