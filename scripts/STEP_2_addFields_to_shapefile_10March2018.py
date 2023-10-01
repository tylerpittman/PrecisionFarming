#!/usr/bin/env python
# This adds prescription farming fields to every shapefile in a specified directory and shortens filename
# 10 March 2018, Tyler Pittman
# python -tt /Users/tylerpittman/Farm/STEP_2_addFields_to_shapefile_10March2018.py
import re, glob, os
from osgeo import ogr, gdal

driver = ogr.GetDriverByName('ESRI Shapefile')

def renamer(files, pattern, replacement):
	for pathname in glob.glob(files):
        	basename= os.path.basename(pathname)
        	new_filename= re.sub(pattern, replacement, basename)
	if new_filename != basename:
            	os.rename(pathname, os.path.join(os.path.dirname(pathname), new_filename))
			  
##########################################################################################

ws = "/Users/tylerpittman/Farm/contour_yields_ll"
os.chdir(ws)

# batch rename all files in directory to truncate after _2017
renamer("*", '(?<=_2018)\w+',  r"")

for fn in os.listdir(ws):
	if fn.endswith("2018.shp") and os.path.isfile(fn):
        	print fn

for i in os.listdir(ws):
	if i.endswith("2018.shp")  and os.path.isfile(i): 
		dataSource = driver.Open(i, 1) #1 is read/write

		#define floating point field:
		fldDef = ogr.FieldDefn('seedTyp', ogr.OFTString)
		fldDef1 = ogr.FieldDefn('seedRt', ogr.OFTReal)
		fldDef2 = ogr.FieldDefn('inocTyp', ogr.OFTString)
		fldDef3 = ogr.FieldDefn('inocRt', ogr.OFTReal)
		fldDef4 = ogr.FieldDefn('fertPTyp', ogr.OFTString)
		fldDef5 = ogr.FieldDefn('fertPRt', ogr.OFTReal)
		fldDef6 = ogr.FieldDefn('fertMTyp', ogr.OFTString)
		fldDef7 = ogr.FieldDefn('fertMRt', ogr.OFTReal)
		fldDef8 = ogr.FieldDefn('fertSTyp', ogr.OFTString)
		fldDef9 = ogr.FieldDefn('fertSRt', ogr.OFTReal)
		#fldDef.SetWidth(16) #16 char string width

		#get layer and add these fields:
		layer = dataSource.GetLayer()
		layer.CreateField(fldDef)
		layer.CreateField(fldDef1)
		layer.CreateField(fldDef2)
		layer.CreateField(fldDef3)
		layer.CreateField(fldDef4)
		layer.CreateField(fldDef5)
		layer.CreateField(fldDef6)
		layer.CreateField(fldDef7)
		layer.CreateField(fldDef8)
		layer.CreateField(fldDef9)

		'''
		for j in range(len(layer)):
			print j	
		'''
		for k in range(len(layer)):
			feature = layer.GetFeature(k)
			
			if feature.GetField('CROP') == 'Wheat':
				seedRt = 80
				inocRt = None
				fertPRt = 75
				fertMRt = 175
				fertSRt = None
				inc1Fert = 5
				seedTyp = 'AC Brigade'
				inocTyp = None	
				fertPTyp = '8-38-16-0'
				fertMTyp = '46-0-0-0'
				fertSTyp = None			
			elif feature.GetField('CROP') == 'Lentils':
				seedRt = 100
				inocRt = 2.8
				fertPRt = 50
				fertMRt = 40
				fertSRt = None
				inc1Fert = 5
				seedTyp = 'CDC Greenstar'
				inocTyp = 'TagTeam'
				fertPTyp = '8-38-16-0'
				fertMTyp = '8-38-16-0'
				fertSTyp = None	
			elif feature.GetField('CROP') == 'Barley':
				seedRt = 75
				inocRt = None
				fertPRt = 75
				fertMRt = 110
				fertSRt = None
				inc1Fert = 5
				seedTyp = 'CDC Copeland'
				inocTyp = None	
				fertPTyp = '8-38-16-0'
				fertMTyp = '46-0-0-0'
				fertSTyp = None	
			elif feature.GetField('CROP') == 'Canary':
				seedRt = 28
				inocRt = None
				fertPRt = 75
				fertMRt = 125
				fertSRt = None
				inc1Fert = 5
				seedTyp = 'CDC Calvi'
				inocTyp = None
				fertPTyp = '8-38-16-0'
				fertMTyp = '46-0-0-0'
				fertSTyp = None	
			else:
				seedRt = 4.6
				inocRt = None
				fertPRt = 50
				fertMRt = 180
				fertSRt = 70
				inc1Fert = 5
				seedTyp = 'InVigor L140P'
				inocTyp = None
				fertPTyp = '15-30-0-0'
				fertMTyp = '46-0-0-0'
				fertSTyp = '21-0-0-4'

			if feature.GetField('COACH') == 's1':
				fMRt = fertMRt + 0.5*inc1Fert
			elif feature.GetField('COACH') == 's2':
				fMRt = fertMRt + inc1Fert
			elif feature.GetField('COACH') == 'a2':
				fMRt = fertMRt - inc1Fert
			elif feature.GetField('COACH') == 'a1':
				fMRt = fertMRt - 0.5*inc1Fert
			else:
				fMRt = fertMRt 
			
			feature.SetField('seedTyp', seedTyp)
			feature.SetField('seedRt', seedRt)
			feature.SetField('inocTyp', inocTyp)
			feature.SetField('inocRt', inocRt)
			feature.SetField('fertPTyp', fertPTyp)
			feature.SetField('fertPRt', fertPRt)
			feature.SetField('fertMTyp', fertMTyp)
			feature.SetField('fertMRt', fMRt)
			feature.SetField('fertSTyp', fertSTyp)
			feature.SetField('fertSRt', fertSRt)

			layer.SetFeature(feature)

		# Close the Shapefile
		source = None

##########################################################################################
