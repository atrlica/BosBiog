# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "H:/BosBiog/BosBiogWorking.gdb/"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("H:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")
    print("found 30m EVI for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR
csRes = arcpy.GetRasterProperties_management(in_raster="H:/BosBiog/processed/bos.cov.V5-canisa+lulc.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

# Run the aggregate on all the 1m cover data

### Forest
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.Forestperv.tif", "H:/BosBiog/processed/bos.cov.Forestperv30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.Forestperv 30m")
##
### Dev
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.Dev-grass.tif", "H:/BosBiog/processed/bos.cov.Dev-grass30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.Dev-grass 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.Dev-barr.tif", "H:/BosBiog/processed/bos.cov.Dev-barr30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.Dev-barr 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.Dev-canPerv.tif", "H:/BosBiog/processed/bos.cov.Dev-canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.Dev-canPerv 30m")
##
### HDRes
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.HDRes-grass.tif", "H:/BosBiog/processed/bos.cov.HDRes-grass30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.HDRes-grass 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.HDRes-barr.tif", "H:/BosBiog/processed/bos.cov.HDRes-barr30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.HDRes-barr 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.HDRes-canPerv.tif", "H:/BosBiog/processed/bos.cov.HDRes-canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.HDRes-canPerv 30m")
##
### LDRes
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.LDRes-grass.tif", "H:/BosBiog/processed/bos.cov.LDRes-grass30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.LDRes-grass 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.LDRes-barr.tif", "H:/BosBiog/processed/bos.cov.LDRes-barr30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.LDRes-barr 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.LDRes-canPerv.tif", "H:/BosBiog/processed/bos.cov.LDRes-canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.LDRes-canPerv 30m")
##
### OVeg
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.OVeg-grass.tif", "H:/BosBiog/processed/bos.cov.OVeg-grass30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.OVeg-grass 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.OVeg-barr.tif", "H:/BosBiog/processed/bos.cov.OVeg-barr30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.OVeg-barr 30m")
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.OVeg-canPerv.tif", "H:/BosBiog/processed/bos.cov.OVeg-canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished cov.OVeg-canPerv 30m")

# Water
arcpy.gp.Aggregate_sa("H:/FragEVI/processed/boston/bos.water_only.tif", "H:/BosBiog/processed/bos.cov.water30m.tif", "30", "SUM", "EXPAND", "DATA")
print("finished cov.water 30m")

