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
csRes = arcpy.GetRasterProperties_management(in_raster="H:/BosBiog/processed/bos.cov.V7-canisa+lulc.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

# Run the aggregate on all the 1m cover data

### Forest
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.Forest.perv.tif", "H:/BosBiog/processed/bos.cov.Forest.perv30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.Forest.perv 30m")
arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.5.Forest.grass.tif", "H:/BosBiog/processed/bos.cov.Forest.grass30m.tif", "30", "SUM", "EXPAND", "DATA")
print("finished cov.Forest.grass 30m")
arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.5.Forest.barr.tif", "H:/BosBiog/processed/bos.cov.Forest.barr30m.tif", "30", "SUM", "EXPAND", "DATA")
print("finished cov.Forest.barr 30m")
arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.5.Forest.canPerv.tif", "H:/BosBiog/processed/bos.cov.Forest.canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
print("finished cov.Forest.canPerv 30m")


### Dev
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.Dev.grass.tif", "H:/BosBiog/processed/bos.cov.Dev.grass30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.Dev.grass 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.Dev.barr.tif", "H:/BosBiog/processed/bos.cov.Dev.barr30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.Dev.barr 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.Dev.canPerv.tif", "H:/BosBiog/processed/bos.cov.Dev.canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.Dev.canPerv 30m")

### HDRes
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.HDRes.grass.tif", "H:/BosBiog/processed/bos.cov.HDRes.grass30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.HDRes.grass 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.HDRes.barr.tif", "H:/BosBiog/processed/bos.cov.HDRes.barr30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.HDRes.barr 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.HDRes.canPerv.tif", "H:/BosBiog/processed/bos.cov.HDRes.canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.HDRes.canPerv 30m")

### LDRes
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.LDRes.grass.tif", "H:/BosBiog/processed/bos.cov.LDRes.grass30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.LDRes.grass 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.LDRes.barr.tif", "H:/BosBiog/processed/bos.cov.LDRes.barr30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.LDRes.barr 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.LDRes.canPerv.tif", "H:/BosBiog/processed/bos.cov.LDRes.canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.LDRes.canPerv 30m")

### OVeg
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.OVeg.grass.tif", "H:/BosBiog/processed/bos.cov.OVeg.grass30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.OVeg.grass 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.OVeg.barr.tif", "H:/BosBiog/processed/bos.cov.OVeg.barr30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.OVeg.barr 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.OVeg.canPerv.tif", "H:/BosBiog/processed/bos.cov.OVeg.canPerv30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.OVeg.canPerv 30m")

### Water
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.water.perv.tif", "H:/BosBiog/processed/bos.cov.water.perv30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.water.perv 30m")
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.cov.V7.water.open.tif", "H:/BosBiog/processed/bos.cov.water.open30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished cov.water.open 30m")

# AOI2
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.aoi2.tif", "H:/BosBiog/processed/bos.aoi230m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished aoi2 30m")

# ISA RR2
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.isa.RR2.tif", "H:/BosBiog/processed/bos.isa30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished isa 30m")

# # lulc.lumped
# arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.lulc.lumped.tif", "H:/BosBiog/processed/bos.lulc.lumped30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished lulc.lumped 30m")

# grass, total
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.turfgrass.tif", "H:/BosBiog/processed/bos.turfgrass30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished turf grass 30m")

# golf grass
##arcpy.gp.Aggregate_sa("H:/BosBiog/processed/bos.golfgrass.tif", "H:/BosBiog/processed/bos.golfgrass30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished golf grass 30m")


