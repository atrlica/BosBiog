# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "H:/BosBiog/BosBiog.gdb"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("H:/BosBiog/processed/bos.isa.RR2.tif")
    print("found 1m raster for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR
csRes = arcpy.GetRasterProperties_management(in_raster="H:/BosBiog/processed/bos.isa.RR2.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

arcpy.Select_analysis(in_features="H:/BosBiog/data/LULC/LU_polys_NAD83UTM19N/LU_polys_NAD83UTM19N.shp", out_feature_class="H:/BosBiog/BosBiogWorking.gdb/LU_BOSTON", where_clause=""""TOWN" = 'BOSTON'""")
print("selected Boston LULC polys") 
arcpy.PolygonToRaster_conversion(in_features="H:/BosBiog/BosBiogWorking.gdb/LU_BOSTON", value_field="LUCODE", out_rasterdataset="H:/BosBiog/processed/bos.lulc.tif", cell_assignment="CELL_CENTER", priority_field="NONE", cellsize=str(cs))
print("rasterized LULC shp to grid "+str(cs)+"m")

