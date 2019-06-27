# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "H:/BosBiog/BosBiogworking.gdb"
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

### This approach is different from FragEVI -- we pre-lump the polygon LULC categories and THEN rasterize, rather than rasterize and lump
# arcpy.Select_analysis(in_features="E:/FragEVI/data/LULC/LU_polys_NAD83UTM19N/LU_polys_NAD83UTM19N.shp", out_feature_class="E:/FragEVI/processed/boston/LU_polys_bos.shp", where_clause=""""TOWN" = 'BOSTON'""")
# print("selected Boston LULC polys")
# arcpy.PolygonToRaster_conversion(in_features="H:/BosBiog/data/LULC/LU_polys_NAD83UTM19N/LU_polys_BOSTON_lump.shp", value_field="LU_LUMP", out_rasterdataset="H:/BosBiog/processed/bos.lulc.lumped.tif", cell_assignment="CELL_CENTER", priority_field="NONE", cellsize=str(cs))
# print("rasterized LULC_LUMP shp to ISARR2 grid "+str(cs)+"m")

### in case you need the raw 1m LULC category raster
arcpy.PolygonToRaster_conversion(in_features="H:/BosBiog/data/LULC/LU_polys_NAD83UTM19N/LU_polys_BOSTON_lump.shp", value_field="LUCODE", out_rasterdataset="H:/BosBiog/processed/bos.lulc.tif", cell_assignment="CELL_CENTER", priority_field="NONE", cellsize=str(cs))
print("rasterized LULC shp to ISARR2 grid "+str(cs)+"m")


### now do it again with the 30m grid guide
try:
    snapR = arcpy.Raster("H:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")
    print("found 30m EVI for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR
csRes = arcpy.GetRasterProperties_management(in_raster="H:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

# arcpy.PolygonToRaster_conversion(in_features="H:/BosBiog/data/LULC/LU_polys_NAD83UTM19N/LU_polys_BOSTON_lump.shp", value_field="LU_LUMP", out_rasterdataset="H:/BosBiog/processed/bos.lulc.lumped30m.tif", cell_assignment="MAXIMUM_COMBINED_AREA", priority_field="NONE", cellsize=str(cs))
# print("rasterized LULC_LUMP shp to EVI grid "+str(cs)+"m")

