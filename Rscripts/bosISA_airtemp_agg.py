# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/BosBiog/BosBiog.gdb"
arcpy.env.overwriteOutput = True

## put 1m ISA into the local UTM
print("working on putting 1m ISA to local UTM")
arcpy.ProjectRaster_management(in_raster="bos.isa1m.nomask.tif", out_raster="E:/BosBiog/BosBiog.gdb/isa1mUTM19N", out_coor_system="PROJCS['NAD_1983_UTM_Zone_19N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-69.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", resampling_type="NEAREST", cell_size="1 1", geographic_transform="", Registration_Point="", in_coor_system="PROJCS['Lambert_Conformal_Conic_2SP',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['false_easting',200000.0],PARAMETER['false_northing',750000.0],PARAMETER['central_meridian',-71.5],PARAMETER['standard_parallel_1',41.71666666666667],PARAMETER['standard_parallel_2',42.68333333333333],PARAMETER['latitude_of_origin',41.0],UNIT['Meter',1.0]]")
## put the 30m AOI grid into local UTM from whatever it is
print("converting 30m AOI grid to local UTM")
arcpy.ProjectRaster_management(in_raster="bos.aoi30m.tif", out_raster="E:/BosBiog/BosBiog.gdb/aoi30mUTM19N", out_coor_system="PROJCS['NAD_1983_UTM_Zone_19N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-69.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", resampling_type="NEAREST", cell_size="30 30", geographic_transform="", Registration_Point="", in_coor_system="PROJCS['UTM_Zone_19_Northern_Hemisphere',GEOGCS['GCS_GRS_1980_IUGG_1980',DATUM['D_unknown',SPHEROID['GRS80',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',-69.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]")
## aggregate the 1m ISA to the local UTM AOI 30m grid
snapR = arcpy.Raster("E:/BosBiog/BosBiog.gdb/aoi30mUTM19N")
arcpy.env.snapRaster = snapR
print("aggregating 1m ISA to local AOI UTM grid")
arcpy.gp.Aggregate_sa("isa1mUTM19N", "E:/BosBiog/BosBiog.gdb/isa30mUTM19N", "30", "MEAN", "EXPAND", "DATA")

