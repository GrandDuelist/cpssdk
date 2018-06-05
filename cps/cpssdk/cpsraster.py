import os 
import json
import shutil
import logging

class CPSPop():
    """
    functions: extract detailed population from worldpop dataset
    """
    def __init__(self,tempDirectory="temp"):
        """
        tempDirectory: location of directory to store generate tiff files for polygons
        """
        if tempDirectory is not None:
            self.tempDirectory = tempDirectory
            self.tempJsonDir  = os.path.join(self.tempDirectory,"json")
            self.tempTiffDir  = os.path.join(self.tempDirectory,'tiff')
            if not os.path.exists(self.tempDirectory): os.makedirs(self.tempDirectory)
            if not os.path.exists(self.tempJsonDir): os.makedirs(self.tempJsonDir)
            if not os.path.exists(self.tempTiffDir): os.makedirs(self.tempTiffDir)
        self.currentSimpleJson = None

    def popInPolygon(self,geo_array,tiff): 
        geo_json = self.geoArrayToGeoJson(geo_array) 

    def popInMultiPolygonFromSimpleJson(self,rasterFilePath,geoJsonFilePath):
        """
        func:
            map population to polygon
        param:
            rasterFilePath: file path information is extracted from  
            geoJsonFilePath:  file path to simple json with multiple polygon
            return: simpleJson with pop attributes
        return:
            population in polygon
        """
        self.tiffInMultiPolygonFromSimpleJson(rasterFilePath,geoJsonFilePath)
        self.tiffToPopSimpleJson(geoJsonFilePath)
        self.cleanTemp()

    def tiffInMultiPolygonFromSimpleJson(self,rasterFilePath,geoJsonFilePath):
        """ 
        param:
            rasterFilePath: file path information is extracted from  
            geoJsonFilePath:  file path to simple json with multiple polygon
        """
        self.currentSimpleJson = geoJsonFilePath
        simplejson = json.load(open(geoJsonFilePath))
        if type(simplejson) == dict:
            outEdge = simplejson['out_edge']
        else:
            outEdge = simplejson
        self.tiffInMultiPolygonFromGeoArray(rasterFilePath,outEdge)

    def tiffInPolygonFromGeoArray(self,rasterFilePath,oneRegion):
        """
        param:
            rasterFilePath:  location of raster file
            outEdge: out edge is a collection of geo arrays and it is a key of simple json file
        """
        geoid = oneRegion['geo_id'] 
        geoarray = oneRegion['geo_array']
        tempjson = self.geoArrayToGeoJson(geoarray)
        tempfile = os.path.join(self.tempJsonDir,"{geoid}.json".format(geoid=geoid))
        json.dump(tempjson,open(tempfile,'wb'))
        outfile = os.path.join(self.tempTiffDir,"{geoid}.tiff".format(geoid=geoid))
        command_line = "gdalwarp -cutline {tempfile} -crop_to_cutline -dstalpha {rasterfile} {outputfile}"
        commandline = command_line.format(tempfile=tempfile,rasterfile=rasterFilePath,outputfile=outfile)
        os.system(commandline)

    def tiffInMultiPolygonFromGeoArray(self,rasterFilePath,outEdge):
        """
        param:
            rasterFilePath:  location of raster file
            outEdge: out edge is a collection of geo arrays and it is a key of simple json file
        """
        for oneRegion in outEdge:
            self.tiffInPolygonFromGeoArray(rasterFilePath,oneRegion)

    def tiffToPopSimpleJson(self,simpleJsonPath=None):
        '''
        func:
            calculate the population of 
        param:
        '''
        from osgeo import gdal
        if simpleJsonPath is None:
            simpleJsonPath = self.currentSimpleJson
        simplejson = json.load(open(simpleJsonPath))
        if type(simplejson) is not dict:
            temp = {'out_edge':simplejson}
            simplejson = temp
        outEdge = simplejson['out_edge']; allPop = 0
        for oneRegion in outEdge:
            geoid = oneRegion['geo_id']
            tifffile = os.path.join(self.tempTiffDir,"{geoid}.tiff".format(geoid=geoid))
            if not os.path.exists(tifffile):
                oneRegion['pop']=0
                continue
            logging.info("count pop in geoid = %s" % str(geoid))
            srcDs = gdal.Open(tifffile)
            band = srcDs.GetRasterBand(1)
            values = band.ReadAsArray()
            sumValue = 0
            for valueArray in values:
                sumValue += sum([value for value in valueArray if value > 0])
            oneRegion['pop'] = sumValue
            allPop += sumValue
        simplejson['pop'] = allPop    
        return(simplejson)

    def geoArrayToGeoJson(self,geoArray):
        """
        geoArray: array of point to present polygon
        """
        oneRegionJson  = {"type":"Polygon", "coordinates": [geoArray]}
        return(oneRegionJson)

    def cleanTemp(self):
        if os.path.exists(self.tempDirectory):
            shutil.rmtree(self.tempDirectory)
# if __name__ == "__main__":
    # cpspop = CPSPop()
    # cpspop.popInPolygon()
