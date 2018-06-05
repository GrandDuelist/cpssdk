from cpsspatio import *
import pandas as pd
#p1 = (21.0122287,52.2296756)
#p2 = (16.9251681,52.406374)
#cpsdistance = CPSDistance()
#r = cpsdistance.GPSDist(p1,p2)
def testCPSCrop():
    cpscrop = CPSCrop()
    cpscrop.setShenzhenPolygonBoundary()
    b = cpscrop.isInPolygon(114.0508539, 22.5463145)

def testGeoJsonMultiplePolygon():
    cpsspatio = CPSSpatio()
    aa = simplejson.load(open(r'C:\Users\zhiha\Downloads\china.json'))
    cpsspatio.setGeoJsonMultiplePolygon(aa)
    cpsspatio.findPointInPolygonJson((118.4187,28.2968))

def testGridCount():
    cpsspatio = CPSSpatio()
    cpsspatio.initSpatioRectBoundary(113.7463515,114.6237079,22.4415225,22.8644043)
    result = cpsspatio.countInGrid([114.23,113.94,114.312,114.11],[22.55,22.65,22.71,22.48],[1,3,4,5])
    print(cpsspatio.grid_location)
    #print(result)
def testSimpleJsonToGeoPandas():
    out_edge=json.load(open(r"E:\drive\W-WorkingOn\1-coding\1-research-projects\7-NYC\spark.sql\json\new-york-block-region.simpleJson"))['out_edge']
    cpsspatio = CPSSpatio()
    test = cpsspatio.simpleJsonToGeoPandas(out_edge)
    print(test)

def testGenerateGuangdong():
    centersPath = r"E:\drive\W-WorkingOn\1-coding\2-visual\1-processing\6-plot-ETC-density\process-data\data\station_volume.csv"
    boundaryPath  = r"E:\drive\W-WorkingOn\1-coding\2-visual\1-processing\6-plot-ETC-density\process-data\data\guangdong-boundary.json"
    cpsspatio = CPSSpatio()
    data = pd.read_csv(centersPath,header=None)
    centers = list(data.apply(lambda x: [x[0],x[1]], axis=1))
    boundary = json.load(open(boundaryPath))['coordinates'][0]
    simpleJson  = cpsspatio.generateVoronoiInBoundary(centers=centers,boundary=boundary)
    print(simpleJson)
testGenerateGuangdong()


