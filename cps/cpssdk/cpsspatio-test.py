from cpsspatio import *
import simplejson
#p1 = (21.0122287,52.2296756)
#p2 = (16.9251681,52.406374)
#cpsdistance = CPSDistance()
#r = cpsdistance.GPSDist(p1,p2)
#cpscrop = CPSCrop()
#cpscrop.setShenzhenPolygonBoundary()
#b = cpscrop.isInPolygon(114.0508539, 22.5463145)
cpsspatio = CPSSpatio()
aa = simplejson.load(open('C:\Users\zhiha\Downloads\china.json'))
cpsspatio.setGeoJsonMultiplePolygon(aa)
cpsspatio.findPointInPolygonJson((118.4187,28.2968))
#%%


