from cpsspatio import *
p1 = (21.0122287,52.2296756)
p2 = (16.9251681,52.406374)
cpsdistance = CPSDistance()
r = cpsdistance.GPSDist(p1,p2)
cpscrop = CPSCrop()
cpscrop.setShenzhenPolygonBoundary()
b = cpscrop.isInPolygon(114.0508539, 22.5463145)
print(b)


