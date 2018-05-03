from math import sin, cos, sqrt, atan2, radians
import json

class CPSSpatio():
    def __init__(self,grid_shape=None):
        self.grid_shape = (50,50)
        if grid_shape:
            self.grid_shape = grid_shape
        self.regions = {}
        self.grids = []
        # self.init()

    def init(self,out_edge_regions=None):
        if out_edge_regions:
            self.initSpatioRegion(out_edge_regions)
            (minx,miny,maxx,maxy) = self.getRectBoundaryFromRegions(self.regions)
            self.initSpatioRectBoundary(minx,maxx,miny,maxy)
            self.initRectToPolygonMapping()

    
    def initSpatioRegion(self,out_edge):
        for one_region in out_edge:
            geo_id = one_region['geo_id']; geo_array = one_region['geo_array']
            self.regions[geo_id] = geo_array

    def getRectBoundaryFromRegions(self,regions):
        minx,miny,maxx,maxy = float('inf'),float('inf'),0,0
        for k,region in regions.items():
            for point in region:
                [x,y] = point
                if x > maxx: maxx = x
                if x < minx: minx = x
                if y > maxy: maxy = y
                if y < miny: miny = y
        return(minx,miny,maxx,maxy)
        
    def initSpatioRectBoundary(self,minx,maxx,miny,maxy):
        self.minmax = (minx,maxx,miny,maxy)
        self.minx = minx; self.maxx=maxx; self.miny=miny; self.maxy = maxy
        (xshape,yshape) = self.grid_shape
        (xstep,ystep) = ((maxx - minx)/float(xshape),(maxy-miny)/float(yshape))
        self.xstep = xstep; self.ystep=ystep
        self.step = (xstep,ystep)
        self.locationOfGrids()

    def initRectToPolygonMapping(self):
        (minx,maxx,miny,maxy) = self.minmax
        (xshape,yshape) = self.grid_shape
        self.grid_regions = [[[]]*xshape for ii in xrange(yshape)]
        for k,v in self.regions.items():
            for point in v:
                (x,y) = (int((point[0]-minx)/self.xstep),int((point[1]-miny)/self.ystep))
                try:
                    if x == xshape: x-=1
                    if y == yshape: y-=1
                    self.grid_regions[x][y].append(k)
                except:
                    print("x=%d" % x)
                    print('y=%d' % y)

    def pointToGridIndex(self,point):
        (x,y) = (int((point[0]-self.minx)/self.xstep),int((point[1]-self.miny)/self.ystep))
        if x == self.grid_shape[0]:
            x = self.grid_shape[0]-1
        if y == self.grid_shape[1]:
            y = self.grid_shape[1]-1
        return(x,y)
     
    def findCandidatesInGrids(self,point):
        (x,y) = self.pointToGridIndex(point)
        candidates = self.grid_regions[x][y]
        return(candidates)

    def pnpoly(self,polygon,point):
        n=len(polygon);
        i,j,c = 0,n-1,False;
        (testx,testy) = (float(point[0]),float(point[1]))
        while(i < n):
            (currentxi,currentyi) = (float(polygon[i][0]),float(polygon[i][1]))
            (currentxj,currentyj) = (float(polygon[j][0]),float(polygon[j][1]))

            if ((currentyi > testy) != (currentyj > testy)) and (testx < (currentxj - currentxi) * (testy-currentyi) / (currentyj - currentyi) + currentxi):
                c = not c
            j = i
            i += 1
        return(c)

    def findPointRegionID(self,point):
        candidates = self.findCandidatesInGrids(point)
        regionid = self.searchPointInRegions(point,candidates)
        if regionid:
            # print("find in candidates")
            return(regionid)
        else:
            # print("find out candidates")
            return(self.searchPointInRegions(point,self.regions.keys()))

    def searchPointInRegions(self,point,candidates):
        for onekey in candidates:
            polygon = self.regions[onekey]
            if self.pnpoly(polygon,point):
                return(onekey)
        return(None)
    
    def setGeoJsonMultiplePolygon(self,geopolygonarray):
        self.polygons = geopolygonarray
    def findPointInPolygonJson(self,point):
        N = len(self.polygons)
        for ii in xrange(N):
            polygon = self.polygons[ii]
            if type(polygon["geometry"]["coordinates"][0][0][0]) is not list:
                temp_array = polygon["geometry"]["coordinates"][0]
                if self.pnpoly(temp_array,point):
                    return(polygon)
            else:
                for one_array in polygon["geometry"]["coordinates"]:
                    temp_array = one_array[0]
                    if self.pnpoly(temp_array,point):
                        return(polygon)
        return(None)

    def locationOfGrids(self):
        (xshape,yshape)= self.grid_shape
        (xstep,ystep) = self.step
        grid_location = [[[0,0]]*yshape for ii in xrange(xshape)]
        x = xstep/2.0+self.minx;  y = ystep/2.0 + self.miny
        for ii in xrange(xshape):
            y = ystep/2.0 + self.miny
            for jj in xrange(yshape):
                grid_location[ii][jj][0] = x
                grid_location[ii][jj][1] = y
                y = y+ystep 
            x = x + xstep; 
        self.grid_location = grid_location
        return(grid_location)

    def countInRegion(self,X,Y,Z=None):
        '''
        X: longitude or x 
        Y: latitude or y
        Z: None or number on locatio X and Y
        return: the density in each grid
        ''' 
        if Z is None: Z = [0] * len(X)
        grid_count = [[0] * self.grid_shape[1] for ii in xrange(self.grid_shape[0])]
        for ii in xrange(len(X)):
            x = X[ii]; y = Y[ii]; p = [x,y]
            grid = self.pointToGridIndex(p)
            grid_x,grid_y = grid
            if grid_x < 0 or grid_x > self.grid_shape[0]: continue
            if grid_y < 0 or grid_y > self.grid_shape[1]: continue
            grid_count[grid_x][grid_y] += Z[ii]
        return(grid_count)
    
def minmaxGeoJson(geojson_path):
    file_path = geojson_path
    # file_path = 'shenzhen_boundary_gps.geoJson'
    data = json.load(open(file_path))
    coordinates = data['coordinates'][0]
    minx=miny=100000; maxx=maxy = 0;
    for point in coordinates:
        (x,y) = point
        minx = min(minx,x);maxx = max(maxx,x);miny=min(miny,y);maxy=max(maxy,y)
    print("minx = "+str(minx)+" maxx = " + str(maxx) + " miny = " + str(miny) + " maxy = "+str(maxy))        

class CPSCrop():
    def __init__(self):
        self.cpsspatio = CPSSpatio()
    def setRectangle(self,minx,maxx,miny,maxy):
        self.minx = minx; self.miny = miny; self.maxx = maxx; self.maxy = maxy
    def isInRectangle(self,x,y):
        isx = (x >= self.minx and x <= self.maxx)
        isy = (y >= self.miny and y <= self.maxy)
        return(isx and isy)
    def setShenzhenRectangle(self):
        self.setRectangle(113.7463515,114.6237079,22.4415225,22.8644043)
    def setPolygonBoundary(self,polygon_list):
        self.polygon_boundary = polygon_list
        
    def setShenzhenPolygonBoundary(self):
        data = json.load(open("data/boundary/shenzhen_boundary_gps.geoJson"))
        polygon_list = data['coordinates'][0]
        self.setPolygonBoundary(polygon_list)
    
    def isInPolygon(self,x,y):
        return(self.cpsspatio.pnpoly(self.polygon_boundary,[x,y]))

class CPSDistance():
    def GPSDist(self,p1,p2):
        #km
        R=  6373.0
        lon1 = radians(p1[0]); lat1 = radians(p1[1])
        lon2 = radians(p2[0]); lat2 = radians(p2[1])
        dlon=lon2-lon1; dlat = lat2-lat1
        a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        distance = R * c
        return(distance)
