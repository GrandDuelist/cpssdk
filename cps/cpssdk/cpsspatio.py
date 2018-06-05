from math import sin, cos, sqrt, atan2, radians
import json
import collections
from cpsspatio_interval import *


class CPSSpatio():
    def __init__(self,grid_shape=None):
        self.grid_shape = (50,50)
        if grid_shape:
            self.grid_shape = grid_shape
        self.regions = {}
        self.grids = []
        self.spatioInterval = CPSSpatioInterval()
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
        minx,miny,maxx,maxy = float('inf'),float('inf'),-float('inf'),-float('inf')
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
        self.grid_regions = [[[]]*yshape for ii in xrange(xshape)]
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
     
    def findCandidatesInGrids(self,point,window,previous):
        (x,y) = self.pointToGridIndex(point)
        candidates = []
        candidates.extend(self.grid_regions[x][y])
        lx = len(self.grid_regions)
        ly = len(self.grid_regions[0])
        for ii in xrange(x-window,x+window+1):
            for jj in xrange(y-window,y+window+1):
                if ii < 0 or jj < 0 or ii > lx-1 or jj > ly-1 or (ii==x and jj==y):
                    continue
                candidates.extend(self.grid_regions[ii][jj]) 
        candidates = list(set(candidates))
        mm = {}; new_cand = []
        for one in previous: 
            mm[one] = True
        for one in candidates:
            if not mm.get(one,False):
                new_cand.append(one)
        return(new_cand,candidates)

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

    def findPointRegionID(self,point,window=0,previous=[],limit=20):
        (candidates,previous) = self.findCandidatesInGrids(point,window,previous)
        regionid = self.searchPointInRegions(point,candidates)
        if regionid:
            # print("find in candidates")
            return(regionid)
        else:
            if window > limit:
                # return(self.searchPointInRegions(point,self.regions.keys()))
                return(None)
            return(self.findPointRegionID(point,window+1,previous,limit))
            # print("find out candidates")
            # return(self.searchPointInRegions(point,self.regions.keys()))

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
        grid_location = [[[0,0] for jj in xrange(yshape)] for ii in xrange(xshape)]
        x = xstep/2.0+self.minx;  y = ystep/2.0 + self.miny
        xs = [0] * xshape
        ys = [0] * yshape 
        for ii in xrange(xshape):
            xs[ii] = x
            x += xstep
        for ii in xrange(yshape):
            ys[ii] = y
            y += ystep
        x = xstep/2.0+self.minx;  y = ystep/2.0 + self.miny
        for ii in xrange(xshape):
            y = ystep/2.0 + self.miny
            for jj in xrange(yshape):
                grid_location[ii][jj][0] = x
                grid_location[ii][jj][1] = y
                y += ystep 
            x += xstep; 
        self.grid_location = grid_location
        self.grid_xs = xs
        self.grid_ys = ys
        return(grid_location,xs,ys)

    def countInGrid(self,X,Y,Z=None):
        """ 
        X: longitude or x 
        Y: latitude or y
        Z: None or number on locatio X and Y
        return: the density in each grid
        """ 
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

    def cutPointByBoundary(self,locations,values,boundaryList,defaultvalue=-999):
        for ii in xrange(len(locations)): 
            for jj in xrange(len(locations[0])):
                p = locations[ii][jj]
                value = values[ii][jj]
                if not self.pnpoly(boundaryList,p):
                    values[ii][jj] = defaultvalue
        return(values)

    def simpleJsonToGeoPandas(self,out_edge):
        import geopandas
        from shapely.geometry import Polygon, Point
        l = len(out_edge[0].keys()) 
        result = collections.defaultdict(list) 
        keys  = out_edge[0].keys()
        for one in out_edge:
            for ii in xrange(l):
                one_key = keys[ii]
                if one_key == 'geo_array':
                    p = Polygon(one[one_key])
                    result[one_key].append(p)
                    continue
                result[one_key].append(one.get(one_key,None))
        data = geopandas.GeoDataFrame(data=result)
        return(data)
    
    def minmaxGeoJson(self,geojson_path):
        file_path = geojson_path
        # file_path = 'shenzhen_boundary_gps.geoJson'
        data = json.load(open(file_path))
        coordinates = data['coordinates'][0]
        minx=miny=100000; maxx=maxy = -10000;
        for point in coordinates:
            (x,y) = point
            minx = min(minx,x);maxx = max(maxx,x);miny=min(miny,y);maxy=max(maxy,y)
        print("minx = "+str(minx)+" maxx = " + str(maxx) + " miny = " + str(miny) + " maxy = "+str(maxy))     

    def gpsPolygonArea(self,polygon):
        '''
        polygon: array contains vertex of the polygon
        return: area in km^2
        '''
        import pyproj    
        import shapely
        import shapely.ops as ops
        from shapely.geometry.polygon import Polygon
        from functools import partial
        
        geom = Polygon(polygon)
        geom_area = ops.transform(
            partial(
                pyproj.transform,
                pyproj.Proj(init='EPSG:4326'),
                pyproj.Proj(
                    proj='aea',
                    lat1=geom.bounds[1],
                    lat2=geom.bounds[3])),
            geom)
        return(geom_area.area/1000000)

    def simpleJsonToIdJson(self,out_edge,key='geo_id',value='geo_array'):
        '''
        out_edge: array of geoid and geo array for polygon 
        return: idJson mapping id to geoarray
        '''
        mm ={}
        for oneRegion in out_edge:
            geoID = oneRegion['geo_id']
            geoArray = oneRegion['geo_array']
            mm[geoID] = geoArray
        return(mm)


    def geoArrayFromPolygonString(self,polygonString):
        pointString = polygonString.replace("POLYGON ((","").replace("))","")
        points = pointString.split(",")
        geo_array = [];
        for one in points:
            attrs = one.strip(" ").split(" ")
            lon = float(attrs[0]); lat = float(attrs[1])
            geo_array.append([lon,lat])
        return(geo_array)

    def generateVoronoiInBoundary(self,centers,boundary):
        '''
        Generate voronoi polygons based on centers and boundary of the city
        :param centers: list of locations of voronoi centers
        :param boundary: boundary of the voronoi partition
        :return: a simple json file with voronoi polygons
        '''
        from scipy.spatial import Voronoi
        from shapely.geometry import Polygon
        import numpy as np
        points = np.array(centers); vor = Voronoi(centers);  mask = Polygon(boundary); id=0
        regions,vertices = self.spatioInterval.voronoi_finite_polygons_2d(vor); allRegions = []
        for ii in range(0,len(regions)):
            region = regions[ii]
            polygon = vertices[region]
            shape = list(polygon.shape)
            shape[0] += 1
            p = Polygon(np.append(polygon, polygon[0]).reshape(*shape)).intersection(mask)
            try:
                poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
                allRegions.append({"geo_array": poly.tolist()+[poly.tolist()[0]], "geo_center":centers[ii], "geo_id":id})
                id += 1
            except:
                pp = p
                for p in pp:
                    poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
                    allRegions.append({"geo_array": poly.tolist()+[poly.tolist()[0]],"geo_center":centers[ii], "geo_id":id})
                    id += 1
        return({"out_edge": allRegions})

    # def generateVoronoiInBoundary(self,simpleJson,centerKeys,boundary):
    #     '''
    #     Generate voronoi polygons based on centers and boundary of the city
    #     :param centers: list of locations of voronoi centers
    #     :param boundary: boundary of the voronoi partition
    #     :return: a simple json file with voronoi polygons
    #     '''
    #     from scipy.spatial import Voronoi
    #     from shapely.geometry import Polygon
    #     import numpy as np
    #     points = np.array(centers); vor = Voronoi(centers);  mask = Polygon(boundary)
    #     regions,vertices = self.spatioInterval.voronoi_finite_polygons_2d(); allRegions = []
    #     for ii in range(0,len(regions)):
    #         region = regions[ii]
    #         polygon = vertices[region]
    #         shape = list(polygon.shape)
    #         shape[0] += 1
    #         p = Polygon(np.append(polygon, polygon[0]).reshape(*shape)).intersection(mask)
    #         try:
    #             poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
    #             allRegions.append({"geo_array": [poly.tolist()+[poly.tolist()[0]]], "geo_center":centers[ii]})
    #         except:
    #             pp = p
    #             for p in pp:
    #                 poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
    #                 allRegions.append({"geo_array": [poly.tolist()+[poly.tolist()[0]]],"geo_center":centers[ii]})
    #     return all_regions

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
