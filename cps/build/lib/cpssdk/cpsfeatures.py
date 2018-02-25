# -*- coding: utf-8 -*-
from .cpsspatio import *
class CPSFeatures():
    def __init__(self):
        self.cpsdist = CPSDistance()
        
    def GPSSpeed(self,start_point,end_point,start_time,end_time):
        dist = self.cpsdist.GPSDist(start_point,end_point)*1000
        time_dif = (end_time-start_time).total_seconds()
        return(float(dist)/float(time_dif))
        
