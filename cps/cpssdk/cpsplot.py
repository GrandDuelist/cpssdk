import os

class CPSPlot():
    def __init__(self):
        pass

    def smooth(self,ys,widthsize):
        ysnew = [0]*len(ys)
        start = (widthsize-1)/2
        end = len(ys)-start
        ysnew[0:start] = ys[0:start]
        ysnew[end:len(ys)] = ys[end:len(ys)]
        for ii in xrange(start,end):
            s = sum(ys[ii-start:ii+start+1])/float(widthsize)
            ysnew[ii] = s
        return(ysnew)
    
    def styles(self):
        pass
