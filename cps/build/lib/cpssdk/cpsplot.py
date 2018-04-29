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
            s = sum(ys[ii-start:ii+start+1])
            s -= ys[ii]
            ysnew[ii] = s/float(widthsize-1) * 0.5 + ys[ii] * 0.5
        return(ysnew)
    
    def styles(self):
        pass
# cpsplot = CPSPlot()
# print(cpsplot.smooth([3,2,5,32,3,5],3))
