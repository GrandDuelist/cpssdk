from cpsplot import CPSPlot

cpsplot = CPSPlot()
def testCDF():
    records = [1,2,2,2,3,3,3,4,5,5,5]
    (cdf,xaxis) = cpsplot.binCount(records,9)
    print(len(cdf))
    print(len(xaxis))

if __name__ == "__main__":
    testCDF()
