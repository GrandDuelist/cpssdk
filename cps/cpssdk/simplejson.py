import json
def loadOutEdge(fh):
    data = json.load(fh)
    return(data['out_edge'])

def dumpRegions(regions,ids=None,fh):
    pass
