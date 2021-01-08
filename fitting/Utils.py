def get_palette(mode):

 palette = {}
 palette['gv'] = [] 
 
 colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

 for c in colors:
  palette['gv'].append(c)
 
 return palette[mode]
 
def getBinning(binsMVV,minx,maxx,bins):
    l=[]
    if binsMVV=="":
        print " do this"
        print binsMVV
        for i in range(0,bins+1):
            l.append(minx + i* (maxx - minx)/bins)
    else:
        print "dot that"
        print binsMVV
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l

def truncate(binning,mmin,mmax):
    res=[]
    for b in binning:
        if b >= mmin and b <= mmax:
            res.append(b)
    return res
