

def geombuilder(charge, mult, geometry, fraglist=None):
    if(fraglist==None):
        fraglist=list(range(0,len(geometry)))
    else:
        fraglist = fraglist
    geometry_frag = "\n%d %d\n" %(charge, mult)
    for j in range(len(fraglist)):
        line = geometry[fraglist[j]]
        geometry_frag += line.lstrip()
    return geometry_frag
