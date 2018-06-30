

def geombuilder(charge, mult, geometry, fraglist=None):
    if(fraglist==None):
        fraglist=list(range(0,len(geometry)))
    else:
        fraglist = fraglist
    geometry_frag = "\n%d %d\n" %(charge, mult)
    for j in range(len(fraglist)):
        line = geometry[fraglist[j]]
        geometry_frag += line
    return geometry_frag
 
def geombuilder_array(natoms, charge_mult, geometries, fraglist_A=None, fraglist_B=None):
    #NOTE function should return array of tuple containing all necessary geometries
    print("pyREX: Building geometries")
    #for i in range(len(geometries)):
    #    for j in range(natoms):
    #        line = geometry[j] 
