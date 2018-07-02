

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

def saptbuilder(natoms, charge_mult, geometries, fraglist_A=None, fraglist_B=None):
    sapt_geometries = []
    print("pyREX: Building SAPT geometries")
    for i in range(len(geometries)):
        sapt_geom = "\n%d %d\n" %(charge_mult[2], charge_mult[3])
        for j in range(len(fraglist_A)):
            line = geometries[i][fraglist_A[j]]
            sapt_geom += line
        sapt_geom += "--\n"
        sapt_geom += "%d %d\n" %(charge_mult[4], charge_mult[5])
        for j in range(len(fraglist_B)):
            line = geometries[i][fraglist_B[j]]
            sapt_geom += line
        print(sapt_geom)
        sapt_geometries.append(sapt_geom)
    return sapt_geometries

def geombuilder_array(natoms, charge_mult, geometries, fraglist_A=None, fraglist_B=None):
    #NOTE function should return array of tuple containing all necessary geometries
    print("pyREX: Building geometries")
    output_geometries = []
    for i in range(len(geometries)):
        dimer_geom = "\n%d %d\n" %(charge_mult[0], charge_mult[1])
        for j in range(natoms):
            line = geometries[i][j]
            dimer_geom += line
        if(fraglist_A):
            frag_A_geom = "\n%d %d\n" %(charge_mult[2], charge_mult[3])
            for j in range(len(fraglist_A)):
                line = geometries[i][fraglist_A[j]]
                frag_A_geom += line
            for j in range(len(fraglist_B)):
                line = geometries[i][fraglist_B[j]]
                frag_A_geom += "@"
                frag_A_geom += line
            frag_B_geom = "\n%d %d\n" %(charge_mult[4], charge_mult[5])
            for j in range(len(fraglist_B)):
                line = geometries[i][fraglist_B[j]]
                frag_B_geom += line
            for j in range(len(fraglist_A)):
                line = geometries[i][fraglist_A[j]]
                frag_B_geom += "@"
                frag_B_geom += line
        else:
            frag_A_geom = None
            frag_B_geom = None
        output_geometries.append((dimer_geom,frag_A_geom,frag_B_geom))
    return output_geometries 
