

start Modifier 0

def testForModification(self,pLList):

    pLList[:] = [l.replace("LATTICE_TEMPERATURE,Value,double","Lattice Temperature,Value,double") for l in pLList]

    return(pLList)

end Modifier 0


start Modifier 1

def testForModification(self,pLList):

    pLList[:] = [l.replace("LATTICE_TEMPERATURE,Value,string","Lattice Temperature,Value,string") for l in pLList]

    return(pLList)

end Modifier 1
