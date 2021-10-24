

start Modifier 0

def testForModification(self,pLList):

    for lineNumber, line in enumerate(pLList):
        if line.find("PHASESHIFTS") >= 0:
            # the line to replace will look like this:
            # 'PHASESHIFTS1e1   3e1   2e2PHASESHIFTS'
            # and it should become
            # '1e1,3e1,2e2'
            [begin, replacement, end] = pLList[lineNumber].split('PHASESHIFTS')
            replacement = ','.join(replacement.split())
            pLList[lineNumber] = begin+replacement+end

    return(pLList)

end Modifier 0
