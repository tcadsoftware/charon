import sys

sweepGoldFile = sys.argv[1]
sweepResultFile = sys.argv[2]

with open(sweepGoldFile,'r') as f:
    sweepGold = f.readlines()
voltages = [p.split()[0] for p in sweepGold[1:]]
dataGold = [float(p.split()[5]) for p in sweepGold[1:]]

with open(sweepResultFile,'r') as f:
    sweepResult = f.readlines()
dataResult = [float(p.split()[5]) for p in sweepResult[1:]]


diffNormSqr, goldNormSqr = 0.0, 0.0
for x in zip(dataGold,dataResult):
    diffNormSqr += (x[0]-x[1])**2
for x in dataGold:
    goldNormSqr += x**2
relErr = (diffNormSqr/goldNormSqr)**0.5

print("relError = "+str(relErr))
if relErr < 0.05:
    print("PASS")
else:
    print("FAIL")
