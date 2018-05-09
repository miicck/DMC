import matplotlib.pyplot as plt
import sys
file = open(sys.argv[1])
allData = file.read()

split = allData.split('\n')

ys = []

for i in range(0,len(split)):
    entry = split[i].strip()
    try:
        ys.append(float(entry))
    except:
        continue

plt.plot(ys)
#plt.ylim(-0.6,-0.2)
#plt.axes().set_yscale("log")
#plt.xlim(0,10)
#plt.xlabel("Bond length (Å)", fontsize=20)
#plt.ylabel("$E_V$ (Hartree)", fontsize=20)
#plt.tight_layout()
#plt.axes().set_aspect('equal', 'datalim')
plt.show()