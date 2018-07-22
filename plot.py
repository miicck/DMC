import matplotlib.pyplot as plt
import numpy as np
import sys

data = open(sys.argv[1]).read()
equil = int(sys.argv[2])

es = []
ps = []
rs = []

count = 0
equil_es = []

for line in data.split('\n'):
	count = count + 1;
	spl = line.split()
	try:
		es.append(float(spl[0]))
		ps.append(float(spl[1]))		
		rs.append(float(spl[2]))
		if (count > equil):
			equil_es.append(float(spl[0]))
	except:
		continue

equil_es = np.array(equil_es)
print "Energy after equilibriation: ", equil_es.mean()

plt.subplot(221)
plt.plot(es)

plt.subplot(222)
plt.plot(ps)

plt.subplot(223)
plt.plot(rs)

plt.show()
