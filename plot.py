import matplotlib.pyplot as plt
import numpy as np
import sys

data = open(sys.argv[1]).read()
ys = []

for line in data.split('\n'):
	try: 
		ys.append(float(line))
	except:
		continue

ys = np.array(ys)
print ys.mean()
plt.plot(ys)
plt.show()
