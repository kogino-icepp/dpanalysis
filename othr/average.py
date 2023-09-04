import numpy as np
f = open('sigma.txt','r')
datalist = f.readlines()
numeric_data = [float(value) for value in datalist]
average_value = np.average(numeric_data)
print(average_value)