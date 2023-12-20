
from statistics import mean, median,variance,stdev

data = [20.9684,17.799,17.9775,11.4138,11.5752,28.7541,44.0505,20.0367,29.7451,28.4959,25.5998,18.8128,42.0191,11.7837,12.2605,137.155,17.1531,16.6207,23.8473,34.1845]

m = mean(data)
median = median(data)
variance = variance(data)
stdev = stdev(data)
print('平均: {0:.2f}'.format(m))
print('中央値: {0:.2f}'.format(median))
print('分散: {0:.2f}'.format(variance))
print('標準偏差: {0:.2f}'.format(stdev))
