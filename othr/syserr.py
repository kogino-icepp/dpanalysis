#測定における諸エラーをここにまとめていつでも計算できるようにしておく
import numpy as np
import math
c = 3*10**8
rho = 4.857*10**(-8)
condc = 2.06*10**(7)#A5052の導電率、変更するかも
diel = 8.85*10**(-12)#誘電率(真空)
perme = 1.26*10**(-6)#透磁率
def epieff(f):
    return np.sqrt(2/(2*3.14*f*perme*condc))
def epiratio(f):
    return 2*epieff(f)/(c/f)
print(epieff(240*10**9))
print(epiratio(240*10**9))