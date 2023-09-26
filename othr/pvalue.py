import numpy as np
from math import gamma
from scipy.special import gammainc
import matplotlib.pyplot as plt
#確率変数と自由度を渡す
def pchi(x,k):
    return 1-(gammainc(k/2,x/2)/gamma(k/2))
def pchindf(x,k):
    return 1-(gammainc(k/2,k*x/(2)))
def psigma(x,sigma):
    return 0
def main():
    #随時p-valueを知りたい関数を追加していく
    s = input()
    if(s=="pchi"):
        print("x : ",end='')
        x = float(input())
        print("k : ",end='')
        k = int(input())
        print("p_value : %f" % pchi(x,k))
    elif(s=="pchindf"):
        print("x : ",end='')
        x = float(input())
        print("k : ",end='')
        k = int(input())
        print("p_value : %f" % pchindf(x,k))
    elif(s=="psigma"):
        x = float(input())
        sigma = float(input())
        print(psigma(x,sigma))
    else:
        print("Can't Find such case!")
if __name__ == "__main__":
    main()