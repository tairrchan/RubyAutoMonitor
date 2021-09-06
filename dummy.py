import os
#import pandas_datareader as pdr
import scipy 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re

       
def dum(**kwargs):
    a = kwargs.get('keyword')
    print(a != None)
    if a is not None:
        print('yes')
       
       
       
        
if __name__ == "__main__":
    xdata = np.linspace(0, 4, 5)
    ydata = np.linspace(10, 14, 5)
    zdata = np.linspace(10, 14, 20)
    

    target_folderdir = r'\\S4\Datenpool\Yuk Tai\Python\Ruby_automaticTracker'
    
    dir_list = os.listdir(target_folderdir)
    dir_list = [x for x in dir_list if re.search('.*(_HRD).*(\.txt)$', x)]
    print(len(dir_list))

            
    dum(keyword='ok')

    tbl = np.empty([0,4])
    a = ['a','b','c','d']
    b = [23,575,2345,1236]
    c = [1,2,3,4]
    tbl = a
    tbl = np.vstack((tbl,b))
    tbl = np.vstack((tbl,c))
    
    print(tbl)
    print(type(tbl))
    print(tbl)