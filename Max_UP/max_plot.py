import numpy as np
from scipy import constants as sp
from matplotlib import pyplot as plt

def dataplot(plot=False):
    delay_file=open('delays.dat', 'r')
    delay_strings=delay_file.readlines()
    delays=[float(d.replace('\n','')) for d in delay_strings]       #list of delays in ps

    te_file=open('electron_temp_map.dat', 'r')
    te_strings_raw=te_file.readlines()
    te_strings=[t.replace('\n','').split() for t in te_strings_raw]
    tes=np.array([[float(layertemp) for layertemp in layer] for layer in te_strings]) #tes[delay][layer] in K
    print(len(tes[0]))
    if plot:
        for i in range(len(tes[0])-1):
            plt.plot(delays, tes[:,i])
        plt.xlabel(r'delays [ps]', fontsize=16)
        plt.ylabel(r'Temperature [K]', fontsize=16)
        plt.show()

dataplot()