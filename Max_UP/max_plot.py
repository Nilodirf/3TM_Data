import numpy as np
from scipy import constants as sp
from matplotlib import pyplot as plt
import os
from npy_append_array import NpyAppendArray

def dataplot(plot=False):
    delay_file=open('delays.dat', 'r')
    delay_strings=delay_file.readlines()
    delays=[float(d.replace('\n','')) for d in delay_strings]       #list of delays in ps

    te_file=open('electron_temp_map.dat', 'r')
    te_strings_raw=te_file.readlines()
    te_strings=[t.replace('\n','').split() for t in te_strings_raw]
    tes=np.array([[float(layertemp) for layertemp in layer] for layer in te_strings]) #tes[delay][layer] in K
    if plot:
        for i in range(len(tes[0])-1):
            plt.plot(delays, tes[:,i])
        plt.xlabel(r'delays [ps]', fontsize=16)
        plt.ylabel(r'Temperature [K]', fontsize=16)
        plt.show()
    return delays, tes

def simplot(sim, offset):
    parentdir='sims'
    simdir=sim

    directory=os.path.join(parentdir, sim)
    files=os.listdir(directory)

    #extract delays
    delayfile=os.path.join(directory,'delays.npy')
    delays=np.load(delayfile, mmap_mode='r')

    #find all the electron temperatures
    tefiles=[f for f in files if str(f).startswith('tes')]
    for tef in tefiles:
        file=os.path.join(directory, tef)
        tes=np.load(file, mmap_mode='r')
        print(tes.shape)
        plt.plot(delays, tes)

    plt.show()


dat_times, dat_tes=dataplot()
simplot('test1', 0.1)


#print(sim_te)
#plt.plot(sim_times, sim_te[0])
#plt.show()