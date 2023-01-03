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
    if plot:
        for i in range(len(tes[0])-1):
            plt.plot(delays, tes[:,i])
        plt.xlabel(r'delays [ps]', fontsize=16)
        plt.ylabel(r'Temperature [K]', fontsize=16)
        plt.show()
    return delays, tes

def simplot(file, offset):
    ### open(file, mode) is an in-built function. 1st argument is the filepath, second argument is the mode. Choose mode 'r' to read the desired file
    input = open(file, 'r')
    content = input.read().split('\n')
    toplot = [line for line in content if not str(line).startswith('#')]
    lines = [line.split('\t') for line in toplot]
    del lines[-1]

    times = np.array([-offset + float(line[0]) for line in lines])

    mags = [line[1] for line in lines]
    mag1bae = [mag.replace('[0.0]]', '').replace('[[', '') for mag in mags]
    mag1 = [m.replace('[', '').replace(']', '') for m in mag1bae]
    m1 = [line.split() for line in mag1]
    m1s = np.array([[float(m1[i][j].replace(',', '')) for i in range(len(m1))] for j in range(len(m1[0]))])
    ####m1s[i][:] is magnetization for all times of (i-1)th monolayer!###
    ###m1s[:][i] is magnetization for timestep i of all monolayers!###

    mags2 = [line[5] for line in lines]
    mag2bae = [mag.replace('[0.0]', '').replace('[[', '') for mag in mags2]
    mag2 = [m.replace('[', '').replace(']', '') for m in mag2bae]
    m2 = [line.split() for line in mag2]
    m2s = np.array([[float(m2[i][j].replace(',', '')) for i in range(len(m2))] for j in range(len(m2[0]))])
    print(len(m1s[0]))
    print(len(m2s[0]))

    temes = [line[2].split('], [') for line in lines]
    tebae = [[te[i].replace('[', '').replace(']', '') for i in range(len(temes[0]))] for te in temes]
    te1bla = [tb[:][0] for tb in tebae]

    te1n = [b.replace(',', '') for b in te1bla]
    te1m = [te.split() for te in te1n]
    te1 = [[float(i) for i in te] for te in te1m]

    tpmes = [line[3].split('], [') for line in lines]
    tpbae = [[tp[i].replace('[', '').replace(']', '') for i in range(len(tpmes[0]))] for tp in tpmes]
    tp1bla = [tb[:][0] for tb in tpbae]
    tp1n = [b.replace(',', '') for b in tp1bla]
    tp1m = [tp.split() for tp in tp1n]
    tp1 = [[float(i) for i in tp] for tp in tp1m]

    tp2s = [line[4].split('], [') for line in lines]
    tp2bae = [[tp[i].replace('[', '').replace(']', '') for i in range(len(tpmes[0]))] for tp in tp2s]
    tp21bla = [tb[:][0] for tb in tp2bae]
    tp21n = [b.replace(',', '') for b in tp21bla]
    tp21m = [tp.split() for tp in tp21n]
    tp21 = [[float(i) for i in tp] for tp in tp21m]

    return (np.array(times), np.array(m1s), np.array(te1), np.array(tp1), np.array(tp21), np.array(m2s))

dat_times, dat_tes=dataplot()
sim_times, sim_ms, sim_te, sim_tp, sim_tp2, sim_m2s=simplot('test2.dat', )