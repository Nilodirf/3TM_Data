import numpy as np
from scipy import constants as sp
from matplotlib import pyplot as plt
import os

def dataplot():
    delay_file=open('max_dat/delays.dat', 'r')
    delay_strings=delay_file.readlines()
    delays=[float(d.replace('\n','')) for d in delay_strings]       #list of delays in ps

    te_file=open('max_dat/electron_temp_map.dat', 'r')
    te_strings_raw=te_file.readlines()
    te_strings=[t.replace('\n','').split() for t in te_strings_raw]
    tes=np.array([[float(layertemp) for layertemp in layer] for layer in te_strings]) #tes[delay][layer] in K
    tp_file=open('max_dat/phonon_temp_map.dat', 'r')
    tp_strings_raw=tp_file.readlines()
    tp_strings=[t.replace('\n','').split() for t in tp_strings_raw]
    tps=np.array([[float(layertemp) for layertemp in layer] for layer in tp_strings]) #tps[delay][layer] in K

    return delays, tes, tps

def simplot(sim, offset):
    parentdir='sims'
    simdir=sim

    directory=os.path.join(parentdir, sim)
    files=os.listdir(directory)


    #get delays
    delayfile=os.path.join(directory,'delays.npy')
    delays=np.load(delayfile, mmap_mode='r')*1e12-offset


    #get all the electron temperatures
    tefiles=[f for f in files if str(f).startswith('tes')]
    for tef in tefiles:
        mat=str(tef).replace('.npy', '').replace('tes','')
        file=os.path.join(directory, tef)
        tes=np.load(file, mmap_mode='r')
        if mat=='Nickel':
            tesNickel=tes
        if mat=='Gold':
            tesGold=tes
        if mat=='Tantalum':
            tesTa=tes
        if mat=='substrate':
            tessub=tes
        # plt.plot(delays, tes)
        # plt.legend([mat], fontsize=14)
        # plt.xlabel(r'delay [ps]', fontsize=16)
        # plt.ylabel(r'$T_e$ [K]', fontsize=16)
        # plt.show()

    #get all lattice temperatures
    tpfiles=[f for f in files if str(f).startswith('tps')]
    for tpf in tpfiles:
        mat=str(tpf).replace('.npy', '').replace('tps', '')
        file=os.path.join(directory, tpf)
        tps=np.load(file, mmap_mode='r')
        if mat=='Nickel':
            tpsNickel=tps
        if mat=='Gold':
            tpsGold=tps
        if mat=='Tantalum':
            tpsTa=tps
        if mat=='substrate':
            tpssub=tps
        # plt.plot(delays, tps)
        # plt.legend([str(tpf).replace('.npy', '').replace('tps', '')], fontsize=14)
        # plt.xlabel(r'delay [ps]', fontsize=16)
        # plt.ylabel(r'$T_p$ [K]', fontsize=16)
        # plt.show()

    #get magnetization
    Nimag=np.load(os.path.join(directory, 'msNickel.npy'))
    Nim0=Nimag[0,0]
    # plt.plot(delays, Nimag/Nim0)
    # plt.legend(['Nickel'], fontsize=14)
    # plt.xlabel(r'delay [ps]', fontsize=16)
    # plt.ylabel(r'$m/m_0$', fontsize=16)
    # plt.show()
    return delays, tesGold, tesNickel, tesTa, tessub, tpsGold, tpsNickel, tpsTa, tpssub, Nimag

dat_times, dat_tes, dat_tps=dataplot()
sim_times, sim_te_au, sim_te_ni, sim_te_ta, sim_te_sub, sim_tp_au, sim_tp_ni, sim_tp_ta, sim_tp_sub, sim_mag_ni=simplot('testkph3', 0.1)

def plot(times, dat, ylabel, color, show):
    plt.plot(times, dat, color=color)
    plt.xlabel(r'delay [ps]', fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    if show:
        plt.show()
    return

plot(dat_times, dat_tes, 'Temperature', 'blue', True)
plot(dat_times, dat_tps, 'Temperature', 'blue', True)
plot(sim_times, sim_te_ni, 'Temperature', 'orange', False)
plot(sim_times, sim_te_au, 'Temperature', 'red', False)
plot(sim_times, sim_te_ta, 'Temperature', 'yellow', False)
plot(sim_times, sim_tp_sub, 'Temperature', 'green', False)
plt.show()

plot(sim_times, sim_tp_ni, 'Tempertature', 'orange', False)
plot(sim_times, sim_tp_au, 'Tempertature', 'red', False)
plot(sim_times, sim_tp_ta, 'Tempertature', 'yellow', False)
plot(sim_times, sim_tp_sub, 'Tempertature', 'green', False)
plt.show()

plot(sim_times, sim_mag_ni, r'm/m_0', 'orange', True)
plt.show()