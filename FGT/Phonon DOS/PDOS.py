import numpy as np
import math
import os
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import constants as sp


#define grids for cp computation
E_new=np.arange(0, 50+1e-4, 1e-4)
Ts=np.arange(1,1000,1)


#prepare some variables for the DOS computation from the band spectra
deltE=5e-1
deltEJ=deltE*sp.electron_volt
deltk=1e-4
V_unit_cell=3.9*3.9*16.3*1e-30

def prep_Es(dE):
    E_range=np.arange(0, 50+dE, dE)
    return(E_range)

#get the grids for DOS computation:
E_range= prep_Es(deltE)

def k_discrete(dk):
    return np.arange(0.,16., dk)

k_disc=k_discrete(deltk)

def E_discrete(E, kspace, E_round_decimal):
    E_discrete=np.round_(E(kspace), E_round_decimal)
    return(E_discrete)


def Bose_Einstein(T, E):
    bef=1/(np.exp(E*sp.eV/sp.k/T[:, np.newaxis])-1)
    return(bef)


def get_bandstructure(method):
    band_files = os.listdir('C:/Users/tgrie/OneDrive/Dokumente/GitHub/3TM_Data/FGT/Phonon DOS/bands')
    band_list=[]
    ac_band_list = []
    op_band_list = []
    for bf in band_files:
        b_s=str(bf)
        dat=open('bands/'+ str(bf), 'r').readlines()[:-1]
        dat_no_lineskips=[line.replace('\n', '') for line in dat]
        dat_singled=[line.split(';') for line in dat_no_lineskips]
        raw_ks=[float(line[0].replace(',', '.')) for line in dat_singled]
        ks=math.pi/3.9e-10*np.array(raw_ks)*1e-9    #rad/nm  #time math.pi because band data is saved on scale [0,2]
        raw_es=[float(line[1].replace(',', '.')) for line in dat_singled]
        es=sp.hbar*2*math.pi*sp.c*np.array(raw_es)/sp.electron_volt*100*1000   #fist from J to eV, then from 1/cm to 1/m, then from eV to meV
        if ks[-1]<16:
            ks[-1]=16
        if ks[0]>0:
            ks[0]=0
        f=interpolate.interp1d(ks, es)
        E_disc = E_discrete(f, k_disc, 4)
        if method== 'opac':
            if b_s.endswith('17.txt') or b_s.endswith('18.txt') or b_s.endswith('16.txt'):
                ac_band_list.append(E_disc)
            else:
                op_band_list.append(E_disc)
        elif method == 'cutoff':
            band_list.append(E_disc)
    if method=='cutoff':
            return band_list
    if method=='opac':
            return ac_band_list, op_band_list

def get_DOS(bands, dE):
    E_hist=prep_Es(deltE)
    for band in bands:
        for point in band:
            E_hist[int(point*1/dE)]+=1
    DOS=E_hist*deltk/deltE
    return DOS

def cp(dos, e, bef, dT):
    #ep=np.sum(e*dos*bef*deltE,axis=-1)/V_unit_cell
    ep=np.trapz(e*dos*bef, e)
    cp=np.diff(ep)/dT
    return(cp)

def cp_cutoff(DOS, cutoff):
    cps=[]
    E1 = E_range[1:int(cutoff)]*1e-3*sp.electron_volt
    E2 = E_range[int(cutoff + 1):]*1e-3*sp.electron_volt
    DOS1 = DOS[1:int(cutoff)]
    DOS2 = DOS[int(cutoff + 1):]
    Es = [E1, E2]
    DOSs = [DOS1, DOS2]

    for i, E in enumerate(Es):
        cps.append(cp(DOSs[i], E, Bose_Einstein(Ts, E), 1))
    return cps

def get_full_cp(DOS, T, E):
    cps=[]
    cps.append(cp(DOS, E, Bose_Einstein(T, E), 1)) #J/UC/K
    cps=np.array(cps)/V_unit_cell*sp.eV
    return cps

def cp_opac(DOSac, DOSop):
    cph=[]
    DOSs=[DOSac, DOSop]
    E=E_range[1:]*1e-3*sp.electron_volt
    for dos in DOS:
        cph.append(cp(dos, E, Bose_Einstein(Ts, E), 1))
    return cph


def plot_bands_DOS(model, bands, DOS, bands_ac=None, bands_op=None, DOS_ac=None, DOS_op=None, cutoff=None):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    if model=='cutoff':
        for band in bands:
            ax1.plot(k_disc, band, color='blue', lw=2.0)
            ax2.plot(DOS[1:cutoff], E_range[1:cutoff], color='blue')
            ax2.plot(DOS[cutoff + 1:], E_range[cutoff+1:], color='red')

    elif model=='opac':
        for band in bands_ac:
            ax1.plot(k_disc, band, lw=2.0, color='orange')
            ax2.plot(DOS_ac, E_range, color='orange')
        for band in bands_op:
            ax1.plot(k_disc, band, lw=2.0, color='purple')
            ax2.plot(DOS_op, E_range, color='purple')
        ax2.plot(DOS, E_range, color='black', label=r'total')

    ax1.set_xlabel(r'k', fontsize=16)
    ax1.set_ylabel(r'E [meV]', fontsize=16)
    ax1.set_xlim(0, 16)
    ax1.set_ylim(0, 50)
    ax1.xaxis.set_ticks([6.7, 10.2])
    ax1.xaxis.set_ticklabels(['K','M'])
    ax1.vlines(6.7, 0, 50, linestyle='dashed', color='black', alpha=0.5, linewidth=0.8)
    ax1.vlines(10.2, 0, 50, linestyle='dashed', color='black', alpha=0.5, linewidth=0.8)

    ax2.set_xlim(0,)
    ax2.sharey(ax1)
    ax2.set_xlabel(r'DOS [1/eV]', fontsize=16)
    plt.savefig('PDOS_cutoff_2.pdf')
    plt.show()

def plot_cp_opac(cpac, cpop):
    plt.plot(Ts[:-1], cpac, label='accoustic phonons')
    plt.plot(Ts[:-1], cpop, label='optical phonons')
    plt.xlabel(r'Temperature [K]', fontsize=16)
    plt.ylabel(r'Lattice heat capcity [J/m$^3$K]', fontsize=16)
    plt.title(r'Heat capacity of different modes', fontsize=18)
    plt.legend(fontsize=16)
    plt.show()

def plot_cp_cutoff(cphs, cutoff):
    plt.plot(Ts[:-1], cphs[0], color= 'blue', label='lower package')
    plt.plot(Ts[:-1], cphs[1], color='red', label='higher package')
    plt.plot(Ts[:-1], np.sum(cphs, axis=0), color='black', ls='dashed', label='total')
    plt.xlabel(r'Temperature [K]', fontsize=16)
    plt.ylabel(r'Lattice Heat capacity [J/m$^3$K]', fontsize=16)
    plt.title(r'cutoff energy=' + str(cutoff * deltE) + 'meV', fontsize=18)
    plt.legend(fontsize=16)
    plt.xlim(0,1000)
    plt.ylim(-0.01e7,1.65e7)
    #plt.savefig('cp_cutoff_2.pdf')
    plt.show()

def get_ab_in_DOS():
    dat=open('ab_in_uk.txt', 'r').readlines()[:-1]
    dat_no_lineskips = [line.replace('\n', '') for line in dat]
    dat_singled = [line.split(';') for line in dat_no_lineskips]
    E_in_cminv = [float(line[0].replace(',', '.')) for line in dat_singled]
    es =np.array(E_in_cminv)*sp.h*sp.c*100/sp.eV #eV
    dos = np.array([float(line[1].replace(',', '.')) for line in dat_singled])
    norm_factor=np.trapz(dos, es) #states/cm
    dos=dos/norm_factor*3*6 #(states/eV/cm)/(states/cm)=1/eV
    print(np.trapz(dos,es))
    T_range=np.arange(1,1000, 1) #K
    plt.plot(es, dos)
    plt.show()
    cp=get_full_cp(dos, T_range, es)[0]
    plt.plot(T_range[:-1], cp/1e6)
    plt.ylim(0,1.1)
    plt.xlim(0,1000)
    plt.xlabel(r'$T_p$ [K]', fontsize=16)
    plt.ylabel(r'$C_p(T)$ [MJ/m$^3$K]', fontsize=16)
    plt.savefig('cp_abin.pdf')
    plt.show()

    return E_range, DOS

es, dos=get_ab_in_DOS()

# #now run the stuff, get the data:
#
# #first for cutoff method:
# E_cutoff= int((35+deltE)/deltE)
# band_list=get_bandstructure('cutoff')
# DOS = get_DOS(band_list, deltE)
#cph_cutoff=cp_cutoff(DOS, E_cutoff)
#
# #then for optical vs. accoustic modes:
# #bands_ac, bands_op=get_bandstructure('opac')
# #DOS_ac= get_DOS(bands_ac, deltE) * deltk / 2 / math.pi / deltE / 1e-3
# #DOS_op= get_DOS(bands_op, deltE) * deltk / 2 / math.pi / deltE / 1e-3
# #cph_opac=cp_opac(DOS_ac, DOS_op)
# #cp_ac=cph_opac[0]
# #cp_op=cph_opac[1]
#
# #plot what you want:
# #plot_bands_DOS('opac', band_list, DOS, bands_ac, bands_op, DOS_ac, DOS_op)
# plot_bands_DOS('cutoff', band_list, DOS, cutoff=E_cutoff)
# #plot_cp_opac(cp_ac, cp_op)
# plot_cp_cutoff(cph_cutoff, E_cutoff)
