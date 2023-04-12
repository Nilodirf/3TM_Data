import numpy as np
import h5py
import scipy
from scipy import io
from matplotlib import pyplot as plt
import os
import math


def templot():
    f = scipy.io.loadmat('Temperatures.mat')

    delay=f['Delay']/1e3

    f1=f['F0p6']
    f2=f['F1p2']
    f3=f['F1p8']
    f4=f['F2p5']
    f5=f['F3p0']

    df1=f['Err_F0p6']
    df2=f['Err_F1p2']
    df3=f['Err_F1p8']
    df4=f['Err_2p5']
    df5=f['Err_F3p0']

    return(delay, f1, f2, f3, f4, f5, df1, df2, df3, df4, df5)

def magplot():
    f= open('mag.txt', 'r')
    fl=f.readlines()
    fd=[line.replace('\n', '').split(',') for line in fl]
    ad=np.array([[float(dp) for dp in line] for line in fd])
    delay=ad[:,0]
    diff=ad[:,1]
    differr=ad[:,2]

    return(delay, diff/diff[0], differr)

def simplot(dat):
    input=open(dat, 'r')
    content=input.read().split('\n')
    toplot=[line for line in content if not str(line).startswith('#')]
    columns=[line.split() for line in toplot]
    del columns[-1]
    times=np.array([float(line[0]) for line in columns])
    te=[float(line[2]) for line in columns]
    tp=[float(line[3]) for line in columns]
    nmag=np.array([float(line[1]) for line in columns])
    return(times, np.array(nmag), np.array(te), np.array(tp))

def fitplot():
    tem=templot()
    mag=magplot()
    sim5=simplot('test5.dat')
    sim6=multiplot('C:/Users/tgrie/OneDrive/Dokumente/Github/3TM_Data/FGT/multilattice_te.dat', 25)
    sim5mag=simplot('test5mag.dat')
    sim6mag=multiplot('C:/Users/tgrie/OneDrive/Dokumente/Github/3TM_Data/FGT/multilattice_mag.dat', 25)
    koop=simplot('koopmans.dat')


    fig, (ax1, ax2)=plt.subplots(2,1)
    ax1.scatter(tem[0]-0.015, tem[1])
    ax1.plot(tem[0]-0.015, tem[1]-tem[6], color='red', ls='dashed', alpha=0.5)
    ax1.plot(tem[0]-0.015, tem[1]+tem[6], color='red', ls='dashed', alpha=0.5)
    #ax1.plot(koop[0], koop[2])
    #ax1.plot(sim3[0], sim3[2])
    ax1.plot(sim5[0], sim5[2])
    ax1.plot(sim6[0], sim6[2], color='purple', label=r'94 J/cm^3')

    ax2.scatter(mag[0]+0.12, mag[1])
    ax2.plot(mag[0]+0.12, mag[1]-mag[2], color='red', ls='dashed', alpha=0.5)
    ax2.plot(mag[0]+0.12, mag[1]+mag[2], color='red', ls='dashed', alpha=0.5)
    #ax2.plot(koop[0], koop[1]/0.8975, lw=2.0, label='lit. params')
    ax2.plot(sim5mag[0], sim5mag[1]/0.8975, label='new params')
    ax2.plot(sim6mag[0], sim6mag[1][0]/0.931/1.04, lw=2.0, color='purple', label=r'470 J/cm^3')

    ax1.set_ylabel(r'$T_e$ [K]', fontsize=16)
    ax2.set_ylabel(r'$m/m(t=0)$', fontsize=16)
    ax2.set_xlabel(r'delay [ps]', fontsize=16)

    ax1.set_xlim(-0.5, 5)
    ax2.set_xlim(-0.5, 5)
    ax1.legend(fontsize=14, bbox_to_anchor = (1., 1.))
    ax2.legend(fontsize=14, bbox_to_anchor = (0.6, 0.7))
    plt.savefig('comparison_te_mag.pdf')
    plt.show()

def koopplot():
    sim=simplot('koopmans.dat')
    tem=templot()
    mag=magplot()

    fig, (ax1, ax2)=plt.subplots(2,1, sharex=True)
    ax1.plot(sim[0], sim[1]/0.8975, lw=2.0, color='green')
    ax1.scatter(mag[0]+0.15, mag[1])


    ax2.plot(sim[0], sim[2], lw=2.0, color='red')
    #ax2.plot(sim[0], sim[3], lw=2.0, color='blue')
    ax2.scatter(tem[0], tem[1])

    ax1.set_xlim(-0.5,5)
    ax1.set_ylim(0.6, 1.05)

    ax1.set_ylabel(r'$m/m(t=0)$', fontsize=16)
    ax2.set_ylabel(r'T [K]', fontsize=16)
    ax2.set_xlabel(r'delay [ps]', fontsize=16)

    plt.savefig('koopfit.pdf')
    plt.show()

def msdplot():
    f=open('first_attempt_MSD.txt', 'r')
    dat=f.read().split('\n')
    lines=[i.split('\t') for i in dat]
    times=[float(j[0])-0.8 for j in lines]
    msd=[float(j[1]) for j in lines]
    plt.xlabel(r'delay [ps]', fontsize=16)
    plt.ylabel(r'normalized MSD [arb. units]', fontsize=16)
    return(times, msd)

def e_lat(tein, t, cp0):
    cph=cp0*(tein/t)**2*np.exp(tein/t)/(np.exp(tein/t)-1)**2
    return(cph*t)

#multiplot is a function that reads your datafiles and returns arrays of the delays, Nickels' magnetization, electron-/ and lattice temperature of all layers for every timestep. Time Offset in ps
def multiplot(file, offset):
    ### open(file, mode) is an in-built function. 1st argument is the filepath, second argument is the mode. Choose mode 'r' to read the desired file
    input = open(file,'r')
    content = input.read().split('\n')
    toplot = [line for line in content if not str(line).startswith('#')]
    lines = [line.split('\t') for line in toplot]
    del lines[-1]

    times = np.array([-offset+ float(line[0]) for line in lines])

    mags = [line[1] for line in lines]
    mag1bae= [mag.replace('[0.0]]', '').replace('[[', '') for mag in mags]
    mag1=[m.replace('[', '').replace(']', '') for m in mag1bae]
    m1=[line.split() for line in mag1]
    m1s=np.array([[float(m1[i][j].replace(',', '')) for i in range(len(m1))] for j in range(len(m1[0]))])
    ####m1s[i][:] is magnetization for all times of (i-1)th monolayer!###
    ###m1s[:][i] is magnetization for timestep i of all monolayers!###

    temes = [line[2].split('], [') for line in lines]
    tebae=[[te[i].replace('[', '').replace(']', '') for i in range(len(temes[0]))] for te in temes]
    te1bla=[tb[:][0] for tb in tebae]

    te1n=[b.replace(',', '') for b in te1bla]
    te1m=[te.split() for te in te1n]
    te1=[[float(i) for i in te] for te in te1m]


    tpmes = [line[3].split('], [') for line in lines]
    tpbae=[[tp[i].replace('[', '').replace(']', '') for i in range(len(tpmes[0]))] for tp in tpmes]
    tp1bla=[tb[:][0] for tb in tpbae]
    tp1n=[b.replace(',', '') for b in tp1bla]
    tp1m=[tp.split() for tp in tp1n]
    tp1=[[float(i) for i in tp] for tp in tp1m]


    tp2s=[line[4].split('], [') for line in lines]
    tp2bae=[[tp[i].replace('[', '').replace(']', '') for i in range(len(tpmes[0]))] for tp in tp2s]
    tp21bla=[tb[:][0] for tb in tp2bae]
    tp21n=[b.replace(',', '') for b in tp21bla]
    tp21m=[tp.split() for tp in tp21n]
    tp21=[[float(i) for i in tp] for tp in tp21m]
    return(np.array(times), np.array(m1s), np.array(te1), np.array(tp1), np.array(tp21))

#fitplot()

def relaxtimes():
    folderpath='C:/Users/tgrie/OneDrive/Dokumente/GitHub/3TM_Data/FGT/powertests'
    folderpathbla='C:/Users/tgrie/OneDrive/Dokumente/GitHub/3TM_Data/FGT/powertests/'
    files=os.listdir(folderpath)
    taue=[]
    taup=[]

    for file in files:
        tecond=5e-3
        tpcond=5e-4
        t, m, te, tp=simplot(folderpathbla + str(file))
        
        ten=np.roll(te, -1)
        tpn=np.roll(tp,-1)
        tpos= t>0
        t=t[tpos]

        terel=np.where(abs(te[tpos]-ten[tpos])<tecond)
        tprel=np.where(abs(tp[tpos]-tpn[tpos])<tpcond)
        #terel=np.where(abs(te[tpos]-te[-1])<tecond)
        #tprel=np.where(abs(tp[tpos]-tp[-1])<tpcond)
        pp=float(str(file).replace('new.dat', '').replace('pp', ''))*math.sqrt(2*math.pi)*15e-15*1e-6
        taue.append([pp, min(t[terel])])
        taup.append([pp, min(t[tprel])])
        
    taue=np.array(taue)
    taup=np.array(taup)
    pps=[taue[i,0] for i in range(len(taue))]
    taues=[taue[i,1] for i in range(len(taue))]
    taups=[taup[i,1] for i in range(len(taup))]

    #plt.scatter(pps, taues, label=r'$T_e$')
    plt.scatter(pps, taups)
    plt.xlabel(r'abs. energy density [J/cm^3]', fontsize=16)
    plt.ylabel(r'Lattice relaxation times [ps]', fontsize=16)
    plt.savefig('relaxtimes.pdf')
    plt.show()

def teqs():
    folderpath='C:/Users/tgrie/OneDrive/Dokumente/GitHub/3TM_Data/FGT/cps'
    folderpathbla='C:/Users/tgrie/OneDrive/Dokumente/GitHub/3TM_Data/FGT/cps/'
    files=os.listdir(folderpath)
    teq=[]

    for file in files:
        t, m, te, tp1, tp2=multiplot(folderpathbla+str(file), 25)
        tint=abs(t-1.)<1e-3
        print(t[tint])
        cp=float(str(file).replace('cp', '').replace('ppte.dat', '').replace('ppmag.dat', ''))/1e6
        teq.append([cp, te[-1], te[tint]])
        #plt.plot(t, (e_lat(0.75*460, tp1, 0.1)+e_lat(0.75*460, tp2, 0.9)-40)/113.5, label=r'$\Delta(C_pT_p)$', color='purple')
        #plt.show()

    teq=np.array(teq)
    plt.scatter(teq[:,0], teq[:,1], label='delay 20 ps')
    plt.scatter(teq[:,0], teq[:,2], label='delay 1 ps')
    plt.xlabel(r'Maximum lattice specific heat [MJ/m^3]', fontsize=16)
    plt.ylabel(r'Electron temperature at different delays')
    plt.legend(fontsize=14)
    plt.savefig('teofcp.pdf')
    plt.show()

def ks():
    folderpath='C:/Users/tgrie/OneDrive/Dokumente/Github/3TM_Data/FGT/ks'
    folderpathbla='C:/Users/tgrie/OneDrive/Dokumente/Github/3TM_Data/FGT/ks/'
    files=os.listdir(folderpath)

    for file in files:
        t, m, te, tp1, tp2=multiplot(folderpathbla+str(file), 25)
        k=float(str(file).replace('k', '').replace('.dat', ''))
        if 10*k%2==0 or k==0.1 or k==0.05:
            elat=e_lat(0.75*460, tp1, k)+e_lat(0.75*460, tp2, (1-k))
            plt.plot(t, (elat-elat[0])/(elat-elat[0])[-1], label='k=' + str(k))

    plt.scatter(np.array(msdplot()[0]), msdplot()[1])
    plt.xlim(-1, 20)
    plt.xlabel(r'delay [ps]', fontsize=16)
    plt.ylabel(r'normalized lattice energy', fontsize=16)
    plt.legend(fontsize=12)
    plt.savefig('ks.pdf')
    plt.show()
                   

#t, m, te, tp, tp2=multiplot('C:/Users/tgrie/OneDrive/Dokumente/Github/3TM_Data/FGT/multilattice_ep.dat', 25)
#plt.scatter(msdplot()[0], msdplot()[1], color='magenta')
#plt.plot(t, (e_lat(0.75*460, tp, 0.1)+e_lat(0.75*460, tp2, 0.9)-40)/44.2, color='red')
#plt.savefig('multi_latdyn.pdf')
#plt.xlim(-0.5, 20)
#plt.text(r'300 mJ/cm^2')
#plt.show()

#plt.plot(t, te)
#plt.scatter(tel[0], tel[6])
#plt.show()

#teqs()
fitplot()
#ks()
#relaxtimes()

#t, m, te, tp, tp2=multiplot('C:/Users/tgrie/Onedrive/Dokumente/Github/3tm_results/10cp90gpp2.dat', 25)
#plt.plot(t, te)
#plt.plot(t, tp)
#plt.plot(t, tp2)

#plt.xlabel(r'delays [ps]', fontsize=16)
#plt.ylabel(r'T [K]', fontsize=16)
#plt.xlim(-1, 5)
#plt.legend(['282 J/cm^3'], fontsize=14)
#plt.savefig('blashit.pdf')
#plt.show()
