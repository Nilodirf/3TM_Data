import numpy as np
import h5py
import scipy
from scipy import io
from matplotlib import pyplot as plt


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
    return(times, nmag, te, tp)

def fitplot():
    tem=templot()
    mag=magplot()
    sim=simplot('test1.dat')
    sim5=simplot('test5.dat')
    sim6=simplot('test6.dat')
    sim5mag=simplot('test5mag.dat')
    sim6mag=simplot('test6mag.dat')
    koop=simplot('koopmans.dat')


    fig, (ax1, ax2)=plt.subplots(2,1)
    ax1.scatter(tem[0]-0.02, tem[1])
    ax1.plot(tem[0]-0.02, tem[1]-tem[6], color='red', ls='dashed', alpha=0.5)
    ax1.plot(tem[0]-0.02, tem[1]+tem[6], color='red', ls='dashed', alpha=0.5)
    ax1.plot(koop[0], koop[2])
    #ax1.plot(sim3[0], sim3[2])
    ax1.plot(sim5[0], sim5[2])

    ax2.scatter(mag[0], mag[1])
    ax2.plot(mag[0], mag[1]-mag[2], color='red', ls='dashed', alpha=0.5)
    ax2.plot(mag[0], mag[1]+mag[2], color='red', ls='dashed', alpha=0.5)
    ax2.plot(koop[0], koop[1]/0.8975, lw=2.0, label='lit. params')
    ax2.plot(sim5mag[0], sim5mag[1]/0.8975, label='new params')
    #ax2.plot(sim6mag[0], sim6mag[1]/0.8975, lw=2.0, color='green', label='higher mag. rate')

    ax1.set_ylabel(r'$T_e$ [K]', fontsize=16)
    ax2.set_ylabel(r'$m/m(t=0)$', fontsize=16)
    ax2.set_xlabel(r'delay [ps]', fontsize=16)

    ax1.set_xlim(-1, 5)
    ax2.set_xlim(-1, 10)
    plt.legend()
    plt.savefig('attempt1.pdf')
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
    plt.xlabel(r'delay [ps]')
    plt.ylabel(r'normalized MSD [arb. units]')
    return(times, msd)

tptimes, msd= msdplot()
te=templot()

plt.scatter(tptimes, msd)
plt.show()
