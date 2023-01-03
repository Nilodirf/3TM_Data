import numpy as np
from matplotlib import pyplot as plt

#ownplot is a function that reads your datafiles and returns arrays of the delays, Nickels' magnetization, electron-/ and lattice temperature of all layers for every timestep. Time Offset in ps
def ownplot(file, offset):
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

    mags2=[line[5] for line in lines]
    mag2bae= [mag.replace('[0.0]', '').replace('[[', '') for mag in mags2]
    mag2=[m.replace('[', '').replace(']', '') for m in mag2bae]
    m2=[line.split() for line in mag2]
    m2s=np.array([[float(m2[i][j].replace(',', '')) for i in range(len(m2))] for j in range(len(m2[0]))])
    print(len(m1s[0]))
    print(len(m2s[0]))

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

    return(np.array(times), np.array(m1s), np.array(te1), np.array(tp1), np.array(tp21), np.array(m2s))

def plot(file, pump_delay): 
    
    fig, (ax1, ax2)=plt.subplots(2,1)

    delays, m, te1, tp1, tp21, m2 = ownplot(folder + file, delay)

    norm_factor=tp1[-1]-tp1[0]                                  #to normalize the phonon temperature (1 at final delay time)
    
    ax1.plot(delays, (te1-te1[0])/norm_factor, label=r'$T_e$')  #this plots the temperature difference, normalized by the factor defined above
    ax1.plot(delays, (tp1-tp1[0])/norm_factor, label=r'$T_{po}$')
    ax1.plot(delays, (tp21-tp21[0])/norm_factor, label=r'$T_{pa}$')


    ax2.plot(delays, m[0], label='$m_{loc}$')                   #plots magnetization
    ax2.plot(delays, m2[0,0]-m2[0], label='$m_{it}$')
    
    ax1.set_ylabel(r'$\Delta T$ normalized', fontsize=14)
    ax1.set_xlim(-1,15)                                         #limits of the x-axis in ps, here it goes until 15 like in the paper
    ax1.legend()

    ax2.set_xlabel(r'delay [ps]', fontsize=14)
    ax2.set_ylabel(r'$m/m_0$', fontsize=14)
    ax2.sharex(ax1)                                             #sets the timescale of the magnetization plot equal to that of temperature plot. comment this with '#' if you want to see the whole time
    ax2.legend()
    plt.show()
    
folder='put path here'

file='put filename here'
plot(file, ' put pump delay time in ps here')
