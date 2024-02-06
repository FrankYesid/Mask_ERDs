################### Functions  ###################
## 1. ERD
################ Database De Bci #################

###############  ERD   #################
# Inputs:
# ERD,R = ERD_ERS(subjects,y,Xtemp,fs,Class,T1,T2,temporal)
# subjects = X["X"][0,0]
# fs = frecuency X["fs"][0,0]
# Class = type of class. 
#   Class 1: left hand
#   Class 2: right hand
#   Class 3: tongue
#   Class 4: both feet
# T1 = 0 -> Default = time of 0s
# T2 = 2 -> Default = time of 2s
# method = 'dps', spectrogram
# --------------------------------------------------
# Outputs:
# ERD = filtered and mean data
# t = Times
# f = Frecuency
##################################################
# Copyright (C) 2018 Signal Processing and Recognition Group
# F.Y. Zapata Castano
##################################################
#################### library #####################
import numpy as np
import numpy.matlib as npmat
import scipy.io as sio
from scipy import signal
import scipy.linalg
import argparse
import time 
import multiprocessing
import sys
import matplotlib.pyplot as plt
from matplotlib import mlab
###################################################
###################################################
# function of mean 'Et{Pnc(t,f):te[Ta,Tb]}\n'
def R_NC(data,r_nc,ul,up): 
    for y in range(9):
        print 'Sujeto', (y+1), 'de 9'
        r_nc[y,0] = np.mean(data[y][0][:,:,:,ul:up],axis=3)
    return r_nc

###################################################
# function of mean 'En{r_nc(f)}\n'
def R_C(data,r_c):
    for y in range(9):
        print 'Sujeto', (y+1), 'de 9'
        r_c[y,0] = np.mean(data[y][0],axis=0)
    return r_c

###################################################
# function of mean 'En{Pnc(t,f)}\n'
def M_C(data,m_c):
    for y in range(9):
        print 'Sujeto', (y+1), 'de 9'
        m_c[y,0] = np.mean(data[y][0],axis=0)
    return m_c

###################################################
# function of ERDs='ERD = m_c/r_c\n'
def m_erd(data2,data1,ERDs):
    ERtm = np.copy(ERDs)
    for y in range(9):
        print 'Sujeto', (y+1), 'de 9'  
        tem1 = np.copy(data2[y][0])
        tem2 = np.copy(np.tile(data1[y][0].T,(88,1,1)))
        for r in range(22):
            ERDs[y,0][r,:,:] = ((tem1[r,:,:]*(1/((tem2[:,:,r]).T)))-1)
        # ERDs[y,0] = ERtm
    return ERDs

###################################################
def ERD(subjects,y,Xtemp,fs,Class,T1,T2,temporal):
    cop = np.copy(Xtemp)
    lon = len(subjects)
    Nw = np.floor(0.7*fs)
    win = np.hamming(Nw)
    NFT = len(win)
    nover = np.floor(0.9*Nw)
    Subjects_class = Xtemp
    for su in range(9):
        Subjects_class = subjects[su,0][y[su,0] == np.int(Class)]
        N_trial = len(Subjects_class)
        channels = np.array(range(22))
        N_channels = len(channels)
        (datf,datt,tem) = signal.spectrogram(Subjects_class[0][:,0],fs=fs,window=win,noverlap=nover)
        # print tem.shape
        tem = np.int(tem[1].shape[0])
        nn = np.int(np.floor((Nw/2)+1))
        Xtemp[su] = [np.zeros((N_trial,N_channels,nn,tem))] # Trials x canal x frecuencia x tiempo
        for tri in range(N_trial):
            print 'Sujeto', (su+1), 'de 9 ...trail : ', (tri+1), 'de', (N_trial)
            for cnl in range(N_channels):
                # Signal
                sig = Subjects_class[tri][:,channels[cnl]]
                # Calcular STFT
                (f,t,X_class) = signal.spectrogram(sig,fs=fs,window=win,noverlap=nover)
                X_class1 = np.abs(X_class)
                Xtemp[su,0][tri,cnl,:,:] = X_class1
    
    # temporal["pot"] = np.copy(Xtemp)
    temporal["t"] = np.copy(t)
    temporal["f"] = np.copy(f)
    # abc = np.mean(Xtemp[0][0][:,:,:,0:23],axis=3)
    # print abc.shape
    
    subj = np.copy(Xtemp)
    subj1= np.copy(Xtemp)
    # print subj[0][0].shape
    temp_1 = np.abs(t-np.int(T1))
    temp_2 = np.abs(t-np.int(T2))
    min_1 = np.min(temp_1)
    min_2 = np.min(temp_2)
    ul = np.nonzero(temp_1 == min_1)[0][0]
    up = np.nonzero(temp_2 == min_2)[0][0]

    r_nc = np.copy(cop)
    print '\n'
    print 'Et{Pnc(t,f):te[Ta,Tb]}\n'
    r_nc = R_NC(subj,r_nc,ul,up)
    print r_nc[0][0].shape

    r_c = np.copy(cop)
    print '\n'
    print 'En{r_nc(f)}\n'
    r_c = R_C(r_nc,r_c)
    print r_c[0][0].shape
    # temporal["rc"] = r_c
    data1 = np.copy(r_c)

    m_c = np.copy(cop)
    print '\n'
    print 'En{Pnc(t,f)}\n'
    m_c = M_C(subj1,m_c)
    print m_c[0][0].shape
    # temporal["mc"] = m_c
    data2 = np.copy(m_c)

    ERDs = np.copy(m_c)
    print '\n'
    print 'ERD = m_c/r_c\n'
    ERDs = m_erd(data2,data1,ERDs*0)

    temporal["ERD"] = np.copy(ERDs)
    return temporal
###################################################
###################### Main #######################
if __name__ == '__main__':
    ## Entradas de variables
    parser = argparse.ArgumentParser(description = "Program")
    parser.add_argument("Arch", help = "The potencials")
    parser.add_argument("erd or ers", help = "", nargs = "?")
    parser.add_argument("Class", help = "", nargs = "?")   
    parser.add_argument("t1", help = "time 1 range", nargs = "?")
    parser.add_argument("t2",  help = "time 2 range", nargs = "?")
    args = parser.parse_args()
    Database = args.Arch              ## Carga de base de datos.
    database = sio.loadmat(Database)  ## 
    temporal = database
    fs = database["fs"][0,0]
    ## Manejo de la entrada de funcion solicitada.
    program_name = sys.argv[0]
    arguments = sys.argv[1:]
    count = len(arguments)
    name = 'n'
    for parametros in sys.argv:
        if parametros == 'erd':
            Class = arguments[arguments.index("erd")+1]
            T1 = arguments[arguments.index("erd")+2]
            T2 = arguments[arguments.index("erd")+3]
            subjects = temporal["X"]
            y = temporal["y"]
            Xtemp = subjects*0
            arch = ERD(subjects,y,Xtemp,fs,Class,T1,T2,temporal)
            print('Done ERDs!')
            name = 'BCI_ERDs.mat'

    if name == 'n':
        print "No se realizo ningun cambio"
    else:
        print ''
        sio.savemat(name,arch) # Save database.
        print 'Nombre del archivo es: %s' % (name)
    print '\nDone Program!\n'    