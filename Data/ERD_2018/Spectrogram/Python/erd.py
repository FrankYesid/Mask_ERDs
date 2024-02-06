##################  Functions  ###################
## 1. ERD-ERS
################ Database De Bci #################
 
#############  Function ERD_ERS  #################
# Inputs:
# ERD,R = ERD_ERS(subject,fs,ro,k,l,method)
# subject = X["X"][0,0]
# fs = frecuency X["fs"][0,0]
# ro = 0 -> time of 0s
# k = 2  -> time of 2s
# method = "power" or "var"
# --------------------------------------------------
# Outputs:
# ERD = filtered and mean data
# R =  mean and sum of channel
##################################################

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

################################################## Filter of ERD:ERS
def filter2(Freq,n,X,fs):
    if n<5:
        n = 5
    ## Filter_desing
    wn_freq = [(Freq[0]/(0.5*fs)), (Freq[1]/(0.5*fs))]
    (b,a) = signal.butter(n,wn_freq,'bandpass')
    tr = len(X)
    Xfreq = X*0
    for k in range(0,tr): ## 273
        (s,c) = X[k,0].shape
        for ch in range(0,c): ## 22
            x = X[k,0][:,ch]
            Xfreq[k,0][:,ch] = signal.filtfilt(b,a,x,axis=-1,padtype='odd',padlen=None,method='pad')
    return Xfreq
###################################################

################################################### Function ERD:ERS
def ERD_ERS(subjects,fs,ro,k,l,method,A,B):
    print 'inicio ERD_ERS'
    for r in range(0,len(A)):
        print 'subject %s de 9' %(r+1)
        subject = subjects[r,0]
        trl = len(subject)
        (s,c) = subject[0,0].shape
        ## reference - time 7s
        t = np.array(range(0,7*fs,1))/float(fs)
        ul = np.where(t == ro)
        up = np.where(t == k)
        Freq = []
        erd = subject[0:4,0]
        R = []
        ## filters

        ## alpha 8-14
        Freq.append(8.5)
        Freq.append(13.5)
        alpha = filter2(Freq,5,subject,fs)
        # print "filter alpha"
        alpha_tmp = np.array(alpha[0,0])
        for a in range(1,trl):
            alpha_tmp = np.vstack((alpha_tmp,alpha[a,0]))
        al_t = np.reshape(alpha_tmp,(trl,s,c))

        ## beta 14-40
        Freq[0] = 14.5
        Freq[1] = 40
        beta = filter2(Freq,5,subject,fs)
        # print "filter beta"
        beta_tmp = np.array(beta[0,0])
        for a in range(1,trl):
            beta_tmp = np.vstack((beta_tmp,beta[a,0]))
        bt = np.reshape(beta_tmp,(trl,s,c))

        ## delta 0.5-4
        Freq[0] = 1
        Freq[1] = 4
        delta = filter2(Freq,5,subject,fs)
        # print "filter delta"
        delta_tmp = np.array(delta[0,0])
        for a in range(1,trl):
            delta_tmp = np.vstack((delta_tmp,delta[a,0]))
        dl = np.reshape(delta_tmp,(trl,s,c))

        ## theta 4-8
        Freq[0] = 4.5
        Freq[1] = 7.5
        theta = filter2(Freq,5,subject,fs)
        # print "filter theta"
        theta_tmp = np.array(theta[0,0])
        for a in range(1,trl):
            theta_tmp = np.vstack((theta_tmp,theta[a,0]))
        tt = np.reshape(theta_tmp,(trl,s,c))

        # metodo aplicado
        if method == 'power':
            m_al = np.mean(al_t**2,axis=0)
            m_bt = np.mean(bt**2,axis=0)
            m_dl = np.mean(dl**2,axis=0)
            m_tt = np.mean(tt**2,axis=0)
            
            R_al = np.sum(m_al[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            R_bt = np.sum(m_bt[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            R_dl = np.sum(m_dl[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            R_tt = np.sum(m_tt[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            
            aa = np.array(m_al - npmat.repmat(R_al,s,1))
            bb = np.array(m_bt - npmat.repmat(R_bt,s,1))
            cc = np.array(m_dl - npmat.repmat(R_dl,s,1))
            dd = np.array(m_tt - npmat.repmat(R_tt,s,1))

            erd_al = 100*(np.divide((aa),(npmat.repmat(R_al,s,1))))
            erd_bt = 100*(np.divide((bb),(npmat.repmat(R_bt,s,1))))
            erd_dl = 100*(np.divide((cc),(npmat.repmat(R_dl,s,1))))
            erd_tt = 100*(np.divide((dd),(npmat.repmat(R_tt,s,1))))

            ## vector de matrices ##
            erd[0] = erd_al
            erd[1] = erd_bt 
            erd[2] = erd_dl 
            erd[3] = erd_tt 

            ## vector de vectores ##
            R.append(R_al)
            R.append(R_bt)
            R.append(R_dl)
            R.append(R_tt)

            print "Done! method power"

        elif method == 'var':
            s_al = np.squeeze(np.std(al_t,[],2))
            s_bt = np.squeeze(np.std(bt,[],2))
            s_dl = np.squeeze(np.std(dl,[],2))
            s_tt = np.squeeze(np.std(tt,[],2))
            
            R_al = np.sum(s_al[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            R_bt = np.sum(s_bt[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            R_dl = np.sum(s_dl[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            R_tt = np.sum(s_tt[ul[0][0]:up[0][0],:],0)/(up[0][0]-ul[0][0]+1)
            
            erd_al = 100*((s_al - npmat.repmat(R_al,s,1))/(npmat.repmat(R_al,s,1)))
            erd_bt = 100*((s_bt - npmat.repmat(R_bt,s,1))/(npmat.repmat(R_bt,s,1)))
            erd_dl = 100*((s_dl - npmat.repmat(R_dl,s,1))/(npmat.repmat(R_dl,s,1)))
            erd_tt = 100*((s_tt - npmat.repmat(R_tt,s,1))/(npmat.repmat(R_tt,s,1)))
            
            ## vector de matrices ##
            erd[0] = erd_al 
            erd[1] = erd_bt 
            erd[2] = erd_dl 
            erd[3] = erd_tt 
            ## vector de vectores ##
            R[0] = R_al   
            R[1] = R_bt
            R[2] = R_dl
            R[3] = R_tt    

            print "Done! method var"

        else:
            print('this method is not defined')

        A[r,0] = erd
        print A.shape
        B[r,0] = R
        # if r == 8:        
        #     break
    temporal["ERD"] = A
    temporal["R"] = B
    temporal["t"] = t
    return temporal
###################################################

###################################################
def calculo_erd_ers(subjects,ro,k,l,fs,method,temporal):
    temp1 = temporal["X"]*0
    temp2 = temporal["X"]*0
    (temporal) = ERD_ERS(subjects,fs,ro,k,l,method,temp1,temp2)   
    return temporal
###################################################

###################### Main #######################
if __name__ == '__main__':

 ## Entradas de variables
    parser = argparse.ArgumentParser(description = "Program_ERD_ERS")
    parser.add_argument("Arch", help = "The potencials")
    parser.add_argument("erd", help = "Event-related (desynchronazation(ERD), synchronazation(ERS)), ", nargs = "?")
    parser.add_argument("ro", help = "initial time", nargs = "?")
    parser.add_argument("k", help = "Final time", nargs = "?")
    args = parser.parse_args()
    Database = args.Arch              ## Carga de base de datos.
    database = sio.loadmat(Database)  ## 
    temporal = database
    subjects = temporal["X"]
    N=len(subjects[0,0][0,0])
    fs = database["fs"][0,0]
    ## Manejo de la entrada de funcion solicitada.
    program_name = sys.argv[0]
    arguments = sys.argv[1:]
    temporal.keys()
    count = len(arguments)
    for parametros in sys.argv:
        if parametros == 'erd':
            Ro = arguments[arguments.index("erd")+1]
            k = arguments[arguments.index("erd")+2]
            method = 'power'
            l=N/fs
            temporal = calculo_erd_ers(subjects,float(Ro),float(k),l,fs,method,temporal)
            print('Done ERD_ERS!')
    sio.savemat('BCI_ERD_ERS.mat',temporal)
    print "Nombre del archivo: BCI_ERD_ERS.mat"
    print "Done Program!"