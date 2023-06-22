#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 16:12:11 2019 
Modif   on Wed Jul  30 02:21:15 2019 
@author: Frank Yesid Zapata Castano.
"""

"""
Código que calcula la sincronización para la base de datos Giga_science

Ejemplo: 
    
Pendientes:
    - organizar los trials malos para cada una de las clases.    
"""

# Modules 
import numpy as np
import numpy.matlib as npm
from scipy.io import loadmat, savemat
import scipy.signal as signal
import scipy.special as special
#import mne
import os

###############################################################################
#                       Functions (ERD/ERS) patterns                          #
###############################################################################


def trigg(sch, TRIG, pre, post, gap):
    #% TRIGG cuts continous sequence into segments.
    #% Missing values (in case s is to short) are substituted by NaN's.
    #%
    #% [X,sz] = trigg(s, TRIG, PRE, PST [, GAP])
    #%
    #% INPUT:
    #%  S   is the continous data sequence (1 channel per column)
    #%  TRIG    defines the trigger points
    #%  PRE     offset of the start of each segment (relative to trigger)
    #%  PST     offset of the end of each segment (relative to trigger)
    #%  GAP number of NaN's to separate trials (default=0)
    #%      TRIG, pre, post and gap are counted in samples
    #%
    #% OUTPUT:
    #%  X   is a matrix of size [sz(1), sz(2)*sz(3)]
    #%  sz  [size(s,2), post-pre+1+gap, length(TRIG)]
    #%  sz(1) is the number of channels NS
    #%  sz(2) is the number of samples per trial
    #%  sz(3) is the number of trials i.e. length(TRIG)
    #%
    #% X3D = reshape(X,sz) returns a 3-dimensional matrix
    #% X2D = reshape(X(K,:),sz(2:3)) returns channel K in a 2-dimensional matrix
    post = round(post)
    nr, nc = sch.shape

    #include leading nan's
#    off = np.min([np.min(np.hstack([TRIG,np.inf]))+pre-1,0])
#    # include following nan's
#    off2 = np.max([np.max(np.hstack([TRIG,-np.inf]))+post-(sch.shape[0]),0])
#    if ((off!=0) or (off2!=0)):
#        sch = np.vstack([npm.repmat(np.NAN,-off,nc),
#                         sch,
#                         npm.repmat(np.NAN,off2,nc)])
#        TRIG = TRIG -off
    TRIG = TRIG.reshape((TRIG.shape[0], 1))
    N = post - pre + gap
    sz = [nc, N, (TRIG.shape[0] * TRIG.shape[1])]
    x = npm.repmat(np.NAN, int(sz[0]), int(sz[1] * sz[2]))
    for m in range(TRIG.shape[0]):
        indsch = ((TRIG[m]) + np.arange(pre + 1, post + 1, 1)).astype(int)
        indxx = ((m + 1) * N + np.arange(-N, -gap, 1)).astype(int)
        x[:, indxx] = sch[indsch, :].T
    return x, sz

###############################################################################

class ERDS:
    def __init__(self, erds, pre_erds, cl, cu):
        self.erds = erds
        self.pre_erds = pre_erds
        self.cl = cl
        self.cu = cu

###############################################################################

def calcErdsMapBP(indx, s1, h2, t, f_borders, f_bandwidths, f_steps,
                  clas, ref, submean, sig, lambd, alpha, refmethod, boxcox, fs, hTRIG, hClasslabel):
    # fs = h.SampleRate;
    # fs = 250.0

    n_butter = 5  # Filter order of Butterworth filter
    triallen = np.round((t[2] - t[0]) * fs) # Trial length (in samples)

    if t[1] != 0:
        # Time vector, reduced resolution
        t_vec_r = np.arange(t(0), t(2) + t(1), t(1))

    # Determine indices of reference interval
    #    ref[0] = np.round(ref[0]*fs)+1-np.round(t[0]*fs)
    #    ref[1] = np.round(ref[1]*fs)+1-np.round(t[0]*fs)

    f_plot = []  # Center frequencies in the plot
    f_low = []   # Lower frequency border for each center frequency
    f_up = []    # Upper frequency border for each center frequency

    for k in range(len(f_borders) - 1):
        f_plot.append(np.arange(f_borders[k], f_borders[k + 1], f_steps[k]))
        f_low.append(np.arange(
            f_borders[k] - f_bandwidths[k] / 2, f_borders[k + 1] - f_bandwidths[k] / 2, f_steps[k]))
        f_up.append(np.arange(f_borders[k] + f_bandwidths[k] / 2,
                              f_borders[k + 1] + f_bandwidths[k] / 2, f_steps[k]))
    f_plot.append(f_borders[-1])
    f_low.append(f_borders[-1] - f_bandwidths[-1] / 2)
    f_up.append(f_borders[-1] + f_bandwidths[-1] / 2)

    f_plot = np.hstack(f_plot)
    f_low = np.hstack(f_low)
    f_up = np.hstack(f_up)
    fn = fs / 2

    hh = {}
    r = {}
    chnList = []
    for chn in range(s1.shape[1]): # Loop over all channels   s1.shape[1]
        if submean: #Subtract evoked components
            # N�mero de canales x (muestras x trials)
            # hh.TRIG = h.TRIG;
            # hh.Classlabel = h.Classlabel;
            # Quitamos los artefactos y solo los canales
            hh['TRIG'] = hTRIG
            hh['Classlabel'] = hClasslabel
            auxIndx = np.logical_and(
                np.in1d(hh['Classlabel'], clas), np.squeeze(indx))
            s_t, _ = trigg(np.vstack(s1[:, chn]), hh['TRIG'][auxIndx], np.round(
                t[0] * fs), np.round(t[2] * fs), 0)

            # Reshape to samples x trials
            temp = np.reshape(
                s_t, (int(s_t.shape[1] / triallen), int(triallen))).T
            average = np.mean(temp, axis=1)
            temp = np.zeros((s1.shape[0], 1))
            hh['trig'] = hh['TRIG'][auxIndx]

            for kk in range(len(hh['trig'])):
                tempx = int(hh['trig'][kk] + np.round(t[0] * fs) + 1)
                tempy = int(hh['trig'][kk] + np.round(t[0]
                                                      * fs) + average.shape[0]) + 1
                temp[tempx:tempy, 0] = average
            s1[:, chn] = s1[:, chn] - temp.reshape((temp.shape[0],))

        # 1 reference per frequency band
        if 'classic' == refmethod:
            refp = np.zeros((1, f_plot.shape[0]))
        # References for each trial per frequency band
        else:
            if submean:
                refp = np.zeros((s_t.shape[1] / triallen, f_low.shape[0]))
            else:
                auxIndx = np.logical_and(
                    np.in1d(hh['Classlabel'], clas), indx)
                s_t, _ = trigg(np.vstack(s1[:, chn]), hh['TRIG'][auxIndx], np.round(
                    t[0] * fs), np.round(t[2] * fs), 0)
                refp = np.zeros((s_t.shape[1] / triallen, f_low.shape[0]))

        erds = np.zeros((int(triallen), f_plot.shape[0]))
        pre_erds = [None] * f_plot.shape[0]
        s_f = []
        for kk in range(f_plot.shape[0]):
            wn = [f_low[kk] / fn, f_up[kk] / fn]
            b_f, a_f = signal.butter(n_butter, wn, 'bandpass')
            smooth_length = np.ceil(2 * fs / f_low[kk])
            # , padlen=3*(max(len(b_f),len(a_f))-1)
            yAux = signal.filtfilt(b_f, a_f, s1[:, chn], padtype='odd')**2
            onesAux = np.ones((1, int(smooth_length))) / smooth_length
            s_f = signal.filtfilt(np.hstack(onesAux), [1], yAux)
            s_f = np.reshape(s_f, (s_f.shape[0], 1))
            s_t, _ = trigg(s_f, hh['TRIG'][auxIndx], np.round(
                t[0] * fs), np.round(t[2] * fs), 0)

            if refmethod == 'classic':
                pre_erds[kk] = np.reshape(
                    s_t, (int((s_t.shape[1]) / triallen), int(triallen))).T
                #% Average activity power over all trials
                activity = np.mean(pre_erds[kk], axis=1)
                refp[0, kk] = np.mean(
                    pre_erds[kk][int(ref[0] * fs + 1):int(ref[1] * fs + 1), :])
                erds[:, kk] = activity / refp[0, kk] - 1
                pre_erds[kk] = pre_erds[kk] / refp[0, kk] - 1

        # erds2 = [None] *s.shaep[1]
        erds2 = []
        refp2 = []
        if t[1] != 0:
            idx = np.round(t_vec_r * fs) + 1
            #r.ERDS{chn}.erds = erds(idx,:);
            erds2.append(erds[idx, :])
            for kk in range(f_plot.shape[0]):
                pre_erds[kk] = pre_erds[kk][idx, :]
        else:
            erds2.append(erds)
        if not(refmethod == 'absolute'):
            refp2.append(refp)
        # r['ERDSerds'] = erds2
        # r['ERDSrefq'] = refp2

        if boxcox == 'boxcox':
            #% FIXME: A normal distribution with known sigma is assumed, which is not
            #% quite correct. It should be replaced with the Student's t-distribution.
            #% This is implemented in boxcox2, but it has to be tested
            #% first. In the meantime, it is safe to use boxcox.
            cl_z = []
            cu_z = []
            z = special.erfinv(1 - alpha) * np.sqrt(2)
            for kk in range(f_plot.shape[0]):
                offset = np.floor(np.min(pre_erds[kk]))
                if offset < 0:
                    pre_erds[kk] = pre_erds[kk] + np.abs(offset)
                if lambd == 0:
                    erds_t = np.log(pre_erds[kk])
                    erds_t_m = np.mean(erds_t, axis=1)
                    erds_t_s = np.std(erds_t, axis=1, ddof=1)
                    ser = erds_t_s / np.sqrt(erds_t.shape[1])
                    cl_t = erds_t_m - z * ser
                    cu_t = erds_t_m + z * ser
                    cl_z.append(np.exp(cl_t) - abs(offset))
                    cu_z.append(np.exp(cu_t) - abs(offset))
                else:
                    erds_t = (pre_erds[kk]**lambd - 1) / lambd
                    erds_t_m = np.mean(erds_t, axis=1)
                    erds_t_s = np.std(erds_t, axis=1, ddof=1)
                    ser = erds_t_s / np.sqrt(erds_t.shape[1])
                    cl_t = erds_t_m - z * ser
                    cu_t = erds_t_m + z * ser
                    cl_z.append((cl_t * lambd + 1) **
                                (1 / lambd) - abs(offset))
                    cu_z.append((cu_t * lambd + 1) **
                                (1 / lambd) - abs(offset))

        # r['ERDScl'] = np.vstack(cl).T
        # r['ERDScu'] = np.vstack(cu).T
        print('End Channel ', chn + 1, ' of ', s1.shape[1])
        chnList.append(ERDS(np.squeeze(erds2), np.squeeze(
            pre_erds), np.vstack(cl_z).T, np.vstack(cu_z).T))
    r['ERDS'] = chnList
    # r['f_plot'] = f_plot
    # r['f_low'] = f_low
    # r['f_up'] = f_up
    # r['n_trials'] = np.sum(np.logical_and(np.in1d(hh['Classlabel'],clas),indx))
    # r['n_trials_class'] = hh['TRIG'][(np.in1d(hh['Classlabel'],clas))].shape[0]
    # r['refmethod'] = refmethod
    # retorna la lista con toda la información del ERD.
    return r

###############################################################################
# ----------------            Parametros de entrada          ---------------- #
###############################################################################
# Sujetos seleccionados de la base de datos.
SS = np.asarray([8])
t = [-2, 0, 5]           # tiempo de la señal.
f_borders = [6, 38]      # Frecuencias centrales de los filtros.
f_bandwidths = [4]       # tamaño del filtro en Hz.
f_steps = [2]            # Translape de entre los filtros.
ref = [0.5, 1.5]         # referencia de la señal.
submean = [0]            # Si realiza promediado de la señal.
submean = True           # .
sig = 'sig'              # Que tipo de significancia realiza.
lambd = 0                # Para realizar la transformación: si es 0  is log, if is different to 0 the Box-Cox transform <1x1>.
# Significancia utilizada default 1%, puede usarse tambien el 5% o el 10%.
alpha = 0.01
refmethod = 'classic'    # Método del calculo de la señal.

# Clases de la base de datos
# clase = 1                     # Clase 1 MI: mano derecha.
# clase = 2                     # Clase 2 MI: mano izquierda.
# clase = 3                     # Clase 3 MOV: mano derecha.
# clase = 4                     # Clase 4 MOV: mano izquierda.
Class_ = np.asarray([1, 2])     # Clases 1 y 2 para la imaginación motora.

############################# Importante mirar las pariciones
# Carga las particiones de los sujetos.
c = loadmat('Particiones/Cv_new.mat')
# Número de folds, para otras pruebas se utilizan 10 folds
Nfolds = 30     
# Número de canales.
Nchans = 64
# lugar donde se encuantra la base de datos.
Carpeta_database = 'F:/Giga_science/data/'
# lugar donde almacena los ERDs.
car_save = 'F:/Giga_science/'
# Carga el nombre de los archivos que comprenden la base de datos.
sub_dats = os.listdir(Carpeta_database)

# Main
if __name__ == '__main__':
    # Ciclo de los sujetos seleccionados.
    for sub_ in range(len(SS)):
        # Sujeto que se utiliza.
        sub = SS[sub_]   # ubicación del sujeto en el vector.
        sa = 2     # Ubicación de los sujetos.
        # sa = sub-1     # Ubicación de los sujetos.
        data_path = Carpeta_database + sub_dats[sa]
        data_s = loadmat(data_path)
        # Ensayos con artefactos por vol
        bad_trials_indx_vol = data_s['eeg']['bad_trial_indices'][0, 0][0, 0][0]
        # Ensayos con artefactos en canales EMG.
        bad_trials_indx_mi = data_s['eeg']['bad_trial_indices'][0, 0][0, 0][1]
        # Al sujeto que se esta realizando el análisis del ERDs.
        print('Load: Subject', sub)
        # Ciclo de la clases
        for class_s in range(len(Class_)):
            # Clase que se usa para el calculo de los ERDs..
            clas = Class_[class_s]
            ERDs = {}
            if class_s == 0:
                # Señal del sujeto imaginación mano izquierda.
                sdata = data_s['eeg']['imagery_left'][0][0][:Nchans, :].T
                bad_trials = bad_trials_indx_mi[0][0].reshape(
                    bad_trials_indx_mi[0][0].shape[0]).astype('int')
            else:
                # Señal del sujeto imaginación mano derecha.
                sdata = data_s['eeg']['imagery_right'][0][0][:Nchans, :].T
                bad_trials = bad_trials_indx_mi[0][1].reshape(
                    bad_trials_indx_mi[0][1].shape[0]).astype('int')
            # Ciclo de los folds.
            for folds in range(Nfolds):
                # Frecuencia de los datos.
                fs = data_s['eeg']['srate'][0][0][0][0]
                # Indices de los estímulos.
                indx_ = data_s['eeg']['imagery_event'][0][0]
                # posiciones de los eventos.
                indx2_ = np.array([i for i, x in enumerate(
                    np.ndarray.tolist(indx_[0])) if x == 1])
                # agregar los trials que se utlizan y revizar los cv #
                indx = np.ones(indx2_.shape[0])
                if class_s == 0:
                    # Carga las particiones de los sujetos en cada los folds.
                    cv = c['Cv_t'][19,folds][0:len(indx) - len(bad_trials)].T
                else:
                    # Carga las particiones de los sujetos en cada los folds.
                    cv = c['Cv_t'][19,folds][len(
                        indx) - len(bad_trials) + 1:-1]
                # Convierte los indices en bool.
                indx2 = cv.astype(bool)
                ind = bad_trials-1
                indx[ind.tolist()] = 0
                indx = indx.astype(bool)
                i_ = 0
                indx_t = np.ones(indx.shape[0]).astype('int').tolist()
                for i in range(indx.shape[0]):
                    if indx[i] == 1:
                        indx_t[i] = indx2[0][i_]
                        i_ += 1
                    else:
                        indx_t[i] = False
                            
                # Comparación de los trials que se van a utilizar.
                indx3 = np.array(indx_t) & indx
                hTRIG = indx2_
                hClasslabel = indx.astype(int) * clas
                if clas == 1:
                    cll = 1
                elif clas == 2:
                    cll = 2
                # Llama la función encargada del calculo del ERDs.
                r = calcErdsMapBP(indx3, sdata, 0, t, f_borders, f_bandwidths, f_steps,
                                  clas, ref, submean, sig, lambd, alpha, refmethod, 'boxcox', fs, hTRIG, hClasslabel)
                ERDs['ERD' + str(folds + 1)] = r
                print('ERD Sub: ', sub, ' Class ', cll, ' fold ', folds + 1)
        # Almacena los folds de un sujeto.
        savemat(car_save + 'GIGASCIENCE_' + sub_dats[sa][1:3] + '/ERD_Sub' + str(sub) + '_CL' + str(
            cll) + '_Nf' + str(folds + 1) + '_.mat', ERDs, format='5', do_compression=True)
