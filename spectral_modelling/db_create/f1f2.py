import spectral_modelling.utils.config as cfg
import numpy as np

# Define finf fsup from smoothed spectrum 

def f1f2(signal_smooth, noise_smooth, freqs):
    minstop = cfg.FMINSTOP
    maxstop = cfg.FMAXSTOP
    fmin = cfg.FREQLIMS[0]
    fmax = cfg.FREQLIMS[1]

    spt_ratio = signal_smooth/noise_smooth
    spt_ratio_dB = 10 * np.log10(spt_ratio)

    # define F min
    freqmin=1./cfg.TDURA
    if(cfg.FREQLIMS[0] < freqmin):
        print("WARNING! Noise window too short ( " + 
              str('{:.3f}'.format(cfg.FREQLIMS[0])) + 
              " seconds) for the selected hp frequency ( " +
              str('{:.3f}'.format(freqmin)) + " Hz)!\n")

    for i in range(len(freqs)):
        iflag=1
        freq = freqs[i]
        for j in range(i, i+minstop):
            if(j < len(freqs)):
                if (spt_ratio_dB[j] < cfg.SNRMIN):
                    iflag=0
        if (iflag == 1):
            fmin=freq
            break

    if (iflag == 0):
        print("NO freqmin found! signal REJECTED!\n")
        return None
   

    # Fmax
    for l in range(len(freqs))[::-1]:
        iflag=1
        freq = freqs[l]
        for m in range(l-maxstop+1,l+1)[::-1]:
            if(m > 0):
                if(spt_ratio_dB[m] < cfg.SNRMIN):
                    iflag=0
        if(iflag == 1):
            fmax=freq
            break

    if (iflag == 0):
        print("NO freqmax found! signal REJECTED!\n")
        return None

    fbmin = np.maximum(cfg.FREQLIMS[0],fmin)
    fbmax = np.minimum(cfg.FREQLIMS[1],fmax)

    if (fbmin > fbmax):
        print("fmin > fmax - signal REJECTED!\n")
        return None
   
    if (fbmin > cfg.F0_SNR or fbmax < cfg.F1_SNR):
        print('fmin ' + str('{:.3f}'.format(fbmin)) + ' > ' + 
              str('{:.3f}'.format(cfg.F0_SNR)) + ' or fmax ' + 
              str('{:.3f}'.format(fbmax)) + ' < ' + 
              str('{:.3f}'.format(cfg.F1_SNR)) + ' - signal REJECTED!\n')

        return None
   
    print("f1f2: band pass selected filter: " + str('{:.3f}'.format(fbmin)) + 
          " Hz - " + str('{:.3f}'.format(fbmax)) + " Hz\n")

    return(np.array([fbmin, fbmax]))
