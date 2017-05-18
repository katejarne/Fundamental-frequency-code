######################################################################################################
#                                                                                                    # 
# 1)Reads wav 1 or 2 ch and plot signal.                                                             #
# 2)Estimate sonogram and fundamental frequency                                                      #
# 3)Make consecutive plots of "step_time" size in file with current wav name for each file in r_dir  #
#                                                                                                    #
#                                         C. Jarne 30-03-2017 V1.0                                   #                           
######################################################################################################

# libraries 
import numpy as np
import scipy
import os
import scipy.stats as stats
import time as tiempos

from scipy.io import wavfile
import wave, struct
import matplotlib.pyplot as pp

from pylab import *
import scipy.signal.signaltools as sigtool
from scipy.io import wavfile
import wave, struct
import scipy.signal as signal
from scipy.fftpack import fft

# Here directory (put the name and path). Directory only with .wav files

#r_dir='/home/.../sound_files'
r_dir ='/home/kathy/Escritorio/Enlace hacia 2017_signal_analysis/sound_files/canary'
#r_dir ='/home/kathy/Escritorio/Enlace hacia 2017_signal_analysis/sound_files/canary_2'
#r_dir='/home/kathy/Escritorio/Enlace hacia 2017_signal_analysis/sound_files/whale'
#r_dir='/home/kathy/Escritorio/Enlace hacia 2017_signal_analysis/sound_files/word_english'
#r_dir='/home/kathy/Escritorio/Enlace hacia 2017_signal_analysis/sound_files/piano'
#r_dir='/home/kathy/Escritorio/Enlace hacia 2017_signal_analysis/sound_files/zebra_finch'

# Parameters

Fmax         = 10000 #maximum frequency for the sonogram [Hz]
step_time    = 1.25  #len for the time serie segment  [Seconds]->>>>>>>>> Change it to zoom in the signal time!!
threshold    = 700

amp_min =-25000
amp_max =25000
start_time = tiempos.time()
for root, sub, files in os.walk(r_dir):
    files = sorted(files)
    for f in files:       
        w     = scipy.io.wavfile.read(os.path.join(root, f))
        print (r_dir)
        base=os.path.basename(f)
        print (base)
        dir = os.path.dirname(base)
        if not os.path.exists(dir):
           os.mkdir(base)        
        print('-------------------------')
        
        a=w[1]
        print('sound vector: ')#, w

        i=w[1].size

        print ('vector size in Frames: ',i)

        x     = w[1]
        x_size= x.size
        tt    = w[1]

        #Comment for stero or not 
        v1    = np.arange(float (i)/float(2)) #stereo
        #v1    = np.arange(float (i))#/float(2)) #not stereo
        c     = np.c_[v1,x]

        print ('vector c:\n' , c)
        print ('vector c1:\n',c[0])
       
        cc=c.T #transpose

        x = cc[0]
        x1= cc[1]
        x2= cc[1]#2

        aa       = scipy.signal.medfilt(scipy.signal.detrend(x2, axis=-1, type='linear'))

        print ('First cc comp:\n ', cc[0])
        print ('Second cc comp:\n', cc[1])
        print ('Third cc comp: \n', cc[1])# cc[2] if stereo


        stop      = (float(i)/float(2))#if stereo
        #stop      = i #if not stereo
        step      = int(step_time*w[0])
        intervalos= np.arange(0, int(stop),step)

        print('intervalos', intervalos)
        print('-------------------')
        print('The step: ',step)
        print('-------------------')

        time1=x*float(1)/w[0]

        ##chop time serie##
        for delta_t in intervalos:
                  
            aa_part                   = aa[delta_t:delta_t+step]
            x1_part                   = x2[delta_t:delta_t+step]#or x1
            x2_part                   = x2[delta_t:delta_t+step]


           ###################################################

            #Figure definition
            
            espectro             = pp.figure(figsize=(11.5,10))
            pp.title('Sound Signal')
            pp.subplot(4,1,1)

            #grid(True)

            #Signal
            pp.plot(x*float(1)/w[0],x1, color='c',label='Time Signal')
            #pp.plot(aa_part*float(1)/w[0],x1_part, color='c',label='Time Signal')
            pp.ylabel('Amplitude [Arbitrary units]')        
            pp.xlim([delta_t*float(1)/w[0],(delta_t+step)*float(1)/w[0]+0.001])
            pp.xticks(np.arange(delta_t*float(1)/w[0],(delta_t+step)*float(1)/w[0]+0.001,0.1),fontsize = 12)
            #pp.yticks(np.arange(-15000,15000+5000,5000),fontsize = 12)
            pp.ylim(amp_min,amp_max)
            #pp.tick_params( axis='x', labelbottom='off')
            #pp.tick_params( axis='y', labelleft='off')
            
            pp.legend(fontsize= 'small',loc=1)


            #Sonogram

            pp.subplot(4,1,2)

            #grid(True)
            nfft_=int(w[0]*0.010)
         
            Pxx, freqs, bins, im = pp.specgram(x1_part, NFFT=int(w[0]*0.008), Fs=w[0], noverlap =int(w[0]*0.005))
            pp.xlim(0, step_time+0.001)
            pp.yticks(np.arange(0,Fmax,1000),fontsize = 12)            
            pp.ylim(0, Fmax)
            
            pp.ylabel('Frequency [Hz]')
            pp.tick_params( axis='x', labelbottom='off')
                                        
            pp.subplot(4,1,3)

            #grid(True)
            nfft_=int(w[0]*0.010)
         
            
            Pxx, freqs, bins, im = pp.specgram(x1_part, NFFT=int(w[0]*0.008), Fs=w[0], noverlap =int(w[0]*0.005))


            #espectro.colorbar(im).set_label('Intensity [dB]')
            #grid(True)
            #pp.xticks(np.arange(0,8+0.001,0.1),fontsize = 12)
            pp.xlim(0, step_time+0.001)   
            #pp.xlim([delta_t*float(1)/w[0],(delta_t+step)*float(1)/w[0]+0.001])
            #pp.xticks(np.arange(delta_t*float(1)/w[0],(delta_t+step)*float(1)/w[0]+0.001,0.1),fontsize = 12)
            pp.yticks(np.arange(0,Fmax,1000),fontsize = 12)            
            pp.ylim(0, Fmax)
            
            pp.ylabel('Frequency [Hz]')             
     
            pp.tick_params( axis='x', labelbottom='off')
            #pp.tick_params( axis='y', labelleft='off')
            
            
            #Filtering freq in the sonogram
 
	    Pxx   = Pxx[(freqs >= 100) & (freqs <= 8000)]
            freqs = freqs[(freqs >= 100) & (freqs <= 8000)]
            #bins  = bins[(freqs >= 3000) & (freqs <= 8000)]
            #debugging printing

            print 'cada componente del specgram:\n'
            print 'Pxx shape:\n',Pxx.shape
            print 'size Pxx:\n',Pxx.size
            print 'maximo de Pxx:\n',np.argmax(Pxx[1])
            print 'Pxx[1] size:\n',Pxx[1].size
	
            Pxx_transpuesto=Pxx.T
            
            print'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
            print 'frequency values:\n', freqs
            print 'frequencies in each temporal bin ', freqs.size
            print "frequency resolution",freqs[1]-freqs[0]
            print 'temporal bins:\n', bins
            print 'size bins: ', bins.size
            if len(bins)>2:
                print 'time bin size', bins[1]-bins[0],' [Sec]'
            intencidad_maxima=[]
            frec_max=[]

            #to find fundamental frequency
            for j, element in enumerate(Pxx_transpuesto):
		intencidad_maxima_indice=np.argmax(element)
		intencidad_maxima.append(element[intencidad_maxima_indice])
		frec_max.append(freqs[intencidad_maxima_indice])

            frequ_size   = len(frec_max)

                  
            #transform into array
            intencidad_maxima_ = np.asarray(intencidad_maxima)
            bins_              = np.asarray(bins)
            frec_max_          = np.asarray(frec_max)
            

            #sonogram intensity threshold             
            bins_              = bins_[intencidad_maxima_>threshold]
            frec_max_          = frec_max_[intencidad_maxima_>threshold]
            intencidad_maxima_ = intencidad_maxima_[intencidad_maxima_>threshold]

            pp.scatter(bins_,frec_max_, color='k',label='Fundamental Frequency')

            pp.subplot(4,1,4)
            #grid(True)
            pp.yticks(np.arange(0,Fmax,1000),fontsize = 12)
            pp.scatter(bins_,frec_max_, color='b',label='Fundamental Frequency')
            pp.tick_params( axis='x', labelbottom='off')
            pp.ylabel('Frequency [Hz]') 
            pp.xlabel('Time [Sec]')
            pp.xlim(0, step_time+0.001)
            pp.ylim(0, Fmax)
            pp.legend(fontsize= 'small',loc=1)

            X=bins_
            Y=frec_max_
            total_bins = 50
            bins = np.linspace(X.min(),X.max(), total_bins)
            delta = bins[1]-bins[0]
            idx  = np.digitize(X,bins)
            running_median = [np.median(Y[idx==k]) for k in range(total_bins)]
            #pp.plot(bins-delta/2,running_median,'r--',lw=4,alpha=.8)

            #pp.tick_params( axis='x', labelbottom='off')        
            figname = "%s.jpg" %(str(base)+'_signal_zoom_'+str(delta_t*float(1)/w[0]))    
            pp.savefig(os.path.join(base,figname),dpi=200)
            pp.close('all')

            ###############################################################
            #save in plot file txt with data if necesary

            f_out     = open('%s.txt' %(str(base)+'_'+str(delta_t*float(1)/float(w[0]))), 'w')                
            xxx       = np.c_[bins_,frec_max_]
            np.savetxt(f_out,xxx,fmt='%f %f',delimiter='\t',header="time   #freq [Hz]") 
                                              
        print ('.---All rigth!!!----.')
        print("--- %s seconds ---" % (tiempos.time() - start_time))







