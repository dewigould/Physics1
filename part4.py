see_results = False
validation_routine= False


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as s

if see_results == True:
    
    #Define sampling rate
    #Range of Convolution
    T= 20.00
    #Set-up ~symmetric interval about zero for nice visualisation of results
    t_start = -T/2.0
    #Time Increment of Sampling
    delta_t = 0.01
    #Number of samples taken
    N = int(T/delta_t)
    #Create array of all sampling times - from t_start to t_start + T
    t_range_input = np.array([t_start + (n*delta_t) for n in range(N)])


#------------------------------------------------
    #Define "Top Hat" Signal Function
    #Normalise the top_hat to avoid floating-point errors later on
    NormFactor = 10.0 #keep normalisation factor separate to use later for plotting
    def h(x):
        if 2.0 <= x <= 4.0:
            return 5.0/NormFactor + 0.0j
        else:
            return 0.0 + 0.0j
    
    #Define Gaussian Response Function       
    def g(x):
        pre_factor = (1.0/np.sqrt(2*(np.pi)))
        result = pre_factor * ((np.e)**(-1.0*((x**2)/2.0)))
        return result + 0.0j



#Function: Convolve two functions together
#InPut:
#       f - signal function
#       g - response function
#       t_range - sampling range of times to perform convolution over
#OutPut:
#       f_sam_padded - updated array of sampled f values including padding
#       g_sam_padded - updated array of sampled g values including padding
#       DFT_f - array of values of DFT of function f
#       DFT_g - array of values of DFT of function g
#       w_range - range of sampling frequencies
#       t_range - updated range of sampling times including padding times
#       conv - normalised result of convolution
#       conv_raw - un-normalised result of convolution
#       test_conv - result of scipy convolution for comparison
def convolve(f,g,N,t_range):
    
    delta_t = max(t_range)/N
    
    #Take N samples of both functions
    f_sam = [f(t) for t in t_range]
    g_sam = [g(t) for t in t_range]

    #Zero Pad the signals to prevent aliasing effects 
    #Avoids issues with cut-offs at end of convolution
    #By Convention, add half as many zeros as there are sampling points to both ends of signal    
    N_padding = int(0.5*N)
    t_before = []
    t_after = []

    f_sam_padded = f_sam
    g_sam_padded = g_sam
  
    #Implement the padding by adding zeros to h and g samples, also extend time array range to account for new points
    for i in range(N_padding):
        
        t_before.append(t_range[0] - (float(i+1))*delta_t)
        t_after.append(t_range[-1] + (float(i+1))*delta_t)
        f_sam_padded.append(0.0+0.0j)
        g_sam_padded.append(0.0+0.0j)
       
    #add zeros to beginning of function
    f_sam_padded = [0.0+0.0j]*N_padding + f_sam_padded
    g_sam_padded = [0.0+0.0j]*N_padding + g_sam_padded
    
    #extend time range before and after previous range
    t_before = sorted(t_before,key=lambda a:a)
    t_range = np.append(t_before,t_range)
    t_range = np.append(t_range,t_after)
    
    #Take discrete fourier transforms of each function
    DFT_f = np.fft.fft(f_sam_padded)
    DFT_g = np.fft.fft(g_sam_padded)
    
    #Generate array of sampling frequencies for visualisation purposes
    w_range = np.fft.fftfreq(t_range.shape[-1])
    
    #Multiply fft's together to get DFT of convolution
    DFT_fg = []
    for i in range(N+(2*N_padding)):
        DFT_fg.append(DFT_f[i]*DFT_g[i])
        
    #Take inverse of multiplication, account for normalisation factor
    conv_raw = np.fft.ifft(DFT_fg)
    conv = conv_raw*delta_t*NormFactor
    conv = np.fft.fftshift(conv)
    conv_raw = np.fft.fftshift(conv_raw)
    
    #Perform convolution using standard library for comparison
    test_conv = s.convolve(f_sam_padded,g_sam_padded,mode='same')
    
    return f_sam_padded,g_sam_padded,DFT_f,DFT_g,w_range,t_range,conv,conv_raw,test_conv


if see_results == True:        
    #Gather results for plotting
    h_sam_padded, g_sam_padded, DFT_h, DFT_g, w_range,t_range,conv, conv_raw,test_conv = convolve(h,g,N,t_range_input)

    #Plot
    f, ax = plt.subplots(2, 2)
    ax[0, 0].plot(t_range,[i*NormFactor for i in h_sam_padded],color="red",label="Signal Function")
    ax[0,0].plot(t_range,g_sam_padded,color="blue",label="Response Function")
    ax[0,0].plot(t_range,conv.real,color = "green",label="Convolution")
    ax[0,0].legend()
    ax[0, 0].set_title("Functions to be Convolved")
    ax[0, 1].plot(w_range,DFT_h.real,label="real")
    ax[0,1].legend()
    ax[0, 1].set_title("DFT of Signal Function")
    ax[1, 0].plot(w_range, DFT_g.real,label = "real")
    ax[1,0].legend()
    ax[1, 0].set_title("DFT of Response Function")
    ax[1, 1].plot(t_range,conv_raw.real,label = "Result from np.fft")
    ax[1,1].plot(t_range,test_conv.real,label = "Comparison result from scipy.signal.convolve")
    ax[1, 1].set_title("Result of Convolution")
    ax[1,1].legend(fontsize = 'x-small')
    


#Validate by convolving a simple Gaussian signal with itself
if validation_routine == True:
    
    #Define sampling rate
    #Range of Convolution
    T= 20.00
    #Set-up ~symmetric interval about zero for nice visualisation of results
    t_start = -T/2.0
    #Time Increment of Sampling
    delta_t = 0.01
    #Number of samples taken
    N = int(T/delta_t)
    #Create array of all sampling times - from t_start to t_start + T
    t_range_input = np.array([t_start + (n*delta_t) for n in range(N)])


    def g(x):
        pre_factor = (1.0/np.sqrt(2*(np.pi)))
        result = pre_factor * ((np.e)**(-1.0*((x**2)/2.0)))
        return result + 0.0j

    #Gather results for plotting
    h_sam_padded, g_sam_padded, DFT_h, DFT_g, w_range,t_range,conv, conv_raw,test_conv = convolve(g,g,N,t_range_input)

    #Plot
    f, ax = plt.subplots(2, 2)
    ax[0, 0].plot(t_range,h_sam_padded,color="red",label="Signal Function")
    ax[0,0].plot(t_range,g_sam_padded,color="blue",label="Response Function")
    ax[0,0].plot(t_range,conv.real,color = "green",label="Convolution")
    ax[0,0].legend()
    ax[0, 0].set_title("Functions to be Convolved")



