import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spintegr
import scipy.interpolate as spinterpol
from scipy.optimize import curve_fit
from scipy.stats import norm
import mat73
from scipy.fft import fft, ifft
import sys

def dur(arr,lim = 0.5):
    lowfound = False
    low = 0
    high = len(arr)-1
    lim = np.max(arr)/2
    for i in range(len(arr)):
        if lowfound and arr[i] <= lim:
            high = i
            return low,high
        if not lowfound and arr[i] >=lim:
            lowfound = True
            low = i

    return low,high



# gauß function
def gauss(x,a,x0,sigma,c):
    return c+a*np.exp(-(x-x0)**2/(2*sigma**2))

def sech2(x,a,x0,b,c):
    return c+ a/((np.cosh((x-x0)/b))**2)

# envelopes an oscillating function
def envelope(x,y,di=2):
    x_new = []
    y_new = []

    for i in range(di,len(y)-di):

        undertest = y[i-2:i+2].copy()
        max1_i = np.nanargmax(undertest)
        undertest[max1_i] = -1
        max2_i = np.nanargmax(undertest)

        undertest = y[i - 2:i + 2].copy()
        avgval = 0.5*(undertest[max1_i]+undertest[max2_i])
        x_new.append(x[i])
        y_new.append(avgval)

    return x_new, y_new

def saveFile(content, path_and_filename,dim):

    if dim == 2:
        with open(path_and_filename,'wb') as f:
            for line in content:
                np.savetxt(f, line)
            f.close()
        #print("successfully saved matrix in "+path_and_filename)
    if dim == 1:
        with open(path_and_filename,'wb') as f:
            np.savetxt(f, content)
        #print("successfully saved array in "+path_and_filename)

def loadFile(path_and_filename):

        with open(path_and_filename,'r') as f:
            mat = np.loadtxt(path_and_filename,dtype=float)
            #print("successfully loaded matrix in "+path_and_filename)
            return mat

def getVmin(arr,start=0.1,end=0.9):
    startindex = int(len(arr)*start)
    endindex = int(len(arr)*end)

    index_min = np.nanargmin(arr[startindex:endindex])
    return arr[startindex+index_min]

# paths
pathtosourcematfiles = "C:/Users/ChrisAcer/Documents/University/MSc Physik/_LABs/QOMP LAB/femtosecond_pulses/"
res_folders = ["results1/","results2/","results3/"]
filenames = ["autocorrelation_0001.mat","autocorrelation_0002.mat","autocorrelation_0003.mat",
             "autocorrelation_0004.mat","autocorrelation_0005.mat"]
filenames_short = ["measure1/","measure2/","measure3/","measure4/","measure5/"]


# ---- load all data----------------
def dataextraction(folderindex,fileindex):

    folder = res_folders[folderindex]
    i = fileindex

    print(" -- FILE ",filenames[i]," ---")
    print(" -- FOLDER ", folder, " ---")
    contents = mat73.loadmat(pathtosourcematfiles + folder + filenames[i])

    # print stage positions
    pos = list(contents["pos"])
    # NEW

    # print & extract data
    print("print & extract")
    v_trace = list(contents["data"])                               # normal voltage traces
    len1 = len(v_trace)    # 500
    len2 = len(v_trace[0]) # 367xxx

    # trace data manipulation - remove offset etc
    print("data manipulation - remove offset")
    v_trace2 = np.zeros((len1, len2))     # without minimum
    for k in range(len1):
        min2 = getVmin(v_trace[k],0,0.1)
        for j in range(len2):
            v_trace2[k][j] = v_trace[k][j] - min2

    # test plot voltage trace
    trace_to_plot_y = v_trace2[75][:500]
    trace_to_plot_x = np.arange(len(trace_to_plot_y))*272 # 272ns 1 data point
    print("plot V trace")
    plt.plot(trace_to_plot_x, trace_to_plot_y)
    plt.grid()
    plt.title("Voltage Trace Signal")
    plt.ylabel("V(t)")
    plt.xlabel("t / ns")
    plt.legend(["raw","without offset"])
    plt.savefig(folder+filenames_short[i]+"trace_signal.jpg")
    plt.close()

    # calculate A(x) by integration (cumulative trapezius)
    print("calc acf -- trapezoidal integration")
    acf2 = [spintegr.cumtrapz(v_trace2[k])[-1] for k in range(len1)]
    acf2_max = np.max(acf2)
    for k in range(len(acf2)):
        acf2[k] /= acf2_max

    saveFile(acf2,folder+filenames_short[i]+"acf_y.txt",1)

    # x-->dT
    print("x --> dTau")
    midpos_index = np.nanargmax(acf2)
    midpos = pos[midpos_index]
    pos2 = np.subtract(pos, midpos)
    pos2 = np.multiply(pos2, 2)
    dTau = np.multiply(pos2, 3e3)
    #divide : ms : 3e8 , us: 3e5,  ns: 3e2 , ps: 3e-1 , fs: 3e-4

    # plot A(t)
    #plt.plot(pos, acf)
    plt.plot(dTau, acf2)
    plt.grid()
    plt.title("A(T) with cumulative trapezoidal method")
    plt.xlabel("T / fs")
    plt.ylabel("A(T) arb. unit")
    plt.legend(["raw","no offset"])
    plt.savefig(folder+filenames_short[i]+"acf_tau.jpg")
    plt.close()

    # save files

    saveFile(pos, folder+filenames_short[i]+"acf_x.txt",1)
    saveFile(dTau ,folder+filenames_short[i]+"acf_t.txt",1)

    saveFile(trace_to_plot_y, folder + filenames_short[i] + "trace_y.txt", 1)
    saveFile(trace_to_plot_x, folder + filenames_short[i] + "trace_x.txt", 1)

if __name__ == '__main__':

    #from matplotlib import rc
    #rc('text', usetex=True)
    for j in range(5):
        for i in range(3):

            # --- which file ---
            fold_index = i  # 0,1,2   robert chris shorouk
            file_index = j # 0,1,2,3,4  measure 1.. 5

            # extract or not
            extract_data = False

            # --- do what ---
            field_routine = False
            interferometric_routine = False
            intensity_routine = False

            if file_index == 0:
                field_routine = True
            elif file_index == 1 or file_index == 2:
                interferometric_routine = True
            elif file_index == 3 or file_index == 4:
                intensity_routine = True
            else:
                continue

            # start ....................

            if extract_data:
                dataextraction(fold_index,file_index)

            # fileindex 0
            if field_routine:

                # load E-Field ACF
                print("loading ACF")

                acf_x = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_x.txt")
                acf_t = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_t.txt")
                acf_y = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_y.txt")

                vt_x = loadFile(res_folders[fold_index]+filenames_short[file_index]+ "trace_x.txt")
                vt_y = loadFile(res_folders[fold_index]+filenames_short[file_index]+ "trace_y.txt")

                #calc delta X and delta T
                dx = np.average(np.diff(acf_x)) / 1e3  # dx in nm
                dT = np.average(np.diff(acf_t))/1e15     # dT in fs
                len_acf = len(acf_x) # 500

                # plot acf(T)
                plt.plot(acf_t, acf_y)
                plt.title("ACF(T); E-Field ACF")
                plt.ylabel("A(T)")
                plt.grid()
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "acf.jpg")
                #plt.show()
                plt.close()

                # FFT of E-Field ACF --> spectrum, I_f
                print("F(acf(T)) -- calculating spectrum")
                I_f = fft(acf_y)
                I_f_abs = []

                for k in range(len(I_f)):   #complex --> real
                    I_f_abs.append((np.abs(I_f[k].real ** 2) + np.abs(I_f[k].imag ** 2)))
                fft_acf_max = np.max(I_f_abs) # get biggest value
                I_f_abs[0] = 0.0

                #saveFile(I_f_abs, res_folders[fold_index] + filenames_short[file_index] + "I_f.txt", 1)

                # I(w) --> I(lambda)
                # k*n/N = f_k * t_n
                X = len_acf*dx
                Tau = len_acf*dT
                print("tau: ",Tau)
                deltaF = 1/Tau
                deltaK = 1/X
                print("dF: ",deltaF)
                fk = [deltaF * k for k in range(len_acf)]
                kk = [deltaK * k for k in range(len_acf)]

                I_f_abs = np.divide(I_f_abs,np.max(I_f_abs))

                # get I_lambda with jacobi determinant
                # I_Lambda = I_f * c / lambda^2
                lambdaK = [3e8 / fk[k] for k in range(1, len(fk))]
                Ilambda = [I_f_abs[k + 1] * 2 * np.pi / lambdaK[k] ** 2 for k in range(len(lambdaK))]
                Ilambda = np.divide(Ilambda,np.max(Ilambda))
                saveFile(np.multiply(lambdaK, 1e9)[:int(len(Ilambda) / 2)],
                         res_folders[fold_index] + filenames_short[file_index] + "I_la_x.txt", 1)
                saveFile(Ilambda[:int(len(Ilambda) / 2)],
                         res_folders[fold_index] + filenames_short[file_index] + "I_la_y.txt", 1)
                plt.plot(np.multiply(lambdaK, 1e9)[:int(len(Ilambda)/2)], Ilambda[:int(len(Ilambda)/2)])
                plt.grid()
                plt.title("Spectrum of Field- ACF (Wavelength)")
                plt.xlabel("lambda / nm")
                plt.ylabel("I / arb. u.")
                plt.xlim(100,1200)
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "I_lambda.jpg")
                #plt.show()
                plt.close()

                # FFT of TRACE
                # FFT of voltage traces
                print("F(V(t))")
                I_ft = fft(vt_y)
                I_ft_abs = [np.abs(I_ft[k]) for k in range(len(I_ft))]

                len_trace = len(vt_y)
                Tau = len_trace * 270 * 1e-9
                print("tau: ", Tau)
                deltaF = 1 / Tau
                print("dF: ", deltaF)
                fk = [deltaF * k for k in range(len_trace)]

                I_ft_abs = np.divide(I_ft_abs, np.max(I_ft_abs))
                I_ft_abs[0] = 0
                plt.plot(np.divide(fk,1e3)[:int(len(I_f_abs)/2)], I_ft_abs[:int(len(I_f_abs)/2)])
                plt.grid()
                plt.title("Spectrum of Voltage-Pulses (Frequency)\nDC-part eliminated")
                plt.xlabel("f / kHz")
                plt.ylabel("I / arb. u.")
                plt.xscale("log")
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "I_f_voltagetrace.jpg")
                #plt.show()
                plt.close()
            # file index 1,2
            if interferometric_routine:

                # load interferometric ACF
                print("loading ACF")
                acf_t = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_t.txt")
                acf_y = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_y.txt")
                plt.plot(acf_t, acf_y)
                plt.title("ACF (T) IFM-ACF")
                plt.grid()
                plt.ylabel("A(T) /")
                #plt.show()
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "acf.jpg")
                plt.close()
                AVG_ACF = np.average(acf_y)  # needed for fit

                # FFT of E-Field ACF --> spectrum, I_f
                print("F(acf(T)) -- calculating spectrum")
                I_f = fft(acf_y)
                I_f_abs = []

                for k in range(len(I_f)):  # complex --> real
                    I_f_abs.append((np.abs(I_f[k].real ** 2) + np.abs(I_f[k].imag ** 2)))
                fft_acf_max = np.max(I_f_abs)  # get biggest value
                I_f_abs[0] = 0.0

                # I(w) --> I(lambda)
                # k*n/N = f_k * t_n
                len_acf = len(acf_y)
                dT = np.average(np.diff(acf_t))*1e-15

                Tau = len_acf * dT
                print("tau: ", Tau)
                deltaF = 1 / Tau
                print("dF: ", deltaF)
                fk = [deltaF * k for k in range(len_acf)]

                I_f_abs = np.divide(I_f_abs, np.max(I_f_abs))

                # get I_lambda with jacobi determinant
                # I_Lambda = I_f * c / lambda^2
                lambdaK = [3e8 / fk[k] for k in range(1, len(fk))]
                Ilambda = [I_f_abs[k + 1] * 2 * np.pi / lambdaK[k] ** 2 for k in range(len(lambdaK))]
                Ilambda = np.divide(Ilambda, np.max(Ilambda))

                saveFile(np.multiply(lambdaK, 1e9)[:int(len(Ilambda) / 2)],
                         res_folders[fold_index] + filenames_short[file_index] + "I_la_x.txt", 1)
                saveFile(Ilambda[:int(len(Ilambda) / 2)],
                         res_folders[fold_index] + filenames_short[file_index] + "I_la_y.txt", 1)

                plt.plot(np.multiply(lambdaK, 1e9)[:int(len(Ilambda) / 2)], Ilambda[:int(len(Ilambda) / 2)])
                plt.grid()
                plt.title("Spectrum of IFM-ACF (Wavelength)\n")
                plt.xlabel("lambda / nm")
                plt.ylabel("I / arb. u.")
                plt.xlim(100, 1200)
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "I_lambda.jpg")
                #plt.show()
                plt.close()
            # file index 3,4
            if intensity_routine:

                # load ACF !!
                acf_x = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_x.txt")
                acf_y = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_y.txt")
                acf_t = loadFile(res_folders[fold_index] + filenames_short[file_index] + "acf_t.txt")

                plt.plot(acf_t, acf_y)
                plt.title("ACF (T) Intensity-ACF")
                plt.grid()
                plt.xlabel("T / fs")
                plt.ylabel("A(T) / arb. u.")
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "acf.jpg")
                #plt.show()
                plt.close()

                # envelope of ACF ( for fit ) !!
                env_x, env_y = envelope(acf_t,acf_y,3)
                plt.plot(env_x,env_y)
                plt.plot(acf_t,acf_y)
                plt.title("Envelope of Intensity-ACF")
                plt.grid()
                plt.xlabel("T / fs")
                plt.ylabel("A(T) / arb. u.")
                #plt.show()
                plt.close()

                # fit Gauß params !!
                popt, pcov = curve_fit(gauss, env_x, env_y, p0=[1,0,10,0.5])
                pcov_arr = np.sqrt(np.diag(pcov))
                f_gauß = gauss(env_x, *popt)

                # fit Sech2 params !!
                popt2, pcov2 = curve_fit(sech2, env_x, env_y, p0=[1, 0, 10, 1])
                pcov2_arr = np.sqrt(np.diag(pcov2))
                f_sech = sech2(env_x, *popt2)

                plt.plot(env_x, env_y, label='smoothed data')
                plt.plot(env_x, f_gauß, label='Gauss fit')
                plt.plot(env_x, f_sech, label='Sech2 fit')

                plt.legend()
                plt.title(f'Fit curves for Intensity Autocorreleation')
                plt.grid()
                plt.xlabel("T / fs")
                plt.ylabel("A(T)")
                plt.savefig(res_folders[fold_index] + filenames_short[file_index] + "fits.jpg")
                #plt.show()
                plt.close()

                l, h = dur(f_gauß)
                pulsduration_gauss = env_x[h]-env_x[l]

                l, h = dur(f_sech)
                pulsduration_sech = env_x[h] - env_x[l]

                saveFile(popt, res_folders[fold_index] + filenames_short[file_index] + "gaussfit.txt", 1)
                saveFile(popt2, res_folders[fold_index] + filenames_short[file_index] + "sech2fit.txt",1)
                saveFile(pcov_arr,res_folders[fold_index] + filenames_short[file_index] + "gaussfit_errors.txt", 1)
                saveFile(pcov2_arr, res_folders[fold_index] + filenames_short[file_index] + "sech2fit_errors.txt", 1)
                saveFile(np.array([pulsduration_gauss,pulsduration_sech]),res_folders[fold_index] + filenames_short[file_index] + "pulse_dur.txt", 1)
                print("PULS DURATIONS: ",pulsduration_gauss,pulsduration_sech)

    # e - field acf - spectrum of lambda

    i_la_x11 = loadFile("results1/measure1/I_la_x.txt")
    i_la_x12 = loadFile("results2/measure1/I_la_x.txt")
    i_la_x13 = loadFile("results3/measure1/I_la_x.txt")

    i_la_y11 = loadFile("results1/measure1/I_la_y.txt")
    i_la_y12 = loadFile("results2/measure1/I_la_y.txt")
    i_la_y13 = loadFile("results3/measure1/I_la_y.txt")

    plt.plot(i_la_x11,i_la_y11)
    plt.plot(i_la_x12, i_la_y12)
    plt.plot(i_la_x13, i_la_y13)
    plt.title("Spectrum of E-Field Autocorrelation")
    plt.legend(["Robert","Christoph","Shorouk"])
    plt.xlim(100,1200)
    plt.grid()
    plt.show()

    # ifm acf - spectrum of lambda without glass

    i_la_x21 = loadFile("results1/measure2/I_la_x.txt")
    i_la_x22 = loadFile("results2/measure2/I_la_x.txt")
    i_la_x23 = loadFile("results3/measure2/I_la_x.txt")

    i_la_y21 = loadFile("results1/measure2/I_la_y.txt")
    i_la_y22 = loadFile("results2/measure2/I_la_y.txt")
    i_la_y23 = loadFile("results3/measure2/I_la_y.txt")

    plt.plot(i_la_x21, i_la_y21)
    plt.plot(i_la_x22, i_la_y22)
    plt.plot(i_la_x23, i_la_y23)
    plt.title("Spectrum of Interferometric Autocorrelation\nGlass Plate removed")
    plt.legend(["Robert", "Christoph", "Shorouk"])
    plt.xlim(100, 1200)
    plt.grid()
    plt.show()

    # ifm acf - spectrum of lambda with glass

    i_la_x31 = loadFile("results1/measure3/I_la_x.txt")
    i_la_x32 = loadFile("results2/measure3/I_la_x.txt")
    i_la_x33 = loadFile("results3/measure3/I_la_x.txt")

    i_la_y31 = loadFile("results1/measure3/I_la_y.txt")
    i_la_y32 = loadFile("results2/measure3/I_la_y.txt")
    i_la_y33 = loadFile("results3/measure3/I_la_y.txt")

    plt.plot(i_la_x31, i_la_y31)
    plt.plot(i_la_x32, i_la_y32)
    plt.plot(i_la_x33, i_la_y33)
    plt.title("Spectrum of Interferometric Autocorrelation\nGlass Plate inserted")
    plt.legend(["Robert", "Christoph", "Shorouk"])
    plt.xlim(100, 1200)
    plt.grid()
    plt.show()
