# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:07:51 2018
@author: Jakub
"""
import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib.pyplot import cm

def getLogHist(a, z, histBins):
    #t0 = time.time()
    a = np.asarray(a)
    hist, _ = np.histogram(a, histBins)
    hist = np.asarray(hist,dtype=np.float)
    for i in range (hist.size):
       hist[i] /= (z**i * (z - 1))
    hist = np.divide(hist,hist.size)
    #print("logHist time: ", time.time() - t0)
    return hist

def getLogHistBins(a, z, maxNodeDegree):
    a = np.asarray(a)
    bins = np.log(maxNodeDegree)/np.log(z)
    bins = int(np.ceil(bins))
    histBins = [z**i for i in range (bins+1)]
    return histBins    

def fitLine(a, a_bins): 
    #t0 = time.time()
    hist_log = []
    hist_bins_log = []
    for i in range (a.size):
        if a[i] > 0:
            hist_log.append(np.log10(a[i]))
            hist_bins_log.append(np.log10(a_bins[i]))
    fit = np.polyfit(hist_bins_log, hist_log, 1)
    #print("fit_line time: ", time.time() - t0)
    return fit  
    
def readArray(input_file):
    t0 = time.time()
    f = open(input_file, 'r')
    lines = f.readlines()
    lines = [x.strip() for x in lines] 
    lines = [float(x) for x in lines]  
    f.close()  
    #print(np.asarray(lines))
    print('reading', str(input_file), "time= " , time.time() - t0)
    return np.asarray(lines)

def readParams(input_file):
    f = open(input_file, 'r')
    lines = f.readlines()    
    lines = [x.replace('\n','')for x in lines]
    lines = [x.replace('\t','') for x in lines]
    lines = [float(x) for x in lines]  
    f.close()  
    NBlocks = lines[0]
    BlockSize = lines[1]
    BlockSize = int(BlockSize)
    NBlocks = int(NBlocks)
    BlockXMin = []
    BlockExponents = []
    for i in range (NBlocks):
        BlockXMin.append(lines[2+i])
        BlockExponents.append(lines[2+NBlocks+i])
    return NBlocks, BlockSize, BlockXMin, BlockExponents

def plotAndFit(nodeDegrees, NBlocks, BlockSize):
    z = 2
    color=iter(cm.rainbow(np.linspace(0,1,NBlocks))) 
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(7,6))
    plt.locator_params(axis='y', nticks=6)
    plt.locator_params(axis='x', nticks=10)
    plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.7)
    for i in range(NBlocks):
        maxNodeDegree = np.max(nodeDegrees[i*BlockSize : (i+1)*BlockSize - 1])
        histBins = getLogHistBins(nodeDegrees[i*BlockSize : (i+1)*BlockSize - 1],z,maxNodeDegree) 
        hist = getLogHist(nodeDegrees[i*BlockSize : (i+1)*BlockSize - 1],z,histBins)
        del histBins[-1]
        fit = fitLine(hist, histBins)
        print("Fit", i, " y = ax + b, [a, b]", fit)
        c = next(color)
        ax1.plot(histBins, hist, color=c, marker="o", linestyle='-', label=i)
        ax1.legend(loc='upper right')
        plt.grid(True)
        plt.title('node degrees histogram - linear')
        ax2.plot(histBins, hist, color=c, marker="o", linestyle='-', label=i)
        plt.yscale('log')
        plt.xscale('log')
        ax2.legend(loc='upper right')
        plt.grid(True)
        plt.title('node degrees histogram - log-log')
    plt.show()

def plotAndFitExpectedSimulation(nodeDegrees, nodeDegreesSym, NBlocks, BlockSize):
    z = 2
    f= open("outputBlockNodeDegreesFit.dat","w+")
    for i in range(NBlocks):
        fig, _ = plt.subplots(nrows=1, ncols=1, figsize=(7,6))
        plt.locator_params(axis='y', nticks=6)
        plt.locator_params(axis='x', nticks=10)
        plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.7)
        maxNodeDegree = np.max(nodeDegrees[i*BlockSize : (i+1)*BlockSize - 1])
        histBins = getLogHistBins(nodeDegrees[i*BlockSize : (i+1)*BlockSize - 1],z,maxNodeDegree) 
        hist = getLogHist(nodeDegrees[i*BlockSize : (i+1)*BlockSize - 1],z,histBins)
        del histBins[-1]
        fit = fitLine(hist, histBins)    
        #print("Fit expected, block ", i+1, " y = ax + b, [a, b]", fit)
        f.write("Fit expected, block " + str(i+1) + " y = ax + b, [a, b] " +  str(fit) +'\n')
        maxNodeDegreeSym = np.max(nodeDegreesSym[i*BlockSize : (i+1)*BlockSize - 1])
        histBinsSym = getLogHistBins(nodeDegreesSym[i*BlockSize : (i+1)*BlockSize - 1],z,maxNodeDegreeSym) 
        histSym = getLogHist(nodeDegreesSym[i*BlockSize : (i+1)*BlockSize - 1],z,histBinsSym)
        del histBinsSym[-1]
        fit = fitLine(histSym, histBinsSym)
        #print("Fit simulation, block ", i+1, " y = ax + b, [a, b]", fit)
        f.write("Fit simulation, block " + str(i+1) + " y = ax + b, [a, b] " + str(fit) +'\n \n')
        plt.plot(histBins, hist, '.b', label="expected")
        plt.plot(histBinsSym, histSym, '.r', label="simulation")
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('k')
        plt.ylabel('P(k)')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.title('block ' + str(i+1))
        plt.savefig("block" + str(i+1) + ".png", dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)
        print(i+1, "done")
        #plt.show()

def plotHHistory(file):
    H = readArray(file)
    plt.plot(H)
    plt.xlabel("% done")
    plt.ylabel("E")
    plt.savefig("hamiltonian.png", dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)     
    #plt.show()

NBlocks, BlockSize, BlockXMin, BlockExponents = readParams('outputParams.dat')
print("NBlocks: ", NBlocks)
print("BlockSize: ", BlockSize)
print("BlockXMin:", BlockXMin)
print("BlocksExponents: ", BlockExponents)
plotHHistory('outputH.dat')
nodeDegrees = readArray('outputInternalND.dat')
nodeDegreesSym = readArray('outputInternalNDSym.dat')
t0 = time.time()
plotAndFitExpectedSimulation(nodeDegrees,nodeDegreesSym,NBlocks,BlockSize)
elapsedTime = time.time() - t0
print("Elapsed time:",time.strftime("%H:%M:%S", time.gmtime(elapsedTime)))