#!/usr/bin/env python
import sys, getopt, os
import array
import re
import ROOT
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex,TRandom3
# from decimal import Decimal
from decimal import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy import optimize
from matplotlib.colors import LinearSegmentedColormap



N=1
# a=120/N
# x0=7

#443x3072,443x2750,443x2500,443x2250,443x2000,443x1750,443x1500,443x1250,443x1000,443x800,443x600,400x1000,300x1000,200x1000,100x1000 cols x rows. Same measurement (27Nov)
MUU=[0.672e-3,0.711e-3,0.785e-3,0.752e-3,0.756e-3,0.802e-3,0.813e-3,1.013e-3,1.120e-3,1.327e-3,1.573e-3,1.237e-3,1.381e-3,1.705e-3,1.966e-3]
MUU=[i*2 for i in MUU] #multiply by 2 because these are 12 hours Amplifier induced events charges per pixel 
E_MUU=[0.042e-3,0.044e-3,0.048e-3,0.049e-3,0.052e-3,0.056e-3,0.060e-3,0.071e-3,0.085e-3,0.096e-3,0.121e-3,0.094e-3,0.116e-3,0.164e-3,0.269e-3]
E_MUU=[i*2 for i in E_MUU] #multiply by 2 because these are 12 hours Amplifier induced events charges per pixel


mu443_600,mu443_800,mu443_1000,mu443_1250,mu443_1500,mu443_1750,mu443_2000,mu443_2250,mu443_2500,mu443_2750,mu443_3072,mu100_1000,mu200_1000,mu300_1000,mu400_1000,a1,x1,c1,summ=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];
A=[90,95,100,105,110,115,120]
X=[8.5,8.3,8,1,7.9,7.7,7.5,7.3,7.1]
C=[0.000775,0.0008,0.000825,0.00085,0.000875,0.0009]
output='output'
counter=0
for a in A:
    for x0 in X:
        for c in C:

            # w, h = 443, 600;
            # w, h = 10, 10;
            # w, h = 443, 1000;
            w, h = 443, 3072;

            SPACE=np.zeros([w,h])
            MU=np.zeros([w,h])
            for x in range(w):
                for y in range(h):
                    MU[x][y]=a/(float(x+x0+.001)**2+float(y+0+.001)**2)+c
            # SPACE=[[0 for x in range(w)] for y in range(h)]
            # MU = [[a/(float(x+7+.001)**2+float(y+7+.001)**2) for x in range(w)] for y in range(h)]


            MU[0][0]=max(MU[1][0],MU[0][1])*1.2 #esto es muy falopa pero es lo que hay


            gRanEng=TRandom3(0)

            events443_600=[[],[]]
            events443_800=[[],[]]
            events443_1000=[[],[]]
            events443_1250=[[],[]]
            events443_1500=[[],[]]
            events443_1750=[[],[]]
            events443_2000=[[],[]]
            events443_2250=[[],[]]
            events443_2500=[[],[]]
            events443_2750=[[],[]]
            events443_3072=[[],[]]
            events100_1000=[[],[]]
            events200_1000=[[],[]]
            events300_1000=[[],[]]
            events400_1000=[[],[]]
            EVENTS=0
            ZEROS=0
            ### finishies initializing for this loop. Simulation starts.

            for ene in range(0,N):
                for n in range(0,len(MU[0,:])):
                    for m in range(0,len(MU[:,0])):
                        # print(str(n))
                        x_poisson = gRanEng.Poisson(MU[m][n])
                        SPACE[m][n]=SPACE[m][n]+x_poisson
                        if ene==N-1:
                            # print("hola"+str(m))
                            if SPACE[m][n]>0.99:
                                EVENTS=EVENTS+SPACE[m][n]
                            else:
                                ZEROS=ZEROS+1.0
                            for i in range(0,int(round(SPACE[m][n]))):
                                # print(str(int(SPACE[m][n])))
                                if n<600:
                                    events443_600[0].append(m)
                                    events443_600[1].append(n)
                                if n<800:
                                    events443_800[0].append(m)
                                    events443_800[1].append(n)
                                if n<1000:
                                    events443_1000[0].append(m)
                                    events443_1000[1].append(n)
                                if n<1000 and m<100:
                                    events100_1000[0].append(m)
                                    events100_1000[1].append(n)
                                if n<1000 and m<200:
                                    events200_1000[0].append(m)
                                    events200_1000[1].append(n)
                                if n<1000 and m<300:
                                    events300_1000[0].append(m)
                                    events300_1000[1].append(n)
                                if n<1000 and m<400:
                                    events400_1000[0].append(m)
                                    events400_1000[1].append(n)
                                if n<1250:
                                    events443_1250[0].append(m)
                                    events443_1250[1].append(n)
                                if n<1500:
                                    events443_1500[0].append(m)
                                    events443_1500[1].append(n)
                                if n<1750:
                                    events443_1750[0].append(m)
                                    events443_1750[1].append(n)
                                if n<2000:
                                    events443_2000[0].append(m)
                                    events443_2000[1].append(n)
                                if n<2250:
                                    events443_2250[0].append(m)
                                    events443_2250[1].append(n)
                                if n<2500:
                                    events443_2500[0].append(m)
                                    events443_2500[1].append(n)
                                if n<2750:
                                    events443_2750[0].append(m)
                                    events443_2750[1].append(n)
                                if n<3072:
                                    events443_3072[0].append(m)
                                    events443_3072[1].append(n)
                            # for i in range(0,int(SPACE[m][n])): #create list of data
                            #     events[0].append(m)
                            #     events[1].append(n)
            
            #file = open(output+'.txt','a')
            #file.write(str(EVENTS/(EVENTS+ZEROS))+' '+str(a)+' '+str(x0)+' \n')
            #file.close() 
            # print(str(EVENTS/(EVENTS+ZEROS)))
            # print(str(float(len(events443_600[0]))/(443*600.0)))
            # print(str(float(len(events443_1000[0]))/(443*1000.0)))
            # print(str(float(len(events443_1500[0]))/(443*1500.0)))
            # print(str(float(len(events443_2000[0]))/(443*2000.0)))
            # print(str(float(len(events443_2500[0]))/(443*2500.0)))
            # print(str(float(len(events443_3072[0]))/(443*3072.0)))

            mu443_600.append(float(len(events443_600[0]))/(443*600.0))
            mu443_800.append(float(len(events443_800[0]))/(443*800.0))
            mu443_1000.append(float(len(events443_1000[0]))/(443*1000.0))
            mu100_1000.append(float(len(events100_1000[0]))/(100*1000.0))
            mu200_1000.append(float(len(events200_1000[0]))/(200*1000.0))
            mu300_1000.append(float(len(events300_1000[0]))/(300*1000.0))
            mu400_1000.append(float(len(events400_1000[0]))/(400*1000.0))
            mu443_1250.append(float(len(events443_1250[0]))/(443*1250.0))
            mu443_1500.append(float(len(events443_1500[0]))/(443*1500.0))
            mu443_1750.append(float(len(events443_1750[0]))/(443*1750.0))
            mu443_2000.append(float(len(events443_2000[0]))/(443*2000.0))
            mu443_2250.append(float(len(events443_2250[0]))/(443*2250.0))
            mu443_2500.append(float(len(events443_2500[0]))/(443*2500.0))
            mu443_2750.append(float(len(events443_2750[0]))/(443*2750.0))
            mu443_3072.append(float(len(events443_3072[0]))/(443*3072.0))
            c1.append(c) #stores every c for every simulation
            x1.append(x0) #stores every x0 for every simulation
            a1.append(a) #stores every a for every simulation
            mu_temp=[mu443_3072[-1],mu443_2750[-1],mu443_2500[-1],mu443_2250[-1],mu443_2000[-1],mu443_1750[-1],mu443_1500[-1],mu443_1250[-1],mu443_1000[-1],mu443_800[-1],mu443_600[-1],mu400_1000[-1],mu300_1000[-1],mu200_1000[-1],mu100_1000[-1]]
            print(mu_temp)
            auxsumm=[np.power(MUU[p]-mu_temp[p],2)/np.power(E_MUU[p],2) for p in range(0,len(MUU))]
            summ.append(abs(np.sum(auxsumm)/len(MUU)-1))

            counter=counter+1
            print(counter)

ind,cmin,xmin,amin,summmin=[],[],[],[],[]

for i in range(0,int(counter/2.0)):
    ind.append(np.argmin(summ))

    cmin.append(c1[ind[i]]) #add them and store them
    amin.append(a1[ind[i]])
    xmin.append(x1[ind[i]])
    summmin.append(summ[ind[i]])

    c1.pop(ind[i]) #delete them so we can look for the second minimum and so on.
    x1.pop(ind[i])
    a1.pop(ind[i])
    summ.pop(ind[i])


with open('summmin.txt', 'w') as f:
    for item in summmin:
        f.write("%s\n" % item)

with open('c1.txt', 'w') as f:
    for item in cmin:
        f.write("%s\n" % item)

with open('x1.txt', 'w') as f:
    for item in xmin:
        f.write("%s\n" % item)

with open('a1.txt', 'w') as f:
    for item in amin:
        f.write("%s\n" % item)


# # # events=np.asarray(events)
# # # events=np.transpose(events)
# # # colors = [(0, 0, 0), (0, 1, 0), (0, 0, 1)]  # R -> G -> B
# # # n_bins = [100,100,100,100]
# # # cmap_name = 'my_list'
# # # fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(6, 9))
# # # fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
# # # for n_bin, ax in zip(n_bins, axs.ravel()):
# # #     cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin) # Fewer bins will result in "coarser" colomap interpolation
# # #     im = ax.hist2d(events[:, 0], events[:, 1], bins=(w,h),cmap=cm)
# # #     fig.tight_layout()
# # #     fig.colorbar(im[3], ax=ax)
# # # plt.show()