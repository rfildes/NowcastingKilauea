  # -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:02:59 2018

@BeckyFildes
"""

#This script has a class of functions that include the 
#2 GR methods, next 3 steps for Nowcasting and the RMSE calculation

#For each subcase:
#Select the correct data input file (lines 34/35 and 41/42)
#Select the correct time window (lines 295-297)
#Set the M_lambda large magnitude value as 3.5 or 4.7 (lines 295-297)

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats
import scipy.optimize as optimization
from datetime import datetime
from sklearn.metrics import mean_squared_error

plt.rcParams["font.size"] = "20"
plt.rcParams["font.family"] = "Times New Roman"

class KilaueaEQs:
    def __init__(self, large_mag, small_mag, start, end):
        self.large_mag = large_mag
        self.small_mag = small_mag
        self.start = start
        self.end = end
        
        #data_file = open('Kilauea4.txt', 'r') #for S subcases, only includes data from within 5 km of the summit caldera
        data_file = open('HawaiiIsland.txt', 'r') #for I subcases, includes data from full island
        total = 0
        for line in data_file:
            total +=1
        data_file.close()
            
        #data_file = open('Kilauea4.txt', 'r') #for S subcases, only includes data from within 5 km of the summit caldera
        data_file = open('HawaiiIsland.txt', 'r') #for I subcases, includes data from full island
        quakes = np.zeros((total, ), dtype=[('year', 'int'), ('month', 'int'), ('day', 'int'), ('hour', 'int'), ('minute', 'int'), ('second', 'int'),('mag', 'float')])

        i =0
        for line in data_file:
            
            line = line.strip()
            
            data = line.split()
            
            date_string = data[0]
            date_array = date_string.replace('-', ' ').split()
            time_string = data[1]
            time_array = time_string.replace(':', ' ').split()
            quakes['hour'][i] = int(time_array[0])
            quakes['minute'][i] = int(time_array[1])
            quakes['second'][i] = int(float(time_array[2]))
            quakes['year'][i] = int(date_array[0])
            quakes['month'][i] = int(date_array[1])
            quakes['day'][i] = int(date_array[2])
            quakes['mag'][i] = float(data[2])
            
            i +=1
                
        data_file.close()
        
        self.quakes = quakes
        mag= []
        small_mag_times = []
        large_mag_times = []
        red_small_count = []
        
        i=0
        j=0
        d0 = datetime(2018, 5, 4, 22, 32, 54) #date and time of M6.9 mainshock 
        for item in self.quakes:
            d1 = datetime(item['year'], item['month'], item['day'], item['hour'], item['minute'], item['second'])
            time_after_main = (d1-d0).days + (d1-d0).seconds/86400
            if (time_after_main >= self.start) and (time_after_main <= self.end): #time delay you want for start of your sequence
                if item['mag'] >= 1:
                    mag.append(item['mag']) 
                if item['mag'] >= self.small_mag:
                    i +=1
                    small_mag_times.append(time_after_main)
                if item['mag'] >= self.large_mag:
                    j +=1
                    large_mag_times.append(time_after_main)
                    red_small_count.append(i)
     
        self.Mags=np.asarray(mag)    #list of all events greater than or equal to Magnitude = 1
        self.small_count = np.arange(i)
        self.large_count = np.arange(j)
        
        self.small_times = np.asarray(small_mag_times)
        self.large_times = np.asarray(large_mag_times)
        self.red_small_count = np.asarray(red_small_count)
        self.plot_times = self.large_times - self.start
        self.small_plot_times = self.small_times - self.start
        
################################################################################################################    
    def plot_GR(self):  
        #uses linear least squares for GR relationship          
        #bin the data to bins of size 0.1 for rounding, the first bin are those that round to 2.5, so the binned data is mags: [2.45, 2.55), 
        fullcount, fullbins = np.histogram(self.Mags, bins=36, range=(1.95, 5.55))  
        #name of bins are the midpoint of each bin, 2.5, 2.6, ...5.5
        Bins = np.arange(2,5.6,0.1)

        self.fullcount=np.flip(fullcount)        #flips order or array so largest mag bins are first
        self.fullBins = np.flip(Bins)            #flips order or array so largest mag bins are first
        self.fullsum = np.cumsum(self.fullcount) #takes cumulative sum of each bin to build next bin
        self.FullCount = np.log10(self.fullsum)  #log10 of cumulative count of bins

        #for Max Mag 4.7
        fullbins47=self.fullBins[8:31] #[) notation
        fullcount47=self.FullCount[8:31]
        slope, intercept, r_value, p_value, std_err = stats.linregress(fullbins47, fullcount47)
        
        #for Max Mag 3.5
        fullbins35=self.fullBins[20:31]
        fullcount35=self.FullCount[20:31]
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(fullbins35, fullcount35)
        
        
        #for Max Mag 4.0
        fullbins40=self.fullBins[15:31]
        fullcount40=self.FullCount[15:31]
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(fullbins40, fullcount40)
      
        #best fit labels
        string1 ='{:.2f}'.format(slope) 
        string2='M + '
        string3='{:.2f}'.format(intercept)
        string4='Log(N$_{c}\geq$M) = '
        label1=string4 + string1 + string2 + string3

        string5 ='{:.2f}'.format(slope2) 
        string6='{:.2f}'.format(intercept2)
        label2= string4 + string5 + string2 + string6
        
        string7 ='{:.2f}'.format(slope3) 
        string8='{:.2f}'.format(intercept3)
        label3= string4 + string7 + string2 + string8
        
        #Plotting the actual data
        plt.figure(figsize=(8,5.5))
        plt.plot(self.fullBins, self.FullCount, 'go')
        
        #Plotting the best fit line
        #plt.plot(fullbins47, slope*fullbins47+intercept, 'r', label=label1, linewidth = 2.2)
        #plt.plot(fullbins35, slope2*fullbins35+intercept2, 'g', label=label2, linewidth = 2.2)
        plt.plot(fullbins40, slope3*fullbins40+intercept3, 'b', label=label3, linewidth = 2.2)
        
        plt.xlabel('Magnitude')
        plt.ylabel('Log(N$_{c}\geq$M)')
        plt.xlim([1.8,7])
        plt.xticks([2,2.5,3,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0])
        plt.ylim([0,4.6])
        plt.yticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5])
        plt.legend(loc='upper right')
        #plt.title('Gutenberg-Richter: days ' + f'{self.start:.2f}' + '-' + f'{self.end:.2f}' + ' Island')
        #pdffile='GR' +  '_' +'bin'+ ' summit_Manuscript' + '.pdf'
        plt.tight_layout()
        #plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
    def maxliklihood(self):
        #maximum liklihood method for GR relationshiop
        #bin the data to bins of size 0.1 for rounding, the first bin are those that round to 2.5, so the binned data is mags: [2.45, 2.55), 
        self.count, bins = np.histogram(self.Mags, bins=16, range=(2.45, 4.05))
        #name of bins are the midpoint of each bin, 2.5, 2.6, ...4.0
        self.Bins = np.arange(2.5,4.1,0.1) #for M = 4.0 
     
        M1= 2.5 #catalog completeness)
        dM= 0.1 # bin size of magnitudes

        #mean of binned mags
        M_avg = sum(self.Bins*self.count)/sum(self.count)
       
        #Maximum Liklihood using Eq 4 from Nava etal. 2017
        b = np.log10(math.e)/(M_avg - (M1-dM/2) )
        print(b)    
   
    def plot_naturaltime(self):
        #defining the nautral time relationship for the data        
        def f(x, a):
            return a*x #fit of equation 4 where a here equals the base 10 term
        pfit, pcov = optimization.curve_fit(f, self.red_small_count, self.large_count)
        perr = np.sqrt(pcov)
        self.perr = perr
        self.pfit = pfit
        self.bval = math.log10(pfit[0])/(self.small_mag-self.large_mag) #natural time b-value calcuated from equation 4
       
        #label for best fit
        string1 = '$\mathrm{N_{c\lambda}} = ($'
        string2 = '{:.4f}'.format(self.pfit[0])
        string3 = ' $\pm$ '
        string4 = '{:.4f}'.format(self.perr[0][0])
        string5 = '$)\mathrm{N_{c\sigma}}$ '
        self.label_string = string1 +string2 + string3 + string4 + string5
       
        #plotting
        fig = plt.figure(figsize=(8,6))        
        plt.plot(self.red_small_count, self.large_count, 'bo')
        plt.plot(self.red_small_count, pfit*self.red_small_count, 'r-', label = self.label_string )        
        #plt.yticks([0,10,30,50,70])
        plt.xticks([0,2500,5000,7500,10000,12500])
        plt.xlabel('N$\mathrm{_{c\sigma}}$' + ' (M$\geq$' + f'{self.small_mag:.1f})', fontsize=20, fontweight=700)
        plt.ylabel('N$\mathrm{_{c\lambda}}$' + ' (M$\geq$' + f'{self.large_mag:.1f})',fontsize=20, fontweight='bold')
        plt.legend(loc = 'upper left', fontsize = 20)
        
        #plt.title('May 17th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        #plt.title('May 29th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        plt.title('June 14th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
         
        pdffile= 'ManKilauea_naturaltime' + f'{self.start:.2f}' + '_' + f'{self.end:.2f}'+ 'sized' + '.pdf'
        plt.tight_layout()
        #fig.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
         
    def plot_small_cumulative(self):
        #plotting the small magnitude earthquakes as a function of clock time
        fig = plt.figure(figsize=(8.7,6))
        
        plt.plot(self.small_plot_times, self.small_count, 'b')
        plt.xlabel('t (days)', fontsize=20, fontweight='bold')
        plt.ylabel('N$\mathrm{_{c\sigma}}$' + ' (M$\geq$' + f'{self.small_mag:.1f})',fontsize=20, fontweight='bold')
        #plt.xticks([40,50,60,70,80,90])
        plt.yticks([0,2500,5000,7500,10000,12500,15000])
        
        #plt.title('May 17th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        #plt.title('May 29th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        plt.title('June 14th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        
        pdffile= 'ManKilauea_small_cumulative' + f'{self.start:.2f}' + '_' + f'{self.end:.2f}' + 'sized' + '.pdf'
        plt.tight_layout()
        #fig.savefig(pdffile, bbox_inches='tight')
        plt.show()
         
    def plot_large_cumulative(self):
        #plotting the large magnitude observed collapse events and nowcasted collapse events as fucntions of clock time
        self.obs = self.large_count #observed
        self.forc = self.red_small_count * self.pfit #nowcasted
        
        #plotting
        fig = plt.figure(figsize=(8,6))
        plt.plot(self.plot_times, self.forc, 'r.', label = 'Nowcasted')
        plt.plot(self.plot_times, self.large_count, '.b', label = 'Actual')
        plt.ylabel('N$\mathrm{_{c\lambda}}$' + ' (M$\geq$' + f'{self.large_mag:.1f})',fontsize=20, fontweight='bold')
        plt.xlabel('t (days)',fontsize=20, fontweight='bold')
        plt.legend(loc='lower right')

        #plt.title('May 17th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        #plt.title('May 29th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        plt.title('June 14th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        
        pdffile= 'ManKilauea_NowcastResult' + f'{self.start:.2f}' + '_' + f'{self.end:.2f}' +'.pdf'
        plt.tight_layout()
        #plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
    def RMS(self):
        #calculating and plotting CDFs of theobserved and nowcasted collapse events
        #CDF observed self.obs
        obspdf=self.obs/sum(self.obs)
        obscdf=np.cumsum(obspdf)
        #CDF of forcasted self.forc
        forcpdf=self.forc/sum(self.forc)
        forccdf=np.cumsum(forcpdf)
        
        #plotting
        fig = plt.figure(figsize=(8,6))
        plt.plot(self.obs, obscdf, 'b.', label='Actual')
        plt.plot(self.forc, forccdf, 'r.', label = 'Nowcasted')
        plt.ylim([-0.02,1.01])
        plt.ylabel('CDF',fontsize=20, fontweight='bold')
        plt.xlabel('N$\mathrm{_{c\lambda}}$' + ' (M$\geq$' + f'{self.large_mag:.1f})',fontsize=20, fontweight='bold')
        plt.legend()
        
        #plt.title('May 17th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        #plt.title('May 29th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
        plt.title('June 14th, 2018 - August 2nd, 2018', fontsize=20, fontweight='bold')
              
        pdffile= 'ManKilauea_CDF' + f'{self.start:.2f}' + '_' + f'{self.end:.2f}' +'.pdf'
        plt.tight_layout()
        #plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
               
        #RMSE Calculation
        self.RMSCDF=math.sqrt(np.average(np.square(forccdf-obscdf))) #equation 7
        #print(self.RMSCDF)
        
#Three different time windows
#myhawaii = KilaueaEQs(4.7, 2.5, 12.2379, 89.9739)  #May 17-Aug 2 subcases a
#myhawaii = KilaueaEQs(4.7, 2.5, 24.5578, 89.9739)  #May 29-Aug 2 subcases b
myhawaii = KilaueaEQs(4.7, 2.5, 40.6158, 89.9739)   #Jun 14- Aug 2 subcases c

myhawaii.plot_GR()
myhawaii.maxliklihood()
myhawaii.plot_naturaltime() 
myhawaii.plot_small_cumulative()
myhawaii.plot_large_cumulative()
myhawaii.RMS()


