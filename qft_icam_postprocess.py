#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculation of absorption spectra by total particles, non-algal particles and
phytoplankton from QFT-ICAM .raw files.

                                                                     
About QFT-ICAM sample labeling: 
1) Use uniform labeling method for all samples. For example, the 2 measurements 
    of the sample "UW5" should be labeled as "UW5#1" and "UW5#2" (here separator
    is '#' and should be used throughout all samples) or simply both as "UW5". 
2) Labels for blank filters measurements for contamination control of each box 
    of GF/F filters (they are NOT reference measurements!) need to contain 
    keywords "batch", "MQ" or "SW" (case insensitive), e.g."*batch_MQ#11". 
    Normal sample labels should not contain any aforementioned keywords.
                                                                     

Before using this script, please 
1) install:
- Python 3.6x or higher version    
- Essential Python packages: numpy, scipy, pandas, matplotlib.
2) prepare files (see examples as well as the main program) of:
- config_PostProc.txt
- rmSpectraIndex.txt
- qft_icam_matched_labels_filtration_volumn.txt
- labels_latlon_datetime_qft.txt

Analysis method detailed in:
1) Röttgers et al. (2016). Quantitative filter technique measurements of spectral 
light absorption by aquatic particles using a portable integrating cavity 
absorption meter (QFT-ICAM). Optics express, 24(2), A1-A20.
2) IOCCG Protocol Series (2018). Inherent Optical Property Measurements and 
Protocols:Absorption Coefficient, Neeley, A. R. and Mannino, A. (eds.), IOCCG 
Ocean Optics and Biogeochemistry Protocols for Satellite Ocean Colour Sensor 
Validation, Volume 1.0, IOCCG, Dartmouth, NS, Canada.

@author: 
    Rüdiger Röttgers (ruediger.roettgers@hzg.de) (raw data process);
    Yangyang Liu (yangyang.liu@awi.de) (postprocess), January 2020.
"""

import glob, os, shutil, math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy import interpolate
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def readConfig(config):
    paramDict = {}
    try:
        with open(config) as f:
            for line in f:
                if line.startswith('#')==False and len(line)>1:
                    key = line.strip().split('=')[0]
                    try:
                        value = line.strip().split('=')[1].strip()
                    except:
                        value = None
                    
                    paramDict[key] = value
    except Exception as e:
        print(e)
    finally:
        return paramDict
    
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------- 
class qftPostProc():
    
    def __init__(self, config='config_PostProc.txt'):
        paramDict = readConfig(config)
        xi = paramDict['xi']
        wlstart, increment, wlend = map(float, xi.split(':'))
        self.xi = np.arange(wlstart, wlend+increment, increment)
        fluorkorr = paramDict['fluorkorr']
        self.fl_wlstart, self.fl_wlend = map(int, fluorkorr.split(':'))

        self.separator = paramDict['separator']
        self.diameter = float(paramDict['diameter'])
        self.beta = float(paramDict['beta'])
        self.externalvolumn = paramDict['externalvolumn']
        self.concurrentInfo = paramDict['concurrentInfo']
        
        self.ODpFile = os.path.join('ODfiles','qft_icam_merged_median_OD_totparticle.txt')
        self.ODpFile_sd = os.path.join('ODfiles','qft_icam_merged_sd_OD_totparticle.txt')
        self.ODdFile = os.path.join('ODfiles','qft_icam_merged_median_OD_NAP.txt')
        self.ODdFile_sd = os.path.join('ODfiles','qft_icam_merged_sd_OD_NAP.txt')
        self.apFile = 'qft_icam_merged_median_abs_totparticle.txt'
        self.apFile_sd = 'qft_icam_merged_sd_abs_totparticle.txt'
        self.adFile = 'qft_icam_merged_median_abs_NAP.txt'
        self.adFile_adjust = 'qft_icam_merged_median_abs_NAP_adjust.txt'
        self.adFile_sd = 'qft_icam_merged_sd_abs_NAP.txt'
        self.aphFile = 'qft_icam_merged_median_abs_phytoplankton.txt'
    
    #--------------------------------------------------------------------------
    def mk_spec(self, line, line_dark):
        """ 
        This function subtracts a dark measurement from either a reference or 
        sample measurement, and interpolate the data in the range of "xi".
        """
        data = line.split('[ms]=')[1].split()
    #    inttime = float(data[0])
        spec = np.array(data[1:], dtype='float')        
        data_dark = line_dark.split('[ms]=')[1].split()
    #    inttime = float(data_dark[0])    
        spec_dark = np.array(data_dark[1:], dtype='float')       
        spec = spec - spec_dark
        #simple stray light correction, subtracting first 20 values
        spec = spec - np.mean(spec[0:20])        
        if len(self.xi) == len(spec) :
            return spec
        elif  len(self.wl_raw) == len(spec):
            f = interpolate.interp1d(self.wl_raw, spec)
            return f(self.xi)
        else:
            print ('>>mk_spec ERROR! Lengths of measured spectrum and xi differ!')
    
    #--------------------------------------------------------------------------    
    def output(self, fw, line_spec, spectra):
        '''
        This function writes spectral data in a file fw.
        '''
        labels = line_spec.split()
        fw.write('%s %s %s %s'%tuple(labels[:4]))
        for spec in spectra:
            fw.write(' %.6f'%spec)
        fw.write('\n')
    
    #--------------------------------------------------------------------------
    def processRaw(self, filename):       
        base = os.path.basename(filename)  
        print('Processing ', base)
        filename_l0 = base.replace('.qft.raw','.l0a')          
        with open(filename, 'r') as f:
            lines = f.readlines()             
        with open(filename_l0, 'w') as fw:  
            for i, line in enumerate(lines):
                wl_old=[]
                if 'wavelen'in line:
                    try:
                        self.wl_raw = np.array(line.split(':')[1].split(), 
                                                   dtype='float')
                    except:
                        self.wl_raw = np.array(line.split(':')[1].split()[:-1],
                                                   dtype='float')   
                    if len(self.wl_raw)==len(wl_old):
                        #avoiding repeated printing of wavelengths in output
                        print(wl_old)
                    else:
                        fw.write('%wavelength: 0 0 0 ')
                        for wl in self.xi:
                            fw.write(' %i'%wl)
                        fw.write('\n')
                    wl_old = self.wl_raw    
                elif line.startswith('%'):
                    fw.write(line)
                elif 'samT_dark' in line:
                    idarkTsam, isampleT = i, i+1
                    try:
                        nexttwolines = lines[i+2]
                        if 'samT_filter' in nexttwolines:
                            #samples
                            isamTf,isampleD, isamDf = i+2, i+3, i+4 
                            try:
                                samT_spec = self.mk_spec(lines[isampleT], lines[idarkTsam])
                                samD_spec = self.mk_spec(lines[isampleD], lines[idarkTsam])
                                samTf_spec = self.mk_spec(lines[isamTf], lines[idarkTsam])
                                samDf_spec = self.mk_spec(lines[isamDf], lines[idarkTsam])
                            except:
                                pass
                            fluor_filt = 1
                        else:
                            #bleached samples
                            isampleD, isamDf = i+2, i+3
                            try:
                                samT_spec = self.mk_spec(lines[isampleT], lines[idarkTsam])
                                samD_spec = self.mk_spec(lines[isampleD], lines[idarkTsam])
                            except:
                                pass
                            fluor_filt = 0 
                    except:
                        pass
                    ##OD calculated with reference before sample
                    k=1
                    while not 'refT_dark' in lines[i-k]:
                        k=k+1 
                    idarkTref, irefT, irefTf, irefD, irefDf = i-k, i-k+1, \
                    i-k+2, i-k+3, i-k+4
                    try:
                        refT_spec = self.mk_spec(lines[irefT], lines[idarkTref])
                        refD_spec = self.mk_spec(lines[irefD], lines[idarkTref])
                        refTf_spec = self.mk_spec(lines[irefTf], lines[idarkTref])
                        refDf_spec = self.mk_spec(lines[irefDf], lines[idarkTref])
                    except:
                        pass
                    
                    #fluorescence correction for samples
                    if fluor_filt == 1:
                        _transm = samT_spec / refT_spec
                        abs_photon = np.sum(refT_spec - samT_spec)
                        abs_photon_fluor = np.sum(refTf_spec - samTf_spec)    
                        ratio = abs_photon / abs_photon_fluor + 0.1
                        Tf_corr = (samTf_spec - refTf_spec*_transm) * ratio
                        from_index = int(np.argwhere(self.xi==self.fl_wlstart)[0])  #nm
                        Tf_corr[:from_index] = 0   #set fluor_corr<670nm to 0
                        to_index = int(np.argwhere(self.xi==self.fl_wlend)[0] ) #nm
                        Tf_corr[to_index:] = 0  #set fluor_corr>800nm to 0
                        samT_spec1 = samT_spec - Tf_corr
            
                        _transm = samD_spec / refD_spec
                        abs_photon = np.sum(refD_spec - samD_spec)
                        abs_photon_fluor = np.sum(refDf_spec - samDf_spec)    
                        ratio = abs_photon/abs_photon_fluor + 0.1
                        Df_corr = (samDf_spec - refDf_spec*_transm)*ratio
                        Df_corr[:from_index] = 0   #set fluor_corr<670nm to 0
                        Df_corr[to_index:] = 0  #set fluor_corr>800nm to 0
                        samD_spec1 = samD_spec - Df_corr
                        ODT = -np.log(samT_spec1/refT_spec)/2.303
                        ODD = -np.log(samD_spec1/refD_spec)/2.303
                    else: #no fluorescence correction for bleached samples
                        ODT = -np.log(samT_spec/refT_spec)/2.303
                        ODD = -np.log(samD_spec/refD_spec)/2.303
                        
                    OD = ODT - ODD
                    self.output(fw, line, OD) #output of OD calculated with reference before sample

                    ##OD calculated with next reference after sample
                    try:
                        k=1
                        while not 'refT_dark' in lines[i+k]:
                            k=k+1
                        
                        idarkTref, irefT, irefTf, irefD, irefDf = i+k, i+k+1, \
                        i+k+2, i+k+3, i+k+4                                    
                        refT_spec = self.mk_spec(lines[irefT], lines[idarkTref])
                        refD_spec = self.mk_spec(lines[irefD], lines[idarkTref])
                        refTf_spec = self.mk_spec(lines[irefTf], lines[idarkTref])
                        refDf_spec = self.mk_spec(lines[irefDf], lines[idarkTref])

                        #fluorescence correction for samples
                        if fluor_filt == 1:
                            _transm = samT_spec / refT_spec
                            abs_photon = np.sum(refT_spec-samT_spec)
                            abs_photon_fluor = np.sum(refTf_spec - samTf_spec)    
                            ratio = abs_photon / abs_photon_fluor + 0.1
                            Tf_corr = (samTf_spec - refTf_spec*_transm) * ratio
                            from_index = int(np.argwhere(self.xi==self.fl_wlstart)[0])  #nm
                            Tf_corr[:from_index] = 0   #set fluor_corr<670nm to 0
                            to_index = int(np.argwhere(self.xi==self.fl_wlend)[0]) #nm
                            Tf_corr[to_index:] = 0  #set fluor_corr>800nm to 0
                            samT_spec2 = samT_spec-Tf_corr
                
                            _transm = samD_spec/refD_spec
                            abs_photon = np.sum(refD_spec-samD_spec)
                            abs_photon_fluor = np.sum(refDf_spec - samDf_spec)    
                            ratio = abs_photon/abs_photon_fluor + 0.1
                            Df_corr = (samDf_spec - refDf_spec*_transm)*ratio
                            Df_corr[:from_index]=0   #set fluor_corr<670nm to 0
                            Df_corr[to_index:]=0  #set fluor_corr>800nm to 0
                            samD_spec2 = samD_spec - Df_corr
                            ODT = -np.log(samT_spec2/refT_spec)/2.303
                            ODD = -np.log(samD_spec2/refD_spec)/2.303
                        else: #no fluorescence correction for bleached samples
                            ODT = -np.log(samT_spec/refT_spec)/2.303 
                            ODD = -np.log(samD_spec/refD_spec)/2.303
                
                        OD = ODT - ODD
                        self.output(fw, line, OD) #output of OD calculated with next reference after sample 
                    except:
                        pass
        
    #--------------------------------------------------------------------------    
    def getWL(self, filename):
        with open(filename, 'r') as f:
            for line in f:
                if 'wavelength' in line:
                   break        
        wl = pd.Series(line.split(' ')[5:]).astype('float')
        
        wl_selected = [350, 400, 650, 676, 700, 710, 715, 750, 800]
        wlpos = {}
        for i, wvl in enumerate(wl_selected):
            key = 'pos'+str(wvl)
            value = np.where(wl>=wvl)[0][0]
            wlpos[key] = value
            
        return wl, line, wlpos
    
    #--------------------------------------------------------------------------
    def extractOD(self, filename):
        #data
        data = pd.read_csv(filename, header=None, comment='%', sep=' ') 
        #label
        for i in range(len(data)):
            if data.iloc[i,1] == 'samT_dark':
                data.iloc[i,2] = 'sam_' + data.iloc[i,2]
            else:
                data.iloc[i,2] = 'BL_' + data.iloc[i,2]
        label = data.iloc[:,2]
        label2 = label.copy()
        for i,lb in enumerate(label):
            if (self.separator !='') & (self.separator in lb):
                label2.iloc[i] = lb[:lb.find(self.separator)]
            else:
                label2.iloc[i] = lb
        label = label2 
        #filtration volumn
        filt_vol = data.iloc[:,3]
        filt_vol = [float(filt_vol[i].replace('filt_vol=','')) for i in 
                        range(len(filt_vol))]
        return data, label, filt_vol
    
    #--------------------------------------------------------------------------
    def plotSpectrum(self, dfx, dfy, figname, ylabel, xlabel='Wavelength [$nm$]', 
                     legend=None):
        fig = plt.figure(figsize=(8.5, 6.5))
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.tick_params(axis="both", labelsize=10)
        ax.plot(dfx, dfy, linewidth=1)
        if legend is not None:
            ax.legend(legend, loc='upper left', bbox_to_anchor=(1, 1),
                      fontsize='small')
        ax.set_xlabel(xlabel,fontsize=12)
        ax.set_ylabel(ylabel,fontsize=12)
        fig.savefig(figname, dpi=200)
        plt.close(fig)
        
    #--------------------------------------------------------------------------
    def saveData(self, filename, dfvarname, headline=None, index=True, 
                 header=False, sep='\t'): 
        with open(filename, 'w') as f:
            if headline is not None:
                f.write(headline)
            dfvarname.to_csv(f, index=index, header=header, sep=sep, 
                             encoding='utf-8') 
    
    #--------------------------------------------------------------------------
    def createDir(self, dirname, overwrite=False):
        if overwrite:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            else:
                shutil.rmtree(dirname) #removes all the existing directories!
                os.makedirs(dirname)
        else:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
  
    #--------------------------------------------------------------------------
    def plotOD_l0a(self, filename):
        '''
        This function plots QFT-ICAM optical density (OD) of each sample in a .l0a file.
        Input:
        filename - str. Path of a QFT-ICAM .l0a file.
        Output:
        A folder named as the base name of the .l0a file in the working direcoty 
        containing plots from individual sample measurements. The legend of each 
        absorption spectrum in every plot gives the row index of each data point 
        when this .l0a data file is imported into a pandas dataframe via 
        "pd.read_csv". These row indices can be used to prepare the file 
        "rmSpectraIndex.txt". Note: if the number of samples (or bleached samples) 
        is greater than 4 (or 2), check the labelling in .raw data!
        '''
        foo = filename.split('/')
        print(f'Importing {foo[-1]}...')
        
        #get wavelengths
        wl,line, wlpos = self.getWL(filename)       
        try:
            #Extract OD spectra
            data, label, filt_vol = self.extractOD(filename)
            label_uni = np.unique(label,return_index=True, return_inverse=True, 
                               return_counts=True)
            #Create target Directory 
            dir_fig = foo[-1][:-4]
            self.createDir(dir_fig, overwrite=True)

            #plot OD spectra
            for i,lb in enumerate(label_uni[0]):
                pos = np.where(label_uni[2] == i)[0]            
                data_for_plot = data.iloc[pos, 4:].transpose()
                lgd = [str(idx) for idx in pos]
                figname = os.path.join(dir_fig, 'check_'+ lb + '.png')
                self.plotSpectrum(wl, data_for_plot, figname, legend=lgd,
                                  ylabel='$Level\ 0\ uncorrected\ OD$', ) 
            print('Plotting Done!') 
        except:
            print(f'Warning: No data found in: {foo[-1]}!!! \n')

    #--------------------------------------------------------------------------    
    def rmOD_median(self, filename, row_index=None, plot=True, smooth=True):
        '''
        This function removes suspicious measurements in QFT-ICAM .l0a data indexed 
        by "row_index" and take the median values and standard deviation of 
        repeatedly measured data.
        Input:
        1) filename - str. Path of a QFT-ICAM .l0a file.
        2) row_index - list, default None. List of numbers as row indices of 
        suspicious data points when this .l0a data file is imported into a pandas 
        dataframe via "pd.read_csv". The indices of each data point are annotated 
        in the legend of the plots produced by the function "plotOD_l0a". 
        3) plot - bool, default True. Plot the output OD spectra.
        4) smooth - bool, default True. Method: Savitzky Golay first order 
        polynomial filter (window size, 5).
        Output:
        1) Tab delimited files in the direcoty "ODfiles" named as 
        "*_OD_median.txt" and "*_OD_sd.txt", respectively, with "*" standing for the 
        base name of .l0a files.
        2) Plots of the spectra in the "*_OD_median.txt" and "*_OD_sd.txt" named as 
        "*_OD_median.png" and "*_OD_sd.png", respectively, in the directory produced 
        by the function "plotOD_l0a".
        '''
        
        filename = filename.strip()
        foo = filename.split('/')
        dir_fig = foo[-1][:-4]
        print(f'Importing {foo[-1]}...')
        wl, line, wlpos = self.getWL(filename)
        #Extract ODf spectra data
        data, label, filt_vol = self.extractOD(filename)
        data.iloc[:,2] = np.array(label)
        data.iloc[:,3] = filt_vol         
        
        if not math.isnan(row_index[0]):
            data.loc[row_index, :] = np.nan
        
        #take the median of repeated measurements in each file
        data_median = data.groupby(pd.Grouper(key=2)).median() 
        data_median.iloc[:,0] = data_median.index
        
        #smooth
        if smooth:
            figname = os.path.join(dir_fig, dir_fig + '_smooth_OD_median.png')
            outname_median = filename.replace('.l0a','_smooth_OD_median.txt') 
            data_median = data_median.interpolate(limit_direction ='both')
            data_median_smooth = savgol_filter(data_median.iloc[:,2:], 5, 1)
            data_median_smooth = pd.DataFrame(data=data_median_smooth, index=
                                        data_median.index)
            data_median = pd.concat([data_median.iloc[:,0], data_median.iloc[:,1], 
                                    data_median_smooth], axis=1)
        else:
            figname = os.path.join(dir_fig, dir_fig + '_OD_median.png')
            outname_median = filename.replace('.l0a','_OD_median.txt')  
        
        #Create target Directory 
        self.createDir('ODfiles')
        self.saveData(os.path.join('ODfiles',outname_median), data_median, 
                      headline=line, index=False) 
       
        #take the standard deviation of repeated measurements in each file
        data_sd = data.groupby(pd.Grouper(key=2)).agg(np.std, ddof=1)
        data_sd.iloc[:,0] = data_sd.index
        data_sd.loc[:,3] = data_median.iloc[:,1]
        data_sd = data_sd.drop(columns=[1])

        outname_sd = filename.replace('.l0a','_OD_sd.txt')  
        self.saveData(os.path.join('ODfiles',outname_sd), data_sd, 
                      headline=line, index=False) 
        
        if plot:
            #plot acdom spectra
            OD = data_median.iloc[:,1:]
            wl_for_plot = wl.iloc[wlpos['pos400']:wlpos['pos750']+1]
            OD_for_plot = OD.iloc[:,wlpos['pos400']:wlpos['pos750']+1].transpose()
            lgd = data_median.iloc[:,0].tolist()
            if len(OD) != 0:
                self.plotSpectrum(wl_for_plot, OD_for_plot, figname, ylabel='$OD$',
                                  legend=lgd)
            else:
                print(f'No plot: after removing suspicious data, no data were left \
                      in {foo}!')
        
        print(f'Removal of suspicious OD spectra and averaging for {foo[-1]} Finished! \n') 
        
    #--------------------------------------------------------------------------     
    def mergeOD(self, keyword='median', plot=True): 
        '''
        This function merge all '*_OD_median.txt' and '*_OD_sd.txt' files, 
        separates total particles (sample), non-algal particles (bleached sample) 
        and contamination control blank filter measurements, and merges data from 
        each catalog, respectively.
        
        The prerequisit for data separation is that labels from blank filter 
        measurements contain the keyword "batch" or "bt" or "MQ" or "SW" (case 
        insensitive), and normal sample labels do not contain any aforementioned keywords.
        Input:
        1) keyword - str, default 'median'. 'meidan' and 'sd' stand for the files to be 
        merged are median or standard deviation of repeterd OD measurements, respectively.
        2) plot - bool, default True. Plot all sample spectra.
        Output:
        1) Eight tab delimited files in the directory "ODfiles" named as 
        "qft_icam_merged_median_OD_all.txt" (OD spectra from all measurements), 
        "qft_icam_merged_sd_OD_all.txt" (standard deviation of OD spectra from all measurements), 
        "qft_icam_merged_median_OD_totparticle.txt" (OD spectra from samples), 
        "qft_icam_merged_sd_OD_totparticle.txt" (standard deviation of OD spectra from samples), 
        "qft_icam_merged_median_OD_NAP.txt" (OD spectra from bleached samples),
        "qft_icam_merged_sd_OD_NAP.txt" (standard deviation of OD spectra from bleached samples),
        "qft_icam_merged_median_OD_batch.txt" (OD spectra from blank filters), and
        "qft_icam_merged_sd_OD_batch.txt" (standard deviation of OD spectra from blank filters), respectively.
        2) Four plots in the directory "ODfiles" named as 
        "qft_icam_merged_median_OD_totparticle.png" (OD spectra from samples).
        "qft_icam_merged_sd_OD_totparticle.png" (standard deviation of OD spectra from samples).
        "qft_icam_merged_median_OD_NAP.png" (OD spectra from bleached samples), 
        "qft_icam_merged_sd_OD_NAP.png" (standard deviation of OD spectra from bleached samples), respectively.
        3) A file named as "qft_icam_merged_labels.txt" containing all the data 
        lables in the working directory. This can be used to construct the files
        "qft_icam_matched_labels_filtration_volumn.txt" and "labels_latlon_datetime_qft.txt". Note that the 
        sample and bleached sample labels are change to "sam_" or "BL_" + original labels.
        '''         
        #merge OD median data files
        print ('Merging OD data files...')
        fnames = sorted(glob.glob(os.path.join('ODfiles', '*OD_'+keyword+'.txt'))) 
        wl, line, wlpos = self.getWL(fnames[0])    
        ODf = pd.DataFrame()
        for f in fnames:
            #Read OD median data file
            print(f'Reading {f} \n')     
            try:
                data = pd.read_csv(f, header=None, comment='%', sep='\t')
                ODf = ODf.append(data)
                del data
            except:
                print(f'Warning: No data found in: {f}!!!')      
        
        #take the median of repeated measurements (if any)
        ODf = ODf.groupby(pd.Grouper(key=0)).median() 
        ODf.sort_index(inplace=True)
        ODf.insert(0,'label',ODf.index, True) 
        ODf.index = range(len(ODf))
        label = ODf.iloc[:,0]
        
        #batch OD
        ODf_batch = ODf[label.str.contains('batch|bt|MQ|SW', case=False)]
        #bleached sample OD
        ODf_bl = ODf[(label.str.contains('BL', case=False) == True) & 
                   (label.str.contains('batch|bt|MQ|SW|PK|PM', case=False) == False)]
        #sample OD
        ODf_sam = ODf[(label.str.contains('sam', case=False) == True) & 
                   (label.str.contains('batch|bt|MQ|SW|bl|PK|PM', case=False) == False)]
        #all except batch
        ODf_all = ODf[~ODf['label'].str.contains('batch|bt|MQ|SW', case=False)]
    
        outnames = ['qft_icam_merged_'+keyword+'_OD_batch.txt', 
                    'qft_icam_merged_'+keyword+'_OD_totparticle.txt',
                   'qft_icam_merged_'+keyword+'_OD_NAP.txt', 
                   'qft_icam_merged_'+keyword+'_OD_all.txt']
        varnames = ['ODf_batch', 'ODf_sam', 'ODf_bl', 'ODf_all']
    
        for i in range(len(outnames)):  
            varname = eval(varnames[i])
            self.saveData(os.path.join('ODfiles',outnames[i]), varname, 
                      headline=line, index=False)
    
        #output the label list        
        label = ODf_all.iloc[:,0]
        label.to_csv('qft_icam_merged_labels.txt', index=False, header=False, 
                     encoding='utf-8')
   
        if plot:
            OD_sam = ODf_sam.iloc[:,1:]
            OD_bl = ODf_bl.iloc[:,1:]
            wl_for_plot = wl.iloc[wlpos['pos400']:wlpos['pos750']+1]
            OD_sam_for_plot = OD_sam.iloc[:,wlpos['pos400']:wlpos['pos750']+1].transpose()
            OD_bl_for_plot = OD_bl.iloc[:,wlpos['pos400']:wlpos['pos750']+1].transpose() 
            #lgd = ODf.iloc[:,0].tolist()
            #plot sample OD
            figname = os.path.join('ODfiles',outnames[1].replace('txt','png'))
            self.plotSpectrum(wl_for_plot, OD_sam_for_plot, figname, ylabel='sample OD')
            
            #plot bleached sample OD
            figname = os.path.join('ODfiles',outnames[2].replace('txt','png'))
            self.plotSpectrum(wl_for_plot, OD_bl_for_plot, figname, ylabel='bleached sample OD')
        
        print ('Merging QFT-ICAM OD data files finished!')
               
    #--------------------------------------------------------------------------
    def calc_abs(self, filename, filename_sd, externalvolumn=None, 
                 lineheight=True, plot=True): 
        '''
        This function calculates absorption coefficients of either total particulate 
        matters or/and non-algal particles, specified by the input file.
        
        Input:
        1) filename - str. Path of OD data file (tab delimited).
        2) externalvolumn - str, default None. If None, filtration volumn will be 
        taken directly from the QFT-ICAM data file; if not, it specifies the path 
        of the external filtration volumn file (tab delimited), in which the first, 
        second and third columns are QFT-ICAM sample labels, matched bleached 
        sample labels and filtration volumn in liter, respectively. If one label is
        missing, e.g. "sam_PS107_UW281" exists but "BL_PS107_UW281_bl" is missing, 
        then the missing label is left blank in the file.
        3) lineheight - bool, default True. If True, calculate absorption line height 
        at 676 nm, vice versa. 
        4) plot - bool, default True. Plot absorption spectra.
        Output:
        1) Light absorption coefficient data and their standard deviation saved in 
        the working directory as "qft_icam_merged_median_abs_*.txt" and 
        "qft_icam_merged_sd_abs*.txt", respectively. The first 2 columns are sample 
        labels and filtration volumn.
        2) Light absorption coefficient data plotted in the working directory as 
        "qft_icam_merged_median_abs_*.png".
        3) Light absorption line height at 676 nm saved in the working directory as 
        "qft_icam_merged_median_lineheight676_*.txt".    
        '''
        wl, line, wlpos = self.getWL(filename)
        data = pd.read_csv(filename, header=None, comment='%', sep='\t')
        data_sd = pd.read_csv(filename_sd, header=None, comment='%', sep='\t')
        
        if externalvolumn is not None:
            data_vol = pd.read_csv(externalvolumn, sep='\t')
            data_vol.columns = range(data_vol.shape[1])
            data_vol_sam = data_vol.drop([1], axis=1).dropna()
            data_vol_bl = data_vol.drop([0], axis=1).dropna()
            data_vol = pd.DataFrame(np.vstack([data_vol_sam.values, data_vol_bl.values]))
            fil_vol = pd.Series(index=data.iloc[:,0])
            for i, lb in enumerate(fil_vol.index):
                pos = np.where(data_vol.iloc[:,0] == lb)[0]
                if len(pos) != 0:
                    fil_vol.iloc[i]  = float(data_vol.iloc[pos[0],1])
                else:
                    fil_vol.iloc[i]  = np.nan
            fil_vol = fil_vol[:,np.newaxis]   
        else:
            fil_vol = data.iloc[:,1][:,np.newaxis]
        
        area_in_centimetersquare = 0.25 * np.pi * self.diameter**2
        absorption = data.iloc[:,2:] * 2.303 * 0.1 * area_in_centimetersquare / self.beta 
        absorption = pd.DataFrame(absorption.values / fil_vol)
        absorption.index = data.iloc[:,0]
        absorption.insert(loc=0, column='filt_vol', value=fil_vol)
        
        absorption_sd = data_sd.iloc[:,2:] * 2.303 * 0.1 * area_in_centimetersquare / self.beta 
        absorption_sd = pd.DataFrame(absorption_sd.values / fil_vol)
        absorption_sd.index = data_sd.iloc[:,0]
        absorption_sd.insert(loc=0, column='filt_vol', value=fil_vol)
     
    #    #remove spectra with negative values
    #    for i in range(len(absorption)):
    #        if any(absorption.iloc[i,pos_400+1:pos_700-100] < 0):
    #            absorption.iloc[i,1:] = np.nan
    #            absorption_sd.iloc[i,1:] = np.nan

        outname = filename.split('/')[-1].replace('OD','abs') 
        outname_sd = filename_sd.split('/')[-1].replace('OD','abs') 
        self.saveData(outname, absorption, headline=line)
        self.saveData(outname_sd, absorption_sd, headline=line)

        print ('Absorption coefficients and their standard deviation calculated!')
        
        #line height
        if lineheight:
            aLH = absorption.iloc[:,wlpos['pos676']+1] - (absorption.iloc[:,wlpos['pos715']+1] - 
                                 absorption.iloc[:,wlpos['pos650']+1]) * (676-650)/(715-
                                                650) - absorption.iloc[:,wlpos['pos650']+1]
            aLH.to_csv(outname.replace('abs','lineheight676'), index=True, 
                       header=False, sep='\t', encoding='utf-8')
            print ('Absorption line height at 676nm calculated!')
            
        #plot
        if plot:
            wl_for_plot = wl.iloc[wlpos['pos400']:wlpos['pos750']+1]
            absorption_for_plot = absorption.iloc[:,wlpos['pos400']+1:wlpos['pos750']+2].transpose()
            absorption_sd_for_plot = absorption_sd.iloc[:,wlpos['pos400']+1:wlpos['pos750']+2].transpose()
            figname = outname.replace('txt','png')
            self.plotSpectrum(wl_for_plot, absorption_for_plot, figname, ylabel='Absorption [$m^{-1}$]')

            figname = outname_sd.replace('txt','png')
            self.plotSpectrum(wl_for_plot, absorption_sd_for_plot, figname, 
                              ylabel='Standard deviation of absorption [$m^{-1}$]')

    #--------------------------------------------------------------------------    
    def calc_ap_ad_aph(self, plot=True):    
        '''
        This function calculates absorption coefficients of total particulate 
        matters, non-algal particles and phytoplankton.
        
        Input:
        plot - bool, default True. Plot phytoplankton absorption spectra.
        Output:
        1) Absorption coefficients of total particulate matters and NAP, and 
        their standard deviations saved in the working directory as 
        "qft_icam_merged_median_abs_totparticle.txt",
        "qft_icam_merged_median_abs_NAP.txt",
        "qft_icam_merged_sd_abs_totparticle.txt", 
        "qft_icam_merged_sd_abs_NAP.txt", respectively. 
        The first 2 columns are sample labels and filtration volumn. 
        2) Adjusted NAP absorption coefficient data saved in the working directory 
        as "qft_icam_merged_median_abs_NAP_adjust.txt". "Adjusted" means bringing 
        NAP absorption back to ap in the NIR so that aph in NIR is 0 (offset the 
        average values in 720-750nm) (IOCCG Protocol Series, 2018).
        3) Adjusted NAP absorption coefficient data plotted in the working directory
        as "qft_icam_merged_median_abs_NAP_adjust.png".
        4) Phytoplankton absorption coefficient data calculated using adjusted NAP 
        absorption and saved in the working directory as 
        "qft_icam_merged_median_abs_phytoplankton.txt". 
        5) Phytoplankton absorption coefficient data plotted in the working 
        directory as "qft_icam_merged_median_abs_phytoplankton.png".
        6) Each ap, a_NAP and aph spectrum plotted in the directories of 
        "abs_totparticle", "abs_NAP", "abs_phytoplankton", respectively. 
        '''
        #Calculate ap and absorption line height at 676nm
        self.calc_abs(self.ODpFile, self.ODpFile_sd, externalvolumn=self.externalvolumn)
        
        #Calculate a_NAP
        self.calc_abs(self.ODdFile, self.ODdFile_sd, externalvolumn=self.externalvolumn, 
                      lineheight=False)
    
        #Calculate aph
        wl,line, wlpos = self.getWL(self.apFile)  
        ap = pd.read_csv(self.apFile, header=None, comment='%', sep='\t')
        anap = pd.read_csv(self.adFile, header=None, comment='%', sep='\t')  
        labels = pd.read_csv(self.externalvolumn, sep='\t').dropna()
        
        aph = pd.DataFrame(index=labels.iloc[:,0], columns=range(len(wl)))
        aph.insert(0,'filtration_volumn_in_liter',labels.values[:,2])
        anap_adjust = pd.DataFrame(index=labels.iloc[:,0].str.replace('sam','BL'),
                                   columns=range(len(wl)))
        anap_adjust.insert(0,'filtration_volumn_in_liter',labels.values[:,2])
        
        for i, lb_sam in enumerate(labels.iloc[:,0]):
            pos_sam = np.where(ap.iloc[:,0]==lb_sam)[0]
            pos_bl = np.where(anap.iloc[:,0]==labels.iloc[i,1])[0]       
            if len(pos_sam)!=0 and len(pos_bl)!=0:
                tmp_ap_NIR = np.nanmedian(ap.iloc[pos_sam,wlpos['pos710']+2:wlpos['pos750']+2])
                tmp_anap_NIR = np.nanmedian(anap.iloc[pos_bl,wlpos['pos710']+2:wlpos['pos750']+2])
                tmp_offset = tmp_anap_NIR - tmp_ap_NIR
                anap_adjust.iloc[i,1:] = anap.values[pos_bl,2:] - tmp_offset 
                aph.iloc[i,1:] = ap.values[pos_sam,2:] - anap_adjust.values[i,1:] 
      
        #remove spectra with negative values         
        for i in range(len(aph)):
            if any(anap_adjust.iloc[i,wlpos['pos400']+1:wlpos['pos700']-5] < 0):
                anap_adjust.iloc[i,:] = np.nan
                aph.iloc[i,:] = np.nan
            if any(aph.iloc[i,wlpos['pos400']+1:wlpos['pos700']-5] < 0):
                aph.iloc[i,:] = np.nan
        
        outname_aph = self.adFile.replace('NAP','phytoplankton') 
        self.saveData(outname_aph, aph, headline=line)
        print ('Phytoplankton absorption coefficients calculated!')
            
        outname_anap_adjust = self.adFile.replace('NAP','NAP_adjust') 
        self.saveData(outname_anap_adjust, anap_adjust, headline=line)
        print ('NAP absorption adjusted!')
            
        #plot
        if plot:
            wl_for_plot = wl.iloc[wlpos['pos400']:wlpos['pos750']+1]
            aph_for_plot = aph.iloc[:,wlpos['pos400']+1:wlpos['pos750']+2].transpose()
            anap_for_plot = anap_adjust.iloc[:,wlpos['pos400']+1:wlpos['pos750']+2].transpose()
            
            self.plotSpectrum(wl_for_plot, aph_for_plot, outname_aph.replace('txt','png'), 
                                      ylabel='$a_{ph}(\lambda)$ [$m^{-1}$]')
            self.plotSpectrum(wl_for_plot, anap_for_plot, outname_anap_adjust.replace('txt','png'), 
                                      ylabel='$a_{NAP}(\lambda)$ [$m^{-1}$]')
            
            #Create target directory for individual spectral plots
            dir_fig_ap = 'abs_' + os.path.basename(self.apFile).split('_')[-1][:-4]
            dir_fig_anap = 'abs_' + os.path.basename(self.adFile).split('_')[-1][:-4]
            dir_fig_aph = dir_fig_anap.replace('NAP','phytoplankton')
            dir_fig = [dir_fig_ap, dir_fig_anap, dir_fig_aph]
            ap.index = ap.iloc[:,0]
            ap = ap.drop([0, 1], axis=1)
            
            for i, dirf in enumerate(dir_fig):
                self.createDir(dirf, overwrite=True)
                #plot abs spectra
                var_fig = ['ap', 'anap_adjust', 'aph']
                print('Plotting each '+ var_fig[i]+ ' spectrum...')
                var = eval(var_fig[i])
                for j in range(len(var)):           
                    data_for_plot = var.iloc[j,wlpos['pos400']+1:wlpos['pos750']+2].transpose()
                    namestr = os.path.join(dirf, var.index[j] + '.png')
                    self.plotSpectrum(wl_for_plot, data_for_plot, namestr, 
                                      ylabel='Absorption [$m^{-1}$]')
    
    #--------------------------------------------------------------------------   
    def preparePangaea(self, filename, NAP=False, adjustNAP=False):
        '''
        This function restructs QFT-ICAM data for uploading to Pangae.
        Input:
        1) filename - str. Path of QFT-ICAM absorption, OD or standard deviation data file (tab delimited).
        2) NAP - bool, default False. If True, input variable "filename" contains NAP data.
        3) adjustNAP - bool, default False. If True, input variable "filename" contains 
        adjusted NAP data (details on the adjustment see self-defined function "qft_calc_aph").   
        Output:
        QFT-ICAM data suitable to be uploaded to Pangaea saved as:
        "*_pangaea.txt", with "*" standing for the base name of input "filename".
        '''
        foo = filename.split('/')
        print(f'Importing {foo[-1]}...')
        
        wl, line, wlpos = self.getWL(filename)
        data = pd.read_csv(filename, header=None, comment='%', sep='\t')
        info = pd.read_csv(self.concurrentInfo, comment='%', sep='\t')
        
        if NAP:
            if adjustNAP:
                label_info = info.iloc[:,0].str.replace('sam','BL')
            else:
                label_info = info.iloc[:,1]
        else:
            label_info = info.iloc[:,0]
        
        pos = pd.Series(np.nan, index=data.index)
        for i,lb in enumerate(data.iloc[:,0]):
            tmp = np.where(label_info==lb)[0]
            if len(tmp) > 0:
                pos.iloc[i] = tmp[0]

        data_pangaea = info.iloc[pos,2:]
        data_pangaea['Date/Time'] = pd.to_datetime(data_pangaea['Date'] + ' ' + 
                                  data_pangaea['Time'])
        data_pangaea = data_pangaea.drop(['Date','Time'], axis=1)
        
        wl_for_data = wl.iloc[wlpos['pos350']:wlpos['pos800']+1]
        data_col = ['wl' + str(wvl) for wvl in wl_for_data]
        new_data = pd.DataFrame(data=data.values[:,2+wlpos['pos350']:wlpos['pos800']+3], 
                                columns=data_col, index=data_pangaea.index)
        data_pangaea[data_col] = new_data
        
        #sort data by 'Date/Time'
        data_pangaea = data_pangaea.sort_values(by='Date/Time')
        
        #sort CTD data by Depth
        label_data_uni = np.unique(data_pangaea['Sample_ID'].tolist(),
                                   return_index=True, return_inverse=True, 
                                   return_counts=True)
        pos_ctd = np.where(label_data_uni[3] > 1)[0]
        for i in pos_ctd:
            pos = np.where(data_pangaea['Sample_ID']==label_data_uni[0][i])[0]
            tmp = data_pangaea.iloc[pos,:].sort_values(
                    by=['Depth_in_meter'])
            tmp.index = data_pangaea.iloc[pos,:].index
            data_pangaea.iloc[pos,:] = tmp
        
        new_col = ['Event', 'Date/Time', 'Sample_ID', 'Latitude', 'Longitude', 
                   'Depth_in_meter', 'filtration_volumn_in_liter']
        data_pangaea = data_pangaea[new_col + [col for col in data_pangaea
                                                       if col not in new_col]]
        
        data_pangaea.to_csv(filename.replace('.txt','_pangaea.txt'), index=False, 
                           header=True, sep='\t', encoding='utf-8')
        print (f'Data from {foo[-1]} ready for Pangaea!')
        
    #--------------------------------------------------------------------------   
    def uploadPangaea(self):
        '''
        Output:
        1) Five tab delimited absorption (and standard deviation) files in the 
            working directory named as
            "qft_icam_merged_median_abs_totparticle_pangaea.txt",
            "qft_icam_merged_sd_abs_totparticle_pangaea.txt",
            "qft_icam_merged_median_abs_NAP_adjust_pangaea.txt",
            "qft_icam_merged_sd_abs_NAP_pangaea.txt",
            "qft_icam_merged_median_abs_phytoplankton_pangaea.txt".
        2) Four tab delimited OD (and standard deviation) files in the directory 
            "ODfiles" named as
            "qft_icam_merged_median_OD_totparticle_pangaea.txt",
            "qft_icam_merged_sd_OD_totparticle_pangaea.txt",
            "qft_icam_merged_median_OD_NAP_pangaea.txt",
            "qft_icam_merged_sd_OD_NAP_pangaea.txt".
        '''
        self.preparePangaea(self.apFile)
        self.preparePangaea(self.apFile_sd)
        self.preparePangaea(self.adFile_adjust, NAP=True, adjustNAP=True)
        self.preparePangaea(self.adFile_sd, NAP=True)
        self.preparePangaea(self.aphFile)
        self.preparePangaea(self.ODpFile)
        self.preparePangaea(self.ODpFile_sd)
        self.preparePangaea(self.ODdFile, NAP=True)
        self.preparePangaea(self.ODdFile_sd, NAP=True)              
    
#-----------------------------------------------------------------------------
# Main Program. 
#-----------------------------------------------------------------------------
if __name__ == '__main__': 

    #Set working directory as the one storing .qft.raw data files.
    wd = '/work/ollie/yliu/Data/cruises/qftlwcc_postproc/qft_icam'
    os.chdir(wd)

    #Prepare your own configuration file 'config_PostProc.txt' in the working directory.
    process = qftPostProc()
    
    
    #Step1:--------------------------------------------------------------------
    #process raw data.
    filenames_raw = sorted(glob.glob(os.path.join(wd,'*.qft.raw')))
    for filename in filenames_raw:
        process.processRaw(filename)


    #Step2:--------------------------------------------------------------------
    #plot QFT-ICAM optical density data from each sample file by file.  
    filenames_l0a = sorted(glob.glob(os.path.join(wd,'*.l0a')))
    for filename in filenames_l0a:
        process.plotOD_l0a(filename)
    
    
    #Step3:--------------------------------------------------------------------
    #1)remove spectra
    #Prepare your own 'rmSpectraIndex.txt' file (tab delimited) based on the 
    #plots generated from Step2 in the working directory!
    SpectraIndex = pd.read_csv('rmSpectraIndex.txt', comment='#', sep='\t')
    for i, filename in enumerate(SpectraIndex.iloc[:,0]):
        row_index = str(SpectraIndex.iloc[i,1])
        row_index = [float(idx) for idx in row_index.split(',')]
        process.rmOD_median(filename, row_index)
    
    #2)merge OD and standard deviation data files and get all data lables.
    process.mergeOD()
    process.mergeOD(keyword='sd')


    #Step4:-------------------------------------------------------------------- 
    #calculate absorption coefficients ap, ad, aph 
    #Prepare your own 'qft_icam_matched_labels_filtration_volumn.txt' file (tab 
    #delimited) based on the file 'qft_icam_merged_labels.txt' generated from 
    #Step3(2) in the working directory!
    process.calc_ap_ad_aph(plot=True)

    
    #Step5:--------------------------------------------------------------------
    #prepare data for upload to Pangaea
    #Prepare your own 'labels_latlon_datetime_qft.txt' file (tab delimited) in 
    #the working directory!
    process.uploadPangaea()
