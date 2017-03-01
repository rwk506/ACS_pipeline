# -*- coding: utf-8 -*-

from pylab import *
from numpy import *
from random import *
from scipy import *
from scipy import optimize
from matplotlib import *
import csv
import sys
import re,pyfits,math
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from uncertainties import *
import uncertainties.unumpy as unumpy
import astLib
import astLib.astStats
import subprocess
import glob

from plotresY import plotresY
from plotres_noY import plotres_noY
import random

rcdefaults()
matplotlib.rc('font',family='Bitstream Vera Serif')


def reformat_ACSfile(ids,f606,f814,f606err,f814err,outputfile):
    with open(outputfile, 'w') as f:
        f.write( 'Star F606W F814W sigF606W sigF814W mass1 massRatio stage Cmprior useDBI' + '\n')
        for i in arange(0,len(f606)):
            f.write( makestr([int(ids[i]), round(f606[i],3),round(f814[i],3),round(f606err[i],3),round(f814err[i],3),'0','0','1',0.95,'1']) + '\n')
    return

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
    return

def makestr(mags):
    mags=np.array(mags)
    string=str(mags[0])
    for i in arange(1,len(mags)-1):
        string=string+' '+str(mags[i])
    string=string+' '+str(mags[-1])
    return string



##########################################################################
######################### BEGIN CLASS DEFINITION #########################
##########################################################################


########################### Initialize Pipeline ##########################

class ACSpipeline:
    '''ACSpipeline(clusterval, basefolder) - initiates processing/analysis pipeline for HST/ACS photometry of Galactic globular clusters.
        
        Some notes on directory structure:
        Top folder: 'self.basefolder/'  
            -->  contains folders for each cluster
            -->  contains 'base9_+self.clusterval+.yaml
            -->  contains other BASE-9 functionalities
        Next level: 'self.basefolder/NGC'+self.clusterval+'/'
            -->  contains original HST/ACS photometry files from survey publications (e.g.: NGC0104R.RDVIQ.cal.adj.zpt)
            -->  will contain photometry files, BASE-9 results files, generated plots etc.
            
        Additional notes/suggestions:
        - Either place this python code in the parent directory to 'self.basefolder'; else, give full pathname when initiating pipeline - realize that changes to or deviations from directory structure may lead to errors!
        - When not using them, keep base9.yaml files named as base9_clusterval.yaml; temporarily rename to base9.yaml to initiate BASE-9 calculations, then return to labeled name
        
        More information about the HST/ACS survey can be found at: http://www.astro.ufl.edu/~ata/public_hstgc/
        BASE-9 code and downloads are available here: https://github.com/argiopetech/BASE
        BASE-9 manuals are available here: https://arxiv.org/abs/1411.3786
        and here: http://base-9.readthedocs.io/en/latest/
        '''
    
    clusterval=str(); basefolder=str();med_values=[]

    def __init__(self,clustername,basefoldername):
        if type(self.clusterval) != str or type(self.basefolder) != str:
            raise ValueError('Inputs need to be in string form.')
        self.basefolder=basefoldername
        self.clusterval=clustername

        ###### load in magnitude file from GP
        self.id, self.x, self.y, self.Vvega, self.Verr, self.VIvega, self.VIerr, self.Ivega, self.Ierr, self.Vground, self.Iground, self.Nv,  self.Ni, self.wV, self.wI, self.DX, self.DY, self.othv, self.othi, self.qfitV, self.qfitI, self.RA, self.Dec = loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'R.RDVIQ.cal.adj.zpt',skiprows=3,unpack=True)
    
        #ALLred[where(ALLnames==float(clusterval1))
        

########################### Look at HB to make cuts ###########################
###  Use this function to visually examine HB area of the CMD
###  (this will help determine where to make cuts)
###
    def checkHB(self,showplt='Yes'):
        '''Use to visually examine CMD and horizontal branch in particular (useful in conjunction with make_cuts() to clean up/remove HB stars).'''
        self.F606Wmag = numpy.array(self.Vvega)
        self.F814Wmag = numpy.array(self.Ivega)
    
        ###### plot to see where to cut HB
        ind=where((self.Ni>=1.) & (self.Nv>=1.))
        
        if showplt=='Yes' or showplt=='yes' or showplt=='y':
            plot(self.F606Wmag[ind]-self.F814Wmag[ind],self.F606Wmag[ind],'k.',ms=1,alpha=0.5)
            xlabel('F606W-F814W');ylabel('F606W')
            plt.gca().invert_yaxis()
            plt.show()
        return


########################### Make HB and proper motion cuts ###########################
### This function is used to make quality cuts to photometry (based on location and errors)
### and to remove HB stars from region (identified from checkHB() function)
### hb1mag and hb1col - make one cut to remove HB stars brighter and bluer than color/mag combo
### hb2mag and hb2col - make another cut to remove HB stars brighter and bluer than color/mag combo
###
    def make_cuts(self,hb1mag,hb1col,hb2mag,hb2col,showplt='Yes',strict='No',TOP=0.):
        ''' Function to remove poorly determined stars from photometry, as well as HB
            Stars with high in the 2.5% tails of the delta-pixel changes from frame to frame are removed (DX, DY)
            Stars with high errors in the 5% tail are removed (Verr, Ierr)
            
            inputs: hb1mag - faint magnitude limit 2 (brighter stars removed)
                    hb1col - color limit 2 (bluer stars removed)
                    hb2mag - faint magnitude limit 2 (brighter stars removed)
                    hb2col - color limit 2 (bluer stars removed)
                    showplt - toggle plot on/off
                    strict - toggle strict rules on for stars with one short exposure (NGC's 5986, 6397, 6779)
                    TOP - turn off point for clusters - used only when strict = 'Yes' to pinpoint bright, short exp. stars
            hbNmag and hbNcol give magnitude/color brighter/bluer limits for removing HB stars'''
        
        ###### quality cuts to data
        self.dx=(abs(percentile(self.DX,2.5))+abs(percentile(self.DX,97.5)))/2.
        self.dy=(abs(percentile(self.DY,2.5))+abs(percentile(self.DY,97.5)))/2.
        self.verr=abs(percentile(self.Verr,95))
        self.ierr=abs(percentile(self.Ierr,95))
        self.indx1=where((self.Ni>=1.) & (self.Nv>=1.) & (abs(self.DX)<=self.dx) & (abs(self.DY)<=self.dy) & (abs(self.Ierr)<=self.ierr) & (abs(self.Verr)<=self.verr))  ### stars to keep in "cleaned" dataset
        
        ###### for clusters with only one short exposure (NGC's 5986, 6397, 6779),
        ###### need to relax proper motion cuts for the brighter stars
        if strict=='Yes' or strict=='yes' or strict=='Y' or strict=='y':
            indxA=where((self.F606Wmag>TOP) & (self.Ni>=1.) & (self.Nv>=1.) & (abs(self.DX)<=self.dx) & (abs(self.DY)<=self.dy) & (abs(self.Ierr)<=self.ierr) & (abs(self.Verr)<=self.verr))
            indxB=where((self.F606Wmag<=TOP) & (self.Ni>=1.) & (self.Nv>=1.) & (abs(self.Ierr)<=self.ierr) & (abs(self.Verr)<=self.verr))
            self.indx1=[concatenate((indxA[0],indxB[0]))]  ### combine "cleaned" short and long exposure stars
    
        ###### HB stars -- remove from sample
        indx2=where((self.F606Wmag-self.F814Wmag<=hb1col) & (self.F606Wmag<=hb1mag))
        indx3=where((self.F606Wmag<=hb2mag)&(self.F606Wmag-self.F814Wmag<=hb2col))
        indx4=np.union1d(np.ravel(indx2),np.ravel(indx3));indx4=[indx4]

        ###### combine sets of indices and plot resulting cleaned CMD
        self.indx = [x for x in self.indx1[0] if x not in indx4[0]]
        if showplt=='Yes' or showplt=='yes' or showplt=='y':
            plot((self.F606Wmag[self.indx]-self.F814Wmag[self.indx]),self.F606Wmag[self.indx],'k.',ms=1,alpha=0.5)
            xlabel('F606W-F814W');ylabel('F606W')
            plt.gca().invert_yaxis()
            plt.show()
        return


########################### Reformat files ###########################
#### reformats photometry in BASE-9 form for ease of use
    def reformat_phot(self):
        '''This function reformats the photometry into "full" and "cleaned" photometry files formatted for use with the BASE-9 program.'''
        
        ### reformat all stars into "Cluster.all.ACS.phot" in BASE-9 input format
        reformat_ACSfile(arange(0,len(self.F606Wmag)),self.F606Wmag,self.F814Wmag,self.Verr,self.Ierr,self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.allACS.phot')
        
        ### reformat "cleaned" photometry into "Cluster.cleaned.ACS.phot" in BASE-9 input format
        reformat_ACSfile(self.indx,self.F606Wmag[self.indx],self.F814Wmag[self.indx],self.Verr[self.indx],self.Ierr[self.indx],self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.cleanedACS.phot')
        
        return 'Complete'
    
    
########################### Take subsample ###########################
#### samples from cleaned photometry to make a smaller, randomly selected sample of stars to run with BASE-9
    def subsample_phot(self,magcut606,magcutTOP=24,seedval=1823,n=3000):
        '''Randomly samples n stars from cleaned photometry and writes a photometry file for the randomly sampled stars. 
            
            magcut606 - n/2 stars are randomly sampled from the cluster stars fainter than the MSTOP and n/2 stars are randomly sampled from the stars in the cluster brighter than the MSTOP
            magcutTOP - a faint limit imposed to avoid the lower main sequence dwarfs (which are not well modeled and have larger photometric uncertainties); default = 24, however, in practice, generally distance modulus + 7 magnitudes is used for analysis purposes
            seedval - the seed value for the randomized selection; default seedval = 1823
            n - number of stars to be randomly chosen, with n/2 stars brighter than the MSTOP and n/2 stars fainter than the MSTOP; default n = 3000'''
        import random
        seed(seedval)
        ### load in "cleaned" photometry to use for subsample
        ids_clean,m606,m814,e606,e814 = loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.cleanedACS.phot',skiprows=1,usecols=(0,1,2,3,4),unpack=True)
        
        ### open new file to write
        with open(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.sampleACS.phot', 'w') as f:
            f.write('Star F606W F814W sigF606W sigF814W mass1 massRatio stage Cmprior useDBI'+'\n')
            #### upper CMD, above TOP, select n/2 unique stars
            indx_upper=where(m606<=magcut606)
            if len(indx_upper[0]) <= floor(n/2.):
                for i in indx_upper[0]:
                    f.write(makestr([ids_clean[i], m606[i], m814[i], e606[i], e814[i], '0', '0', '1', .95, '1'])+'\n')
                samp_upper=indx_upper[0]
            else:
                samp_upper=random.sample(indx_upper[0],int(floor(n/2.)))
                for i in samp_upper:
                    f.write(makestr([ids_clean[i], m606[i], m814[i], e606[i], e814[i], '0','0','1',.95,'1'])+'\n')

            #### lower MS, below TOP, select n/2 fainter than MSTOP but brighter than faint limit cutoff (distmod+9)
            indx_lower=where((m606>=magcut606)&(m606<=magcutTOP))
            samp_lower=random.sample(indx_lower[0],len(samp_upper))#int(floor(n/2.)))
            
            #### write new file
            for i in samp_lower:
                f.write(makestr([ids_clean[i], m606[i], m814[i], e606[i], e814[i], '0','0','1',.95,'1'])+'\n')
        print 'N* above MSTOP ', len(samp_upper)
        print 'N* below MSTOP ', len(samp_lower)
        print 'N* total sample', len(samp_upper)+len(samp_lower)
        return



########################### Check sample ###########################
#### check process to make sure sample looks good/reasonable (sanity check!)
    def check_sample(self,titleon='Yes'):
        '''This functions allows the user to double-check the results of the pre-processing of the ACS photometry, with cleaned photometry in black, and the randomly selected sample of n stars in green.'''
        ids_clean,m606,m814 = loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.cleanedACS.phot',skiprows=1,usecols=(0,1,2),unpack=True)
        plot(m606-m814,m606,'k.',ms=1,alpha=1)
        ids_clean,m606,m814 = loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.sampleACS.phot',skiprows=1,usecols=(0,1,2),unpack=True)
        plot(m606-m814,m606,'g.',ms=3,alpha=1)
        gca().invert_yaxis()
        xlabel('F606W-F814W');ylabel('F606W')
        if titleon=='Yes' or titleon=='yes' or titleon=='y':
            title('NGC '+str(self.clusterval))
        plt.show()
        
        return





#####################################################################
########################### POST-MCMC RUN ###########################
#####################################################################


########################### Plot MCMC sampling ###########################
#### following the BASE-9 calculations, visualize Adapative Metropolis MCMC chain and the resulting posteriors
####
    def plot_MCMCresults(self,vnum,meansline=27,saveYN='Yes',plotpriors='No',printtable='No',startn=1001,noY='False'):
        '''Following the BASE-9 calculations, this function enables the user to examine the MCMC chains for each parameter and the resulting posterior distributions. An option also will provide a printed LaTeX table of the results.
            
            vnum - file number of the results file to plot (e.g.: for a BASE-9 results file NGC0104.v2.res, vnum='2')
            meansline - default value is 27 - DO NOT CHANGE unless the .yaml file format has been changed in recent BASE-9 versions
            plotpriors - default is no, but allows the user to plot horizontal lines indicating the prior distribution over the MCMC chain and over the posterior distributions
            printtable - default is no, user may elect to print the results in LaTeX format
            startn - the line to start in the MCMC results file; the default is 1001 to skip the burnin, however it may be desirable or necessary in certain cases to change this
            noY - if Y is included as a parameter, leave as False. If helium is set to a specific value during the calculations, then toggle to "True".'''
        
        ##### load in results from BASE-9 calculations
        res=loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.v'+str(vnum)+'.res',skiprows=startn)
        
        ##### get priors used in BASE-9 from base9.yaml file to plot, if desired
        if noY=='False':
            if plotpriors=='Yes' or plotpriors=='yes' or plotpriors=='y':
                priors=[]
                base9 = open(self.basefolder+'/base9_'+self.clusterval+'.yaml','r')
                for i,line in enumerate(base9.read().split('\n')):
                    if i in arange(meansline+1,meansline+7):  ### limit intake to range where priors are, take values as floats, and store
                        priors.append(float(line[28:]))
                base9.close()
                priors.insert(0, priors.pop())  ### moves ages from last (in yaml file) to first (for plotting tool)
                priors=priors[:4]  ### use only age, metallicity, distance, and Av
                priors.append(0)  ### for Y prior; unknown, really
                plotresY(res,priors=priors)   ### plot using outside function
            else:
                plotresY(res)
            
            ##### print table in LaTeX format, if desired
            if printtable=='Yes' or printtable=='yes' or printtable=='y':
                labels_mult=['log(Age) ','Y$_{A}$  ','$[Fe/H]$ ','$\mu_{V}$','A$_{V}$  ']
                psigmas=[' ','$\pm$0.01',' -- ','$\pm$0.02','$\pm$0.003']
                #age,feh,distmod,av, ya, yb,lamd
                print ''
                for i in [0,2,3,4,1]:
                    ### second column should be published values; third column should be prior used in .yaml file
                    print labels_mult[i], ' & ', str(priors[i])+psigmas[i] ,' & ', str(round(median(res[:,i]),3))+'$^{+'+str(round(percentile((res[:,i]-median(res[:,i])),95),3))+'}_{'+str(round(percentile((res[:,i]-median(res[:,i])),5),3))+'}$'

        ###### if there is no Y sampling.... essentially repeat above but under scenario of a fixed Y value
        if noY=='True':
            if plotpriors=='Yes' or plotpriors=='yes' or plotpriors=='y':
                priors=[]
                base9 = open(self.basefolder+'/base9_'+self.clusterval+'.yaml','r')
                for i,line in enumerate(base9.read().split('\n')):
                    if i in arange(meansline+1,meansline+7):  ### limit intake to range where priors are, take values as floats, and store
                        priors.append(float(line[28:]))
                base9.close()
                priors.insert(0, priors.pop())  ### moves ages from last (in yaml file) to first (for plotting tool)
                priors=priors[:4]  ### use only age, metallicity, distance, and Av
                plotres_noY(res,priors=priors)  ### plot using outside function
            else:
                plotres_noY(res)
            
            ##### print table in LaTeX format, if desired
            if printtable=='Yes' or printtable=='yes' or printtable=='y':
                labels_mult=['log(Age) ','$[Fe/H]$ ','$\mu_{V}$','A$_{V}$  ']
                psigmas=[' ','$\pm$0.05','$\pm$0.05','$\pm$0.02']
                #age,feh,distmod,av, ya, yb,lamd
                print ''
                for i in [0,1,2,3]:
                    ### second column should be published values; third column should be prior used in .yaml file
                    print labels_mult[i], ' & ', str(priors[i])+psigmas[i] ,' & ', str(round(median(res[:,i]),3))+'$^{+'+str(round(percentile((res[:,i]-median(res[:,i])),95),3))+'}_{'+str(round(percentile((res[:,i]-median(res[:,i])),5),3))+'}$'

        ##### save if desired
        if saveYN=='Yes' or saveYN=='yes' or saveYN=='y':
            savefig(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.v'+str(vnum)+'.samplingACS.png', bbox_inches='tight')

        plt.show()
        return



########################### Re-write .yaml and run makeCMD ###########################
#### use results from BASE-9 to re-write base9.yaml file to prep running ./makeCMD with BASE-9, which generates isochrone
####
    def run_makeCMD(self,vnum,startingline=18,meansline=27,saveYN='Yes',startn=1001, noY='False', Yval=0.):
        '''A function designed to prep a base9_cluster.yaml file from the results file with median values from the posterior distributions. The base9.yaml file that results from this can be used to generate an isochrone with makeCMD.
            
            vnum - file number of the results file to plot (e.g.: for a BASE-9 results file NGC0104.v2.res, vnum='2')
            startingline - default value is 18 - DO NOT CHANGE unless the .yaml file format has been changed in current BASE-9 versions
            meansline - default value is 27 - DO NOT CHANGE unless the .yaml file format has been changed in current BASE-9 versions
            startn - the line to start in the MCMC results file; the default is 1001 to skip the burnin, however it may be desirable or necessary in certain cases to change this
            noY - if Y is included as a parameter, leave as False. If helium is set to a specific value during the calculations, then toggle to "True"
            Yval - if noY is set to true, specify the Y value that has been used'''
        
        self.med_values=[]  #### set up array for medians
        res=loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.v'+str(vnum)+'.res',skiprows=startn)  #### results file
        labels=['                Fe_H:        ','                distMod:     ','                Av:          ','                Y:           ','                logAge:      ']  #### prefixes for replacing file lines in base9.yaml
        
        if noY=='False':
            ### get median values
            #med_values=[]
            for i in [2,3,4,1,0]:  ## metallicity, distance, Av, YA, age
                self.med_values.append(median(res[:,i]))
    
        if noY=='True':
            for i in [1, 2, 3, 10, 0]:
                if i!=10:
                    self.med_values.append(median(res[:,i]))
                if i==10:
                    self.med_values.append(Yval)
        
        print self.med_values  ### print median posterior distribution values obtained from results file
        

        ### rename base9_cluster.yaml to base9.yaml
        command='cp '+self.basefolder+'/base9_'+self.clusterval+'.yaml '+self.basefolder+'/base9.yaml'
        os.system(str(command))

        ### take base9.yaml, insert median values for priors and starting values
        base9 = open(self.basefolder+'/base9.yaml','r')
        base9path = self.basefolder+'/base9.yaml'
        n=0;m=0
        for i,line in enumerate(base9.read()):
            if i==5:
                replace_line(base9path,i,'        outputFileBase: "NGC'+self.clusterval+'/NGC'+self.clusterval+'.v'+vnum+'"'+'\n')
            if i in [1,2,3,4,6]:  ### limit intake to range where priors are, take values as floats, and store
                linenum=startingline+i
                replace_line(base9path,linenum,labels[m][4:]+str(self.med_values[n])+'\n')
                m+=1
            if i in [1,2,3,4,6]:  ### limit intake to range where priors are, take values as floats, and store
                linenum=meansline+i
                replace_line(base9path,linenum,labels[n]+str(self.med_values[n])+'\n')
                n+=1
        base9.close()
        return



########################### Plot results on CMD ###########################
#### Plot cleaned and sample photometry with resulting isochrone from BASE-9 analysis
####
    def plot_CMDresults(self,vnum,saveYN='Yes',titleon='No'):
        '''A function to plot the photometry from HST/ACS with the resulting isochrone from the BASE-9 results.
            vnum - file number of the results file to plot (e.g.: for a BASE-9 results file NGC0104.v2.res, vnum='2')'''
        
        #### load model from BASE-9 results (generated using run_makeCMD() and makeCMD in BASE-9)
        model606,model814=loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.v'+str(vnum)+'.ms',skiprows=1,usecols=(14,17),unpack=True)
        
        #### plot photometry and results
        ids_clean,m606,m814 = loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.cleanedACS.phot',skiprows=1,usecols=(0,1,2),unpack=True)
        plot(m606-m814,m606,'.',ms=1.5,alpha=.6,color='gray')
        ids_clean,m606,m814 = loadtxt(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.sampleACS.phot',skiprows=1,usecols=(0,1,2),unpack=True)
        plot(m606-m814,m606,'k.',ms=3,alpha=.6)
        plot(model606-model814,model606,'c-')
        
        gca().invert_yaxis()
        xlabel('F606W-F814W');ylabel('F606W')

        if titleon=='Yes' or titleon=='yes' or titleon=='y':
            title('NGC '+str(self.clusterval))
        
        #### save
        if saveYN=='Yes' or 'yes' or 'y':
            savefig(self.basefolder+'/NGC'+self.clusterval+'/NGC'+self.clusterval+'.v'+str(vnum)+'ACS.png', bbox_inches='tight')
        plt.show()
        
        return




