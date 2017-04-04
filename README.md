# ACS_pipeline
A pre-processing and post-processing pipeline for formatting and analyzing Hubble Space Telescope Galactic Globular Cluster Treasury photometry with BASE-9

<br />

**Contents**

[Introduction](#Intro)<br />
[Summary of Data and Code](#Summary)<br />
[Directory Structure](#Structure)<br />
[Repository Contents](#RepoStuff)<br />
[ACS Pipeline Functions](#Functions)<br />
[Pre-Processing](#Process) <br />
[Post-Computation](#Compute)<br />
[Test Case Example](#Code)<br />
[Additional Information](#Other)<br />


<br />
<a name="Intro"/>
<h4>Introduction</h4>

The package presented here is a data processing and analysis pipeline to be used in conjunction with and complement to the [BASE-9 software](https://github.com/argiopetech/BASE). Bayesian Analysis for Stellar Evolution with Nine Variables (BASE-9) is a software suite that recovers star cluster and stellar parameters from photometry via hierarchical Bayesian modeling with physics-based theoretical models. BASE-9 uses a Markov chain Monte Carlo (MCMC) technique along with brute force numerical integration to estimate the posterior probability distribution for parameters of interest to the user.

The ACS_pipeline repository is designed to complement the BASE-9 software to ease the use and analysis of the Bayesian results. The pipeline was created for photometry derived from Hubble Space Telescope (HST) Advanced Camera for Survey (ACS) images of Galactic globular clusters -- specifically, the HST ACS Galactic Globular Cluster Treasury Program (GO Cycle 14 Proposal 10775; [Sarajedini et al. 2007](https://arxiv.org/abs/astro-ph/0612598)). The tabulated photometry files from this program, derived via the methods described in Sarajedini et al. 2007, are available for download [here](http://www.astro.ufl.edu/~ata/public_hstgc/), under "Databases". Additional information on this HST Treasury Program is available at that page, as well as in a wealth of academic publications.

The code presented in this repository was used in a 2017 publication to analyze 69 globular clusters (available [here](https://arxiv.org/abs/1702.08856) on arXiv). A similar pipeline was used for [these](https://arxiv.org/abs/1604.06074) [two](https://arxiv.org/abs/1609.01527) publications for inclusion of additional filters.


<br/>
<a name="Summary"/>
<h4>Summary of Data and Code</h4>

The excerpt below from [WK et al. 2017](https://arxiv.org/abs/1702.08856) summarizes the primary goals of the pre-processing performed via the ACS pipeline package:

>"The ACS Globular Cluster Treasury program has provided observations of several hundred thousand stars [in two HST filters, F606W and F814W]. To help BASE-9 perform effectively, we make modest quality cuts on the photometry for the clusters before randomly selecting a subsample of stars. To rid the data of any stars with poorly determined photometry, we remove stars for which both filters fall into the outer 5% tail of the photometric error distributions. 

>Additionally, we remove stars in the outer 2.5% tails of the distributions of X and Y pixel location errors from frame to frame, as high pixel location errors may indicate non-cluster members. The exceptions to this are the clusters NGC 5986, NGC 6397, and NGC 6779. For these clusters, we ignore the location error constraint as it removes the majority of bright stars above the main sequence turn-off that were observed in the short exposure images. However, we still remove stars in the outer 5% tail of the photometric error distributions for both filters. 

>With the cleaned photometry, we randomly select a subsample of ~ 3000 stars, with half above the main sequence turnoff point (MSTOP) of the cluster and half below the MSTOP. If there are fewer than 1500 stars above the MSTOP, we match the number of stars above and below the MSTOP. This procedure is adopted to ensure a reasonable sample of stars on the sub-giant and red-giant branches, while balancing their contribution with MS stars, without the computational cost of running 10,000 or more stars per cluster."


The ACS pipeline here also formats the cleaned and sample photometry to output files such that they can be easily used with BASE-9 by putting the filepath into the BASE-9 .yaml file.

Additionally, the ACS pipeline code herein allows a post-BASE-9 analysis, visualizing the results of the BASE-9 computations and calculating the resulting ("median") isochrone with BASE-9.

This package is *not* meant to include the actual BASE-9 computations or code; these *must* be acquired and installed separately from the [BASE-9 repository](https://github.com/argiopetech/BASE).

Manuals for the use and implementation of BASE-9 are available [here](https://arxiv.org/abs/1411.3786) and [here](http://base-9.readthedocs.io/en/latest/). We assume the user is familiar with the setup and use of BASE-9 *prior* to the incorporation of the ACS pipeline.





<br/>
<a name="Structure"/>
<h4>Directory Structure</h4>

When initializing the pipeline for a cluster, there are two initial inputs - the directory filepath and the cluster name: ACSpipeline(clusterval, basefolder)

The cluster name/identified is passed as "clusterval". This should be a string indicating the NGC number for the cluster. For example, if analyzing NGC 1261, the user would pass "1261" as clusterval. For NGC 104 (as in the example test case), clusterval would be '0104', with a leading zero.

The directory filepath ("basefolder") should be the home or parent directory that you are working in, including the BASE-9 code and this pipeline code. Within the "basefolder", there should be a directory for each cluster the user is analyzing, named NGC+clusterval.

Directory structure:<br/>
Top folder: "basefolder"<br/>
    -->  contains folders for each cluster (e.g. "./NGC" + clusterval)<br/>
    -->  contains .yaml files for BASE-9 (e.g.: "base9_" + clusterval + ".yaml")<br/>
    -->  contains other BASE-9 functionalities (e.g. makeCMD, singlePopMcmc, etc.)<br/>
    <br/>
Within each cluster-level folder: 'basefolder/NGC'+clusterval<br/>
    -->  contains original HST/ACS photometry files from survey publications (e.g.: NGC0104R.RDVIQ.cal.adj.zpt)<br/>
    -->  will contain photometry files, BASE-9 results files, generated plots etc.<br/>

Additional notes/suggestions:
- Either place this python code in the parent directory to 'self.basefolder'; else, give full pathname when initiating pipeline - realize that changes to or deviations from directory structure may lead to errors!
- When not using them, keep base9.yaml files named as base9_clusterval.yaml; temporarily rename to base9.yaml to initiate BASE-9 calculations, then return to labeled name
- Soft links should work fine when organizing your directory; see example image from this GitHub repo directory setup [here](https://github.com/rwk506/ACS_pipeline/blob/master/Directory_Structure.png).






<br/>
<a name="RepoStuff"/>
<h4>Repository Contents</h4>

This repo contains the following files:

- **ACS_Pipeline.ipynb**: Example use of the ACS pipeline pre-processing and analysis for test case NGC 104 <br/>
- **ACSpipe.py**: The primary python code, described further below<br/>
- **base9.yaml**: A typical BASE-9 yaml file, formatted for use, may be overwritten to generate isochrones (see additional notes on function calls, below)<br/>
- **base9_0104.yaml**: A backup BASE-9 yaml file for NGC 104, formatted specifically for use with NGC 104<br/>
- **plotresY.py**: See [BASE9_PlottingResults GitHub repo](https://github.com/rwk506/BASE9_PlottingResults) on visualizing Adaptive Metropolis MCMC chains<br/>
- **plotres_noY.py**: Same as previous for scenarios in which helium is an assumed, stationary parameter.<br/>
- **NGC0104** - Directory containing the following files:
  - *NGC0104.allACS.phot*: Photometry file formatted for BASE-9 including all stars from raw HST Treasury data
  - *NGC0104.cleanedACS.phot*: Photometry file formatted for BASE-9 from cleaned photometry
  - *NGC0104.sampleACS.phot*: Photometry file formatted for BASE-9 with randomly selected sample of stars
  - *NGC0104.v2.ms*: BASE-9 output file from makeCMD
  - *NGC0104.v2.res*: BASE-9 output file from singlePopMcmc
  - *NGC0104.v2.samplingACS.png*: Output from visualizing Adaptive Metropolis MCMC sampling chain
  - *NGC0104.v2ACS.png*: Output plot of resulting CMD with BASE-9 results
  - *NGC0104R.RDVIQ.cal.adj.zpt*: Raw photometry file from HST/ACS GCC Treasury Program


<br/>
<a name="Functions"/>
<h4>ACS Pipeline Functions</h4>

<br/>
<a name="Process"/>
<h4>Pre-processing</h4>

- **ACSpipeline(clusterval, basefolder)** - initiates processing/analysis pipeline for HST/ACS photometry of Galactic globular clusters<br/>

- **checkHB(self,showplt='Yes')** - Use to visually examine CMD and horizontal branch in particular (useful in conjunction with make_cuts() to clean up/remove HB stars).<br/>

- **make_cuts(self,hb1mag,hb1col,hb2mag,hb2col,showplt='Yes',strict='No',TOP=0.)** - Function to remove poorly determined stars from photometry, as well as HB. Stars with high in the 2.5% tails of the delta-pixel changes from frame to frame are removed (DX, DY). Stars with high errors in the 5% tail are removed (Verr, Ierr). Inputs: hbNmag and hbNcol give magnitude/color brighter/bluer limits for removing HB stars
    - *hb1mag* - faint magnitude limit 2 (brighter stars removed)
    - *hb1col* - color limit 2 (bluer stars removed)
    - *hb2mag* - faint magnitude limit 2 (brighter stars removed)
    - *hb2col* - color limit 2 (bluer stars removed)
    - *showplt* - toggle plot on/off
    - *strict* - toggle strict rules on for stars with one short exposure (NGC's 5986, 6397, 6779)
    - *TOP* - turn off point for clusters - used only when strict = 'Yes' to pinpoint bright, short exp. stars
            
- **reformat_phot(self)** - This function reformats the photometry into "full" and "cleaned" photometry files formatted for use with the BASE-9 program<br/>

- **subsample_phot(self,magcut606,magcutTOP=24,seedval=1823,n=3000)** - Randomly samples n stars from cleaned photometry and writes a photometry file for the randomly sampled stars. 
    - *magcut606* - n/2 stars are randomly sampled from the cluster stars fainter than the MSTOP and n/2 stars are randomly sampled from the stars in the cluster brighter than the MSTOP
    - *magcutTOP* - a faint limit imposed to avoid the lower main sequence dwarfs (which are not well modeled and have larger photometric uncertainties); default = 24, however, in practice, generally distance modulus + 7 magnitudes is used for analysis purposes
    - *seedval* - the seed value for the randomized selection; default seedval = 1823
    - *n* - number of stars to be randomly chosen, with n/2 stars brighter than the MSTOP and n/2 stars fainter than the MSTOP; default n = 3000

- **check_sample(self,titleon='Yes')** - This functions allows the user to double-check the results of the pre-processing of the ACS photometry, with cleaned photometry in black, and the randomly selected sample of n stars in green<br/>

<br/>
<br/>
<a name="Compute"/>
<h4>Post-Computation</h4>

- **plot_MCMCresults(self,vnum,meansline=27,saveYN='Yes',plotpriors='No',printtable='No',startn=1001,noY='False')** - Following the BASE-9 calculations, this function enables the user to examine the MCMC chains for each parameter and the resulting posterior distributions. An option also will provide a printed LaTeX table of the results.
    - *vnum* - file number of the results file to plot (e.g.: for a BASE-9 results file NGC0104.v2.res, vnum='2')
    - *meansline* - default value is 27 - DO NOT CHANGE unless the .yaml file format has been changed in recent BASE-9 versions
    - *plotpriors* - default is no, but allows the user to plot horizontal lines indicating the prior distribution over the MCMC chain and over the posterior distributions
    - *printtable* - default is no, user may elect to print the results in LaTeX format
    - *startn* - the line to start in the MCMC results file; the default is 1001 to skip the burnin, however it may be desirable or necessary in certain cases to change this
    - *noY* - if Y is included as a parameter, leave as False. If helium is set to a specific value during the calculations, then toggle to "True".

- **run_makeCMD(self,vnum,startingline=18,meansline=27,saveYN='Yes',startn=1001, noY='False', Yval=0.)** - A function designed to prep a base9_cluster.yaml file from the results file with median values from the posterior distributions. The base9.yaml file that results from this can be used to generate an isochrone with makeCMD.
    - *vnum* - file number of the results file to plot (e.g.: for a BASE-9 results file NGC0104.v2.res, vnum='2')
    - *startingline* - default value is 18 - DO NOT CHANGE unless the .yaml file format has been changed in current BASE-9 versions
    - *meansline* - default value is 27 - DO NOT CHANGE unless the .yaml file format has been changed in current BASE-9 versions
    - *startn* - the line to start in the MCMC results file; the default is 1001 to skip the burnin, however it may be desirable or necessary in certain cases to change this
    - *noY* - if Y is included as a parameter, leave as False. If helium is set to a specific value during the calculations, then toggle to "True"
    - *Yval* - if noY is set to true, specify the Y value that has been used
            

- **plot_CMDresults(self,vnum,saveYN='Yes',titleon='No')** - A function to plot the photometry from HST/ACS with the resulting isochrone from the BASE-9 results.
    - *vnum* - file number of the results file to plot (e.g.: for a BASE-9 results file NGC0104.v2.res, vnum='2')




<br/>
<a name="Code"/>
<h4>Test Case Example</h4>

An worked-through example of this pipeline in action is provided for the test cluster NGC 104. The Jupyter notebook of this analysis can be downloaded and edited for your own use or perused [here](https://github.com/rwk506/ACS_pipeline/blob/master/ACS_Pipeline.ipynb).




<br/>
<a name="Other"/>
<h4>Additional Information</h4>

Author: RWK

License: None, free to use and edit as people wish.

Contact: May be made through GitHub

<br /><br />
