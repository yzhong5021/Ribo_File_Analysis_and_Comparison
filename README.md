# Ribo_File_Analysis_and_Comparison
 Analysis and comparison of ribo files.
 ```python
ribo_commands(ribofile1,
               ribofile2 = None,
               experiments=[], experiments2=[]
               savefigs=False,
               plot=True,
               minlength=28, maxlength=32,
               do_regioncounts=True,
               lengthdist=True,
               metagene=True, metagene_radius=50,
               comparetranscriptcounts=True, regions=[],
               individualtranscript=False, annotations = "",  transcripts=[], transcriptregions = [], numpeaks = 3):
```
               
 ribofile1: The path to the first ribo file you are analyzing.
 
 ribofile2 (optional): The path to the second ribo file you are analyzing. If not specified, comparisons will be between experiments in the first ribo file.

 experiments/experiments2: Specify which experiments to analyze.
 
 savefigs: Save the generated figures.
 
 plot: Plot the figures when run.
 
 minlength/maxlength: Minimum and maximum ribosome footprint lengths included in the analysis.
 
 regioncounts: Plot and compare the number of footprints in the UTR5, UTR5 junction, CDS, UTR3 junction, and UTR3 regions.
 
 lengthdist: Plot and compare the distribution of ribosome footprint lengths.
 
 metagene: Plot and compare the metagene coverage around the start and stop sites.
    metagene_radius: The radius from the metagene in which to plot and compare.
    
 comparetranscriptcounts: plot and compare number of transcript counts.
    regions: Specify the regions ("UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3") you would like to compare. If left empty, analyzes the entire transcript.
    
 individualtranscript: plot and compare individual transcript coverages.
    annotations: The path to the annotation file used; required for individualtranscript to work.
    transcripts: Specify transcripts to analyze. If left empty, uses transcript count outliers.
    transcriptregions: Specify regions ("UTR5", "CDS", "UTR3") to analyze. If left empty, analyzes the entire transcript.
    numpeaks: Specify the number of peaks to identify.

 

