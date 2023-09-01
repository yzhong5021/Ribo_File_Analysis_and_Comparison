import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score as r2_score
from sklearn.metrics import mean_squared_error as mse
import pandas as pd
import math
import ribopy

#########PREPARING_RC_ARRAYS#################

TRANSCRIPT_REGIONS = ["UTR5", "UTR5_junction", "CDS",
            "UTR3_junction", "UTR3"]

def parse_regioncounts_file(regions):
    """
    Parses a region counts file and returns a
    list of counts for each region
    """
    num = len(regions[TRANSCRIPT_REGIONS[0]])

    counts = []
    for i in range(num):
        counts.append([regions[j][i] for j in TRANSCRIPT_REGIONS])
    return counts, TRANSCRIPT_REGIONS

#########ALL_RC_GRAPHS_COMBINED###############\

def RC_graphs_combined(counts, regions, experimentlist, plot, savefigs, filename, counts2=[], regions2=[], experimentlist2=[]):
    num_regions = len(regions)
    num_counts = len(counts)
    
    width = 0.8 / (num_counts * 2) 
    
    offsets = np.linspace(-width*(num_counts-1), width*(num_counts-1), num_counts)

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    for i in range(num_counts):
        rects1 = ax.bar(np.arange(num_regions)+offsets[i], counts[i], width, label=experimentlist[i])
        
        for rect in rects1:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')
            
    for i in range(num_counts):
        if(i<len(counts2)):
            rects2 = ax.bar(np.arange(num_regions)+ offsets[i]+ width, counts2[i], width, label=experimentlist2[i]+" (Ribo 2)")

            for rect in rects2:
                height = rect.get_height()
                ax.annotate('{}'.format(height),
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),
                            textcoords="offset points",
                            ha='center', va='bottom')

    ax.set_ylabel('Abundance')
    ax.set_title('Region Counts')
    ax.set_xticks(np.arange(num_regions))
    ax.set_xticklabels(regions)
    ax.legend(bbox_to_anchor=(1, 1), title='Counts')
    ax.set_ylim(0, max(sum(abundance) for abundance in counts) * 1.1)
    
    if savefigs:
        plt.savefig(filename + ".png")
    if plot:
        plt.show()

##############STACKED_GRAPHS###################
      
def stacked_RC_graphs(counts, experimentlist, savefigs, plot, filename, counts2 = [], experimentlist2 = []):
    percentcounts = []
    for count_list in counts:
        percentcounts.append([(count*100)/sum(count_list) for count in count_list])
    percentcounts2 = []
    for count_list in counts2:
        percentcounts2.append([(count*100)/sum(count_list) for count in count_list])
    percents = {}
    percents2 = {}
    for index, percent_list in enumerate(percentcounts):
        percents[index] = np.array([percent_list[i] for i in range(len(percent_list))])
    for index, percent_list in enumerate(percentcounts2):
        percents2[index] = np.array([percent_list[i] for i in range(len(percent_list))])
    sections = ["UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3"]
    categories = list(percents.keys())
    categories2 = list(percents2.keys())
    x = list(range(len(categories)))
    x2 = list(range(len(categories), len(categories) + len(categories2)))

    width = 0.3
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    bottom = np.zeros(len(categories))
    bottom2 = np.zeros(len(categories2))
    for i, section in enumerate(sections):
        percent_data = []
        percent_data2 = []
        for j in percents:
            percent_data.append(percents[j][i])
        rects = ax.bar(x, percent_data, width, label=section, bottom=bottom)
        bottom += percent_data
        if(len(counts2) > 0):
            for j in percents2:
                percent_data2.append(percents2[j][i])
            rects = ax.bar(x2, percent_data2, width, label=section + " (Ribo 2)", bottom=bottom2)
            bottom2 += percent_data2
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Region Counts')
    ax.set_title('Region Counts by Percentage')
    ax.set_xticks(x+x2)
    for experiment in experimentlist2:
        experiment = experiment + "_2"
    ax.set_xticklabels(experimentlist+experimentlist2)
    ax.legend(bbox_to_anchor=(1, 1), title='Sections')
    ax.set_ylim(0, 100)
    if(savefigs):
        plt.savefig(filename+".png")
    if(plot):
        plt.show()


#####PLOT_LD_GRAPHS#####
   
def plot_LD_graphs(lengths_list, abundances_list, experiments, file_prefix, plot, savefigs, logscale=False,
                   abundances_list2=[], experiments2=[]):
    compare = True
    if(len(experiments2)==0):
        compare = False
    experiments = list(experiments)
    experiments2 = list(experiments2)
    all_abundances = abundances_list+abundances_list2
    all_experiments = experiments+experiments2
    index = len(abundances_list)-1
    colors = ["Red", "Black", "Green", "Yellow", "Magenta", "Cyan"]
    colors2 = ["Blue", "Magenta", "Cyan", "Red", "Black", "Green"]
    linestyles = ["solid", "dotted", "dashed", "dashdot"]
    regions = ["UTR5 ", "UTR5 Junction ", "CDS ", "UTR3 Junction ", "UTR3 "]
    ####MULTIPLE_REGIONS####
    plt.figure(figsize=(15, 15))

    plt.title(file_prefix, fontsize=20)

    plt.xlabel('Length', fontsize=15)
    plt.ylabel('Abundance', fontsize=15)
    if(compare):
        if(isinstance(all_abundances[0][0], list)):
            for i, abundances in enumerate(abundances_list):
                for j, experiment in enumerate(experiments):
                    plt.plot(lengths_list, abundances[j], color=colors[i%6], linestyle=linestyles[i%4], label= regions[i%5] + experiment)
            for i, abundances in enumerate(abundances_list2):
                for j, experiment in enumerate(experiments2):
                    plt.plot(lengths_list, abundances[j], color=colors2[i%6], linestyle=linestyles[i%4], label= regions[i%5] + experiment+" (Ribo 2)")
        else:
            for i, experiment in enumerate(experiments):
                plt.plot(lengths_list, abundances_list[i], color=colors[i%6], linestyle=linestyles[i%4], label= experiment)
            for i, experiment in enumerate(experiments2):
                plt.plot(lengths_list, abundances_list2[i], color=colors2[i%6], linestyle=linestyles[i%4], label= experiment+" (Ribo 2)")
            
    else:
        for i, experiment in enumerate(experiments):
            plt.plot(lengths_list, abundances_data, color=colors[(i % 6)], linestyle=linestyles[i % 4], label=experiment+title)

    plt.legend()
    plt.xticks(lengths_list)
    plt.locator_params(axis='y', nbins=20)
#     if(isinstance(all_abundances[0][0], list)):
#         max_y = max(max(max(all_abundances)))
#     else:
#         max_y = max(max(all_abundances))
#     plt.yticks([i * (max_y/20) for i in range(0, 21)])
    plt.legend()
    if(logscale):
        plt.yscale("log")
    if(savefigs):
        plt.savefig(file_prefix + ".png")
    if(plot):
            plt.show()
    
####PLOT_METAGENE####
def plot_metagene_graphs(meta_radius, start_abunlist, stop_abunlist, indices, experiments, plot, savefigs, start_abunlist2 = [], stop_abunlist2 = [], experiments2 = []):
    colors = ["Red", "Black", "Green", "Yellow", "Magenta", "Cyan"]
    colors2 = ["Blue", "Magenta", "Cyan", "Red", "Black", "Green"]
    compare = True
    if(len(start_abunlist2) == 0):
        compare = False
    # COMBINED_METAGENE_GRAPHS
    plt.figure(figsize=(30, 15))
    plt.title("All Metagene Positions", fontsize=20)
    plt.xlabel("Position", fontsize=15)
    plt.ylabel("Abundance", fontsize=15)
    if(compare):
        for i, experiment in enumerate(experiments):
            plt.plot(indices, start_abunlist[i], color=colors[i % 6], linestyle="solid", label=experiment + " start")
            plt.plot(indices, stop_abunlist[i], color=colors[i % 6], linestyle="dotted", label=experiment + " stop")
        for i, experiment in enumerate(experiments2):
            plt.plot(indices, start_abunlist2[i], color=colors2[i % 6], linestyle="solid", label=experiment + " start (Ribo 2)")
            plt.plot(indices, stop_abunlist2[i], color=colors2[i % 6], linestyle="dotted", label=experiment + " stop (Ribo 2)")
    else:
        for i, experiment in enumerate(experiments):
            plt.plot(indices, start_abunlist[i], color=colors[i % 6], linestyle="solid", label=experiment + " start")
            plt.plot(indices, stop_abunlist[i], color=colors[i % 6], linestyle="dotted", label=experiment + " stop")
            
    plt.xticks(range(-meta_radius, meta_radius + 1))
    # max_y = max(max(max(start_abun) for start_abun in start_abunlist), max(max(stop_abun) for stop_abun in stop_abunlist), max(max(start_abun2) for start_abun2 in start_abunlist2), max(max(stop_abun2) for stop_abun2 in stop_abunlist2))
    # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
    plt.locator_params(axis='y', nbins=20)
    plt.legend()
    if savefigs:
        plt.savefig('combined_metagene.png')
    if plot:
        plt.show()

    if(not compare):
        for i, experiment in enumerate(experiments):
            ####METAGENE_GRAPHS_COMBINED_BY_EXPERIMENT####
            plt.figure(figsize=(30, 15))
            plt.title(experiment + " Metagene Positions", fontsize=20)
            plt.xlabel("Position", fontsize=15)
            plt.ylabel("Abundance", fontsize=15)
            plt.plot(indices, start_abunlist[i], color="red", linestyle="solid", label=experiment + " start")
            plt.plot(indices, stop_abunlist[i], color="blue", linestyle="dotted", label=experiment + " stop")
            plt.xticks(range(-meta_radius, meta_radius + 1))
            max_y = max(max(start_abunlist[i]), max(stop_abunlist[i]))
            plt.yticks([i * (max_y / 20) for i in range(0, 21)])
            plt.legend()
            if savefigs:
                plt.savefig(experiment + '_combined_metagene.png')
            if plot:
                plt.show()

            ####START_METAGENE####
            plt.figure(figsize=(30, 15))
            plt.title(experiment + " Start Metagene Positions", fontsize=20)
            plt.xlabel("Position", fontsize=15)
            plt.ylabel("Abundance", fontsize=15)
            plt.plot(indices, start_abunlist[i], color="red", linestyle="solid", label=experiment + " start")
            plt.xticks(range(-meta_radius, meta_radius + 1))
            # max_y = max(start_abunlist[i])
            # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
            plt.locator_params(axis='y', nbins=20)
            plt.legend()
            if savefigs:
                plt.savefig(experiment + '_start_metagene.png')
            if plot:
                plt.show()

            ####STOP_METAGENE####
            plt.figure(figsize=(30, 15))
            plt.title(experiment + " Stop Metagene Positions", fontsize=20)
            plt.xlabel("Position", fontsize=15)
            plt.ylabel("Abundance", fontsize=15)
            plt.plot(indices, stop_abunlist[i], color="blue", linestyle="solid", label=experiment + " stop")
            plt.xticks(range(-meta_radius, meta_radius + 1))
            # max_y = max(stop_abunlist[i])
            # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
            plt.locator_params(axis='y', nbins=20)
            plt.legend()
            if savefigs:
                plt.savefig(experiment + '_stop_metagene.png')
            if plot:
                plt.show()
    else:
        ####START_METAGENE_COMPARE####
        plt.figure(figsize=(30, 15))
        plt.title("Compared Start Metagene Positions", fontsize=20)
        plt.xlabel("Position", fontsize=15)
        plt.ylabel("Abundance", fontsize=15)
        for i,experiment in enumerate(experiments):
            plt.plot(indices, start_abunlist[i], color=colors[i%6], label=experiment + " start")
        for i,experiment in enumerate(experiments2):
            plt.plot(indices, start_abunlist2[i], color=colors2[i%6], label=experiment + " start (Ribo 2)")

        plt.xticks(range(-meta_radius, meta_radius + 1))
        # max_y = max(max(max(start_abunlist),max(start_abunlist2)))
        # print(max_y)
        # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
        plt.locator_params(axis='y', nbins=20)
        plt.legend()
        if savefigs:
            plt.savefig('compare_start_metagene.png')
        if plot:
            plt.show()


        ####STOP_METAGENE_COMPARE####
        plt.figure(figsize=(30, 15))
        plt.title("Compared Stop Metagene Positions", fontsize=20)
        plt.xlabel("Position", fontsize=15)
        plt.ylabel("Abundance", fontsize=15)
        for i,experiment in enumerate(experiments):
            plt.plot(indices, stop_abunlist[i], color=colors[i%6], label=experiment + " stop")
        for i,experiment in enumerate(experiments2):
            plt.plot(indices, stop_abunlist2[i], color=colors2[i%6], label=experiment + " stop (Ribo 2)")

        plt.xticks(range(-meta_radius, meta_radius + 1))
        # max_y = max(max(max(stop_abunlist),max(stop_abunlist2)))
        # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
        plt.locator_params(axis='y', nbins=20)
        plt.legend()
        if savefigs:
            plt.savefig('compare_stop_metagene.png')
        if plot:
            plt.show()


########TC_SCATTERPLOTS#################

def scatterplot(x, y, labels, xaxis, yaxis, title, filename, plot, savefigs):
    if(plot):
        print(title)

    ####SPEARMAN####
        a, b = spearmanr(x, y)
        print("Spearman correlation coefficient: "+str(a)+"\n")


    #####PLOT_WITH_LABELS######
    plt.figure(figsize=(20, 20))
    plt.xscale("log")
    plt.yscale("log")
    scatter = plt.scatter(x, y, color='blue', alpha=0.5)

    x_low_percentile = np.percentile(x, 0.03)
    x_high_percentile = np.percentile(x, 99.97)
    y_low_percentile = np.percentile(y, 0.03)
    y_high_percentile = np.percentile(y, 99.97)

    ####LABEL_&_IDENTIFY_OUTLIERS
    x_no_outliers = []
    y_no_outliers = []
    labels_outliers = []

    for i, label in enumerate(labels):
        if x[i] >= x_low_percentile and x[i] <= x_high_percentile and y[i] >= y_low_percentile and y[i] <= y_high_percentile:
            x_no_outliers.append(x[i])
            y_no_outliers.append(y[i])
        else:
            if(plot):
                print(label)
                plt.annotate(label.split("|")[0], (x[i], y[i]), textcoords="offset points", xytext=(0, 10), ha='left', fontsize = 15)
            labels_outliers.append(label)

    ####REGRESSION####
    newX = np.logspace(0, math.log(max(max(x,y)),10)+0.5, base=10)
    def myExpFunc(x, a, b):
        return a * np.power(x, b)
    if(plot):
        popt, pcov = curve_fit(myExpFunc, x_no_outliers, y_no_outliers)
        plt.plot(newX, myExpFunc(newX, *popt), 'r-', 
                 label="({0:.3f}*x**{1:.3f})".format(*popt))
        print("Exponential Fit: y = (a*(x**b))")
        print("\ta = popt[0] = {0}\n\tb = popt[1] = {1}".format(*popt))
        x_np = np.array(x)
        r2 = r2_score(y, popt[0]*(x_np**popt[1]))
        sig3 = math.sqrt(mse(y,popt[0]*(x_np**popt[1])))
        print("\nR2 Score: " + str(r2))
        print("\nResidual 3-Sigma value: " + str(int(sig3)))

    ####FORMAT####
    if(plot):
        plt.xlabel(xaxis, fontsize = 15)
        plt.ylabel(yaxis, fontsize = 15)
        plt.title(title, fontsize = 20)
        plt.grid(True)
        if(savefigs):
            plt.savefig(filename+".png")
        plt.show()


####PLOT_RIBO_DISTRIBUTIONS_OF_INDIVIDUAL_TRANSCRIPTS####

def plot_trans_dist(transcript, abundances, numpeaks, fullname, plot, savefigs, start = 1):
    shortname = fullname.split("|")[0]
    
    positions = list(range(start, start + len(abundances)))
    
    plt.figure(figsize = (20,20))
    plt.plot(positions, abundances)
    plt.xlabel("Position", fontsize = 15)
    plt.ylabel("Abundance", fontsize = 15)
    plt.title(shortname, fontsize = 20)

    ####TICKS####
    # maxabunval = max(abundances)
    # maxposval = positions[-1]
    # minposval = positions[0]
    # plt.yticks([i * (maxabunval/20) for i in range(0, 21)])
    # plt.xticks([i * (maxposval/20)+minposval for i in range(0, 21)])
    plt.locator_params(axis='y', nbins=20)
    plt.locator_params(axis='x', nbins=20)  

    ####PEAKS####
    max_values = []
    indexed_abundances = list(enumerate(abundances, start=1))
    indexed_abundances.sort(key=lambda x: x[1], reverse=True)

    for i in range(numpeaks):
        if i < len(indexed_abundances):
            value, index = indexed_abundances[i]
            max_values.append((value, index))
    print("Peaks:")
    for value, index in max_values:
        print(f"Position: {value}, Abundance:{index}")
    plt.grid(True)
    if(savefigs):
        plt.savefig(fullname + ".png")
    if(plot):
        plt.show()

    ####ZOOM_IN_TO_TOP_PEAK####
    if(len(max_values)>0):
        index, value = max_values[0]
        if(index-15 < 0):
            leftbound = 0
            rightbound = index + (15-index) + 15
        elif(index+15 < len(max_values)):
            leftbound = index - (15-index)-15
            rightbound = len(max_values)-1
        else: 
            leftbound = index-15
            rightbound = index+15
        peakabundances = abundances[leftbound:rightbound]
        peakpositions = positions[leftbound:rightbound]
        plt.figure(figsize = (15,15))
        plt.plot(peakpositions, peakabundances)
        plt.title(shortname + " Zoom-in at Highest Peak", fontsize = 15)
        plt.xlabel("Position", fontsize = 11)
        plt.ylabel("Abundance", fontsize = 11)
    
        ####TICKS_FOR_ZOOM####
        
        # max_y = max(peakabundances)
        # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
        plt.locator_params(axis='y', nbins=20)
        plt.xticks(peakpositions)
        plt.grid(True)
        if(savefigs):
            plt.savefig(fullname+"_Zoom.png")
        if(plot):
            plt.show()


def findoutliers(x, y, labels):
    x_low_percentile = np.percentile(x, 0.03)
    x_high_percentile = np.percentile(x, 99.97)
    y_low_percentile = np.percentile(y, 0.03)
    y_high_percentile = np.percentile(y, 99.97)
    labels_outliers = []
    for i, label in enumerate(labels):
        if ((x[i] < x_low_percentile or x[i] > x_high_percentile) and (y[i] < y_low_percentile or y[i] > y_high_percentile)):
            labels_outliers.append(label)
    return(labels_outliers)




###########################################################################
############################ actual commands ##############################
###########################################################################

def do_regioncounts(
                ribo1,
                ribo2=None,
                region_names=["UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3"],
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                regioncounts=True,
                experiments=[], experiments2=[]):
    regions = {"UTR5":[],"UTR5_junction":[],"CDS":[],"UTR3_junction":[],"UTR3":[]}
    regioncounts = {"UTR5":[],"UTR5_junction":[],"CDS":[],"UTR3_junction":[],"UTR3":[]}
    regions2 = {"UTR5":[],"UTR5_junction":[],"CDS":[],"UTR3_junction":[],"UTR3":[]}
    regioncounts2 = {"UTR5":[],"UTR5_junction":[],"CDS":[],"UTR3_junction":[],"UTR3":[]}
    for i in region_names:
        regions[i] = ribo1.get_region_counts(i,
                                            sum_lengths=True,
                                            sum_references=True,
                                            range_lower=minlength,
                                            range_upper=maxlength,
                                            experiments=experiments
                                            )
        for experiment in experiments:
            regioncounts[i].append(int(regions[i].loc[:, experiment]))
        if(compareribos):
            regions2[i] = ribo2.get_region_counts(i,
                                sum_lengths=True,
                                sum_references=True,
                                range_lower=minlength,
                                range_upper=maxlength,
                                experiments=experiments2
                                )
            for experiment in experiments2:
                regioncounts2[i].append(int(regions2[i].loc[:, experiment]))
    
    RCcounts, regions = parse_regioncounts_file(regioncounts)
    if compareribos:
        RCcounts2, regions2 = parse_regioncounts_file(regioncounts2)
    else:
        RCcounts2 = []
        regions2 = []
        
    RC_graphs_combined(RCcounts, regions=regions, counts2 = RCcounts2, regions2 = regions2, experimentlist=experiments,
                        experimentlist2 = experiments2, filename="region_counts", plot=plot, savefigs=savefigs)
    stacked_RC_graphs(RCcounts, experiments, savefigs,
                        plot, filename="region_counts%", counts2 = RCcounts2, experimentlist2 = experiments2)

def do_lengthdist(
                ribo1,
                ribo2=None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                experiments=[], experiments2=[]
    ):
    CDS = ribo1.get_length_dist("CDS", experiments=experiments)
    UTR5j = ribo1.get_length_dist("UTR5_junction", experiments=experiments)
    UTR3j = ribo1.get_length_dist("UTR3_junction", experiments=experiments)
    UTR5 = ribo1.get_length_dist("UTR5", experiments=experiments)
    UTR3 = ribo1.get_length_dist("UTR3", experiments=experiments)
    CDScounts = []
    UTR5jcounts = []
    UTR3jcounts = []
    UTR3counts = []
    UTR5counts = []
    for experiment in experiments:
        CDScounts.append(
            CDS.loc[minlength:maxlength, experiment].astype(int).tolist())
        UTR5jcounts.append(
            UTR5j.loc[minlength:maxlength, experiment].astype(int).tolist())
        UTR3jcounts.append(
            UTR3j.loc[minlength:maxlength, experiment].astype(int).tolist())
        UTR3counts.append(
            UTR3.loc[minlength:maxlength, experiment].astype(int).tolist())
        UTR5counts.append(
            UTR5.loc[minlength:maxlength, experiment].astype(int).tolist())
    allcounts = []
    allcounts.append(UTR5counts)
    allcounts.append(UTR5jcounts)
    allcounts.append(CDScounts)
    allcounts.append(UTR3jcounts)
    allcounts.append(UTR3counts)

    lengths_list = [i+minlength for i in range(maxlength+1-minlength)]
    
    if(not compareribos):
        plot_LD_graphs(lengths_list, UTR5counts,
                        experiments, "UTR5LD", plot = plot, savefigs=savefigs)
        plot_LD_graphs(lengths_list, UTR5jcounts,
                        experiments, "UTR5jLD", plot = plot, savefigs=savefigs)
        plot_LD_graphs(lengths_list, CDScounts,
                        experiments, "CDSLD", plot = plot, savefigs=savefigs)
        plot_LD_graphs(lengths_list, UTR3jcounts,
                        experiments, "UTR3jLD", plot = plot,  savefigs=savefigs)
        plot_LD_graphs(lengths_list, UTR3counts,
                        experiments, "UTR3LD",plot = plot, savefigs=savefigs)
        plot_LD_graphs(lengths_list, allcounts, experiments, "All LD", logscale = True, plot = plot, savefigs = savefigs)

    if(compareribos):
        CDS_2 = ribo2.get_length_dist("CDS", experiments=experiments2)
        UTR5j_2 = ribo2.get_length_dist("UTR5_junction", experiments=experiments2)
        UTR3j_2 = ribo2.get_length_dist("UTR3_junction", experiments=experiments2)
        UTR5_2 = ribo2.get_length_dist("UTR5", experiments=experiments2)
        UTR3_2 = ribo2.get_length_dist("UTR3", experiments=experiments2)
        CDScounts_2 = []
        UTR5jcounts_2 = []
        UTR3jcounts_2 = []
        UTR3counts_2 = []
        UTR5counts_2 = []
        
        for experiment in experiments2:
            CDScounts_2.append(
                CDS_2.loc[minlength:maxlength, experiment].astype(int).tolist())
            UTR5jcounts_2.append(
                UTR5j_2.loc[minlength:maxlength, experiment].astype(int).tolist())
            UTR3jcounts_2.append(
                UTR3j_2.loc[minlength:maxlength, experiment].astype(int).tolist())
            UTR3counts_2.append(
                UTR3_2.loc[minlength:maxlength, experiment].astype(int).tolist())
            UTR5counts_2.append(
                UTR5_2.loc[minlength:maxlength, experiment].astype(int).tolist())
            
        allcounts_2 = []
        for experiment in experiments:
            allcounts_2.append(UTR5counts_2)
            allcounts_2.append(UTR5jcounts_2)
            allcounts_2.append(CDScounts_2)
            allcounts_2.append(UTR3jcounts_2)
            allcounts_2.append(UTR3counts_2)

        plot_LD_graphs(lengths_list, UTR5counts, experiments, "UTR5 Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR5counts_2, experiments2=experiments2)
        plot_LD_graphs(lengths_list, UTR5jcounts, experiments, "UTR5 Junction Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR5jcounts_2, experiments2=experiments2)
        plot_LD_graphs(lengths_list, CDScounts, experiments, "CDS Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=CDScounts_2, experiments2=experiments2)
        plot_LD_graphs(lengths_list, UTR3jcounts, experiments, "UTR3 Junction Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR3jcounts_2, experiments2=experiments2)
        plot_LD_graphs(lengths_list, UTR3counts, experiments, "UTR3 Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR3counts_2, experiments2=experiments2)
        
        allcounts2 = [UTR5counts_2, UTR5jcounts_2,
                    CDScounts_2, UTR3jcounts_2, UTR3counts_2]

        plot_LD_graphs(lengths_list, allcounts, experiments, "Compare All Lengths", abundances_list2 = allcounts2, experiments2 = experiments2,
                        logscale=True, plot = plot, savefigs=savefigs)

def do_metagene(
                ribo1,
                ribo2=None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                metagene=True, metagene_radius=50,
                experiments=[], experiments2=[]
):
    starttest = ribo1.get_metagene(
    "start", experiments=experiments, range_lower=minlength, range_upper=maxlength)
    stoptest = ribo1.get_metagene(
        "stop", experiments=experiments, range_lower=minlength, range_upper=maxlength)
    
    startarray = []
    stoparray = []
    for experiment in starttest.index:
        values = starttest.loc[experiment,
                                (0 - 1 - metagene_radius):(0 + metagene_radius)].tolist()
        startarray.append(values)
        values = stoptest.loc[experiment,
                                (0 - 1 - metagene_radius):(0 + metagene_radius)].tolist()
        stoparray.append(values)
        
    indices = list(range(-metagene_radius, metagene_radius+1))
    if(not compareribos):
        plot_metagene_graphs(meta_radius=metagene_radius,
                                start_abunlist=startarray,
                                stop_abunlist=stoparray,
                                indices=indices,
                                experiments=experiments,
                                plot=plot,
                                savefigs=savefigs)
    if(compareribos):
        starttest2 = ribo2.get_metagene(
            "start", experiments=experiments2, range_lower=minlength, range_upper=maxlength)
        stoptest2 = ribo2.get_metagene(
            "stop", experiments=experiments2, range_lower=minlength, range_upper=maxlength)
        startarray2 = []
        stoparray2 = []
        for experiment in starttest2.index:
            values = starttest2.loc[experiment,
                                    (0 - 1 - metagene_radius):(0 + metagene_radius)].tolist()
            startarray2.append(values)
            values = stoptest2.loc[experiment,
                                    (0 - 1 - metagene_radius):(0 + metagene_radius)].tolist()
            stoparray2.append(values)
        plot_metagene_graphs(meta_radius=metagene_radius,
                                start_abunlist=startarray,
                                stop_abunlist=stoparray,
                                indices=indices,
                                experiments=experiments,
                                plot=plot,
                                savefigs=savefigs,
                                start_abunlist2 = startarray2,
                                stop_abunlist2 = stoparray2,
                                experiments2 = experiments2)

def do_comparetranscriptcounts(
                ribo1,
                ribo2=None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                regions=[],
                experiments=[], experiments2=[]
):
    if not compareribos:
        if(len(experiments)==0):
            print("Please include the experiments you would like to compare.")
            return
        UTR5transcriptcounts = ribo1.get_region_counts(
            "UTR5", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR5jtranscriptcounts = ribo1.get_region_counts(
            "UTR5_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        CDStranscriptcounts = ribo1.get_region_counts(
            "CDS", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR3jtranscriptcounts = ribo1.get_region_counts(
            "UTR3_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR3transcriptcounts = ribo1.get_region_counts(
            "UTR3", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        transcriptnames = UTR5transcriptcounts.index.tolist()
        UTR5tc = []
        UTR5jtc = []
        CDStc = []
        UTR3jtc = []
        UTR3tc = []
        for experiment in experiments:
            UTR5tc.append(UTR5transcriptcounts.loc[:, experiment].tolist())
            UTR5jtc.append(UTR5jtranscriptcounts.loc[:, experiment].tolist())
            CDStc.append(CDStranscriptcounts.loc[:, experiment].tolist())
            UTR3jtc.append(UTR3jtranscriptcounts.loc[:, experiment].tolist())
            UTR3tc.append(UTR3transcriptcounts.loc[:, experiment].tolist())
        if(len(regions) == 0):
            allexperimentscounts = []
            for i in range(len(experiments)):
                allexperimentscounts.append([sum(sublist) for sublist in zip(UTR5tc[i], UTR5jtc[i], CDStc[i], UTR3jtc[i], UTR3tc[i])])
            for i in range(len(experiments)):
                for j in range(i + 1, len(experiments)):
                    scatterplot(allexperimentscounts[i], allexperimentscounts[j], transcriptnames, experiments[i], experiments[j], experiments[i] + " Versus " + experiments[j] +
                                " Transcript Counts",
                                experiments[i] + "_Versus_" + experiments[j] + "_Transcript Counts", plot, savefigs)
        else:
            for i in range(len(experiments)):
                for j in range(i + 1, len(experiments)):
                    if "UTR5" in regions:
                        scatterplot(UTR5tc[i], UTR5tc[j], transcriptnames, experiments[i], experiments[j],
                                    experiments[i] + " Versus " +
                                    experiments[j] + " UTR5 Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments[j] + "_UTR5_Transcript Counts", plot, savefigs)
                    if "UTR5_junction" in regions:
                        scatterplot(UTR5jtc[i], UTR5jtc[j], transcriptnames, experiments[i], experiments[j],
                                    experiments[i] + " Versus " + experiments[j] +
                                    " UTR5 junction Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments[j] + "_UTR5_junction_Transcript Counts", plot, savefigs)
                    if "CDS" in regions:
                        scatterplot(CDStc[i], CDStc[j], transcriptnames, experiments[i], experiments[j],
                                    experiments[i] + " Versus " +
                                    experiments[j] + " CDS Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments[j] + "_UTR5_Transcript Counts", plot, savefigs)
                    if "UTR3" in regions:
                        scatterplot(UTR3jtc[i], UTR3jtc[j], transcriptnames, experiments[i], experiments[j],
                                    experiments[i] + " Versus " + experiments[j] +
                                    " UTR3 junction Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments[j] + "_UTR3_junction_Transcript Counts", plot, savefigs)
                    if "UTR3_junction" in regions:
                        scatterplot(UTR3tc[i], UTR3tc[j], transcriptnames, experiments[i], experiments[j],
                                    experiments[i] + " Versus " +
                                    experiments[j] + " UTR3 Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments[j] + "_UTR3_Transcript Counts", plot, savefigs)
    else:
        UTR5transcriptcounts = ribo1.get_region_counts(
            "UTR5", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR5jtranscriptcounts = ribo1.get_region_counts(
            "UTR5_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        CDStranscriptcounts = ribo1.get_region_counts(
            "CDS", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR3jtranscriptcounts = ribo1.get_region_counts(
            "UTR3_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR3transcriptcounts = ribo1.get_region_counts(
            "UTR3", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        UTR5transcriptcounts_2 = ribo2.get_region_counts(
            "UTR5", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
        UTR5jtranscriptcounts_2 = ribo2.get_region_counts(
            "UTR5_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
        CDStranscriptcounts_2 = ribo2.get_region_counts(
            "CDS", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
        UTR3jtranscriptcounts_2 = ribo2.get_region_counts(
            "UTR3_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
        UTR3transcriptcounts_2 = ribo2.get_region_counts(
            "UTR3", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
        transcriptnames = UTR5transcriptcounts.index.tolist()
        UTR5tc = []
        UTR5jtc = []
        CDStc = []
        UTR3jtc = []
        UTR3tc = []
        UTR5tc_2 = []
        UTR5jtc_2 = []
        CDStc_2 = []
        UTR3jtc_2 = []
        UTR3tc_2 = []
        for experiment in experiments:
            UTR5tc.append(UTR5transcriptcounts.loc[:, experiment].tolist())
            UTR5jtc.append(UTR5jtranscriptcounts.loc[:, experiment].tolist())
            CDStc.append(CDStranscriptcounts.loc[:, experiment].tolist())
            UTR3jtc.append(UTR3jtranscriptcounts.loc[:, experiment].tolist())
            UTR3tc.append(UTR3transcriptcounts.loc[:, experiment].tolist())
        for experiment in experiments2:
            UTR5tc_2.append(UTR5transcriptcounts_2.loc[:, experiment].tolist())
            UTR5jtc_2.append(UTR5jtranscriptcounts_2.loc[:, experiment].tolist())
            CDStc_2.append(CDStranscriptcounts_2.loc[:, experiment].tolist())
            UTR3jtc_2.append(UTR3jtranscriptcounts_2.loc[:, experiment].tolist())
            UTR3tc_2.append(UTR3transcriptcounts_2.loc[:, experiment].tolist())         
        if(len(regions) == 0):
            allexperimentscounts = []
            allexperimentscounts_2 = []
            for i in range(len(experiments)):
                allexperimentscounts.append([sum(sublist) for sublist in zip(UTR5tc[i], UTR5jtc[i], CDStc[i], UTR3jtc[i], UTR3tc[i])])
            for i in range(len(experiments2)):
                allexperimentscounts_2.append([sum(sublist) for sublist in zip(UTR5tc_2[i], UTR5jtc_2[i], CDStc_2[i], UTR3jtc_2[i], UTR3tc_2[i])])
            for i in range(len(experiments)):
                for j in range(len(experiments2)):
                    scatterplot(allexperimentscounts[i], allexperimentscounts_2[j], transcriptnames, experiments[i], experiments2[j] + " (Ribo 2)", experiments[i] + " Versus " + experiments2[j] +
                                " (Ribo 2) Transcript Counts",
                                experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_Transcript Counts", plot, savefigs)
        else:
            for i in range(len(experiments)):
                for j in range(len(experiments2)):
                    if "UTR5" in regions:
                        scatterplot(UTR5tc[i], UTR5tc_2[j], transcriptnames, experiments[i], experiments2[j]+  "(Ribo 2)",
                                    experiments[i] + " Versus " +
                                    experiments2[j] + " (Ribo 2) UTR5 Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR5_Transcript Counts", plot, savefigs)
                    if "UTR5_junction" in regions:
                        scatterplot(UTR5jtc[i], UTR5jtc_2[j], transcriptnames, experiments[i], experiments2[j]+ " (Ribo 2)",
                                    experiments[i] + " Versus " + experiments2[j] +
                                    " (Ribo 2) UTR5 junction Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_UTR5_junction_Transcript Counts", plot, savefigs)
                    if "CDS" in regions:
                        scatterplot(CDStc[i], CDStc_2[j], transcriptnames, experiments[i], experiments2[j]+ " (Ribo 2)",
                                    experiments[i] + " Versus " +
                                    experiments2[j] + " CDS Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR5_Transcript Counts", plot, savefigs)
                    if "UTR3" in regions:
                        scatterplot(UTR3jtc[i], UTR3jtc_2[j], transcriptnames, experiments[i], experiments2[j] + " (Ribo 2)",
                                    experiments[i] + " Versus " + experiments2[j] +
                                    " (Ribo 2) UTR3 junction Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR3_junction_Transcript Counts", plot, savefigs)
                    if "UTR3_junction" in regions:
                        scatterplot(UTR3tc[i], UTR3tc_2[j], transcriptnames, experiments[i], experiments2[j] + " (Ribo 2)",
                                    experiments[i] + " Versus " +
                                    experiments2[j] + " (Ribo 2) UTR3 Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR3_Transcript Counts", plot, savefigs)



def do_individualtranscript(
                ribo1,
                ribo2=None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                transcriptregions=[],
                numpeaks = 3,
                experiments=[], experiments2=[],
                comparetranscriptcounts=False,
                transcripts=[],
                annotations = []
):
    realregions = ["UTR5", "CDS", "UTR3"]
    if(len(transcripts)>=1):
        temptranscripts = transcripts
        transcripts = {}
        for region in realregions:
            transcripts[region] = temptranscripts

    if("UTR5_junction" in transcriptregions or "UTR3_junction" in transcriptregions):
        print("Warning: UTR5_junction and UTR3_junction are not currently supported for individual transcript analysis.")
        return

    all_tc = {}
    # transcriptnames = None
    for region in realregions:
        region_counts = ribo1.get_region_counts(
        region, sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        for experiment in experiments:
            all_tc[region] = [region_counts.loc[:, experiment].tolist()]
  
    transcriptnames = region_counts.index.tolist()

    if (not comparetranscriptcounts) and (len(transcripts) == 0) and compareribos:
        all_tc_2 = {}
        for region in realregions:
            region_counts = ribo2.get_region_counts(
                region, sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
            for experiment in experiments2:
                all_tc_2[region] = [region_counts.loc[:, experiment].tolist()]

    if(not compareribos):
        # We only have one ribo!
        for i in range(len(experiments)):
            for j in range(i + 1, len(experiments)):
                if (len(transcripts)==0):
                    transcripts = {}
                    for region in transcriptregions:
                        transcripts[region] = findoutliers(all_tc[region][i], all_tc[region][j], transcriptnames)
                        if(transcripts[region] == []):
                            print("No significant outlier transcripts found. Please specify transcripts to plot.")

                coverage = ribo1.get_coverage(experiments[i],  range_lower=minlength, range_upper=maxlength)
                coverage2 = ribo1.get_coverage(experiments[j], range_lower=minlength, range_upper=maxlength)

                if(len(transcriptregions) == 0):
                    print("Transcripts used in plotting: ")
                    print(transcripts["CDS"])
                    for transcript in transcripts["CDS"]:
                        plot_trans_dist(transcript, abundances = coverage[transcript], numpeaks = numpeaks,
                                    fullname=experiments[i]+ " "+transcript, plot=plot, savefigs=savefigs)
                        plot_trans_dist(transcript, abundances = coverage2[transcript], numpeaks = numpeaks,
                                    fullname=experiments[j]+ " "+transcript, plot=plot, savefigs=savefigs)

                else:
                    if annotations == "":
                        print("\nTo view transcript coverages in specific regions, please specify an annotation file.")
                        break
                    column_names = ['transcript', 'start', 'end', 'feature', 'score', 'strand']
                    print("Transcripts used in plotting: ")
                    for transcriptregion in transcriptregions:
                        print(transcripts[transcriptregion])
                    if(transcriptregions == []):
                        print(transcripts)
                    df = pd.read_csv(annotations, sep='\t', header=None, names=column_names)
                    for transcriptregion in transcriptregions:
                        for transcript in transcripts[transcriptregion]:  
                            rows = df[(df['transcript'] == transcript) & (df['feature'] == transcriptregion)]
                            indices = rows[['start', 'end']]
                            indices_values = indices.iloc[0]  
                            startvalue = indices_values['start']
                            endvalue = indices_values['end']
                            plot_trans_dist(transcript, abundances=coverage[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[i]+ " " + transcript), plot=plot, savefigs=savefigs)
                            plot_trans_dist(transcript, abundances=coverage2[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[j]+ " " + transcript), plot=plot, savefigs=savefigs)
                    if(len(transcriptregions) == 0):
                        for transcript in transcripts:
                            rows = df[(df['transcript'] == transcript) & (df['feature'] == transcriptregion)]
                            indices = rows[['start', 'end']]
                            indices_values = indices.iloc[0]  
                            startvalue = indices_values['start']
                            endvalue = indices_values['end']
                            plot_trans_dist(transcript, abundances=coverage[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[i]+ " " + transcript), plot=plot, savefigs=savefigs)
                            plot_trans_dist(transcript, abundances=coverage2[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[j]+ " " + transcript), plot=plot, savefigs=savefigs)

    else:
        for i in range(len(experiments)):
            for j in range(len(experiments2)):
                if (len(transcripts)==0):
                    transcripts = {}
                    for region in realregions:
                        transcripts[region] = findoutliers(all_tc[region][i], all_tc_2[region][j], transcriptnames)
                        if(transcripts[region] == []):
                            print("No significant outlier transcripts found. Please specify transcripts to plot.")

                coverage = ribo1.get_coverage(experiments[i],  range_lower = minlength, range_upper = maxlength)
                coverage2 = ribo2.get_coverage(experiments2[j], range_lower = minlength, range_upper = maxlength)

                if(len(transcriptregions) == 0):
                    print("Transcripts used in plotting: ")
                    print(transcripts["CDS"])
                    for transcript in transcripts["CDS"]:
                        plot_trans_dist(transcript, abundances = coverage[transcript], numpeaks = numpeaks,
                                    fullname = experiments[i] + " " + transcript, plot=plot, savefigs=savefigs)
                        plot_trans_dist(transcript, abundances = coverage2[transcript], numpeaks = numpeaks, 
                                    fullname = experiments2[j] + " (Ribo 2) " + transcript, plot = plot, savefigs = savefigs)

                else:
                    if annotations == "":
                        print("\nTo view transcript coverages in specific regions, please specify an annotation file.")
                        return
                    print("Transcripts used in plotting: ")
                    for transcriptregion in transcriptregions:
                        print(transcripts[transcriptregion])
                    if(transcriptregions == []):
                        print(transcripts)
                    column_names = ['transcript', 'start', 'end', 'feature', 'score', 'strand']
                    df = pd.read_csv(annotations, sep='\t', header=None, names=column_names)
                    for transcriptregion in transcriptregions:
                        for transcript in transcripts[transcriptregion]:  
                            rows = df[(df['transcript'] == transcript) & (df['feature'] == transcriptregion)]
                            indices = rows[['start', 'end']]
                            indices_values = indices.iloc[0]  
                            startvalue = indices_values['start']
                            endvalue = indices_values['end']
                            plot_trans_dist(transcript, abundances=coverage[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion + " " + transcript), plot=plot, savefigs=savefigs)
                            plot_trans_dist(transcript, abundances=coverage2[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion + " (Ribo 2) " + transcript), plot=plot, savefigs=savefigs)
        

  

def ribo_commands(ribofile1,
               ribofile2 = None,
               savefigs=False,
               plot=True,
               minlength=28, maxlength=32,
               experiments=[], experiments2=[],
               regioncounts=True,
               lengthdist=True,
               metagene=True, metagene_radius=50,
               comparetranscriptcounts=True, regions=[],
               individualtranscript=False, transcripts=[], transcriptregions = [], annotations = "", numpeaks = 3,
               ):
    """
    This is the main function!

    ribofile1: the first ribo file to be analyzed
    ribofile2: the second ribo file to be analyzed (optional)
    if you want to compare two ribo files, you must specify both ribo files
    savefigs: whether or not to save the figures
    plot: whether or not to plot the figures

    minlength: the minimum length of a ribosome footprint to be included in the analysis
    maxlength: the maximum length of a ribosome footprint to be included in the analysis

    experiments: the experiments to analyze in the first ribo file.
    experiments2: the experiments to analyze in the second ribo file.
    
    regioncounts: whether or not to plot the region counts

    lengthdist: whether or not to plot the length distribution

    metagene: whether or not to plot the metagene coverage
        - metagene_radius: the radius of the metagene coverage

    comparetranscriptcounts: whether or not to compare the transcript counts of the two ribo files
        - regions: the regions to compare the transcript counts of.
        If left empty, sum up all the regions and compare at once.

    individualtranscript: analyze the coverage of individual transcripts.
        - transcripts: the transcripts to analyze. If left empty, analyzes outlying transcripts.
        - transcriptregions: the regions to analyze. If left empty, analyzes all regions.
        - annotations: the annotations file to use. (Mandatory if individualtranscript is True)
        - numpeaks: Zooms in on certain peaks.


    Here are the possible commands:
    
    
    """
    ####SETUP####
    line = "------" * 10

    compareribos = False
    ribo1 = ribopy.Ribo(ribofile1)
    if(ribofile2 is not None):
        ribo2 = ribopy.Ribo(ribofile2)
        # ribos = [ribo1, ribo2]
        compareribos = True
    else:
        ribo2 = None
        
    if (len(experiments) == 0):
        experiments = ribo1.experiments
    if(compareribos):
        if(len(experiments2) == 0):
            experiments2 = ribo2.experiments

    ribo1.print_info()
    
    if(compareribos):
        print("\n" + line + "\n")
        print("Second Ribo File:\n")
        ribo2.print_info()
    print("\n" + line + "\n")
    region_names = ["UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3"]


    #### REGIONCOUNTS####
    if regioncounts:
        do_regioncounts(ribo1, ribo2, region_names,
                compareribos, savefigs,
                plot, minlength, maxlength,
                regioncounts, experiments, experiments2)
    ####LENGTHDIST####
    if lengthdist:
        do_lengthdist(ribo1, ribo2, compareribos,
                   savefigs, plot, minlength,
                   maxlength, experiments, experiments2)

    #### METAGENE####
    if metagene:
        do_metagene(ribo1, ribo2, compareribos,
                 savefigs, plot, minlength,
                 maxlength, metagene, metagene_radius,
                 experiments, experiments2)

    ####TRANSCRIPT_COUNTS####
    if comparetranscriptcounts:
        do_comparetranscriptcounts(ribo1, ribo2, compareribos,
                                savefigs, plot, minlength,
                                maxlength, regions,
                                experiments, experiments2)

    ####INDIVIDUAL_TRANSCRIPT####
    if individualtranscript:
        do_individualtranscript(ribo1, ribo2, compareribos,
                             savefigs, plot, minlength,
                             maxlength, transcriptregions, numpeaks,
                             experiments, experiments2,
                             comparetranscriptcounts, transcripts,
                             annotations)