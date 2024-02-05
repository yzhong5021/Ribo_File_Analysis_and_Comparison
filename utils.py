import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import f_oneway
from scipy.stats import spearmanr, t, ttest_ind
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from sklearn.metrics import r2_score as r2_score
from sklearn.metrics import mean_squared_error as mse
import pandas as pd
import math
import ribopy
from datetime import datetime
import random
from scipy.signal import find_peaks
from math import sin, cos, pi, sqrt
import warnings
from scipy import stats
from scipy import signal
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import csv
from scipy import stats
import statsmodels.api as sm



#########PREPARING_RC_ARRAYS#################

TRANSCRIPT_REGIONS = ["UTR5", "UTR5_junction", "CDS",
            "UTR3_junction", "UTR3"]

def parse_regioncounts_file(regions):

    num = len(regions[TRANSCRIPT_REGIONS[0]])

    counts = []
    for i in range(num):
        counts.append([regions[j][i] for j in TRANSCRIPT_REGIONS])
    return counts, TRANSCRIPT_REGIONS

#########ALL_RC_GRAPHS_COMBINED###############\

def RC_graphs_combined(riboname, riboname2, counts, regions, experimentlist, plot, savefigs, filename, counts2=[], regions2=[], experimentlist2=[]):
    num_regions = len(regions)
    num_counts = len(counts)
    
    width = 0.8 / (num_counts * 2) 
    
    offsets = np.linspace(-width*(num_counts-1), width*(num_counts-1), num_counts)

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    for i in range(num_counts):
        rects1 = ax.bar(np.arange(num_regions)+offsets[i], counts[i], width, label=experimentlist[i] + " " + riboname)
        
        for rect in rects1:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')
            
    for i in range(num_counts):
        if(i<len(counts2)):
            rects2 = ax.bar(np.arange(num_regions)+ offsets[i]+ width, counts2[i], width, label=experimentlist2[i]+" " + riboname2)

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
        hi = random.randint(0,1000000)
        plt.savefig(filename + "_"+str(hi)+"_region_counts_combined.png")
    if plot:
        plt.show()

##############STACKED_GRAPHS###################
      
def stacked_RC_graphs(riboname, riboname2, counts, experimentlist, savefigs, plot, filename, counts2 = [], experimentlist2 = []):
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
        rects = ax.bar(x, percent_data, width, label=section + " " + riboname, bottom=bottom)
        bottom += percent_data
        if(len(counts2) > 0):
            for j in percents2:
                percent_data2.append(percents2[j][i])
            rects = ax.bar(x2, percent_data2, width, label=section + " " + riboname2, bottom=bottom2)
            bottom2 += percent_data2
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Region Counts')
    ax.set_title('Region Counts by Percentage')
    ax.set_xticks(x+x2)
    for experiment in experimentlist2:
        experiment = experiment + "_2"
    experimentlist = list(experimentlist)
    experimentlist2 = list(experimentlist2)
    experimentlist.extend(experimentlist2)
    ax.set_xticklabels(experimentlist)

    ax.legend(bbox_to_anchor=(1, 1), title='Sections')
    ax.set_ylim(0, 100)
    if(savefigs):
        hi = random.randint(0,1000000)
        plt.savefig(filename+"_"+str(hi)+"_region_counts.pdf")
    if(plot):
        plt.show()


#####PLOT_LD_GRAPHS#####
   
def plot_LD_graphs(riboname, riboname2, lengths_list, abundances_list, experiments, file_prefix, plot, savefigs, logscale=False,
                   abundances_list2=[], experiments2=[], filename = ""):
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
    plt.figure(figsize=(2.4, 2.2))

    plt.title(file_prefix, fontsize=8)

    plt.xlabel('Length', fontsize=8)
    plt.ylabel('Abundance', fontsize=8)
    if(compare):
        if(isinstance(all_abundances[0][0], list)):
            for i, abundances in enumerate(abundances_list):
                for j, experiment in enumerate(experiments):
                    plt.plot(lengths_list, abundances[j], color=colors[i%6], linestyle=linestyles[i%4], label= regions[i%5] + experiment+" "+riboname)
            for i, abundances in enumerate(abundances_list2):
                for j, experiment in enumerate(experiments2):
                    plt.plot(lengths_list, abundances[j], color=colors2[i%6], linestyle=linestyles[i%4], label= regions[i%5] + experiment+" "+riboname2)
        else:
            for i, experiment in enumerate(experiments):
                plt.plot(lengths_list, abundances_list[i], color=colors[i%6], linestyle=linestyles[i%4], label= experiment+" "+riboname)
            for i, experiment in enumerate(experiments2):
                plt.plot(lengths_list, abundances_list2[i], color=colors2[i%6], linestyle=linestyles[i%4], label= experiment+" " + riboname2)
            
    else:
        if(isinstance(all_abundances[0][0], list)):
            total_abundances = []
            for h in range(len(experiments)):
                print(len(abundances_list))
                print(len(abundances_list[0]))
                print(len(abundances_list[0][0]))
                total_abundances.append([i+j+k+l+m for i,j,k,l,m in zip(abundances_list[0][h],abundances_list[1][h],abundances_list[2][h],abundances_list[3][h],abundances_list[4][h])])                         
            for j, experiment in enumerate(experiments):
                plt.plot(lengths_list, total_abundances[j], color=colors[j%6], linestyle=linestyles[j%4], label= experiment +" "+ riboname)
        else:
            for i, experiment in enumerate(experiments):
                plt.plot(lengths_list, abundances_list[i], color=colors[i%6], linestyle=linestyles[i%4], label= experiment+" " + riboname)

        # for i, experiment in enumerate(experiments):
        #     plt.plot(lengths_list, abundances_data, color=colors[(i % 6)], linestyle=linestyles[i % 4], label=experiment+title)

    plt.legend()
    plt.xticks(lengths_list)
    plt.locator_params(axis='y', nbins=20)
#     if(isinstance(all_abundances[0][0], list)):
#         max_y = max(max(max(all_abundances)))
#     else:
#         max_y = max(max(all_abundances))
#     plt.yticks([i * (max_y/20) for i in range(0, 21)])
    plt.legend()

    if(savefigs):
        hi = random.randint(0,1000000)
        plt.savefig(filename +"_"+ hi+ "_length_distributions.pdf")
    if(plot):
            plt.show()
    
####PLOT_METAGEN/PERIODICITYE####
def plot_metagene_graphs(riboname, riboname2, meta_radius, start_abunlist, stop_abunlist, indices, experiments, plot, savefigs, site_periodicity, start_abunlist2 = [], stop_abunlist2 = [], experiments2 = [], filename = ""):

    # def codonperiodicity(metagene_coverage, period, name, color):
    #     fft_result = np.fft.fft(metagene_coverage)
    #     n = len(fft_result)
    #     frequencies = np.fft.fftfreq(n)

    #     amplitude_spectrum = np.abs(fft_result) / n

    #     max_amplitude_index = np.argmax(amplitude_spectrum)
    #     max_amplitude_frequency = frequencies[max_amplitude_index]

    #     plt.plot(frequencies[:n // 2], amplitude_spectrum[:n // 2], color = color, label = name)
    #     peak_indices = np.where(frequencies[:n // 2] >= 0.25)[0]
    #     largest_peak_index = peak_indices[np.argmax(amplitude_spectrum[peak_indices])]
    #     largest_peak_frequency = frequencies[largest_peak_index]
    #     largest_peak_amplitude = amplitude_spectrum[largest_peak_index]
    #     plt.annotate(str(round(largest_peak_amplitude, 3)), 
    #                 xy=(largest_peak_frequency, largest_peak_amplitude), 
    #                 xytext=(largest_peak_frequency - 0.03, largest_peak_amplitude + 0.03))

    def codonperiodicity(metagene_coverage, name, color, offset):
        fft_result = np.fft.fft(metagene_coverage)
        n = len(fft_result)
        frequencies = np.fft.fftfreq(n)

        amplitude_spectrum = np.abs(fft_result) / n 
        maxamp = max(amplitude_spectrum)
        normalized_amplitude_spectrum = [i/maxamp for i in amplitude_spectrum]
        periods = [1/i if i != 0 else float("inf") for i in frequencies]
        plt.plot(periods, normalized_amplitude_spectrum, color=color, label=name)

        plt.xlim(0, 5)
        # Find the index of the next highest amplitude
        temp = normalized_amplitude_spectrum.copy()
        temp.sort(reverse = True)
        index = normalized_amplitude_spectrum.index(temp[1])

        # Annotate the next highest amplitude value
        next_highest_amplitude = normalized_amplitude_spectrum[index]
        next_highest_amplitude_period = 3

        annotation_text = str(round(next_highest_amplitude, 4))
        plt.annotate(annotation_text,
                    xy=(next_highest_amplitude_period, next_highest_amplitude),
                    xytext=(next_highest_amplitude_period + 0.17, next_highest_amplitude + offset), color  = color)
         
        plt.xticks([3], fontsize = 10)
        return(amplitude_spectrum)

    colors = ["Red", "Blue", "Green", "Black", "Magenta", "Cyan"]
    compare = True
    if(len(start_abunlist2) == 0):
        compare = False
    # COMBINED_METAGENE_GRAPHS
    # plt.figure(figsize=(2.4, 2.2))
    # plt.title(filename, fontsize=8)
    # plt.xlabel("Position", fontsize=8)
    plt.figure(figsize=(15,9))
    plt.xlabel("Position Relative to Start/Stop Codon")
    # plt.ylabel("Abundance", fontsize=8)
    plt.ylabel("Footprint Abundance")
    if(compare):
        colors2 = ["Blue", "Magenta", "Cyan", "Red", "Black", "Green"]
        for i, experiment in enumerate(experiments):
            plt.plot(indices, start_abunlist[i], color=colors[i % 6], linestyle="solid", label=experiment.split("_")[0] + " " + experiment.split("_")[-1] + " start "+ riboname)
            plt.plot(indices, stop_abunlist[i], color=colors[i % 6], linestyle="dotted", label=experiment.split("_")[0] + " " + experiment.split("_")[-1] + " stop " + riboname)
        for i, experiment in enumerate(experiments2):
            plt.plot(indices, start_abunlist2[i], color=colors2[i % 6], linestyle="solid", label=experiment + " start " + riboname2)
            plt.plot(indices, stop_abunlist2[i], color=colors2[i % 6], linestyle="dotted", label=experiment + " stop " + riboname2)
    else:
        offset = [-0.0015, 0, 0.0015]

        for i, experiment in enumerate(experiments):
            plt.plot(indices, start_abunlist[i], color=colors[i % 6], linestyle="solid", label=experiment.split("_")[0] + " " + experiment.split("_")[-1] + " Start", linewidth = 0.7)
            plt.plot(indices, stop_abunlist[i], color=colors[i % 6], linestyle="dotted", label=experiment.split("_")[0] + " " + experiment.split("_")[-1] + " Stop", linewidth = 0.7)
        # plt.xticks([], fontsize = 4)
        # plt.yticks(fontsize = 4)
        plt.locator_params(axis='y', nbins=6)
        # plt.legend()
        if savefigs:
            plt.savefig(filename+'combined_metagene.pdf')
        plt.xticks([x for x in range(-50, 51) if (x) % 3 == 0], fontsize = 10)
        plt.yticks([])
        plt.title("Metagene Coverage, "+experiments[i].split("_")[0])

        if plot:
            plt.legend(prop={'size': 10})
            plt.show()

        plt.figure(figsize=(15, 9))
        amplitudes = []
        for i, experiment in enumerate(experiments):
            print(experiment)
            temp = start_abunlist[i].index(max(start_abunlist[i]))
            test_start_abunlist = start_abunlist[i][temp:]
            temp2 = stop_abunlist[i].index(max(stop_abunlist[i]))
            test_stop_abunlist = stop_abunlist[i][:temp]


            if("start" in site_periodicity):
                a = codonperiodicity(test_start_abunlist, experiment + " Start", colors[i], offset[i])
                amplitudes.append(a)

            if("stop" in site_periodicity):
                codonperiodicity(test_stop_abunlist, experiment + " Stop", colors[i], offset[i])
        alternative = ["greater", "less", "less"]
        for i in range(len(amplitudes)):
            for j in range(i, len(amplitudes)):
                result = ttest_ind(amplitudes[i], amplitudes[j], alternative = alternative[i])
                p_value = result.pvalue
                print("P-Value for " + experiments[i] + " vs " + experiments[j] + " :", p_value)
        plt.xlabel("Period (Frequency⁻¹)")
        plt.ylabel("Normalized Amplitude")
        plt.xlim((0, 10))
        plt.legend(prop={'size': 10})
        plt.title("Codon Periodicity Strength (Fourier Transform), "+experiments[i].split("_")[0])
        plt.grid(True)
        plt.show()


            
    # plt.xticks([i * 2 - meta_radius for i in range(meta_radius + 1)])
    # max_y = max(max(max(start_abun) for start_abun in start_abunlist), max(max(stop_abun) for stop_abun in stop_abunlist), max(max(start_abun2) for start_abun2 in start_abunlist2), max(max(stop_abun2) for stop_abun2 in stop_abunlist2))
    # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
    # plt.xticks([], fontsize = 4)
    # plt.yticks(fontsize = 4)
    # plt.locator_params(axis='y', nbins=6)
    # plt.legend()
    # if savefigs:
    #     plt.savefig(filename+'combined_metagene.pdf')
    # if plot:
    #     plt.show()


    # if(not compare):
    #     for i, experiment in enumerate(experiments):
    #         ####METAGENE_GRAPHS_COMBINED_BY_EXPERIMENT####
    #         # plt.figure(figsize=(2.4, 2.2))
    #         # plt.title(experiment + " Metagene Positions", fontsize=8)
    #         # plt.xlabel("Position", fontsize=8)
    #         # plt.ylabel("Abundance", fontsize=8)
    #         # plt.plot(indices, start_abunlist[i], color="red", linestyle="solid", label=experiment + " start", linewidth = 2)
    #         # plt.plot(indices, stop_abunlist[i], color="blue", linestyle="dotted", label=experiment + " stop", linewidth = 2)
    #         # plt.xticks(range(-meta_radius, meta_radius + 1), fontsize = 2)
    #         # max_y = max(max(start_abunlist[i]), max(stop_abunlist[i]))
    #         # plt.yticks([i * (max_y / 10) for i in range(0, 11)], fontsize = 5)
    #         # plt.legend(fontsize = 8, posi)
    #         plt.figure(figsize=(24, 22))
    #         plt.title(experiment + " Metagene Positions", fontsize=25)  # Increased title font size
    #         plt.xlabel("Position", fontsize=8)  # Increased label font size
    #         plt.ylabel("Abundance", fontsize=8)  # Increased label font size

    #         # Plot start and stop lines with improved styling
    #         plt.plot(indices, start_abunlist[i], color="red", linestyle="solid", label=experiment + " start")
    #         plt.plot(indices, stop_abunlist[i], color="blue", linestyle="dotted", label=experiment + " stop")

    #         # Adjust x-axis ticks for better visibility
    #         plt.xticks([i * 2 - meta_radius for i in range(meta_radius + 1)])
    #   # Increased tick font size

    #         max_y = max(max(start_abunlist[i]), max(stop_abunlist[i]))
    #         # Adjust y-axis ticks for better visibility and use a larger range
    #         plt.yticks([i * (max_y / 5) for i in range(6)])  # Increased tick font size

    #         # Adjust legend position and font size
    #         plt.legend(fontsize = 20)  # Change "best" to the desired legend location

    #         if savefigs:
    #             plt.savefig(filename+ experiment + '_combined_metagene.png')
    #         if plot:
    #             plt.show()

            # ####START_METAGENE####
            # plt.figure(figsize=(30, 15))
            # plt.title(experiment + " Start Metagene Positions "+riboname, fontsize=20)
            # plt.xlabel("Position", fontsize=15)
            # plt.ylabel("Abundance", fontsize=15)
            # plt.plot(indices, start_abunlist[i], color="red", linestyle="solid", label=experiment + " start " + riboname)
            # plt.xticks(range(-meta_radius, meta_radius + 1))
            # # max_y = max(start_abunlist[i])
            # # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
            # plt.locator_params(axis='y', nbins=20)
            # plt.legend()
            # if savefigs:
            #     plt.savefig(filename + experiment + '_start_metagene.png')
            # if plot:
            #     plt.show()

            # ####STOP_METAGENE####
            # plt.figure(figsize=(30, 15))
            # plt.title(experiment + " Stop Metagene Positions " + riboname, fontsize=20)
            # plt.xlabel("Position", fontsize=15)
            # plt.ylabel("Abundance", fontsize=15)
            # plt.plot(indices, stop_abunlist[i], color="blue", linestyle="solid", label=experiment + " stop " + riboname)
            # plt.xticks(range(-meta_radius, meta_radius + 1))
            # # max_y = max(stop_abunlist[i])
            # # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
            # plt.locator_params(axis='y', nbins=20)
            # plt.legend()
            # if savefigs:
            #     plt.savefig(filename + experiment + '_stop_metagene.png')
            # if plot:
            #     plt.show()
    # else:
    #     ####START_METAGENE_COMPARE####
    #     plt.figure(figsize=(30, 15))
    #     plt.title("Compared Start Metagene Positions", fontsize=20)
    #     plt.xlabel("Position", fontsize=15)
    #     plt.ylabel("Abundance", fontsize=15)
    #     for i,experiment in enumerate(experiments):
    #         plt.plot(indices, start_abunlist[i], color=colors[i%6], label=experiment + " start "+ riboname)
    #     for i,experiment in enumerate(experiments2):
    #         plt.plot(indices, start_abunlist2[i], color=colors2[i%6], label=experiment + " start "+ riboname2)

    #     plt.xticks(range(-meta_radius, meta_radius + 1))
    #     # max_y = max(max(max(start_abunlist),max(start_abunlist2)))
    #     # print(max_y)
    #     # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
    #     plt.locator_params(axis='y', nbins=20)
    #     plt.legend()
    #     if savefigs:
    #         hi = random.randint(0,1000000) 
    #         plt.savefig(filename + hi+ 'compare_start_metagene.png')
    #     if plot:
    #         plt.show()


    #     ####STOP_METAGENE_COMPARE####
    #     plt.figure(figsize=(30, 15))
    #     plt.title("Compared Stop Metagene Positions", fontsize=20)
    #     plt.xlabel("Position", fontsize=15)
    #     plt.ylabel("Abundance", fontsize=15)
    #     for i,experiment in enumerate(experiments):
    #         plt.plot(indices, stop_abunlist[i], color=colors[i%6], label=experiment + " stop " + riboname)
    #     for i,experiment in enumerate(experiments2):
    #         plt.plot(indices, stop_abunlist2[i], color=colors2[i%6], label=experiment + " stop "+ riboname2)

    #     plt.xticks(range(-meta_radius, meta_radius + 1))
    #     plt.locator_params(axis='y', nbins=20)
    #     plt.legend()
    #     if savefigs:
    #         hi = random.randint(0,1000000) 
    #         plt.savefig(filename + hi+ 'compare_stop_metagene.png')
    #     if plot:
    #         plt.show()


########TC_SCATTERPLOTS#################

def scatterplot(riboname, riboname2, x, y, labels, xaxis, yaxis, title, filename, plot, savefigs, regression = False, correctUMI = False, strength = 0.4, allreadsx = None, allreadsy = None, proteinabundance = None, line = None):
    if(plot):
        print(title)

    ####SPEARMAN####
        a, b = spearmanr(x, y)
        print("Spearman correlation coefficient: "+str(a)+"\n")

    #####PLOT_WITH_LABELS######


    # Create a figure for the scatterplot with the same axis limits

    plt.figure(figsize=(15,15))
    # plt.figure(figsize=(24,22))
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1, 31622)
    plt.ylim(1, 31622)
    if(not correctUMI):
        
        plt.scatter(x, y, color='blue', alpha=0.5)

    x_low_percentile = np.percentile(x, 0.03)
    x_high_percentile = np.percentile(x, 99.97)
    y_low_percentile = np.percentile(y, 0.03)
    y_high_percentile = np.percentile(y, 99.97)
    
    ####REGRESSION_OR_CORRECTION####
    newX = np.logspace(0, math.log(max(max(x,y)),10)+0.5, base=10)
    # newX = np.logspace(0, 4, base = 10)
    def myExpFunc(x, a, b):
        return a * np.power(x, b)

    if(correctUMI): 
        def power_weight(residuals, power=2):
            return np.abs(residuals)**power
        def irls(x, y, max_iterations=10, tolerance=1e-4):
            popt, _ = curve_fit(myExpFunc, x, y)  # Initial curve fitting

            for _ in range(max_iterations):
                # Weighted exponential curve fitting using current parameters
                # _, pcov = curve_fit(myExpFunc, x, y, p0=popt)
                _, pcov = curve_fit(myExpFunc, x, y, sigma=power_weight(y - myExpFunc(x, *popt)))

                # Calculate residuals
                residuals = 0.4*(y - myExpFunc(x, *popt))

                # Update data points using residuals and weights
                y -= residuals

                # Re-fit the curve with updated data points
                popt, _ = curve_fit(myExpFunc, x, y, sigma=1/np.abs(residuals), p0=popt)

                # Check for convergence
                if np.max(np.abs(residuals)) < tolerance:
                    break

                return y
        tempx = x.copy() 
        tempy = y.copy()
        otherx = x.copy()
        othery = y.copy()

        indices_to_remove = []
        indices_to_remove_2 = []
        for i, count in enumerate(x):
            if count == 0 or y[i] == 0:
                indices_to_remove.append(i)
        #         if count < 10 or count > x_high_percentile or y[i] < 10 or y[i] > y_high_percentile:
        #             indices_to_remove_2.append(i)

        tempx = [tempx[i] for i in range(len(tempx)) if i not in (indices_to_remove or indices_to_remove_2)]
        tempy = [tempy[i] for i in range(len(tempy)) if i not in (indices_to_remove or indices_to_remove_2)]
        otherx = [otherx[i] for i in range(len(otherx)) if i not in indices_to_remove]
        othery = [othery[i] for i in range(len(othery)) if i not in indices_to_remove]
        # newlabels = [labels[i] for i in range(len(labels)) if i not in indices_to_remove]

        x_log = np.log(tempx)
        y_log = np.log(tempy)
        x_full_log = np.log(otherx)
        y_full_log = np.log(othery)
        newy = irls(x_log, y_log)
        adjusted_y = np.round(np.exp(newy))
        slope, intercept, r_value, p_value, std_err = linregress(np.exp(x_full_log), adjusted_y)
        # popt, pcov = curve_fit(myExpFunc, x_log, y_log)


        # plt.plot(newX, myExpFunc(newX, *popt), 'r-', 
        #             label="({0:.3f}*x**{1:.3f})".format(*popt))

        # popt, _ = curve_fit(myExpFunc, x_full_log, y_full_log)

        # residuals = y_full_log - myExpFunc(x_full_log, *popt)

        # adjusted_y_log = y_full_log - residuals * strength

        # adjusted_y = np.round(np.exp(adjusted_y_log))
        # regression_line = float(slope) * float(x) + intercept
        # plt.plot(x, regression_line, label='Regression Line', color='green')

        print("R^2 Value: ",r_value**2)
        plt.scatter(np.exp(x_full_log), adjusted_y, color = "red", label='Adjusted Data', alpha=0.5)


    if(plot):
        plt.plot(x, x, label='x = y')
        if(regression):

            popt, pcov = curve_fit(myExpFunc, x, y)
            plt.plot(newX, myExpFunc(newX, *popt), 'r-', 
                    label="({0:.3f}*x**{1:.3f})".format(*popt))
            print("Exponential Fit: y = (a*(x**b))")
            print("\ta = popt[0] = {0}\n\tb = popt[1] = {1}".format(*popt))
            x_np = np.array(x)
            r2 = r2_score(y, popt[0]*(x_np**popt[1]))
            # sig3 = math.sqrt(mse(y,popt[0]*(x_np**popt[1])))
            print("\nR2 Score: " + str(r2))
            # plt.annotate("R^2 Value (Exp): " + str(round(r2, 3)), xy = (1010, 300), fontsize = 15)
            # plt.annotate("Regression Slope (Exp): "+str(round((a), 3)), xy = (1010, 180), fontsize = 15)
            # print("\nResidual 3-Sigma value: " + str(int(sig3)))
            x = np.array(x)
            y = np.array(y)
            epsilon = 1e-10
            x = np.maximum(x, epsilon)
            y = np.maximum(y, epsilon)  

            log_X = np.log10(x)
            log_Y = np.log10(y)

            slope, intercept, r_value, p_value, std_err = linregress(log_X, log_Y)

            def predict(x):
                return slope * x + intercept

            log_Y_pred = [predict(log_x) for log_x in log_X]
            Y_pred = 10 ** np.array(log_Y_pred)

            r2 = r2_score(x, Y_pred)
            print("R^2 Value: ",r2)
            print("Regression Slope: ", slope)
            # plt.annotate("R^2 Value: " + str(round(r2, 3)), xy = (1010, 70), fontsize = 15)
            # plt.annotate("Regression Slope: "+str(round(slope, 5)), xy = (1010, 30), fontsize = 15)

            

    ####FORMAT####
    if(plot):
        plt.xlabel(xaxis.split("_")[0]+"_"+xaxis.split("_")[-1] + " (Log Transformed)", fontsize = 16)
        plt.ylabel(yaxis.split("_")[0]+"_"+yaxis.split("_")[-1] + " (Log Transformed)", fontsize = 16)
        plt.yticks()
        plt.xticks([10**0, 10**1, 10**2, 10**3, 10**4], fontsize = 10)
    
        plt.title(title, fontsize = 16)
        # plt.xlabel(xaxis.split("_")[0]+"_"+xaxis.split("_")[-1], fontsize = 20)
        # plt.ylabel(yaxis.split("_")[0]+"_"+yaxis.split("_")[-1], fontsize = 20)
        # plt.yticks(fontsize = 16)
        # plt.xticks(fontsize = 16)
        # plt.title(title, fontsize = 24)
        plt.grid(True)
        hi = random.randint(0,1000000)
        if(savefigs):
            plt.savefig(title+filename+str(hi)+".pdf")
        plt.show()
    
    if(correctUMI):
        return(adjusted_y, newlabels, indices_to_remove)


####PLOT_RIBO_DISTRIBUTIONS_OF_INDIVIDUAL_TRANSCRIPTS####

def plot_trans_dist(riboname, transcript, abundances, numpeaks, fullname, plot, savefigs, start = 1):
    shortname = fullname.split("|")[0]
    
    positions = list(range(start, start + len(abundances)))
    
    plt.figure(figsize = (15,15))
    plt.plot(positions, abundances)
    plt.xlabel("Position", fontsize = 15)
    plt.ylabel("Abundance", fontsize = 15)
    plt.title(shortname, fontsize = 20)

    ####TICKS####
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
    plt.yticks()
    plt.xticks([])
    for value, index in max_values:
        print(f"Position: {value}, Abundance:{index}")
    if(savefigs):
        plt.savefig(riboname+"_"+fullname + ".png")
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
        plt.title(riboname+" "+shortname + " Zoom-in at Highest Peak", fontsize = 15)
        plt.xlabel("Position", fontsize = 11)
        plt.ylabel("Abundance", fontsize = 11)
    
        ####TICKS_FOR_ZOOM####
        
        # max_y = max(peakabundances)
        # plt.yticks([i * (max_y / 20) for i in range(0, 21)])
        plt.locator_params(axis='y', nbins=20)
        plt.xticks()
        plt.yticks([1,2,3,4,5])
        if(savefigs):
            plt.savefig(riboname+" "+fullname+"_Zoom.png")
        if(plot):
            plt.show()

def findoutliers(x, y, labels):
    x_low_percentile = np.percentile(x, 0.03)
    x_high_percentile = np.percentile(x, 99.97)
    y_low_percentile = np.percentile(y, 0.03)
    y_high_percentile = np.percentile(y, 99.97)
    labels_outliers = []
    for i, label in enumerate(labels):

        if ((x[i] > x_high_percentile) or (y[i] > y_high_percentile)):
            labels_outliers.append(label)



    return(labels_outliers)
    
def protein_abundance(allcounts, proteinabundance, line, experiment, transcript_lengths, mRNAabundances, regions, cutoff = 0.2, threshold = 1.5, newlabels = None, indicestoremove = None):
    # Read CSV file using pandas
    df = pd.read_csv(proteinabundance)
    def plotTranscriptVsProtein(transcripts, proteinabundances, experiment):

        transcriptcount = np.array(transcripts)
        allproteinabundance = np.array(proteinabundances)
        indicestoremovee = []
        for i, count in enumerate(transcriptcount):
            if(count == 0 or allproteinabundance[i] == 0):
                indicestoremovee.append(i)
        transcriptcount = [transcriptcount[i] for i in range(len(transcriptcount)) if i not in indicestoremovee]
        allproteinabundance = [allproteinabundance[i] for i in range(len(allproteinabundance)) if i not in indicestoremovee]
        log_transcriptcount = np.log(np.array(transcriptcount))
        data = {
            'log_transcriptcount': log_transcriptcount,
            'allproteinabundance': allproteinabundance,
        }
        

        plt.scatter(log_transcriptcount, allproteinabundance, label=experiment, color = "blue", alpha=0.7)

        slope, intercept, r_value, p_value, std_err = linregress(log_transcriptcount, allproteinabundance)
        trendline = slope * log_transcriptcount + intercept
        plt.plot(log_transcriptcount, trendline, color='red', label=experiment)
        slope = float(slope)
        slope = slope*max(transcriptcount) * max(allproteinabundance)
        print("Slope: ", slope)
        corr, _ = spearmanr(log_transcriptcount, allproteinabundance)

        print(f"Spearman Correlation for {experiment}: {corr}") 

        plt.xlabel('Transcript Counts (Log Transform)')
        plt.ylabel('Protein Abundance')
        if(len(regions) == 0):
            plt.title("Protein Abundance vs Transcript Counts, "+experiment)
            
        else:
            plt.title("Protein Abundance vs Transcript Counts, "+ experiment + " " + regions[0])
        plt.legend()

        plt.show()
    
    proteinabundances = []
    finaltranscriptcount = []
    labels = []
    if(newlabels == None):
        if(len(allcounts) == 1):
            for gene, abundance in zip(df["gene"], df[line]):
                j = allcounts[0].get(gene)
                if(abundance != 0):
                    if j is not None and j != 0:

                        transcript_length = next((int(k.split('\t')[-1]) for k in transcript_lengths if k.split("|")[5] == gene), None)
                        mRNAabundance = next((int(k.split('\t')[2]) for k in mRNAabundances if len(k.split("|")) > 4 and k.split("|")[5] == gene), None)
                        
                        if transcript_length !=0 and mRNAabundance !=0 and transcript_length is not None:
                            finaltranscriptcount.append(float(j)/transcript_length)
                            proteinabundances.append(float(abundance))
                            labels.append(gene)
            plotTranscriptVsProtein(finaltranscriptcount, proteinabundances, experiment[0])
        else:
            transcript_length = [int(k.split('\t')[-1]) for k in transcript_lengths]
            mRNAabundance = [int(k.split('\t')[2]) for k in mRNAabundances]
            firstScreen = []
            PCRcounts = list(allcounts[0].values())
            UMIcounts = list(allcounts[1].values())
            # print("PCRcounts", PCRcounts)
            # print("UMIcounts", UMIcounts)
            for i in range(len(allcounts[0])):
                if(PCRcounts[i]*threshold<UMIcounts[i]):
                    firstScreen.append(i)
            abundancesindex = []
            indicesToCompare = []
            genes = df['gene'].tolist()
            for i in range(len(transcript_length)):
                gene = transcript_lengths[i].split("|")[5]
                if(gene in genes):
                    if(i in firstScreen):
                        indicesToCompare.append(i)
                        abundancesindex.append(genes.index(gene))

            # PCRtranscripts = [(float(PCRcounts[i]) / transcript_length[i] / mRNAabundance[i]) if mRNAabundance[i]>0 else float(PCRcounts[i]/transcript_length[i])for i in indicesToCompare]
            # UMItranscripts = [(float(UMIcounts[i]) / transcript_length[i] / mRNAabundance[i]) if mRNAabundance[i]>0 else float(UMIcounts[i]/transcript_length[i]) for i in indicesToCompare]
            PCRtranscripts = [(float(PCRcounts[i])) for i in indicesToCompare]
            UMItranscripts = [(float(UMIcounts[i])) for i in indicesToCompare]

            proteinabundances = [float(df[line][i]) for i in abundancesindex]

            PCRtrans, PCRabun, correlation1 = plotTranscriptVsProtein(PCRtranscripts, proteinabundances, experiment[0])
            UMItrans, UMIabun, correlation2 = plotTranscriptVsProtein(UMItranscripts, proteinabundances, experiment[1])

    else:
        for gene, abundance in zip(df["gene"], df[line]):
            j = allcounts[0].get(gene)
            if j is not None and j != 0:

                transcript_length = next((int(k.split('\t')[-1]) for i,k in enumerate(transcript_lengths) if k.split("|")[5] == gene and i not in indicestoremove), None)
                mRNAabundance = next((int(k.split('\t')[2]) for i,k in enumerate(mRNAabundances) if len(k.split("|")) > 4 and k.split("|")[5] == gene and i not in indicestoremove), None)
                
                if transcript_length !=0 and mRNAabundance !=0 and transcript_length is not None and mRNAabundance is not None:
                    finaltranscriptcount.append(float(j) / transcript_length / mRNAabundance)
                    proteinabundances.append(float(abundance))
                    labels.append(gene)
        print(allcounts[0])

        plotTranscriptVsProtein(finaltranscriptcount, proteinabundances, experiment[0])

def count_between_datasets(alldatasetcounts, allproteinabundances, additionalexperiments, regions, genes):
    everything = []
    corrpa = []
    corrtc = []
    nocorrpa = []
    nocorrtc = []
    for gene in genes:
        proteinabundances = []
        transcriptcounts = []
        for i, dataset in enumerate(alldatasetcounts):
            transcriptcounts.append(dataset[gene])
            proteinabundances.append(allproteinabundances[i][gene])
        everything.append([transcriptcounts, proteinabundances])
    for i in everything:
        if(not any(0 in sublist for sublist in i)):
            tc1 = i[0][0]
            tc2 = i[0][1]
            pa1 = i[1][0]
            pa2 = i[1][1]
            if(tc1 > tc2 and pa1 > pa2):
                corrpa.append(pa1-pa2)
                corrtc.append(tc1-tc2)
            elif(tc1 < tc2 and pa1 < pa2):
                corrpa.append(pa2-pa1)
                corrtc.append(tc2-tc1)
            elif(tc1 < tc2 and pa1 > pa2):
                nocorrpa.append(pa1-pa2)
                nocorrtc.append(tc2-tc1)
            else:
                nocorrpa.append(pa2-pa1)
                nocorrtc.append(tc1-tc2)
    bar_labels = ['Correlation', 'No Correlation']
    bar_lengths = [len(corrpa), len(nocorrpa)]
    plt.figure(figsize = (10,6))
    plt.bar(bar_labels, bar_lengths, color=['green', 'red'])
    bars = plt.bar(bar_labels, bar_lengths, color=['green', 'red'])

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom')

    plt.title(additionalexperiments[0]+", "+additionalexperiments[1])
    plt.ylabel('Number of Genes')
    plt.show()

def trend_between_datasets(alldatasetcounts, allproteinabundances, additionalexperiments, regions, genes, transcript_lengths):
    
    UMIcounts = [alldatasetcounts[0][gene] for gene in genes]
    PCRcounts = [alldatasetcounts[1][gene] for gene in genes]
    allproteinabundance = np.array([allproteinabundances[0][gene] for gene in genes])
    indicestoremove = []
    tl = {i.split("|")[5]:i.split("|")[-1] for i in transcript_lengths}
    transcript_lengths = [tl[gene] for gene in genes]
    for i, count in enumerate(UMIcounts):
        if(count == 0 or PCRcounts[i] == 0 or allproteinabundance[i] == 0 or transcript_lengths == 0):
            indicestoremove.append(i)
    UMIcounts = [UMIcounts[i] for i in range(len(UMIcounts)) if i not in indicestoremove]
    PCRcounts = [PCRcounts[i] for i in range(len(PCRcounts)) if i not in indicestoremove]
    allproteinabundance = [allproteinabundance[i] for i in range(len(allproteinabundance)) if i not in indicestoremove]
    transcript_lengths = [transcript_lengths[i] for i in range(len(transcript_lengths)) if i not in indicestoremove]

    log_UMIcounts = np.log(np.array(UMIcounts).astype(float) / np.array(transcript_lengths).astype(float))

    log_PCRcounts = np.log(np.array(PCRcounts).astype(float) / np.array(transcript_lengths).astype(float))

    plt.scatter(log_UMIcounts, allproteinabundance, label=additionalexperiments[0], color = "red", alpha=0.3)
    plt.scatter(log_PCRcounts, allproteinabundance, label=additionalexperiments[1], color = "blue", alpha=0.3)

    slope_umi, intercept_umi, r_value_umi, p_value_umi, std_err_umi = linregress(log_UMIcounts, allproteinabundance)
    trendline_umi = slope_umi * log_UMIcounts + intercept_umi
    plt.plot(log_UMIcounts, trendline_umi, color='red', label=additionalexperiments[0])

    slope_pcr, intercept_pcr, r_value_pcr, p_value_pcr, std_err_pcr = linregress(log_PCRcounts, allproteinabundance)
    trendline_pcr = slope_pcr * log_PCRcounts + intercept_pcr
    plt.plot(log_PCRcounts, trendline_pcr, color='blue', label=additionalexperiments[1])

    corr_umi, _ = spearmanr(UMIcounts, allproteinabundance)
    corr_pcr, _ = spearmanr(PCRcounts, allproteinabundance)
    print(f"{additionalexperiments[0]} Spearman correlation coefficient = {corr_umi}")
    print(f"{additionalexperiments[1]} Spearman correlation coefficient = {corr_pcr}")
    print(f"{additionalexperiments[0]} Standard Error = {std_err_umi}")
    print(f"{additionalexperiments[1]} Standard Error = {std_err_pcr}")

    plt.xlabel("Transcript Counts (Log-Transformed, Normalized for Transcript Length)")
    plt.ylabel("Protein Abundance")
    plt.yticks([])
    plt.xticks([])
    plt.title('Transcript Counts v Protein Abundance')
    plt.legend()

    plt.show()

###########################################################################
############################ actual commands ##############################
###########################################################################

def do_regioncounts(
                ribo1,
                riboname,
                ribo2=None,
                riboname2=None,
                region_names=["UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3"],
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                regioncounts=True,
                experiments=[], experiments2=[], 
                filename = ""):
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
        
    RC_graphs_combined(riboname, riboname2, RCcounts, regions=regions, counts2 = RCcounts2, regions2 = regions2, experimentlist=experiments,
                        experimentlist2 = experiments2,  plot=plot, savefigs=savefigs, filename = filename)
    stacked_RC_graphs(riboname, riboname2, RCcounts, experiments, savefigs,
                        plot, counts2 = RCcounts2, experimentlist2 = experiments2, filename = filename)

def do_lengthdist(
                ribo1,
                riboname, 
                ribo2=None,
                riboname2 = None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                experiments=[], experiments2=[], filename = ""
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

        plot_LD_graphs(riboname, riboname2, lengths_list, UTR5counts,
                        experiments, "UTR5 Length Distributions", plot = plot, savefigs=savefigs, filename = filename+"_UTR5")
        plot_LD_graphs(riboname, riboname2, lengths_list, UTR5jcounts,
                        experiments, "UTR5 Junction Length Distributions", plot = plot, savefigs=savefigs, filename = filename+"_UTR5_junction")
        plot_LD_graphs(riboname, riboname2, lengths_list, CDScounts,
                        experiments, "CDS Length Distributions", plot = plot, savefigs=savefigs, filename = filename+"_CDS")
        plot_LD_graphs(riboname, riboname2,lengths_list, UTR3jcounts,
                        experiments, "UTR3 Junction Length Distributions", plot = plot,  savefigs=savefigs, filename = filename+"_UTR3_junction")
        plot_LD_graphs(riboname, riboname2, lengths_list, UTR3counts,
                        experiments, "UTR3 Length Distributions",plot = plot, savefigs=savefigs, filename = filename+"_UTR3")
        plot_LD_graphs(riboname, riboname2, lengths_list, allcounts, 
                       experiments, "Length Distributions", logscale = True, plot = plot, savefigs = savefigs, filename = filename)

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

        plot_LD_graphs(riboname, riboname2, lengths_list, UTR5counts, experiments, "UTR5 Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR5counts_2, experiments2=experiments2, filename = filename+"_UTR5")
        plot_LD_graphs(riboname, riboname2, lengths_list, UTR5jcounts, experiments, "UTR5 Junction Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR5jcounts_2, experiments2=experiments2, filename = filename+"_UTR5_junction")
        plot_LD_graphs(riboname, riboname2, lengths_list, CDScounts, experiments, "CDS Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=CDScounts_2, experiments2=experiments2, filename = filename+"_CDS")
        plot_LD_graphs(riboname, riboname2, lengths_list, UTR3jcounts, experiments, "UTR3 Junction Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR3jcounts_2, experiments2=experiments2, filename = filename+"_UTR3_junction")
        plot_LD_graphs(riboname, riboname2, lengths_list, UTR3counts, experiments, "UTR3 Length Distribution Comparison", plot = plot, savefigs = savefigs, logscale=False,
                    abundances_list2=UTR3counts_2, experiments2=experiments2, filename = filename+"_UTR3")
        
        allcounts2 = [UTR5counts_2, UTR5jcounts_2,
                    CDScounts_2, UTR3jcounts_2, UTR3counts_2]

        plot_LD_graphs(riboname, riboname2, lengths_list, allcounts, experiments, "Compare All Lengths", abundances_list2 = allcounts2, experiments2 = experiments2,
                        logscale=True, plot = plot, savefigs=savefigs, filename = filename)

def do_metagene(ribo1,
                riboname,
                ribo2=None,
                riboname2 = None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                metagene=True, metagene_radius=50,
                experiments=[], experiments2=[],
                filename = "", site_periodicity = ["start", "stop"]
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
        plot_metagene_graphs(riboname, riboname2,
                             meta_radius=metagene_radius,
                                start_abunlist=startarray,
                                stop_abunlist=stoparray,
                                indices=indices,
                                experiments=experiments,
                                plot=plot,
                                savefigs=savefigs,
                                filename = filename,
                                site_periodicity = site_periodicity)
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
        plot_metagene_graphs(riboname = riboname, riboname2= riboname2,
                             meta_radius=metagene_radius,
                                start_abunlist=startarray,
                                stop_abunlist=stoparray,
                                indices=indices,
                                experiments=experiments,
                                plot=plot,
                                savefigs=savefigs,
                                start_abunlist2 = startarray2,
                                stop_abunlist2 = stoparray2,
                                experiments2 = experiments2,
                                filename = filename, 
                                site_periodicity = site_periodicity)

def do_comparetranscriptcounts(
                ribo1,
                riboname,
                ribo2=None,
                riboname2 = None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                regions=[],
                experiments=[], experiments2=[],
                regression = False,
                correctUMI = False,
                annotations = "",
                strength = 0.4,
                filename = ""
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
                    return(scatterplot(riboname, riboname2, allexperimentscounts[i], allexperimentscounts[j], transcriptnames, experiments[i], experiments[j],  experiments[i].split()[-1] + " Versus " + experiments[j].split()[-1] +
                                " Transcript Counts",
                                filename + " " + experiments[i].split("_")[-1] + "_Versus_" + experiments[j].split("_")[-1] + "_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength))
        else:
            for i in range(len(experiments)):
                for j in range(i + 1, len(experiments)):
                    if "UTR5" in regions:
                        return(scatterplot(riboname, riboname2, UTR5tc[i], UTR5tc[j], transcriptnames, experiments[i], experiments[j],
                                    filename + " " + experiments[i].split("_")[-1] + " Versus " +
                                    experiments[j].split("_")[-1] + " UTR5 Transcript Counts",
                                    experiments[i].split("_")[-1] + "_Versus_" + experiments[j].split("_")[-1] + "_UTR5_Transcript Counts", plot, savefigs, regression = regression,correctUMI = correctUMI, strength = strength))
                    if "UTR5_junction" in regions:
                        return(scatterplot(riboname, riboname2, UTR5jtc[i], UTR5jtc[j], transcriptnames, experiments[i], experiments[j],
                                    filename + " " + experiments[i].split("_")[-1] + " Versus " + experiments[j].split("_")[-1] +
                                    " UTR5 junction Transcript Counts",
                                    experiments[i].split("_")[-1] + "_Versus_" + experiments[j].split("_")[-1] + "_UTR5_junction_Transcript Counts", plot, savefigs, regression = regression,correctUMI = correctUMI, strength = strength))
                    if "CDS" in regions:
                        return(scatterplot(riboname, riboname2, CDStc[i], CDStc[j], transcriptnames, experiments[i], experiments[j],
                                    experiments[j].split("_")[0] + " " + experiments[j].split("_")[-1] + " vs " +
                                    experiments[i].split("_")[-1] + " CDS Transcript Counts",
                                    experiments[j].split("_")[-1] + "_Versus_" + experiments[i].split("_")[-1] + "_CDS_Transcript Counts", plot, savefigs, regression = regression,correctUMI = correctUMI, strength = strength))
                    if "UTR3_junction" in regions:
                        return(scatterplot(riboname, riboname2, UTR3jtc[i], UTR3jtc[j], transcriptnames, experiments[i], experiments[j],
                                    filename + " " + experiments[i].split("_")[-1] + " Versus " + experiments[j].split("_")[-1] +
                                    " UTR3 junction Transcript Counts",
                                    experiments[i].split("_")[-1] + "_Versus_" + experiments[j].split("_")[-1] + "_UTR3_junction_Transcript Counts", plot, savefigs, regression = regression,correctUMI = correctUMI, strength = strength))
                    if "UTR3" in regions:
                        return(scatterplot(riboname, riboname2, UTR3tc[i], UTR3tc[j], transcriptnames, experiments[i], experiments[j],
                                    filename + " "+ experiments[i].split("_")[-1] + " Versus " +
                                    experiments[j].split("_")[-1] + " UTR3 Transcript Counts",
                                    experiments[i].split("_")[-1] + "_Versus_" + experiments[j].split("_")[-1] + "_UTR3_Transcript Counts", plot, savefigs, regression = regression,correctUMI = correctUMI, strength = strength))
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
        rawtc = [UTR5transcriptcounts, UTR5jtranscriptcounts, CDStranscriptcounts, UTR3jtranscriptcounts, UTR3transcriptcounts]
        rawtc2 = [UTR5transcriptcounts_2, UTR5jtranscriptcounts_2, CDStranscriptcounts_2, UTR3jtranscriptcounts_2, UTR3transcriptcounts_2]
        tc = [[],[],[],[],[]]
        tc2 = [[],[],[],[],[]]

        for experiment in experiments:
            for i in range(len(rawtc)):
                tc[i].append(rawtc[i].loc[:, experiment].tolist())

        for experiment in experiments2:
            for i in range(len(rawtc2)):
                tc2[i].append(rawtc2[i].loc[:, experiment].tolist())

        if(not len(tc[0][0]) == len(tc2[0][0])):
            
            transcriptnames_2 = UTR5transcriptcounts_2.index.tolist() 
            print("Note: the transcripts are not the same for both ribo files. Only comparing common transcripts.")
            transcriptnames = list(j.split("|")[0] for j in transcriptnames)
            transcriptnames_2 = list(j.split("|")[0] for j in transcriptnames_2)

            temptranscriptnames = []
            temptranscriptnames2 = []
            temptranscriptcounts = [[],[],[],[],[]]
            temptranscriptcounts2 = [[],[],[],[],[]]
            for i in temptranscriptcounts:
                for j in range(len(experiments)):
                    i.append([])
            for i in temptranscriptcounts2:
                for j in range(len(experiments2)):
                    i.append([])

            for i in range(len(transcriptnames_2)):
                if(transcriptnames_2[i] in transcriptnames):
                    temptranscriptnames.append(transcriptnames_2[i])
                    temptranscriptnames2.append(transcriptnames_2[i])
                    index = transcriptnames.index(transcriptnames_2[i])
                    for ind, j in enumerate(temptranscriptcounts):
                        for k in range(len(experiments)):
                            j[k].append(tc[ind][k][index])
                    for ind, j in enumerate(temptranscriptcounts2):
                        for k in range(len(experiments2)):
                            j[k].append(tc2[ind][k][i])

            transcriptnames = temptranscriptnames
            tc2 = temptranscriptcounts2
            tc = temptranscriptcounts

        if(len(regions) == 0):
            allexperimentscounts = []
            allexperimentscounts_2 = []
            for i in range(len(experiments)):
                allexperimentscounts.append([sum(sublist) for sublist in zip(tc[0][i], tc[1][i], tc[2][i], tc[3][i], tc[4][i])])
            for i in range(len(experiments2)):
                allexperimentscounts_2.append([sum(sublist) for sublist in zip(tc2[0][i], tc2[1][i], tc2[2][i], tc2[3][i], tc2[4][i])])
            for i in range(len(experiments)):
                for j in range(len(experiments2)):
                    scatterplot(riboname, riboname2, allexperimentscounts[i], allexperimentscounts_2[j], transcriptnames, experiments[i], experiments2[j] + " (Ribo 2)", experiments[i] + " Versus " + experiments2[j] +
                                " (Ribo 2) Transcript Counts",
                                experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength)
        else:
            for i in range(len(experiments)):
                for j in range(len(experiments2)):
                    if "UTR5" in regions:
                        scatterplot(riboname, riboname2, tc[0][i], tc2[0][j], transcriptnames, experiments[i], experiments2[j]+  "(Ribo 2)",
                                    experiments[i] + " Versus " +
                                    experiments2[j] + " (Ribo 2) UTR5 Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR5_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength)
                    if "UTR5_junction" in regions:
                        scatterplot(riboname, riboname2, tc[1][i], tc2[1][j], transcriptnames, experiments[i], experiments2[j]+ " (Ribo 2)",
                                    experiments[i] + " Versus " + experiments2[j] +
                                    " (Ribo 2) UTR5 junction Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_UTR5_junction_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength)
                    if "CDS" in regions:
                        scatterplot(riboname, riboname2, tc[2][i], tc2[2][j], transcriptnames, experiments[i], experiments2[j]+ " (Ribo 2)",
                                    experiments[i] + " Versus " +
                                    experiments2[j] + " CDS Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR5_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength)
                    if "UTR3" in regions:
                        scatterplot(riboname, riboname2, tc[3][i], tc2[3][j], transcriptnames, experiments[i], experiments2[j] + " (Ribo 2)",
                                    experiments[i] + " Versus " + experiments2[j] +
                                    " (Ribo 2) UTR3 junction Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR3_junction_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength)
                    if "UTR3_junction" in regions:
                        scatterplot(riboname, riboname2, tc[4][i], tc2[4][j], transcriptnames, experiments[i], experiments2[j] + " (Ribo 2)",
                                    experiments[i] + " Versus " +
                                    experiments2[j] + " (Ribo 2) UTR3 Transcript Counts",
                                    experiments[i] + "_Versus_" + experiments2[j] + "_(Ribo 2)_UTR3_Transcript Counts", plot, savefigs, regression = regression, correctUMI = correctUMI, strength = strength)

def do_individualtranscript(
                ribo1,
                riboname,
                ribo2=None,
                riboname2 = None,
                compareribos=False,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                transcriptregions=[],
                numpeaks = 3,
                experiments=[], experiments2=[],
                comparetranscriptcounts=False,
                transcripts=[],
                annotations = [],
                filename = ""
):
    realregions = ["UTR5", "CDS", "UTR3"]
    if(len(transcripts)>0):
        temptranscripts = transcripts
        transcripts = {}
        for region in realregions:
            transcripts[region] = temptranscripts

    if("UTR5_junction" in transcriptregions or "UTR3_junction" in transcriptregions):
        print("UTR5_junction and UTR3_junction are not currently supported for individual transcript analysis.")
        return

    tc = {"UTR5":[], "CDS":[], "UTR3":[]}
    # transcriptnames = None
    for region in realregions:
        region_counts = ribo1.get_region_counts(
        region, sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments)
        for experiment in experiments:
            tc[region].append(region_counts.loc[:, experiment].tolist())
  
    transcriptnames = region_counts.index.tolist()

    if (len(transcripts) == 0 and compareribos):
        tc2 = {"UTR5":[], "CDS":[], "UTR3":[]}
        for region in realregions:
            region_counts = ribo2.get_region_counts(
                region, sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiments2)
            for experiment in experiments2:
                tc2[region].append(region_counts.loc[:, experiment].tolist())
        if(not(len(tc["CDS"][0])==len(tc2["CDS"][0]))):
            print("Note: the transcripts are not the same for both ribo files. Only comparing common transcripts.")
            transcriptnames_2 = region_counts.index.tolist()
            shortenedtranscriptnames = list(j.split("|")[0] for j in transcriptnames)
            shortenedtranscriptnames_2 = list(j.split("|")[0] for j in transcriptnames_2)
            print(transcriptnames_2[0])
            print(transcriptnames[0])
            print(len(transcriptnames_2))
            print(len(transcriptnames))
            temptranscriptnames = []
            realtranscriptnames = []
            realtranscriptnames2 = []
            temptranscriptcounts = {"UTR5":[], "CDS":[], "UTR3":[]}
            temptranscriptcounts2 = {"UTR5":[], "CDS":[], "UTR3":[]}
            for region in realregions:
                for j in range(len(experiments)):
                    temptranscriptcounts[region].append([])
                for j in range(len(experiments2)):
                    temptranscriptcounts2[region].append([])
            for i in range(len(shortenedtranscriptnames_2)):
                if(shortenedtranscriptnames_2[i] in shortenedtranscriptnames):
                    temptranscriptnames.append(shortenedtranscriptnames_2[i])
                    index = shortenedtranscriptnames.index(shortenedtranscriptnames_2[i])
                    for region in realregions:
                        for k in range(len(experiments)):
                            temptranscriptcounts[region][k].append(tc[region][k][index])
                    for region in realregions:
                        for k in range(len(experiments2)):
                            temptranscriptcounts2[region][k].append(tc2[region][k][i])
                    realtranscriptnames.append(transcriptnames[index])
                    realtranscriptnames2.append(transcriptnames_2[i])
            print(len(realtranscriptnames))
            for i, region in enumerate(realregions):
                tc[region] = temptranscriptcounts[region]
                tc2[region] = temptranscriptcounts2[region]
            transcriptnames = realtranscriptnames

    if(not compareribos):
        # We only have one ribo!
        for i in range(len(experiments)):
            for j in range(i + 1, len(experiments)):
                if (len(transcripts)==0):
                    transcripts = {}
                    for region in realregions:
                        transcripts[region] = findoutliers(tc[region][i], tc[region][j], transcriptnames)
                        if(transcripts[region] == []):
                            print("No significant outlier transcripts found. Please specify transcripts to plot.")

                coverage = ribo1.get_coverage(experiments[i],  range_lower=minlength, range_upper=maxlength)
                coverage2 = ribo1.get_coverage(experiments[j], range_lower=minlength, range_upper=maxlength)

                if(len(transcriptregions) == 0):
                    print("Transcripts used in plotting: ")
                    print(transcripts["CDS"])
                    for transcript in transcripts["CDS"]:
                        plot_trans_dist(riboname, transcript, abundances = coverage[transcript], numpeaks = numpeaks,
                                    fullname=experiments[i]+ " "+transcript, plot=plot, savefigs=savefigs)
                        plot_trans_dist(riboname, transcript, abundances = coverage2[transcript], numpeaks = numpeaks,
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
                            plot_trans_dist(riboname, transcript, abundances=coverage[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[i]+ " " + transcript), plot=plot, savefigs=savefigs)
                            plot_trans_dist(riboname, transcript, abundances=coverage2[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[j]+ " " + transcript), plot=plot, savefigs=savefigs)
                    if(len(transcriptregions) == 0):
                        for transcript in transcripts:
                            rows = df[(df['transcript'] == transcript) & (df['feature'] == transcriptregion)]
                            indices = rows[['start', 'end']]
                            indices_values = indices.iloc[0]  
                            startvalue = indices_values['start']
                            endvalue = indices_values['end']
                            plot_trans_dist(riboname, transcript, abundances=coverage[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[i]+ " " + transcript), plot=plot, savefigs=savefigs)
                            plot_trans_dist(riboname, transcript, abundances=coverage2[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion +" "+ experiments[j]+ " " + transcript), plot=plot, savefigs=savefigs)

    else:
        for i in range(len(experiments)):
            for j in range(len(experiments2)):
                if (len(transcripts)==0):
                    transcripts = {}
                    for region in realregions:
                        transcripts[region] = findoutliers(tc[region][i], tc2[region][j], transcriptnames)
                        if(transcripts[region] == []):
                            print("No significant outlier transcripts found. Please specify transcripts to plot.")

                coverage = ribo1.get_coverage(experiments[i],  range_lower = minlength, range_upper = maxlength)
                coverage2 = ribo2.get_coverage(experiments2[j], range_lower = minlength, range_upper = maxlength)

                if(len(transcriptregions) == 0):
                    print("Transcripts used in plotting: ")
                    print(transcripts["CDS"])

                    for transcript in transcripts["CDS"]:
                        plot_trans_dist(riboname, transcript, abundances = coverage[transcript], numpeaks = numpeaks,
                                    fullname = experiments[i] + " " + transcript, plot=plot, savefigs=savefigs)
                        
                        plot_trans_dist(riboname2, realtranscriptnames2[transcriptnames.index(transcript)], abundances = coverage2[realtranscriptnames2[transcriptnames.index(transcript)]], numpeaks = numpeaks, 
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
                            plot_trans_dist(riboname, transcript, abundances=coverage[transcript][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion + " " + transcript), plot=plot, savefigs=savefigs)
                            plot_trans_dist(riboname2, transcript, abundances=coverage2[coverage2[realtranscriptnames2[transcriptnames.index(transcript)]]][startvalue:endvalue], start=startvalue, numpeaks=numpeaks, 
                                            fullname=(transcriptregion + " (Ribo 2) " + transcript), plot=plot, savefigs=savefigs) 

def do_protein_abundance(ribo1, experiments, minlength, maxlength, transcript_lengths, mRNA, regions, proteinabundances, lineA, cutoff, threshold, testdata = [], newlabels = None, indicestoremove = None):
    comparing = []
    with open(transcript_lengths, 'r') as f:
        lines = f.readlines()
    with open(mRNA, 'r') as f:
        mRNAlines = f.readlines()
    if(len(testdata) == 0):
        for experiment in experiments:
        
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
            rawtc = [UTR5transcriptcounts, UTR5jtranscriptcounts, CDStranscriptcounts, UTR3jtranscriptcounts, UTR3transcriptcounts]

            if(regions == []):
                tc = [[],[],[],[],[]]
                allcount = []
                CDStranscriptcounts.to_csv("CDStc_"+experiment+".csv")
                for i in range(len(rawtc)):
                    tc[i].append(rawtc[i].loc[:, experiment].values.tolist())
                allcount = [sum(sublist) for sublist in zip(tc[0][0], tc[1][0], tc[2][0], tc[3][0], tc[4][0])]
            else:
                if "UTR5" in regions:
                    allcount = UTR5transcriptcounts.loc[:, experiment].values.tolist()
                if "UTR5_junction" in regions:
                    allcount = UTR5jtranscriptcounts.loc[:, experiment].values.tolist()
                if "CDS" in regions:
                    allcount = CDStranscriptcounts.loc[:, experiment].values.tolist()
                if "UTR3_junction" in regions:
                    allcount = UTR3jtranscriptcounts.loc[:, experiment].values.tolist()
                if "UTR3" in regions:
                    allcount = UTR3transcriptcounts.loc[:, experiment].values.tolist()
            labels = CDStranscriptcounts.index.tolist()
            allcounts = {}
            for i in range(len(labels)):
                allcounts[labels[i].split("|")[5]] = allcount[i]
            comparing.append(allcounts)
        protein_abundance(comparing, proteinabundances, lineA, experiments, lines, mRNAlines, regions, cutoff, threshold)
    else:
        allcounts = {}
        for i in range(len(newlabels)):
            allcounts[newlabels[i].split("|")[5]] = testdata[i]
        comparing.append(allcounts)
        protein_abundance(comparing, proteinabundances, lineA, experiments, lines, mRNAlines, regions, cutoff, threshold, newlabels = newlabels, indicestoremove = indicestoremove)

def do_trend_between_datasets(ribofile1, additionaldatasets, proteinabundances, experiments, additionalexperiments, minlength, maxlength, regions, lineA, additionallines, transcript_lengths = None):
    alldatasetcounts = []
    allproteinabundances = []
    alllabels = []
    additionaldatasets.append(ribofile1)
    additionalexperiments.append(experiments[0])
    additionallines.append(lineA)
    with open(transcript_lengths, 'r') as f:
        lines = f.readlines()
    for i,dataset in enumerate(additionaldatasets):
        experiment = additionalexperiments[i]
        tempribo = ribopy.Ribo(dataset)
        UTR5transcriptcounts = tempribo.get_region_counts(
            "UTR5", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiment)
        UTR5jtranscriptcounts = tempribo.get_region_counts(
            "UTR5_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiment)
        CDStranscriptcounts = tempribo.get_region_counts(
            "CDS", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiment)
        UTR3jtranscriptcounts = tempribo.get_region_counts(
            "UTR3_junction", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiment)
        UTR3transcriptcounts = tempribo.get_region_counts(
            "UTR3", sum_lengths=True, sum_references=False, range_lower=minlength, range_upper=maxlength, experiments=experiment)
        rawtc = [UTR5transcriptcounts, UTR5jtranscriptcounts, CDStranscriptcounts, UTR3jtranscriptcounts, UTR3transcriptcounts]

        if(regions == []):
            tc = [[],[],[],[],[]]
            for i in range(len(rawtc)):
                tc[i].append(rawtc[i].loc[:, experiment].values.tolist())
            allcount = [sum(sublist) for sublist in zip(tc[0][0], tc[1][0], tc[2][0], tc[3][0], tc[4][0])]

        else:
            if "UTR5" in regions:
                allcount = UTR5transcriptcounts.loc[:, experiment].values.tolist()
            if "UTR5_junction" in regions:
                allcount = UTR5jtranscriptcounts.loc[:, experiment].values.tolist()
            if "CDS" in regions:
                allcount = CDStranscriptcounts.loc[:, experiment].values.tolist()
            if "UTR3_junction" in regions:
                allcount = UTR3jtranscriptcounts.loc[:, experiment].values.tolist()
            if "UTR3" in regions:
                allcount = UTR3transcriptcounts.loc[:, experiment].values.tolist()
        labels = CDStranscriptcounts.index.tolist()
        alllabels.append(labels)
        max_value = max(allcount)
        allcount = [i/max_value for i in allcount]
        allcounts = {}
        for i in range(len(labels)):
            allcounts[labels[i].split("|")[5]] = allcount[i]
        
        
        alldatasetcounts.append(allcounts)
    df = pd.read_csv(proteinabundances)

    genedatabase = set(df["gene"])

    alllabels_set_0 = set(i.split("|")[5] for i in alllabels[0])
    alllabels_set_1 = set(i.split("|")[5] for i in alllabels[1])

    genes = list(genedatabase.intersection(alllabels_set_0, alllabels_set_1))

    if genes:
        allproteinabundances = [{gene: df.loc[df["gene"] == gene, line].values[0] for gene in genes} for line in additionallines]

        trend_between_datasets(alldatasetcounts, allproteinabundances, additionalexperiments, regions, genes, transcript_lengths = lines)
        count_between_datasets(alldatasetcounts, allproteinabundances, additionalexperiments, regions, genes)
    else:
        print("No common genes found.")
    # df = pd.read_csv(proteinabundances)
    # for line in additionallines:
    #     proteinabundances = dict(zip(df["gene"], df[line]))
    #     allproteinabundances.append(proteinabundances)
    
    # genes = []
    # genedatabase = [str(i) for i in df["gene"]]
    # alllabels = [[i.split("|")[5] for i in alllabels[0]], [i.split("|")[5] for i in alllabels[1]]]
    # genes = [gene for gene in genedatabase if gene in alllabels[0] and gene in alllabels[1]]
   
    # trend_between_datasets(alldatasetcounts, allproteinabundances, additionalexperiments, regions, genes)



def ribo_commands(ribofile1,
                ribofile2 = None,
                savefigs=False,
                plot=True,
                minlength=28, maxlength=32,
                experiments=[], experiments2=[],
                regioncounts=True,
                lengthdist=True,
                metagene=True, metagene_radius=50, site_periodicity = ["start", "stop"],
                comparetranscriptcounts=True, regions=[], regression = False,
                individualtranscript=False, transcripts=[], transcriptregions = [], annotations = "", numpeaks = 3,
                correctUMI = False, strength = 0.4,
                proteinabundances = None, transcript_lengths=None, lineA = None, mRNA = None, threshold = 1, cutoff = 1,
                additionaldatasets = [], additionalexperiments = [], additionallines = [],
                filename = "", saveData = False,
                test = False):
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
    riboname = ribo1.info["Reference"]
    riboname2 = ""
    if(ribofile2 is not None):
        ribo2 = ribopy.Ribo(ribofile2)
        riboname2 = ribo2.info["Reference"]
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

    if(test):
        testdata, newlabels, indicestoremove = do_comparetranscriptcounts(ribo1, riboname, ribo2, riboname2, compareribos,
                savefigs, plot, minlength,
                maxlength, regions,
                experiments, experiments2, regression = regression, correctUMI = correctUMI, strength = strength, annotations = annotations, filename = filename)
        
        do_protein_abundance(ribo1, experiments, minlength, maxlength, transcript_lengths, mRNA, regions, proteinabundances, lineA, cutoff, threshold, testdata = testdata, newlabels = newlabels, indicestoremove = indicestoremove)
    else:
        ####PROTEIN_ABUNDANCES####
        if len(additionaldatasets)>0:
            do_trend_between_datasets(ribofile1, additionaldatasets, proteinabundances, experiments, additionalexperiments, minlength, maxlength, regions, lineA, additionallines, transcript_lengths = transcript_lengths)
        else:
            if proteinabundances != None:
                do_protein_abundance(ribo1, experiments, minlength, maxlength, transcript_lengths, mRNA, regions, proteinabundances, lineA, cutoff, threshold)      

    ####LENGTHDIST####
    if lengthdist:
        do_lengthdist(ribo1, riboname, ribo2, riboname2, compareribos,
                   savefigs, plot, minlength,
                   maxlength, experiments, experiments2, filename)
        
    ####REGIONCOUNTS####
    if regioncounts:
        do_regioncounts(ribo1, riboname, ribo2, riboname2, region_names,
                compareribos, savefigs,
                plot, minlength, maxlength,
                regioncounts, experiments, experiments2, filename)

    ####METAGENE####
    if metagene:
        do_metagene(ribo1, riboname, ribo2, riboname2, compareribos,
                 savefigs, plot, minlength,
                 maxlength, metagene, metagene_radius,
                 experiments, experiments2, filename, site_periodicity)

    ####TRANSCRIPT_COUNTS####
    if comparetranscriptcounts:
        do_comparetranscriptcounts(ribo1, riboname, ribo2, riboname2, compareribos,
                                savefigs, plot, minlength,
                                maxlength, regions,
                                experiments, experiments2, regression = regression, correctUMI = correctUMI, strength = strength, annotations = annotations, filename = filename)

    ####INDIVIDUAL_TRANSCRIPT####
    if individualtranscript:
        do_individualtranscript(ribo1, riboname, ribo2, riboname2, compareribos,
                             savefigs, plot, minlength,
                             maxlength, transcriptregions, numpeaks,
                             experiments, experiments2,
                             comparetranscriptcounts, transcripts,
                             annotations, filename)