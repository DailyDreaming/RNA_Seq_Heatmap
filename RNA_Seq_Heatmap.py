###############################################################################
# Filename: RNA_Seq_Heatmap.py                                                #
# Written by: Lon Blauvelt                                                    #
# Data, stylesheet, and guidance Provided by Prof. Christopher Vollmers       #
#                                                                             #
# Creates a heatmap of RNA-seq expression over 24 hours (2hr intervals)       #
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
import math
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import matplotlib.patches as mplpatches

# Location:
# /home/lifeisaboutfishtacos/.local/lib/python3.5/site-packages/matplotlib/mpl-data/stylelib
plt.style.use('vollmers')

fig_width      = 2
fig_height     = 3
panel1_width    = 0.4
panel1_height   = 0.7
panel1_x_margin = 0.2
panel1_y_margin = 0.2
panel_spacing  = 0.025

plt.figure(figsize=(fig_width, fig_height))

# Yellow = (1,1,0)
# Blue   = (0,0,1)
R = np.linspace(255 / 255, 56 / 255, 101)
G = np.linspace(225 / 255, 66 / 255, 101)
B = np.linspace(40 / 255, 157 / 255, 101)

###############################################################################
# Read data and prepare it for plotting.                                      #
# Data File: data.txt                                                         #
#                                                                             #
# Columns contained within the file:                                          #
# 1.  Gene_symbol            (Example: 'Narf')                                #
# 2.  Ensembl_ID             (Example: 'ENSMUSG00000000056')                  #
# 3.  Gene_position          (Example: 'chr11:121098566-121117170')           #
# 4.  Gene_type              (Example: 'protein_coding')                      #
# 5.  FPKM_CT0               (Example: '221313')                              #
# 6.  FPKM_CT3               (Example: '236590')                              #
# 7.  FPKM_CT6               (Example: '434754')                              #
# 8.  FPKM_CT9               (Example: '850649')                              #
# 9.  FPKM_CT12              (Example: '947489')                              #
# 10. FPKM_CT15              (Example: '268738')                              #
# 11. FPKM_CT18              (Example: '180969')                              #
# 12. FPKM_CT21              (Example: '266168')                              #
# 13. Period_length(h)       (Example: '23.0')                                #
# 14. Peak_phase(CT)         (Example: '11.6')                                #
# 15. pMMC-beta              (Example: '0.0096')                              #
# 16. Amplitude (fold)       (Example: '5.2')                                 #
###############################################################################

FPKM_list = []
phaseCT_list = []
first_line = 0

splice_sequence_filename = 'data.txt'

for line in open(splice_sequence_filename):
    if first_line is not 0:
        split_line = line.strip().split('\t')
        FPKM_list.append([float(split_line[13]), \
                          int(split_line[4]), \
                          int(split_line[5]), \
                          int(split_line[6]), \
                          int(split_line[7]), \
                          int(split_line[8]), \
                          int(split_line[9]), \
                          int(split_line[10]), \
                          int(split_line[11])])
        phaseCT_list.append(int(float(split_line[13]) / 2) * 2)
    if first_line is 0:
        first_line = 1

FPKM_list.sort(reverse=True, key=lambda x: x[0])
phaseCT_set = set(phaseCT_list)

###############################################################################
# Plot the Gene "Heatmap"                                                     #
###############################################################################

panel1 = plt.axes([panel1_x_margin, panel1_y_margin, panel1_width, panel1_height], frameon=True)

y = 0
for genes in FPKM_list:
    x = 0
    gene_array = np.array([genes[1], genes[2], genes[3], genes[4], genes[5], genes[6], genes[7], genes[8]])
    normalized = ((gene_array - min(gene_array)) / (max(gene_array) - min(gene_array))) * 100
    for point in normalized:
        rectangle = mplpatches.Rectangle([x, y], 2, 1,facecolor=(R[int(point)], G[int(point)], B[int(point)]), linewidth=0)
        panel1.add_patch(rectangle)
        x = x + 2
    y = y + 1

# Set Axes Limits
panel1.set_xlim([0, 16])
panel1.set_ylim([0, 1265])

# Set Axes Ticks
panel1.set_xticks([1,3,5,7,9,11,13,15])
panel1.set_yticks([0,200,400,600,800,1000,1200])

# Set Axes Tick Labels
panel1.set_xticklabels([0,'',6,'',12,'',18,''])
panel1.set_yticklabels([0,200,400,600,800,1000,1200])

# Set Axes Labels
panel1.set_ylabel("Number of Genes")
panel1.set_xlabel("CT")

# Turn On or Off Plot Axes
panel1.tick_params(axis='both', which='both', \
                   bottom='on', labelbottom='on', \
                   left='on', labelleft='on', \
                   right='off', labelright='off', \
                   top='off', labeltop='off')

plt.savefig('RNA_Seq_Heatmap.pdf')
