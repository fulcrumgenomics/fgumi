#################################################################################
# The MIT License                                                               #
# Copyright (c) 2017 Fulcrum Genomics LLC                                       #
# Permission is hereby granted, free of charge, to any person obtaining a copy  #
# of this software and associated documentation files (the "Software"), to deal #
# in the Software without restriction, including without limitation the rights  #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
# copies of the Software, and to permit persons to whom the Software is         #
# furnished to do so, subject to the following conditions:                      #
# The above copyright notice and this permission notice shall be included in    #
# all copies or substantial portions of the Software.                           #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     #
# THE SOFTWARE.                                                                 #
#################################################################################

#################################################################################
# R script to generate QC plots for simplex sequencing data QC'd with
# CollectSimplexSeqMetrics tools.  All arguments are positional and required:
#   1.  The family size file with data on families by start/stop and ss
#   2.  The simplex yield metrics file
#   3.  The UMI metrics file
#   4.  An output PDF to write
#   5+. One or more strings to use in plot titles to ID the sample/dataset
#################################################################################

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)
library(scales)  # Require by ggplot

# Standard colors
fgblue  = "#2E9DD7"
fggreen = "#155936"
fgred   = "firebrick3"
fgcolors = c(fgblue, fggreen, fgred)

args           = commandArgs(trailingOnly=T)
familyData     = read.table(args[1], sep="\t", header=T)
yieldData      = read.table(args[2], sep="\t", header=T)
umiData        = read.table(args[3], sep="\t", header=T)
outputFile     = args[4]
sampleInfo     = paste(args[5:length(args)])

pdf(outputFile, width=11, height=8.5)

# Plot #1 - Family Size Distributions (CS and SS)
ggplot(familyData) + aes(x=family_size) +
  geom_line(aes(y=ss_count, color="SS Families")) +
  geom_line(aes(y=cs_count, color="By Coord+Strand")) +
  scale_x_continuous(trans="log2", minor_breaks=seq(0,max(familyData$family_size), by=2)) +
  scale_color_manual(values=c("SS Families"=fgblue, "By Coord+Strand"=fggreen)) +
  labs(x="Family Size (log2 scaled)", y="Count of Families", title=paste("Family Size Distributions for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #2 - Cumulative Family Size Distributions
ggplot(familyData) + aes(x=family_size) +
  geom_line(aes(y=ss_fraction_gt_or_eq_size, color="SS Families"), alpha=0.5) +
  geom_line(aes(y=cs_fraction_gt_or_eq_size, color="By Coord+Strand"), alpha=0.5) +
  scale_x_continuous(trans="log2", minor_breaks=seq(0,max(familyData$family_size), by=2)) +
  scale_color_manual(values=c("SS Families"=fgblue, "By Coord+Strand"=fggreen)) +
  labs(x="Family Size (log2 scaled)", y="Fraction of Families at >= Family Size", title=paste("Cumulative Family Size Distributions for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #3 - Simplex Yield
ggplot(yieldData) +
  aes(x=read_pairs) +
  geom_area(aes(y=ss_families, fill="SS Families")) +
  geom_area(aes(y=ss_consensus_families, fill="Consensus-Eligible")) +
  scale_fill_manual(values=c("SS Families"=fgblue, "Consensus-Eligible"=fggreen)) +
  labs(x="Read Pairs", y="Count of SS Families", title=paste("Simplex Yield by Input Read Pairs for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #4 - Mean SS Family Size by Read Pairs
ggplot(yieldData) +
  aes(x=read_pairs, y=mean_ss_family_size) +
  geom_line(color=fgblue) +
  labs(x="Read Pairs", y="Mean SS Family Size", title=paste("Mean SS Family Size by Input Read Pairs for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #5 - Singleton Fraction by Read Pairs
ggplot(yieldData) +
  aes(x=read_pairs, y=ss_singleton_fraction) +
  geom_line(color=fgred) +
  ylim(0, 1) +
  labs(x="Read Pairs", y="Singleton Fraction", title=paste("Singleton Fraction by Input Read Pairs for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #6 - UMI distribution (for UMIs that don't contain any no-calls)
ggplot(subset(umiData, !grepl("N", umiData$umi, fixed=T))) +
  aes(x=raw_observations, y=unique_observations) +
  geom_point(color=fggreen) +
  labs(x="Observations of UMI in Raw Reads", y="Unique Observations (Tag Families w/UMI)", title=paste("UMI Representation in", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #7 - Distribution of reads within tag families
ggplot(familyData) + aes(x=family_size) +
  geom_line(aes(y=cs_count * family_size, color="By Coord+Strand")) +
  geom_line(aes(y=ss_count * family_size, color="SS Families")) +
  scale_x_continuous(trans="log2", minor_breaks=seq(0,max(familyData$family_size), by=2)) +
  scale_color_manual(values=c("By Coord+Strand"=fgblue, "SS Families"=fggreen)) +
  labs(x="Family Size (log2 scaled)", y="Reads Allocated to Families of Size N", title=paste("Read Distribution Among Families for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

dev.off()
