
# This file was produced by plot-vcfstats, the command line was:
#   plot-vcfstats Samtools.stats -p Samtools.plots
#
# Edit as necessary and recreate the plots by running
#   python plot.py
#
# Title abbreviations:
# 	 0 .. Sampl .. Sample.vcf.gz
# 	 1 .. Samto .. Samtools.vcf.gz
# 	 2 .. Sampl + Samto .. Sample.vcf.gz
#

# Set to 1 to plot in PDF instead of PNG
pdf_plots = 1

# Use logarithimic X axis for allele frequency plots
af_xlog = 0

# Plots to generate, set to 0 to disable
plot_venn_snps = 1
plot_venn_indels = 1
plot_tstv_by_sample = 1
plot_hethom_by_sample = 1
plot_snps_by_sample = 1
plot_indels_by_sample = 1
plot_singletons_by_sample = 1 
plot_depth_by_sample = 1
plot_SNP_count_by_af = 1
plot_Indel_count_by_af = 1
plot_SNP_overlap_by_af = 1
plot_Indel_overlap_by_af = 1
plot_dp_dist = 1
plot_hwe = 1
plot_concordance_by_af = 1
plot_r2_by_af = 1
plot_discordance_by_sample = 1
plot_tstv_by_af = 1
plot_indel_dist = 1
plot_tstv_by_qual = 1
plot_substitutions = 1


# Set to 1 to use sample names for xticks instead of numeric sequential IDs
#   and adjust margins and font properties if necessary
sample_names   = 0
sample_margins = {'right':0.98, 'left':0.07, 'bottom':0.2}
sample_font    = {'rotation':45, 'ha':'right', 'fontsize':8}

if sample_names==0: sample_margins=(); sample_font=();


#-------------------------------------------------


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

import numpy
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1: raise ValueError("The function 'smooth' only accepts 1 dimension arrays.")
	if x.size < window_len: return x
	if window_len<3: return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']: raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
	s = numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
	if window == 'flat': # moving average
		w = numpy.ones(window_len,'d')
	else:
		w = eval('numpy.'+window+'(window_len)')
	y = numpy.convolve(w/w.sum(),s,mode='valid')
	return y[(window_len/2-1):-(window_len/2)]



if plot_venn_snps:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.bar([1,2,3],[1000,0,0],align='center',color='orange',width=0.3)
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(0.5,3.5)
	plt.xticks([1,2,3],('Sampl\nts/tv 0.52\n1,000','Sampl + Samto\nts/tv 0.00\n0','Samto\nts/tv 0.00\n0'))
	plt.title('Number of SNPs')
	plt.subplots_adjust(right=0.95,bottom=0.15)
	plt.savefig('venn_bars.snps.png')
	if pdf_plots: plt.savefig('venn_bars.snps.pdf')
	plt.close()


if plot_venn_indels:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.bar([1,2,3],[1261,739,747],align='center',color='red',width=0.3)
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(0.5,3.5)
	plt.xticks([1,2,3],('Sampl\n1,261','Sampl + Samto\n739','Samto\n747'))
	plt.title('Number of indels')
	plt.subplots_adjust(right=0.95,bottom=0.15)
	plt.savefig('venn_bars.indels.png')
	if pdf_plots: plt.savefig('venn_bars.indels.pdf')
	plt.close()



dat = {}
with open('counts_by_af.snps.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] == '#': continue
		id = int(row[0])
		if id not in dat: dat[id] = []
		dat[id].append([float(row[1]),float(row[2])])

if plot_SNP_count_by_af:
	fig = plt.figure(figsize=(2*4.33070866141732,3.93700787401575*0.7))
	ax1 = fig.add_subplot(111)
	ax1.set_ylabel('Number of sites')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_yscale('log')
	if af_xlog: ax1.set_xscale('log')
	ax1.set_xlabel('Non-reference allele frequency')
	ax1.set_xlim(-0.05,1.05)
	has_data = 0

	if 0 in dat and len(dat[0])>2:
		ax1.plot([row[0] for row in dat[0]], [row[1] for row in dat[0]], '-o',markersize=3, color='orange',mec='orange',label='Sampl')
		has_data = 1

	if 1 in dat and len(dat[1])>2:
		ax1.plot([row[0] for row in dat[1]], [row[1] for row in dat[1]], '-o',markersize=3, color='red',mec='red',label='Samto')
		has_data = 1

	if 2 in dat and len(dat[2])>2:
		ax1.plot([row[0] for row in dat[2]], [row[1] for row in dat[2]], '-o',markersize=3, color='darkgreen',mec='darkgreen',label='Sampl + Samto')
		has_data = 1

	if has_data:
		ax1.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)
		plt.title('SNP count by AF')
		plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)
		plt.savefig('counts_by_af.snps.png')
		if pdf_plots: plt.savefig('counts_by_af.snps.pdf')
		plt.close()




dat = {}
with open('counts_by_af.indels.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] == '#': continue
		id = int(row[0])
		if id not in dat: dat[id] = []
		dat[id].append([float(row[1]),float(row[2])])

if plot_Indel_count_by_af:
	fig = plt.figure(figsize=(2*4.33070866141732,3.93700787401575*0.7))
	ax1 = fig.add_subplot(111)
	ax1.set_ylabel('Number of sites')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_yscale('log')
	if af_xlog: ax1.set_xscale('log')
	ax1.set_xlabel('Non-reference allele frequency')
	ax1.set_xlim(-0.05,1.05)
	has_data = 0

	if 0 in dat and len(dat[0])>2:
		ax1.plot([row[0] for row in dat[0]], [row[1] for row in dat[0]], '-o',markersize=3, color='orange',mec='orange',label='Sampl')
		has_data = 1

	if 1 in dat and len(dat[1])>2:
		ax1.plot([row[0] for row in dat[1]], [row[1] for row in dat[1]], '-o',markersize=3, color='red',mec='red',label='Samto')
		has_data = 1

	if 2 in dat and len(dat[2])>2:
		ax1.plot([row[0] for row in dat[2]], [row[1] for row in dat[2]], '-o',markersize=3, color='darkgreen',mec='darkgreen',label='Sampl + Samto')
		has_data = 1

	if has_data:
		ax1.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)
		plt.title('Indel count by AF')
		plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)
		plt.savefig('counts_by_af.indels.png')
		if pdf_plots: plt.savefig('counts_by_af.indels.pdf')
		plt.close()




dat = []

with open('overlap_by_af.snps.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append(row)

if plot_SNP_overlap_by_af and len(dat)>1:
	fig = plt.figure(figsize=(2*4.33070866141732,3.93700787401575*0.7))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[0] for row in dat], [row[1] for row in dat],'-o',markersize=3, color='red',mec='red')
	ax1.set_ylabel('Fraction found in Samto [%]')
	ax1.set_xscale('log')
	ax1.set_xlabel('Non-reference allele frequency in Sampl')
	ax1.set_xlim(0,1.01)
	plt.title('SNP overlap by AF')
	plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)
	plt.savefig('overlap_by_af.snps.png')
	if pdf_plots: plt.savefig('overlap_by_af.snps.pdf')
	plt.close()


dat = []

with open('overlap_by_af.indels.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append(row)

if plot_Indel_overlap_by_af and len(dat)>1:
	fig = plt.figure(figsize=(2*4.33070866141732,3.93700787401575*0.7))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[0] for row in dat], [row[1] for row in dat],'-o',markersize=3, color='red',mec='red')
	ax1.set_ylabel('Fraction found in Samto [%]')
	ax1.set_xscale('log')
	ax1.set_xlabel('Non-reference allele frequency in Sampl')
	ax1.set_xlim(0,1.01)
	plt.title('Indel overlap by AF')
	plt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)
	plt.savefig('overlap_by_af.indels.png')
	if pdf_plots: plt.savefig('overlap_by_af.indels.pdf')
	plt.close()


dat = []
with open('tstv_by_af.0.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])


if plot_tstv_by_af and len(dat)>2:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[0] for row in dat], [row[1] for row in dat], '-o',color='k',mec='k',markersize=3)
	ax1.set_ylabel('Number of sites',color='k')
	ax1.set_yscale('log')
	#ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	for tl in ax1.get_yticklabels(): tl.set_color('k')
	ax1.set_xlabel('Non-ref allele frequency')
	ax2 = ax1.twinx()
	ax2.plot([row[0] for row in dat], [row[2] for row in dat], '-o',color='orange',mec='orange',markersize=3)
	ax2.set_ylabel('Ts/Tv',color='orange')
	ax2.set_ylim(0,0.5+max(3,max(row[2] for row in dat)))
	ax1.set_xlim(0,1)
	for tl in ax2.get_yticklabels(): tl.set_color('orange')
	plt.subplots_adjust(right=0.88,left=0.15,bottom=0.11)
	plt.title('Sampl')
	plt.savefig('tstv_by_af.0.png')
	if pdf_plots: plt.savefig('tstv_by_af.0.pdf')
	plt.close()



dat = []
with open('tstv_by_qual.0.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_tstv_by_qual and len(dat)>2:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[1] for row in dat], [row[2] for row in dat], '^-', ms=3, mec='orange', color='orange')
	ax1.set_ylabel('Ts/Tv',fontsize=10)
	ax1.set_xlabel('Number of sites\n(sorted by QUAL, descending)',fontsize=10)
	ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')
	ax1.set_ylim(min(2,min(row[2] for row in dat))-0.3,0.3+max(2.2,max(row[2] for row in dat)))

	plt.subplots_adjust(right=0.88,left=0.15,bottom=0.15)
	plt.title('Sampl')
	plt.savefig('tstv_by_qual.0.png')
	if pdf_plots: plt.savefig('tstv_by_qual.0.pdf')
	plt.close()



dat = []
with open('indels.0.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_indel_dist and len(dat)>0:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.bar([row[0]-0.5 for row in dat], [row[1] for row in dat], color='orange')# , edgecolor='orange')
	ax1.set_xlabel('InDel Length')
	ax1.set_ylabel('Count')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(-20,20)
	plt.subplots_adjust(bottom=0.17)
	plt.title('Sampl')
	plt.savefig('indels.0.png')
	if pdf_plots: plt.savefig('indels.0.pdf')
	plt.close()

dat = [
	[0,'A>C',93],
	[1,'A>G',107],
	[2,'A>T',87],
	[3,'C>A',81],
	[4,'C>G',74],
	[5,'C>T',67],
	[6,'G>A',63],
	[7,'G>C',75],
	[8,'G>T',68],
	[9,'T>A',93],
	[10,'T>C',104],
	[11,'T>G',88],
]

if plot_substitutions:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	cm  = mpl.cm.get_cmap('autumn')
	n = 12
	col = []
	for i in list(range(n)): col.append(cm(1.*i/n))
	ax1 = fig.add_subplot(111)
	ax1.bar([row[0] for row in dat], [row[2] for row in dat], color=col)
	ax1.set_ylabel('Count')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(-0.5,n+0.5)
	plt.xticks([row[0] for row in dat],[row[1] for row in dat],rotation=45)
	plt.title('Sampl')
	plt.savefig('substitutions.0.png')
	if pdf_plots: plt.savefig('substitutions.0.pdf')
	plt.close()



dat = []
with open('tstv_by_af.1.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])


if plot_tstv_by_af and len(dat)>2:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[0] for row in dat], [row[1] for row in dat], '-o',color='k',mec='k',markersize=3)
	ax1.set_ylabel('Number of sites',color='k')
	ax1.set_yscale('log')
	#ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	for tl in ax1.get_yticklabels(): tl.set_color('k')
	ax1.set_xlabel('Non-ref allele frequency')
	ax2 = ax1.twinx()
	ax2.plot([row[0] for row in dat], [row[2] for row in dat], '-o',color='red',mec='red',markersize=3)
	ax2.set_ylabel('Ts/Tv',color='red')
	ax2.set_ylim(0,0.5+max(3,max(row[2] for row in dat)))
	ax1.set_xlim(0,1)
	for tl in ax2.get_yticklabels(): tl.set_color('red')
	plt.subplots_adjust(right=0.88,left=0.15,bottom=0.11)
	plt.title('Samto')
	plt.savefig('tstv_by_af.1.png')
	if pdf_plots: plt.savefig('tstv_by_af.1.pdf')
	plt.close()



dat = []
with open('tstv_by_qual.1.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_tstv_by_qual and len(dat)>2:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[1] for row in dat], [row[2] for row in dat], '^-', ms=3, mec='red', color='red')
	ax1.set_ylabel('Ts/Tv',fontsize=10)
	ax1.set_xlabel('Number of sites\n(sorted by QUAL, descending)',fontsize=10)
	ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')
	ax1.set_ylim(min(2,min(row[2] for row in dat))-0.3,0.3+max(2.2,max(row[2] for row in dat)))

	plt.subplots_adjust(right=0.88,left=0.15,bottom=0.15)
	plt.title('Samto')
	plt.savefig('tstv_by_qual.1.png')
	if pdf_plots: plt.savefig('tstv_by_qual.1.pdf')
	plt.close()



dat = []
with open('indels.1.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_indel_dist and len(dat)>0:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.bar([row[0]-0.5 for row in dat], [row[1] for row in dat], color='orange')# , edgecolor='orange')
	ax1.set_xlabel('InDel Length')
	ax1.set_ylabel('Count')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(-20,20)
	plt.subplots_adjust(bottom=0.17)
	plt.title('Samto')
	plt.savefig('indels.1.png')
	if pdf_plots: plt.savefig('indels.1.pdf')
	plt.close()

dat = [
	[0,'A>C',0],
	[1,'A>G',0],
	[2,'A>T',0],
	[3,'C>A',0],
	[4,'C>G',0],
	[5,'C>T',0],
	[6,'G>A',0],
	[7,'G>C',0],
	[8,'G>T',0],
	[9,'T>A',0],
	[10,'T>C',0],
	[11,'T>G',0],
]

if plot_substitutions:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	cm  = mpl.cm.get_cmap('autumn')
	n = 12
	col = []
	for i in list(range(n)): col.append(cm(1.*i/n))
	ax1 = fig.add_subplot(111)
	ax1.bar([row[0] for row in dat], [row[2] for row in dat], color=col)
	ax1.set_ylabel('Count')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(-0.5,n+0.5)
	plt.xticks([row[0] for row in dat],[row[1] for row in dat],rotation=45)
	plt.title('Samto')
	plt.savefig('substitutions.1.png')
	if pdf_plots: plt.savefig('substitutions.1.pdf')
	plt.close()



dat = []
with open('depth.1.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_dp_dist:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[0] for row in dat], [row[2] for row in dat], '-^', color='k')
	ax1.set_ylabel('Number of genotypes [%]',color='k')
	ax1.set_xlabel('Depth')
	ax2 = ax1.twinx()
	ax2.plot([row[0] for row in dat], [row[1] for row in dat], '-o', color='red')
	ax2.set_ylabel('Cumulative number of genotypes [%]',color='red')
	for tl in ax2.get_yticklabels(): tl.set_color('red')
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.87)
	plt.title('Samto')
	plt.savefig('depth.1.png')
	if pdf_plots: plt.savefig('depth.1.pdf')
	plt.close()



dat = []
with open('tstv_by_af.2.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])


if plot_tstv_by_af and len(dat)>2:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[0] for row in dat], [row[1] for row in dat], '-o',color='k',mec='k',markersize=3)
	ax1.set_ylabel('Number of sites',color='k')
	ax1.set_yscale('log')
	#ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	for tl in ax1.get_yticklabels(): tl.set_color('k')
	ax1.set_xlabel('Non-ref allele frequency')
	ax2 = ax1.twinx()
	ax2.plot([row[0] for row in dat], [row[2] for row in dat], '-o',color='darkgreen',mec='darkgreen',markersize=3)
	ax2.set_ylabel('Ts/Tv',color='darkgreen')
	ax2.set_ylim(0,0.5+max(3,max(row[2] for row in dat)))
	ax1.set_xlim(0,1)
	for tl in ax2.get_yticklabels(): tl.set_color('darkgreen')
	plt.subplots_adjust(right=0.88,left=0.15,bottom=0.11)
	plt.title('Sampl + Samto')
	plt.savefig('tstv_by_af.2.png')
	if pdf_plots: plt.savefig('tstv_by_af.2.pdf')
	plt.close()



dat = []
with open('tstv_by_qual.2.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_tstv_by_qual and len(dat)>2:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.plot([row[1] for row in dat], [row[2] for row in dat], '^-', ms=3, mec='darkgreen', color='darkgreen')
	ax1.set_ylabel('Ts/Tv',fontsize=10)
	ax1.set_xlabel('Number of sites\n(sorted by QUAL, descending)',fontsize=10)
	ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')
	ax1.set_ylim(min(2,min(row[2] for row in dat))-0.3,0.3+max(2.2,max(row[2] for row in dat)))

	plt.subplots_adjust(right=0.88,left=0.15,bottom=0.15)
	plt.title('Sampl + Samto')
	plt.savefig('tstv_by_qual.2.png')
	if pdf_plots: plt.savefig('tstv_by_qual.2.pdf')
	plt.close()



dat = []
with open('indels.2.dat', 'r') as f:
	reader = csv.reader(f, 'tab')
	for row in reader:
		if row[0][0] != '#': dat.append([float(x) for x in row])

if plot_indel_dist and len(dat)>0:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	ax1 = fig.add_subplot(111)
	ax1.bar([row[0]-0.5 for row in dat], [row[1] for row in dat], color='orange')# , edgecolor='orange')
	ax1.set_xlabel('InDel Length')
	ax1.set_ylabel('Count')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(-20,20)
	plt.subplots_adjust(bottom=0.17)
	plt.title('Sampl + Samto')
	plt.savefig('indels.2.png')
	if pdf_plots: plt.savefig('indels.2.pdf')
	plt.close()

dat = [
	[0,'A>C',0],
	[1,'A>G',0],
	[2,'A>T',0],
	[3,'C>A',0],
	[4,'C>G',0],
	[5,'C>T',0],
	[6,'G>A',0],
	[7,'G>C',0],
	[8,'G>T',0],
	[9,'T>A',0],
	[10,'T>C',0],
	[11,'T>G',0],
]

if plot_substitutions:
	fig = plt.figure(figsize=(4.33070866141732,3.93700787401575))
	cm  = mpl.cm.get_cmap('autumn')
	n = 12
	col = []
	for i in list(range(n)): col.append(cm(1.*i/n))
	ax1 = fig.add_subplot(111)
	ax1.bar([row[0] for row in dat], [row[2] for row in dat], color=col)
	ax1.set_ylabel('Count')
	ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	ax1.set_xlim(-0.5,n+0.5)
	plt.xticks([row[0] for row in dat],[row[1] for row in dat],rotation=45)
	plt.title('Sampl + Samto')
	plt.savefig('substitutions.2.png')
	if pdf_plots: plt.savefig('substitutions.2.pdf')
	plt.close()

