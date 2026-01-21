import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np


def scatter():
	sns.set(style='whitegrid')
	dataframe = pd.read_csv('tabla_positives_negatives_hel.csv', sep=";") #tabla contacts

	sns.displot(x='Helicity', y='Total contacts', hue='Type',data=dataframe)

	def annotate(data, **kws):
	    r, p = sp.stats.pearsonr(data['Total contacts/length'], data['Mass'])
	    ax = plt.gca()
	    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
	            transform=ax.transAxes)

	#g.map_dataframe(annotate) #get r and p values
	plt.xlabel("Helicity")
	plt.ylabel("Average contacts")
	plt.show()
	#plt.savefig("./hel_positives_negatives_displot.png", dpi=300)


def lmplot():
	sns.set(style='whitegrid')
	dataframe = pd.read_csv('pd_dataframe_csv.csv', sep=";") #tabla contacts

	g=sns.lmplot(x='H', y='8a9l_2_avg',data=dataframe)

	def annotate(data, **kws):
	    r, p = sp.stats.pearsonr(data['H'], data['8a9l_2_avg'])
	    ax = plt.gca()
	    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
	            transform=ax.transAxes)

	g.map_dataframe(annotate) #get r and p values
	
	plt.xlabel("H")
	plt.ylabel("Average contacts")
	#plt.show()
	plt.savefig("./properties/H_Helicity.png", dpi=300)


def barplot():
	dataframe = pd.read_csv('pd_dataframe_csv.csv', sep=";") #tabla ALL contacts
	plt.figure(figsize=(15, 7))
	#cols = ['red' if x == 'EXTER001' else 'green' if x == 'ASYNP113' else 'skyblue' for x in dataframe.sort_values('8a9l_2_tot',ascending = False).ID]
	cols = ['skyblue' if x == 'EXTER001' else 'green' if x == 'ASYNP025' else 'skyblue' for x in dataframe.sort_values('8a9l_2_tot',ascending = False).ID]
	sns.barplot(x='ID', y='8a9l_2_tot', data=dataframe, order=dataframe.sort_values('8a9l_2_tot',ascending = False).ID, palette=cols)
	plt.xticks(rotation=90, fontsize=5)
	

	plt.xlabel("Peptide ID")
	plt.ylabel("Average contacts")
	plt.savefig("./barplot/8a9l_2_tot_barplot_peptides_total_ID_ll37.png", dpi=300)

	#plt.show()

def violin():
	tot_contacts_with=np.load("./data_rg/rg_asy_ct_with_peptide.npy")
	tot_contacts_alone=np.load("./data_rg/rg_asy_ct_alone.npy")
	
	g=sns.violinplot(data=[tot_contacts_with, tot_contacts_alone])
	#g.set_xticks(2) # <--- set the ticks first
	g.set_xticklabels(['With peptide','Withoutout peptide'])
	plt.xlabel("Simulation")
	plt.ylabel("Rg aSyn[100:140]")
	#plt.title("Avearge Rg ll37")

	from scipy.stats import mannwhitneyu, normaltest
	from statannotations.Annotator import Annotator
	pvalue_score = mannwhitneyu(tot_contacts_with, tot_contacts_alone, alternative="two-sided")
	print(pvalue_score[1])

	#plt.show()
	plt.savefig("rg_asyn_2.png", dpi=300)

def rg_distr():
	rg_hist_0=np.load("./data_rg/rg_hist[0]_asy_ct.npy")
	rg_hist_1=np.load("./data_rg/rg_hist[1]_asy_ct.npy")
	rg_hist_0_alone=np.load("./data_rg/rg_hist[0]_asy_ct_alone.npy")
	rg_hist_1_alone=np.load("./data_rg/rg_hist[1]_asy_ct_alone.npy")

	plt.plot(rg_hist_0,rg_hist_1, 'b')
	plt.plot(rg_hist_0_alone,rg_hist_1_alone, 'r')
	plt.show()



#scatter()
barplot()
#violin()
#rg_distr()
#lmplot()