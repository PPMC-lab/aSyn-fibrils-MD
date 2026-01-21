import matplotlib.pyplot as plt
import numpy as np
import  sys
from matplotlib import colors


n_chains_asyn=209
n_chains_pep=105
len_ll37=38
n_frames=9000

data={}
#Load single np array
#data_single=np.load("./data_8a9l_180_all_negatives_hel/8a9l_180_all_negatives_hel_tot_contact_pep_chunk_1000.npy")

#Load independent .npy chunks
for i in range(1000, 10000, 500):
	data[i]=np.load("./data_chunk_ll37/8a9l_ll37_1000ns_cubic_contact_map.npy_chunk_"+str(i)+".npy")

#Generate final array (either sum or append)
#data_append=np.concatenate((data[1000], data[2000],data[3000],data[4000],data[5000],data[6000],data[7000],data[8000],data[9000]),axis = 0)
#data_append=np.concatenate((data[1000], data[5500]),axis = 0)
data_sum=np.sum((data[1000], data[1500], data[2000], data[2500], data[3000], data[3500], data[4000], data[4500], data[5000], data[5500],data[6000], data[6500],data[7000], data[7500],data[8000],data[8500],data[9000], data[9500]), axis=0) #summing n chunks chunks
#data_sum=np.sum((data[1000], data[5500]), axis=0) #summing n chunks chunks
#data_min=np.mean((data[1000], data[2000], data[3000], data[4000]), axis=0)

def array_2d_cmap(data):
	#print(data.shape)
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	fig, ax = plt.subplots()

	im=plt.imshow(data, cmap='hot', origin='lower')
	#ax.set_xticklabels(custom_tick_xlabels)
	plt.xlabel("# res aSyn")
	plt.ylabel("# res LL37")

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)

	#Colorbar
	plt.colorbar(im,cax=cax, label="p(cont.)")
	labels = [item.get_text() for item in ax.get_xticklabels()]
	#print(labels)
	#custom_tick_xlabels=['-10','30','40','50','60','70','80','90','100'] #core
	
	#what  to show
	#ax.set_xticklabels(custom_tick_xlabels)

	#plt.title("avg_contact_map_asyn_asyn_1000ns")
	plt.tight_layout()

	#plt.show()
	plt.savefig(f'avg_contact_ll37_8a9l.png', dpi=300)


def array_2d_peptide_contacts(data):
	#print(data.shape)

	#averages_peptides=data.mean(axis=1) #to calculate avg contacts per peptide (all peptides)
	#print(averages_peptides)
	#for avg in averages_peptides:
	#		print(avg.round(2))
	
	#np.save(f'./data_chunk_ll37/tot_contacts_each_ll37_disordered.npy',averages_peptides)
	#sys.exit()
	fig, ax = plt.subplots()
	

	# make a color map of fixed colors
	cmap = colors.ListedColormap(['white', 'red'])
	bounds=[0,40,600]
	norm = colors.BoundaryNorm(bounds, cmap.N)

	im=plt.imshow(data, interpolation='nearest', cmap=cmap, norm=norm,origin='lower', aspect="80")

	#divider = make_axes_locatable(ax)
	#cax = divider.append_axes("right", size="5%", pad=0.1)

	#Colrobar
	cbar=plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds, ticks=[20, 300], shrink=0.25)
	cbar.ax.set_yticklabels(['No contact', 'Contact'])  # vertically oriented colorbar
	labels = [item.get_text() for item in ax.get_xticklabels()]

	custom_tick_xlabels=['-100','100','200','300','400','500','600','700','800','900','1000']
	custom_tick_ylabels_monomer=['','','','']
	
	#what  to show
	#ax.set_yticklabels(custom_tick_ylabels_monomer)
	ax.set_xticklabels(custom_tick_xlabels)
	plt.xlabel("time (ns)")
	plt.ylabel("# n LL-37")
	#plt.title("tot_contacts_6xyo_ll37_cubic_1000ns")
	plt.tight_layout()
	plt.savefig(f'tot_contacts_LL37_8a9l_discrete_legend.png', dpi=300)
	#plt.show()

def array_2d_peptide_contacts_curves(data):
	fig, ax = plt.subplots()
	plt.plot(data[83])
	labels = [item.get_text() for item in ax.get_xticklabels()]
	print(labels)
	custom_tick_labels=['-100','100','300','500','700','900']
	
	#what  to show
	#ax.set_xticklabels(custom_tick_labels)
	plt.axhline(y=40, color='r', linestyle='--')
	plt.xlabel("time (ns)")
	plt.ylabel("totalot contacts ")
	plt.title("LL-37[83]")
	plt.savefig(f'tot_contacts_8a9l_LL37_83_threshold.png', dpi=300)
	#plt.show()

def array_2d_peptide_asyn(data):
	fig, ax = plt.subplots()
	#plt.imshow(data, cmap='hot', origin='lower')
	plt.plot(data)

	#Colrobar
	#plt.colorbar(fraction=0.05, pad=0.04)

	#what  to show
	plt.xlabel("# aSyn chain")
	plt.ylabel("Avg contacts with any LL-37")
	#plt.title("avg_contacts_8a9l_ll37[43]_cubic_1000ns")
	plt.tight_layout()
	plt.savefig(f'avg_contacts_8a9l_LL37_cubic_1000ns.png', dpi=300)
	#plt.show()

def array_1d(data):
	plt.plot(data)

	#what  to show
	plt.xlabel("# asyn chain")
	plt.ylabel("# avg contacts ")
	plt.title("avg_tot_contacts_asyn_ll37_1000ns")
	#plt.savefig(f'avg_tot_contacts_asyn_ll37_1000ns.png', dpi=300)
	plt.show()

#Plot type
array_2d_cmap(data_sum/n_frames/n_chains_pep) #contact map, avg over time and n_peptides
#array_2d_peptide_contacts(data_append.T) #avg peptide contacts over time 
#array_2d_peptide_contacts_curves(data_append.T) #avg peptide X contacts over time 
#array_2d_peptide_asyn(data_sum/n_frames) #specific peptide contacts with each aSyn chain
#array_1d(data_single/n_frames) #avg over time 1d plot