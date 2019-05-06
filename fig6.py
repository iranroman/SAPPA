import matplotlib.pyplot as plt
import numpy as np
import h5py

with h5py.File('parameter_analysis.mat', 'r') as file:    
	X_m_1 = np.asarray(list(file['X_m_1'])).T
	X_m_0_5 = np.asarray(list(file['X_m_0_5'])).T
	X_m_0 = np.asarray(list(file['X_m_0'])).T
	X_0_5 = np.asarray(list(file['X_0_5'])).T
	X_1 = np.asarray(list(file['X_1'])).T
	X_f_1 = np.asarray(list(file['X_f_1'])).T
	X_f_0_28 = np.asarray(list(file['X_f_0_28'])).T
	
vmin = -350
vmax = 0
fsize = 15

titles = ["-1.0","-0.5","0","0.5","1.0"]
for it,X in enumerate([X_m_1,X_m_0_5,X_m_0,X_0_5,X_1]):

	xvalues = np.linspace(0,0.5,10).round(2)
	yvalues = np.linspace(0,1.0,10).round(2)

	fig, ax = plt.subplots(figsize=(7,7))	
	ax.imshow(X,vmin=vmin, vmax=vmax)
	ax.set_xticks(np.arange(0, 10, 1))
	ax.set_yticks(np.arange(0, 10, 1))
	ax.set_xticklabels(xvalues,fontsize=fsize-2)
	ax.set_yticklabels(yvalues,fontsize=fsize-2)
	ax.set_title('A = '+titles[it],fontsize=fsize,fontweight='bold')
	ax.set_ylabel('D',fontsize=fsize)
	ax.set_xlabel(r"$\tau$",fontsize=fsize)

	for i in xrange(X.shape[0]):
	    for j in xrange(X.shape[1]):
	        c = X[j,i]
	        ax.text(i, j, str('%.1f' % c), va='center', ha='center')

	plt.savefig('test.eps', format='eps', dpi=100)
	plt.show()

titles = ["1.00","0.28"]
for it,X in enumerate([X_f_1,X_f_0_28]):

	xvalues = np.linspace(0,0.3,10).round(2)
	yvalues = np.linspace(0,0.4,10).round(2)

	fig, ax = plt.subplots(figsize=(7,7))	
	ax.imshow(X,vmin=vmin, vmax=vmax)
	ax.set_xticks(np.arange(0, 10, 1))
	ax.set_yticks(np.arange(0, 10, 1))
	ax.set_xticklabels(xvalues,fontsize=fsize-2)
	ax.set_yticklabels(yvalues,fontsize=fsize-2)
	ax.set_title('f = '+titles[it],fontsize=fsize,fontweight='bold',style='italic')
	ax.set_ylabel('D',fontsize=fsize)
	ax.set_xlabel(r"$\tau$",fontsize=fsize)

	for i in xrange(X.shape[0]):
	    for j in xrange(X.shape[1]):
	        c = X[j,i]
	        ax.text(i, j, str('%.1f' % c), va='center', ha='center')

	plt.savefig('test.eps', format='eps', dpi=100)
	plt.show()
