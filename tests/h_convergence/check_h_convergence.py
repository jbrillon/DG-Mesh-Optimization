# ---------------------------------------------
# Description: 
#  - Plot h-convergence and check rates of
#    convergence vs expected
# ---------------------------------------------
# Import libraries
import numpy as np # NumPy: contains basic numerical routines
import matplotlib.pyplot as plt # Matlab-like plotting
from matplotlib.lines import Line2D
from matplotlib import rc as matplotlibrc
matplotlibrc('text.latex', preamble='\usepackage{color}')
matplotlibrc('text',usetex=True)
matplotlibrc('font', family='serif')
# ---------------------------------------------
# Common file parameters
# ---------------------------------------------
subdirectories = ['Data/','Figures/']
data_fileType = 'txt'
figure_fileType = 'png'
# ---------------------------------------------
# Common plotting parameters
# ---------------------------------------------
clr = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown']
mrkr = ['o','s','^','d','v','>','<']
# Font sizes
figTitle_FontSize = 18
axisTitle_FontSize = 16
axisTickLabel_FontSize = 14
legend_fontSize = 16
# ---------------------------------------------
# Load data
# ---------------------------------------------
fields = ["nDOF","L2error","CFL","P","Nel","np"]
# ---------------------------------------------
# h-convergence
# ---------------------------------------------
P_min = 2; P_max=5
P_store = range(P_min,P_max+1)
expectedConvRate = -1.0*(np.array(P_store)+1.0)
computedConvRate = np.zeros(len(P_store))
tol = 0.1
print('')
print('=====================================================')
for polyorder in P_store:
	data_fileName = 'Data/h_convergence_P' + str(polyorder)
	ConvData = np.loadtxt(data_fileName+'.'+data_fileType, unpack=False)

	y_label = "L2-error"
	x_label = "$h^{-1}$"
	P = ConvData[:,fields.index("P")][0]
	npts = ConvData[:,fields.index("np")][0]
	x_data = (1.0)/np.float64(ConvData[:,fields.index("Nel")])

	# compute convergence rate
	err = ConvData[:,fields.index("L2error")]
	dx = 1.0/x_data
	ngrids = dx.size
	orders_of_accuracy = np.empty(ngrids-1)
	for igrid in range(1,ngrids):
		last_slope = np.log(err[igrid]/err[igrid-1])/np.log(dx[igrid]/dx[igrid-1])
		before_last_slope = 1.0*last_slope
		if(igrid>1):
			before_last_slope = np.log(err[igrid-1]/err[igrid-2])/np.log(dx[igrid-1]/dx[igrid-2])
		slope_avg = 0.5*(last_slope + before_last_slope)
		computedConvRate[polyorder-P_min] = slope_avg

	# Plot L2error vs 1/h
	figure_title = "$P=%i$, $N_{el}\\in[%i,%i]$" % (P,ConvData[:,fields.index("Nel")][0],ConvData[:,fields.index("Nel")][-1])
	figure_fileName = "h_convergence_P_%i" % (P)
	print('Plotting: ' + figure_fileName)
	fig, ax1 = plt.subplots(figsize=(6,6))
	ax1.grid()
	ax1.set_title(figure_title,fontsize=axisTitle_FontSize)
	ax1.set_xlabel(x_label,fontsize=axisTitle_FontSize)
	ax1.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
	ax1.semilogy(1.0/x_data, ConvData[:,fields.index("L2error")],label="h-extension", color=clr[1], marker=mrkr[0], markersize=6, mfc='None', linestyle='-')
	index_align = -1
	y_align = ConvData[:,fields.index("L2error")][index_align]
	x_align = x_data[index_align]
	n = -(P+1)
	shift = np.log10(1.0/(y_align*((x_align*10.0)**np.float(n))))
	# Confirm the order of accuracy
	dx_n = (x_data**n)*10.0**(shift+float(n))
	ax1.loglog(1.0/x_data, 1.0/dx_n,label="$h^{-(P+1)}$", color='k', mfc='None', linestyle='--',linewidth=1.0)
	leg = ax1.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
	plt.tight_layout()
	print('\t ... Saving figure ...')
	plt.savefig(subdirectories[1] + figure_fileName + '.' + figure_fileType,bbox_inches='tight',format=figure_fileType,dpi=500)
	if(polyorder!=P_max):
		print('-----------------------------------------------------')
print('=====================================================')
print("Convergence Rate Summary:")
print("P\t Computed\t Expected \t tol")
for i in range(0,len(P_store)):
	if(np.abs(computedConvRate[i]-expectedConvRate[i])<=tol):
		pass_or_fail = "PASSED"
	else:
		pass_or_fail = "FAILED"
	print("%i\t %1.4f\t %.1f\t %.2f\t %s" % (P_store[i],computedConvRate[i],expectedConvRate[i],tol,pass_or_fail))
print('=====================================================')