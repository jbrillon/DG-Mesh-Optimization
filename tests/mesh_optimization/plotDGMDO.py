# Import libraries
import numpy as np
from numpy import linalg as LA
import subprocess
#-----------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib import rc as matplotlibrc
matplotlibrc('text.latex', preamble='\usepackage{color}')
matplotlibrc('text', usetex=True)
matplotlibrc('font', family='serif')
#=====================================================
restart_level = 0 # USER INPUT
#=====================================================
fig_directory = "Figures/MDO/"
for i in range(0,restart_level):
	fig_directory += "restart/"
# Create the directory if it doesn't already exist
subprocess.call(["mkdir", fig_directory])
figure_filetype = "pdf"
clr = ['tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:red']
mrkr = ['o','s','^','d','v','>','<']
# Font sizes
axisTitle_FontSize = 16
axisTickLabel_FontSize = 14
legend_fontSize = 16
#=====================================================
# Open the files for storing stuff
# global L2error_file, vert_file, err_file, iter_file
mdo_data_dir = "./Data/MDO/"
nVert, xL, xR = np.loadtxt(mdo_data_dir+"plot_parameters.txt",unpack=True)
nVert = int(nVert)
data_dir = mdo_data_dir
for i in range(0,restart_level):
	data_dir += "restart/"
L2error_store = np.loadtxt(data_dir+"L2error_MDO_store.txt")
sol_store = np.loadtxt(data_dir+"vertices_MDO_store.txt")
err_store = np.loadtxt(data_dir+"err_MDO_store.txt")
iter_store = np.loadtxt(data_dir+"iter_MDO_store.txt")
alpha_restart_store = [0.025,0.025]
iter_restart_store = [0]
data_dir = mdo_data_dir
for i in range(0,restart_level):
	data_dir += "restart/"
	iter_restart = np.loadtxt(data_dir+"iter_restart.txt")
	iter_restart = int(iter_restart)
	iter_restart_store.append(iter_restart)
#=====================================================
# 				   POST PROCESSING
#=====================================================
print("---------------------------------------------")

# midpoint = (nVert-1)/2
# cL = range(0,midpoint+1)
# norm = mpl.colors.Normalize(vmin=0, vmax=midpoint)
# cmapL = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.copper)
# cmapL.set_array([])
# cmapL.set_clim([0, midpoint])
# cR = range(midpoint+1,nVert)
# norm = mpl.colors.Normalize(vmin=midpoint, vmax=(nVert-1))
# cmapR = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.copper)
# cmapR.set_array([])
# cmapR.set_clim([midpoint,nVert-1])


x_label = '$x$'
y_label = 'Iteration'
figure_filename = "path_of_vertices"
figure_title = "Optimization Path of Vertices"
print('Plotting: ' + figure_filename)
fig, ax = plt.subplots(figsize=(9,6))
# plt.grid()
plt.xlim([xL,xR])
ax.set_title(figure_title,fontsize=axisTitle_FontSize)
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

for i in range(0,nVert):
	plt.plot(sol_store[:,i], iter_store, markersize=6, mfc='None', linestyle='-')
	# if(i>midpoint):
	# 	plt.plot(sol_store[:,i], iter_store, markersize=6, mfc='None', linestyle='-',color=cmapR.to_rgba(i))
	# else:
	# 	plt.plot(sol_store[:,i], iter_store, markersize=6, mfc='None', linestyle='-',color=cmapL.to_rgba(midpoint-i))
	# if()

# leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
plt.tight_layout()
# plt.show()
print('\t ... Saving figure ...')
plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
print("---------------------------------------------")
#-----------------------------------------------------
print("---------------------------------------------")
y_label = 'Magnitude'
x_label = 'Iteration'
figure_filename = "convergence"
figure_title = "Convergence of $||\\nabla f||$ and L2-error"
print('Plotting: ' + figure_filename)
fig, ax = plt.subplots(figsize=(6,6))
plt.grid()
ax.set_title(figure_title,fontsize=axisTitle_FontSize)
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

	# plt.plot(sol_store[:,i], iter_store, color=clr[0], markersize=6, mfc='None', linestyle='-',label=name)
plt.semilogy(iter_store, err_store, color=clr[0],markersize=6, mfc='None', linestyle='-', label="$||\\nabla f||$")
plt.semilogy(iter_store, L2error_store, color=clr[1],markersize=6, mfc='None', linestyle='-', label="L2-error")

# if(restart_level>0):
# 	for i, restart_iter in enumerate(iter_restart_store):
# 		name = "$\\alpha = %.3f$" % alpha_restart_store[i]
# 		plt.plot(restart_iter, err_store[restart_iter], marker=mrkr[i],color=clr[0],markersize=6, mfc='None', linestyle='None', label=name)

leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
plt.tight_layout()
# plt.show()
print('\t ... Saving figure ...')
plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
print("---------------------------------------------")
#-----------------------------------------------------
print("---------------------------------------------")

data_fileName = 'u_exact'
u_exact = np.loadtxt("./Data/"+data_fileName+'.'+"txt", unpack=False)

data_fileName = 'x'
x = np.loadtxt("./Data/"+data_fileName+'.'+"txt", unpack=False)

x_label = '$x$'
y_label = 'Iteration'
figure_filename = "path_of_vertices_with_solution"
figure_title = "Optimization Path of Vertices with Exact Solution"
print('Plotting: ' + figure_filename)
fig, ax = plt.subplots(figsize=(9,6))

# plt.grid()
plt.xlim([xL,xR])
ax.set_title(figure_title,fontsize=axisTitle_FontSize)
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel("$u(x)$",rotation=90,fontsize=axisTitle_FontSize,color=clr[-1])  # we already handled the x-label with ax1
for eln in range(0,(nVert+1)):
	ax2.plot(x[eln,:],u_exact[eln,:],linestyle="-",color=clr[-1],alpha=1.0)
ax2.tick_params(axis='y',labelcolor=clr[-1])
plt.setp(ax2.get_yticklabels(),fontsize=axisTickLabel_FontSize,color=clr[-1]);

for i in range(0,nVert):
	ax.plot(sol_store[:,i], iter_store, markersize=6, mfc='None', color='k', linestyle='-',alpha=0.7)
	# if(i>midpoint):
	# 	plt.plot(sol_store[:,i], iter_store, markersize=6, mfc='None', linestyle='-',color=cmapR.to_rgba(i))
	# else:
	# 	plt.plot(sol_store[:,i], iter_store, markersize=6, mfc='None', linestyle='-',color=cmapL.to_rgba(midpoint-i))
	# if()


# leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
plt.tight_layout()
# plt.show()
print('\t ... Saving figure ...')
plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
print("---------------------------------------------")
#-----------------------------------------------------
print("---------------------------------------------")
x_label = '$x$'
y_label = "$u(x)$"
figure_filename = "exact_sol_equispaced_vertices"
figure_title = "Exact Solution with Equispaced Vertices"
print('Plotting: ' + figure_filename)
fig, ax = plt.subplots(figsize=(9,6))

# plt.grid()
plt.xlim([xL,xR])
ax.set_title(figure_title,fontsize=axisTitle_FontSize)
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

for eln in range(0,(nVert+1)):
	ax.plot(x[eln,:],u_exact[eln,:],linestyle="-",color=clr[-1],alpha=1.0)

for i in range(0,nVert):
	ax.plot(sol_store[0,i]*np.ones(2), np.linspace(1.0,2.0,2), markersize=6, mfc='None', color='k', linestyle='-',alpha=0.7)
# ax.plot(sol_store[:,-1], np.linspace(xL,xR,nVert), markersize=6, mfc='None', color='k', linestyle='-',alpha=0.7)

# leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
plt.tight_layout()
# plt.show()
print('\t ... Saving figure ...')
plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
print("---------------------------------------------")
#-----------------------------------------------------
print("---------------------------------------------")
x_label = '$x$'
y_label = "$u(x)$"
figure_filename = "exact_sol_optimal_vertices"
figure_title = "Exact Solution with Optimal Vertices"
print('Plotting: ' + figure_filename)
fig, ax = plt.subplots(figsize=(9,6))

# plt.grid()
plt.xlim([xL,xR])
ax.set_title(figure_title,fontsize=axisTitle_FontSize)
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

for eln in range(0,(nVert+1)):
	ax.plot(x[eln,:],u_exact[eln,:],linestyle="-",color=clr[-1],alpha=1.0)

for i in range(0,nVert):
	ax.plot(sol_store[-1,i]*np.ones(2), np.linspace(1.0,2.0,2), markersize=6, mfc='None', color='k', linestyle='-',alpha=0.7)
# ax.plot(sol_store[:,0], np.linspace(1.0,2.0,nVert), markersize=6, mfc='None', color='k', linestyle='-',alpha=0.7)


# leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
plt.tight_layout()
# plt.show()
print('\t ... Saving figure ...')
plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
print("---------------------------------------------")