import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt('input.txt', unpack = True)
l = data[6]
hh = l/20
t, w = np.loadtxt('w_data.txt', delimiter = ' ', unpack = True)
w2 = np.loadtxt('w2_data.txt', usecols = 1, unpack = True)
theta2 = np.loadtxt('theta2_data.txt', usecols = 1, unpack = True)
M2 = np.loadtxt('M2_data.txt', usecols = 1, unpack = True)
Q2 = np.loadtxt('Q2_data.txt', usecols = 1, unpack = True) 
theta = np.loadtxt('theta_data.txt', usecols = 1, unpack = True)
M = np.loadtxt('M_data.txt', usecols = 1, unpack = True)
Q = np.loadtxt('Q_data.txt', usecols = 1, unpack = True)
wsage = np.loadtxt('sage_w.txt', usecols = 1, unpack = True)
thetasage = np.loadtxt('sage_theta.txt', usecols = 1, unpack = True)
Msage = np.loadtxt('sage_M.txt', usecols = 1, unpack = True)
Qsage = np.loadtxt('sage_Q.txt', usecols = 1, unpack = True)
fig, (ax0,ax4, ax3, ax2, ax1) = plt.subplots(5, figsize = (20,15), sharex = True)
fig.suptitle("Beam solution")
ax1.plot(t,w)
ax1.set_title("w(x)")
ax2.plot(t,theta)
ax2.set_title("theta(x)")
ax3.plot(t,M)
ax3.set_title("M(x)")
ax4.plot(t,Q)
ax4.set_title("Q(x)")
ax1.spines['bottom'].set_position('zero')
ax2.spines['bottom'].set_position('zero')
ax3.spines['bottom'].set_position('zero')
ax4.spines['bottom'].set_position('zero')
ax1.spines['right'].set_color('none')
ax1.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.spines['top'].set_color('none')
ax4.spines['right'].set_color('none')
ax4.spines['top'].set_color('none')
ax0.set_title("Beam:")
ax0.axhline(xmin = 0, xmax = l, linewidth = 10 )
ax0.arrow(4*hh, 0, 0, 2, width = 0.05, head_length = 0.5, label = "R_a", color ='red')
ax0.arrow(11*hh,0, 0, 2, width = 0.05, head_length = 0.5, label = "R_b", color ='green')
ax0.arrow(16*hh,0, 0, 2, width = 0.05, head_length = 0.5, label = "R_c", color ='blue')
ax0.arrow(6*hh,0, 0, 2, width = 0.07, head_length = 0.5, label = "P", color ='orange')
ax0.arrow(2*hh,1, 0.15, 0, width = 0.17, head_length = 0.5, label = "M", color ='black')
ax0.arrow(2*hh,-1, -0.15, 0, width = 0.17, head_length = 0.5, color ='black')
ax0.plot([8*hh,17*hh],[1,1], linewidth = 8, color = 'yellow', label = 'q')
ax0.plot([2*hh,2*hh],[-1,1], linewidth = 3, color = 'black')
ax0.plot([0,0],[-2,2], linewidth = 10, color = 'chocolate', label = 'Z')
ax0.plot([l,l],[-2,0], linewidth = 10, color = 'purple', label = 'K')
ax0.set_ylim([-7,7])
ax0.set_xlim([0,l])
ax0.legend(loc='upper right', ncol = 8)
ax0.grid(True, linestyle = '-.')
ax1.grid(True, linestyle = '-.')
ax2.grid(True, linestyle = '-.')
ax3.grid(True, linestyle = '-.')
ax4.grid(True, linestyle = '-.')
ax0.set_xticks(np.arange(0.0, l + 1.0, 1.0))
ax0.tick_params(axis='x', which='both', labelsize=15, labelbottom=True)
ax1.tick_params(axis='x', which='both', labelsize=15, labelbottom=True)
ax2.tick_params(axis='x', which='both', labelsize=15, labelbottom=True)
ax3.tick_params(axis='x', which='both', labelsize=15, labelbottom=True)
ax4.tick_params(axis='x', which='both', labelsize=15, labelbottom=True)
ax1.plot(t, w2, color = 'red', linestyle = ':', alpha = 0.8, linewidth = 1.5)
ax2.plot(t, theta2, color = 'red', linestyle = ':', alpha = 0.8, linewidth = 1.5)
ax3.plot(t, M2, color = 'red', alpha = 0.8)
ax4.plot(t, Q2, color = 'red', alpha = 0.8)
ax1.plot(t, wsage, color = 'yellow', alpha = 0.2, linewidth = 3)
ax2.plot(t, thetasage, color = 'yellow', alpha = 0.2, linewidth = 3)
ax3.plot(t, Msage, color = 'yellow', alpha = 0.2, linewidth = 3)
ax4.plot(t, Qsage, color = 'yellow', alpha = 0.2, linewidth = 3)
fig.savefig('balka.png')
