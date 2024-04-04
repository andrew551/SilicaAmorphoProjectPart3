import numpy as np
import matplotlib .pyplot as plt
import  matplotlib as mpl

import matplotlib
import matplotlib.pyplot as pl
#pl.switch_backend('agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib import rc
from matplotlib import rcParams




def read_file_table(filename):
    with open(filename) as f:
        read_array=[]
        for line in f: # read rest of lines
            read_array.append([float(x) for x in line.split()])
    return np.array(read_array)


def read_grun(filename):
    
    npinput=read_file_table(filename)

    T=npinput[:,0]
    gamma=npinput[:,1]

    return T, gamma



filenames=[f"gruneisen_T_{x}_.dat" for x in range(-2, 3)]



labels=[str(x) for x in range(-2, 3)]
Ts=[]
gammas=[]

#colors
cmap = cm.hsv
matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'mathtext.fontset' : 'dejavuserif'})
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',
       r'\usepackage{amssymb}',   
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               
] 

n=3
c_cyan_nature=np.array([68/255, 128/255, 244/255,1.])
c_red=np.array([206./255., 30./255., 65./255.,1.0])
dark_red='#9e0b00'
c_green_nature='#63cf8c'#np.array([96./255., 172./255., 63./255.,1.0])
c_blue_nature=np.array([54./255., 79./255., 156./255.,1.0])
ec=[c_blue_nature,c_green_nature,c_red, 'yellow', 'purple']
#markers=["o","v","^","<",">","p","s","P","X"]


fig = pl.figure()
#pl.tick_params(labelsize=12)

plt.axhline(y=0, color='black')


for i in range(len(filenames)):
    temp_freq,temp_gamma=read_grun(filenames[i])
    Ts.append(temp_freq)
    gammas.append(temp_gamma)
    plt.semilogx(temp_freq,temp_gamma,color=ec[i],markersize=5,label=labels[i],marker="o")


plt.xlabel(r'$T(K)$')
plt.ylabel(r'$\gamma$')
plt.legend()
plt.ylim([-40,45])
fig.savefig("gruneisen_small.pdf", dpi=300, bbox_inches="tight", transparent=True)

fig.savefig("gruneisen_small.png", dpi=600, bbox_inches="tight", transparent=False)

plt.show()

ev = np.loadtxt('frequencies.dat')
freq = ev[:, 2]#sorted(ev)
print(freq[:20])
#v = np.sqrt(np.abs(ev)*9.648e27)/18.8e10
xx = np.linspace(-np.max(freq)*1.02, np.max(freq)*1.02, 4000)


sigma = 5
conv = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-xx**2/(2*sigma**2))

histo, _ = np.histogram(freq, bins=np.append(xx, xx[-1]*2 - xx[-2])) # bins have one more element than x as they represent edges
dos = np.convolve(conv, histo, mode = 'same')


IS_gamma = read_file_table('gruneisen_IS.dat')
rr_gamma= read_file_table('gruneisen_rr.dat')
pr_n = read_file_table('pr.dat')

plt.clf()
f, axarr = plt.subplots(4, sharex=True, figsize=(15, 30))
x = IS_gamma[:, 0]
plt.xlabel('$\hbar\omega$ $(cm^{-1})$')
axarr[0].axhline(y=0, color='black', linestyle='--')
axarr[1].axhline(y=0, color='black', linestyle='--')

axarr[0].scatter(x, rr_gamma[:, 1], s=2, color=c_cyan_nature)
axarr[1].scatter(x, IS_gamma[:, 1], s=2, color=c_green_nature)

axarr[0].set_ylim((-20, 20))
axarr[1].set_ylim((-20, 20))

axarr[2].scatter(x, pr_n[3:, 1], s=2, color=c_red)
axarr[3].plot(xx[xx > 0], dos[xx > 0])

axarr[0].set_ylabel("$\gamma_i$")
axarr[1].set_ylabel("$\gamma_i$")
axarr[2].set_ylabel("1 / PR")
axarr[3].set_ylabel("DOS (arb. units)")

axarr[0].set_title("Anharmonic Gruneisen")
axarr[1].set_title("Strain Gruneisen")



plt.savefig("gamma_rr.png", dpi=600, bbox_inches="tight", transparent=False)