from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

params = {'backend': 'agg',
          'axes.labelsize': 24,  
          'axes.titlesize': 24,
          'axes.labelweight': 'heavy',
          'legend.fontsize': 20,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'text.usetex': True,
          'figure.figsize': [8,6],
          'figure.dpi': 300,
          'savefig.dpi': 300,
          'font.family': 'serif',
          'font.serif': ['Times'],
          'font.weight': 'heavy',
          'lines.linewidth': 2
}
plt.rcParams.update(params)

file_path = "datafile3a.txt"  # Replace with the actual file path
table = ascii.read(file_path)

#exclude interacting systems and systems in which the secondary is a compact object
excluded_flags = ["CV", "MT", "TRI", "AM CVn type", "CV/SWD", "CV/Sq", "CV/TRI"]
excluded_type2 = ["WD", "NS", "BH"]
table = table[~np.isin(table["Flag"], excluded_flags) & ~np.isin(table["Type2"], excluded_type2)]

#limit the mass of the core + guaranteeing that M2 and a have numerical values
table = table[(table["M1"] < 0.87) & (table["M1"] > 0.2) & (table["a"] > 0)& (table["M2"] > 0)]


def estimate_giant_properties(M1):
    if M1 <= 0.47:
        M_MS = 1.19
        M_giant = 0.9 * M_MS
        R_giant = 440 * (M_MS)**(-0.47) * ((M1) / 0.6)**(5.1) 
    else:
        M1=M1+0.028
        if 0.47 <= M1 < 0.53:
            M_MS =  7.42 * M1 - 2.93
        elif 0.53 <= M1 < 0.59:
            M_MS =  13.10 * M1 - 5.95
        elif 0.59 <= M1 < 0.61:
            M_MS =  15.41 * M1 - 7.31
        elif 0.61 <= M1 < 0.633:
            M_MS =  36.90 * M1 - 20.36
        elif 0.633 <= M1 < 0.835:
            M_MS =  5.08 * M1 - 0.21
        elif 0.835 <= M1 <= 0.9:
            M_MS =  32.26 * M1 - 22.90
        M_giant = 0.75 * M_MS
        R_giant = 440 * (M_MS)**(-0.47) * ((M1-0.028) / 0.6)**(5.1) #we decrease the 0.028 because we want the core mass during CEE not the WD mass
    return M_giant, R_giant


m_giant_list = []
r_giant_list = []

for M1 in table["M1"]:
    m_giant, r_giant = estimate_giant_properties(M1)
    m_giant_list.append(m_giant)
    r_giant_list.append(r_giant)

table["M_giant"] = m_giant_list
table["R_giant"] = r_giant_list
table["q_before"] = table["M2"]/table["M_giant"] 


#limit the mass ratio to include only q<1
table = table[table["q_before"]<1] 


file_path = "Tab4.dat"  # Replace with the actual file path
table_SG = ascii.read(file_path)

file_path = "Tab5.dat"  # Replace with the actual file path
table_AGB_RGB = ascii.read(file_path)


print(table)



AGB_data = table_AGB_RGB[table_AGB_RGB["Type"] == "AGB"]
RGB_data = table_AGB_RGB[table_AGB_RGB["Type"] == "RGB"]


plt.figure()
sc = plt.scatter(table["q_before"], table["a"]/table["R_giant"], c=table["M1"], cmap='rainbow', s=70, linewidth=0.3)
plt.xscale("log")
plt.yscale("log")
q = np.linspace(0.15, 1.2, 100)
plt.loglog(q, 0.31 * q - 0.009, '-', color='k', label=r'Eq. (13)')
plt.plot(AGB_data["q"], AGB_data["a_f"] / AGB_data["R1"], 'o', color='k', label='AGB Simulations')
plt.plot(RGB_data["q"], RGB_data["a_f"] / RGB_data["R1"], 's', color='k', label='RGB Simulations')
plt.plot(table_SG["q"], table_SG["a_f"] / table_SG["R1"], '*', color='k', label='Supergiant Simulations', markersize=9)
plt.xlabel("q")
plt.ylabel('$\\rm a_{\\rm f}/ \\rm R_{\\rm 1}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("observations_log.png")


plt.figure()
sc = plt.scatter(table["q_before"], table["a"], c=table["M1"], cmap='rainbow', s=70, linewidth=0.3)
plt.plot(AGB_data["q"], AGB_data["a_f"] , '^', color='k', label='AGB Simulations')
plt.plot(RGB_data["q"], RGB_data["a_f"] , 'o', color='k', label='RGB Simulations')
plt.plot(table_SG["q"], table_SG["a_f"], '*', color='k', label='Supergiant Simulations', markersize=9)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('q')
plt.ylabel('$\\rm a_{\\rm f}~(\\rm R_{\odot})$')
cbar = plt.colorbar(sc)
cbar.set_label('$\\rm M_{\\rm c}~(\\rm M_{\odot})$')
plt.grid()
plt.tight_layout()
plt.savefig("observations_af.png")
