from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

params = {'backend': 'agg',
          'axes.labelsize': 24,  
          'axes.titlesize': 24,
          'axes.labelweight': 'heavy',
          'legend.fontsize': 14,
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
excluded_type2 = [ "WD", "NS", "BH"]
table = table[~np.isin(table["Flag"], excluded_flags) & ~np.isin(table["Type2"], excluded_type2)]

#limit the mass of the core + guaranteeing that M2 and a have numerical values
table = table[(table["M1"] < 0.87) & (table["M1"] > 0.2) & (table["a"] > 0)& (table["M2"] > 0) & (table["Name"] !='KOI-3278')]
# KOI-3278 is reintroduced below

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




df1 = pd.read_csv('Tab3.dat', delim_whitespace=True)
df2 = pd.read_csv('Tab4.dat', delim_whitespace=True)

df = pd.concat([df1, df2], ignore_index=True)

df_agb = df[df["Type"] == "AGB"]
df_sg = df[df["Type"] == "SG"]
df_rgb = df[df["Type"] == "RGB"]

##### Intermediate-period observtions ######

# The objects below are in the following order J2117+0332, J1111+5515, J1314+3818, J2034-5037, J0107-2827, KOI-3278. 
# To estimate R1 for KOI-3278, we considered P_orb= 978 days (onset of CEE, Belloni et al., 2024a) to estimate "a" and, consequently, R thorough Eggleton's formula
# For the other objects R_1 is already given in Belloni et al.(2024b).
G=6.67e-8
MS_newobs=np.array([3.37,3.23,2.44,2.88,2.55,1.18]) # in Msun
WD_newobs=np.array([1.24,1.37,1.32,1.42,1.27,0.53]) # in Msun
Sec_newobs=np.array([1.11,1.14,0.71, 0.96, 0.98,0.91])  # in Msun
R1_newobs=np.array([1123,1273,1336,1375,1262,213])  # in Rsun
Porb_newobs=np.array([17.9,32.2,45.5,46.1,49,88.2])*24*3600 # in seconds
a_final_newobs =((G*(WD_newobs+Sec_newobs)* 2e+33 *Porb_newobs**2/4/np.pi**2)**(1/3))/6.957e+10 # in Rsun

plt.figure()
sc = plt.scatter(table["q_before"], table["a"]/table["R_giant"], c=table["M1"], cmap='rainbow', s=70, linewidth=0.3)
plt.xscale("log")
plt.yscale("log")
q = np.linspace(0.15, 1.2, 100)
plt.loglog(q, 0.31 * q -0.01, '-', color='k', label='Eq. (6)')

plt.plot(df_agb["q"], df_agb["a_f"] / df_agb["R1"], 's', color='k', label='AGB Simulations')
plt.plot(df_rgb["q"], df_rgb["a_f"] / df_rgb["R1"], 'P', color='k', label='RGB Simulations')
plt.plot(df_sg["q"], df_sg["a_f"] / df_sg["R1"], '*', color='k', label='Supergiant Simulations', markersize=9)
plt.plot(Sec_newobs/MS_newobs, a_final_newobs / R1_newobs, 'd', color='deeppink', markersize=8, label='Gaia Observations')
plt.plot(Sec_newobs[-1]/MS_newobs[-1], a_final_newobs[-1]/ R1_newobs[-1], 'o', color='hotpink', markersize=11, label=' KOI 3278 ')


plt.xlabel("q")
plt.ylabel('$\\rm a_{\\rm f}/ \\rm R_{\\rm 1}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("observations_log.png")


plt.figure()
sc = plt.scatter(table["q_before"], table["a"], c=table["M1"], cmap='rainbow', s=70, linewidth=0.3)

plt.plot(df_agb["q"], df_agb["a_f"] , 's', color='k', label='AGB Simulations')
plt.plot(df_rgb["q"], df_rgb["a_f"] , 'P', color='k', label='RGB Simulations')
plt.plot(df_sg["q"], df_sg["a_f"], '*', color='k', label='Supergiant Simulations', markersize=9)
plt.plot(Sec_newobs/MS_newobs, a_final_newobs, 'd', color='deeppink', markersize=8, label='Gaia Observations')
plt.plot(Sec_newobs[-1]/MS_newobs[-1], a_final_newobs[-1], 'o', color='hotpink', markersize=11, label=' KOI 3278 ')
plt.xscale("log")
plt.yscale("log")
plt.xlabel('q')
plt.ylabel('$\\rm a_{\\rm f}~(\\rm R_{\odot})$')
cbar = plt.colorbar(sc)
cbar.set_label('$\\rm M_{\\rm c}~(\\rm M_{\odot})$')
plt.grid()
plt.tight_layout()
plt.savefig("observations_af.png")
