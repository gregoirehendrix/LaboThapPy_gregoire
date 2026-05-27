# -*- coding: utf-8 -*-
"""
T-s diagram — Air Brayton (recuperated, IC, RH)
T_hot = 1000°C, eta = 44.20%   —   Author: Grégoire Hendrix
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from CoolProp.CoolProp import PropsSI
import os

plt.rcParams.update({
    'font.family':'serif','font.size':28,'axes.labelsize':28,
    'xtick.labelsize':24,'ytick.labelsize':24,'legend.fontsize':18,
    'figure.dpi':150,'axes.grid':True,'grid.alpha':0.25,'grid.linestyle':'--',
})
SAVE_DIR=(r"C:\Users\gregoire.hendrix@johncockerill.com"
          r"\OneDrive - John Cockerill\Documents\Cockerill"
          r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR,exist_ok=True)

F='Air'
def _s(T_C,P_bar): return PropsSI('S','T',T_C+273.15,'P',P_bar*1e5,F)/1e3
def _ib(Ta,Tb,P,n=300):
    T=np.linspace(Ta,Tb,n)
    return np.array([PropsSI('S','T',t+273.15,'P',P*1e5,F)/1e3 for t in T]),T

ST={1:(13.20,1.013),2:(111.66,2.266),3:(18.13,2.266),4:(118.27,5.066),
    5:(691.75,5.066),6:(969.67,5.066),7:(779.81,2.266),8:(978.24,2.266),
    9:(787.26,1.013),10:(224.56,1.013)}
S={k:_s(*v) for k,v in ST.items()}

print(f"\n{'St':<4}{'T[°C]':>9}{'P[bar]':>8}{'s[kJ/kgK]':>12}")
for k in ST: print(f"  {k:<3} {ST[k][0]:>9.2f} {ST[k][1]:>8.3f} {S[k]:>12.4f}")

C_HEAT='#d62728';C_TURB='#1f77b4';C_COOL='#2ca02c'
C_RECUP='#9467bd';C_ATM='gray';C_ISO='#bbbbbb'

fig,ax=plt.subplots(figsize=(14,10))

# Guide isobars
for P in [1.013,2.266,5.066]:
    T=np.linspace(-20,1060,400)
    sg=np.array([PropsSI('S','T',t+273.15,'P',P*1e5,F)/1e3 for t in T])
    ax.plot(sg,T,color=C_ISO,lw=1.0,zorder=1)

# Segments
ax.plot([S[1],S[2]],[ST[1][0],ST[2][0]],color=C_TURB,lw=2.5,zorder=4)  # Comp1
s,T=_ib(ST[2][0],ST[3][0],ST[2][1]); ax.plot(s,T,color=C_COOL,lw=2.5,zorder=4)  # IC
ax.plot([S[3],S[4]],[ST[3][0],ST[4][0]],color=C_TURB,lw=2.5,zorder=4)  # Comp2
s,T=_ib(ST[4][0],ST[5][0],ST[4][1]); ax.plot(s,T,color=C_RECUP,lw=2.5,ls='--',zorder=4)  # Recup cold
s,T=_ib(ST[5][0],ST[6][0],ST[5][1]); ax.plot(s,T,color=C_HEAT,lw=2.5,zorder=4)  # Heater
ax.plot([S[6],S[7]],[ST[6][0],ST[7][0]],color=C_TURB,lw=2.5,zorder=4)  # Turb1
s,T=_ib(ST[7][0],ST[8][0],ST[7][1]); ax.plot(s,T,color=C_HEAT,lw=2.5,zorder=4)  # Reheater
ax.plot([S[8],S[9]],[ST[8][0],ST[9][0]],color=C_TURB,lw=2.5,zorder=4)  # Turb2
s,T=_ib(ST[9][0],ST[10][0],ST[9][1]); ax.plot(s,T,color=C_RECUP,lw=2.5,zorder=4)  # Recup hot
s,T=_ib(ST[10][0],ST[1][0],ST[10][1]); ax.plot(s,T,color=C_ATM,lw=1.8,ls='--',zorder=3)  # Atm

# Labels — offsets chosen to avoid overlaps
off={1:(-0.04,-70),2:(0.03,15),3:(-0.04,-70),4:(-0.04,-70),
     5:(-0.14,15),6:(0.03,15),7:(-0.14,15),8:(0.03,15),
     9:(0.03,15),10:(0.03,-70)}
for k,(T,P) in ST.items():
    ax.scatter(S[k],T,color='k',s=70,zorder=7)
    ds,dT=off[k]; ax.text(S[k]+ds,T+dT,str(k),fontsize=21,fontweight='bold',zorder=8)

# Isobar labels at top
sall=list(S.values()); Tg_top=1020
for P,lbl in [(1.013,'1.013 bar'),(2.266,'2.266 bar'),(5.066,'5.066 bar')]:
    sg=PropsSI('S','T',Tg_top+273.15,'P',P*1e5,F)/1e3
    ax.text(sg-0.02,Tg_top,lbl,fontsize=15,color='gray',style='italic',ha='right')

leg=[
    Line2D([0],[0],color=C_HEAT,lw=2.5,label='Heat addition (heater + reheater)'),
    Line2D([0],[0],color=C_TURB,lw=2.5,label='Turbines / Compressors'),
    Line2D([0],[0],color=C_COOL,lw=2.5,label='Intercooler'),
    Line2D([0],[0],color=C_RECUP,lw=2.5,ls='--',label='Recuperator — cold side'),
    Line2D([0],[0],color=C_RECUP,lw=2.5,label='Recuperator — hot side'),
    Line2D([0],[0],color=C_ATM,lw=1.8,ls='--',label='Atmospheric cooling (open cycle)'),
]
ax.legend(handles=leg,loc='upper left',fontsize=16)
ax.set_xlim(min(sall)-0.15, max(sall)+0.20)
ax.set_ylim(-50,1060)
ax.set_xlabel(r'Specific entropy  $s$  [kJ kg$^{-1}$ K$^{-1}$]')
ax.set_ylabel(r'Temperature  $T$  [°C]')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
fig.tight_layout()
path=os.path.join(SAVE_DIR,"fig_brayton_Ts.pdf")
fig.savefig(path,bbox_inches='tight'); print(f"Saved → {path}"); plt.show()