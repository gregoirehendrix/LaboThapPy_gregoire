# -*- coding: utf-8 -*-
"""
T-s diagram — Steam Rankine cycle (Siemens SST-PAC DCRH)
T_salt = 565°C, P_high = 160 bar, eta = 46.40%   —   Author: Grégoire Hendrix
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from CoolProp.CoolProp import PropsSI
import os

plt.rcParams.update({
    'font.family':'serif','font.size':28,'axes.labelsize':28,
    'xtick.labelsize':24,'ytick.labelsize':24,'legend.fontsize':20,
    'figure.dpi':150,'axes.grid':True,'grid.alpha':0.25,'grid.linestyle':'--',
})
SAVE_DIR=(r"C:\Users\gregoire.hendrix@johncockerill.com"
          r"\OneDrive - John Cockerill\Documents\Cockerill"
          r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR,exist_ok=True)

F='Water'
def _Ts(h,p): return PropsSI('T','H',h,'P',p,F)-273.15, PropsSI('S','H',h,'P',p,F)/1e3
def _ib(h0,h1,P,n=300):
    hh=np.linspace(h0,h1,n)
    T=np.array([PropsSI('T','H',h,'P',P,F)-273.15 for h in hh])
    s=np.array([PropsSI('S','H',h,'P',P,F)/1e3 for h in hh])
    return s,T
def _expand(h_su,p_su,p_ex,eta=0.93):
    s=PropsSI('S','H',h_su,'P',p_su,F)
    return h_su-eta*(h_su-PropsSI('H','P',p_ex,'S',s,F))

P_high=160e5; P_rh=38.41e5; P_low=0.095e5; P_dea=11.64e5
P_hp_ext=40.50e5; P_fwh5=19.11e5
P_lp4=6.08e5; P_lp3=3.17e5; P_lp2=1.38e5; P_lp1=0.44e5

# State enthalpies [J/kg]
h1=118832.1; h2=1067803.4; h3=3467116.9; h4=3070962.4
h5=3569326.7; h6=2384410.5
# FWH path
h_plp=120287.0; h_fwh1=307249.7; h_fwh2=436260.0; h_fwh3=548742.7
h_fwh4=651333.2; h_dea=792227.8; h_php=813213.2; h_fwh5o=880980.5

# Extraction enthalpies
h_x_hp  =_expand(h3,P_high,P_hp_ext)
h_x_fwh5=_expand(h5,P_rh,P_fwh5)
h_x_dea =_expand(h_x_fwh5,P_fwh5,P_dea)
h_x_lp4 =_expand(h_x_dea,P_dea,P_lp4)
h_x_lp3 =_expand(h_x_lp4,P_lp4,P_lp3)
h_x_lp2 =_expand(h_x_lp3,P_lp3,P_lp2)
h_x_lp1 =_expand(h_x_lp2,P_lp2,P_lp1)

# Print states
print("\n=== Rankine state points ===")
for lbl,h,p in [('1 condensate',h1,P_low),('2 eco inlet',h2,P_high),
                ('3 HP turb in',h3,P_high),('4 HP turb ex',h4,P_rh),
                ('5 IC turb in',h5,P_rh),('6 IC turb ex',h6,P_low)]:
    T,s=_Ts(h,p); print(f"  {lbl:<20} T={T:7.2f}°C  s={s:.4f} kJ/kgK")

# Saturation dome
T_sat=np.linspace(0.1+273.15,373.9+273.15,600)
s_liq=np.array([PropsSI('S','T',T,'Q',0,F) for T in T_sat])/1e3
s_vap=np.array([PropsSI('S','T',T,'Q',1,F) for T in T_sat])/1e3

# Segments
s_boi,T_boi = _ib(h2,h3,P_high)         # boiler
s_rh, T_rh  = _ib(h4,h5,P_rh)           # reheater
s_cd, T_cd  = _ib(h6,h1,P_low)          # condenser
# HP turbine
T_HP=[_Ts(h,p)[0] for h,p in [(h3,P_high),(h_x_hp,P_hp_ext),(h4,P_rh)]]
s_HP=[_Ts(h,p)[1] for h,p in [(h3,P_high),(h_x_hp,P_hp_ext),(h4,P_rh)]]
# IC turbine
h_IC=[h5,h_x_fwh5,h_x_dea,h_x_lp4,h_x_lp3,h_x_lp2,h_x_lp1,h6]
p_IC=[P_rh,P_fwh5,P_dea,P_lp4,P_lp3,P_lp2,P_lp1,P_low]
T_IC=[_Ts(h,p)[0] for h,p in zip(h_IC,p_IC)]
s_IC=[_Ts(h,p)[1] for h,p in zip(h_IC,p_IC)]
# Feedwater
s_fw_lp,T_fw_lp=_ib(h_plp,h_dea,P_dea)   # LP FWH section
s_fw_hp,T_fw_hp=_ib(h_php,h2,P_high)       # HP FWH section
# Pump jumps
T1c,s1c=_Ts(h1,P_low); Tp0,sp0=_Ts(h_plp,P_dea)   # LP pump
Tda,sda=_Ts(h_dea,P_dea); Thp,shp=_Ts(h_php,P_high) # HP pump

C_HEAT='#d62728'; C_TURB='#1f77b4'; C_COND='#2ca02c'
C_FW='#ff7f0e'; C_DOME='k'

fig,ax=plt.subplots(figsize=(14,10))
ax.plot(s_liq,T_sat-273.15,color=C_DOME,lw=2.0,zorder=5)
ax.plot(s_vap,T_sat-273.15,color=C_DOME,lw=2.0,zorder=5)
ax.fill_betweenx(T_sat-273.15,s_liq,s_vap,alpha=0.06,color='steelblue')

ax.plot(s_boi,T_boi,color=C_HEAT,lw=2.5,zorder=4)
ax.plot(s_rh, T_rh, color=C_HEAT,lw=2.5,zorder=4)
ax.plot(s_HP, T_HP, color=C_TURB,lw=2.5,zorder=4)
ax.plot(s_IC, T_IC, color=C_TURB,lw=2.5,zorder=4)
ax.plot(s_cd, T_cd, color=C_COND,lw=2.5,zorder=4)
ax.plot([s1c,sp0],[T1c,Tp0],color=C_FW,lw=1.8,ls='--',zorder=3)  # LP pump
ax.plot(s_fw_lp,T_fw_lp,color=C_FW,lw=1.8,ls='--',zorder=3)
ax.plot([sda,shp],[Tda,Thp],color=C_FW,lw=1.8,ls='--',zorder=3)  # HP pump
ax.plot(s_fw_hp,T_fw_hp,color=C_FW,lw=1.8,ls='--',zorder=3)

# Extraction markers
for h,p in [(h_x_hp,P_hp_ext),(h_x_fwh5,P_fwh5),(h_x_dea,P_dea),
            (h_x_lp4,P_lp4),(h_x_lp3,P_lp3),(h_x_lp2,P_lp2),(h_x_lp1,P_lp1)]:
    T_,s_=_Ts(h,p); ax.scatter(s_,T_,color='#9467bd',s=55,marker='D',zorder=6)

# State labels
off={(h1,P_low,  '1',(  0.05,-35)),(h2,P_high,'2',(  -0.05, 12)),
    (h3,P_high,  '3',( -0.05, 12)),(h4,P_rh,  '4',(  0.05, -20)),
    (h5,P_rh,    '5',( -0.08, 12)),(h6,P_low,  '6',(  0.05,-35))}
for h,p,lbl,(ds,dT) in off:
    T_,s_=_Ts(h,p); ax.scatter(s_,T_,color='k',s=80,zorder=7)
    ax.text(s_+ds,T_+dT,lbl,fontsize=22,fontweight='bold',zorder=8)

# Isobar annotations
T2,s2=_Ts(h2,P_high); ax.text(s2+1.6,T2+70,'160 bar',fontsize=20,color=C_HEAT,style='italic')
T4,s4=_Ts(h4,P_rh);   ax.text(s4+0.15,T4+15,'38.4 bar',fontsize=20,color=C_HEAT,style='italic')
T1c2,s1c2=_Ts(h1,P_low); ax.text(s1c2+4,T1c2+25,'0.095 bar',fontsize=20,color=C_COND,style='italic',ha='right')

leg=[
    Line2D([0],[0],color=C_HEAT,lw=2.5,label='Heat addition (boiler + RH)'),
    Line2D([0],[0],color=C_TURB,lw=2.5,label='Turbine expansion (HP + IC)'),
    Line2D([0],[0],color=C_COND,lw=2.5,label='Condenser'),
    Line2D([0],[0],color=C_FW,lw=1.8,ls='--',label='Feedwater (pumps + FWH 1–6 + DEA)'),
    Line2D([0],[0],color=C_DOME,lw=2.0,label='Saturation dome'),
    plt.scatter([],[],color='#9467bd',marker='D',s=55,label='Extraction points'),
]
ax.legend(handles=leg,loc='upper left',fontsize=20)

# Auto axes
all_s=[_Ts(h,p)[1] for h,p in [(h1,P_low),(h2,P_high),(h3,P_high),
                                 (h4,P_rh),(h5,P_rh),(h6,P_low)]]
ax.set_xlim(min(all_s)-0.3, max(all_s)+0.3)
ax.set_ylim(-15,620)
ax.set_xlabel(r'Specific entropy  $s$  [kJ kg$^{-1}$ K$^{-1}$]')
ax.set_ylabel(r'Temperature  $T$  [°C]')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
fig.tight_layout()
path=os.path.join(SAVE_DIR,"fig_rankine_Ts.pdf")
fig.savefig(path,bbox_inches='tight'); print(f"Saved → {path}")