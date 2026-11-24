"""How much a shell loses in transit, against the three Ceara Rise water columns of ceara_bracket.py, and why.

Prints the calcite and aragonite solubility (Mucci 1983 with the Millero 1995 pressure terms) at surface and
deep T and P; calcite Omega at 4500 m with and without pressure and at +-1 C; and, for tests of 20, 5 and 1 ug
(Takahashi and Be 1984 speeds) and a slow 10 m/day, the days spent below each saturation horizon on the way to
the station's seafloor and the mass lost, with the in-situ rate (Naviaux et al. 2019, Mar. Chem. 215) and with
its laboratory fits (Naviaux et al. 2019, GCA 246), about ten times faster, as an upper bound.  The loss is
summed at constant surface area (BET 2 m2/g); a value over 100% means the shell is gone before the seafloor.

    python scripts/ceara_transit.py            (about a minute)
"""
import math
from reaction import marine
from reaction.kinetics import seawater_rate_law
from reaction.app.engine import sinking_speed
from reaction.app.scenario import Scenario, WATER_PROFILE_CEARA as P198, WATER_PROFILE_CEARA_2015 as P484, WATER_PROFILE_CEARA_1972 as P39
M=100.087; SSA=2.0
pb=lambda z:1.01325+z*1033*9.81/1e5
st={"198 (1994)":(P198,4707),"484 (2015)":(P484,5218),"39 (1972)":(P39,4877)}
print("== Ksp calcite, log10 ==")
for t,p in [(25,1.01325),(2,1.01325),(2,pb(4500)),(25,pb(4500))]:
    c=marine.constants(35,t,p); print(f"T {t} C P {p:.0f} bar: logKsp cal {math.log10(c.ksp_calcite):.3f} arag {math.log10(c.ksp_aragonite):.3f}")
print("== Omega calcite at 4500 m: actual / P=1 atm / T+1C / T-1C ==")
for k,(pr,bot) in st.items():
    t,s,ta,dic=Scenario(water_profile=pr).profile_at(4500.)
    f=lambda T,P: marine.solve(s,T,P,ta,dic).omega_calcite
    print(k, f"{f(t,pb(4500)):.3f} {f(t,1.01325):.3f} {f(t+1,pb(4500)):.3f} {f(t-1,pb(4500)):.3f}  bottom T {t:.2f}")
def transit(pr,bottom,ug,lab=False,mineral="calcite",z0=1000.,dz=1.):
    card=Scenario(water_profile=pr); w=sinking_speed(ug)*86400 if ug>0 else -ug  # m/day
    law=seawater_rate_law(mineral,lab=lab); lost=0; below=0; z=z0
    while z<bottom:
        t,s,ta,dic=card.profile_at(z); cs=marine.solve(s,t,pb(z),ta,dic)
        om=cs.omega_calcite if mineral=="calcite" else cs.omega_aragonite
        r=law.rate_per_m2(t,math.log10(om)); days=dz/w
        lost+=r*SSA*M*days*86400; below+=days if om<1 else 0; z+=dz
    return w,below,min(100*lost,100.0)
print("== transit to the seafloor: speed m/day, days below horizon, % lost ==")
for k,(pr,bot) in st.items():
    for ug,label in [(20,"20 ug"),(5,"5 ug"),(1,"1 ug fragment"),(-10,"10 m/day (aggregate-free slow end)")]:
        for lab in (False,True):
            w,b,l=transit(pr,bot,ug,lab); 
            print(f"{k} bottom {bot} {label} {'lab' if lab else 'insitu'}: w {w:.0f} m/d, {b:.2f} d below, calcite lost {l:.3f}%", end="")
            w,b2,l2=transit(pr,bot,ug,lab,"aragonite"); print(f" | aragonite {b2:.1f} d below, lost {l2:.2f}%")
# Omega needed for 1%/day loss
law=seawater_rate_law("calcite")
for om in (0.95,0.9,0.8,0.7,0.5):
    print(om, f"%/day {law.rate_per_m2(2.0,math.log10(om))*SSA*M*86400*100:.3f}", f"lab {seawater_rate_law('calcite',lab=True).rate_per_m2(2.0,math.log10(om))*SSA*M*86400*100:.3f}")
