#plot Efficient Frontier
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as Inp
import scipy
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches

Portfolio=np.load("PortfolioOptimizationWindWaveOcean(50_0_50).npz",allow_pickle=True)
LCOE_Stage1=Portfolio['RelaxedSolutionsLCOE']
Std_Stage1=Portfolio['RelaxedSolutionsVar']
#in some simulations I forgot to multiply the rated power by 10^-6

LCOE_Stage1=LCOE_Stage1[LCOE_Stage1!=None]
Std_Stage1=Std_Stage1[(Std_Stage1!=None)]

Std_Stage1=(Std_Stage1)**(1/2)


LCOE_Stage2=Portfolio['MINLPSolutionsLCOE']
Std_Stage2=Portfolio['MINLPSolutionsVar']

LCOE_Stage2=LCOE_Stage2[LCOE_Stage2!=None]
Std_Stage2=Std_Stage2[(Std_Stage2!=None)]

Std_Stage2=Std_Stage2**(1/2)


plt.plot(Std_Stage2[LCOE_Stage2>90], LCOE_Stage2[LCOE_Stage2>90], c='r',linestyle='-', linewidth=1, label="200MW Kite + 300MW Wind")

Portfolio=np.load("PortfolioOptimizationWindWaveOcean(0_0_50).npz",allow_pickle=True)
LCOE_Stage1=Portfolio['RelaxedSolutionsLCOE']
Std_Stage1=Portfolio['RelaxedSolutionsVar']
#in some simulations I forgot to multiply the rated power by 10^-6

LCOE_Stage1=LCOE_Stage1[LCOE_Stage1!=None]
Std_Stage1=Std_Stage1[(Std_Stage1!=None)]#*10**12

Std_Stage1=(Std_Stage1)**(1/2)


LCOE_Stage2=Portfolio['MINLPSolutionsLCOE']
Std_Stage2=Portfolio['MINLPSolutionsVar']

LCOE_Stage2=LCOE_Stage2[LCOE_Stage2!=None]
Std_Stage2=Std_Stage2[(Std_Stage2!=None)]#*10**12

Std_Stage2=Std_Stage2**(1/2)


plt.plot(Std_Stage2[LCOE_Stage2>60], LCOE_Stage2[LCOE_Stage2>60], c='k',linestyle='-',linewidth=1, label="200MW Kite")


Portfolio=np.load("PortfolioOptimizationWindWaveOcean(50_0_0).npz",allow_pickle=True)
LCOE_Stage1=Portfolio['RelaxedSolutionsLCOE']
Std_Stage1=Portfolio['RelaxedSolutionsVar']
#in some simulations I forgot to multiply the rated power by 10^-6

LCOE_Stage1=LCOE_Stage1[LCOE_Stage1!=None]
Std_Stage1=Std_Stage1[(Std_Stage1!=None)]

Std_Stage1=(Std_Stage1)**(1/2)


LCOE_Stage2=Portfolio['MINLPSolutionsLCOE']
Std_Stage2=Portfolio['MINLPSolutionsVar']

LCOE_Stage2=LCOE_Stage2[LCOE_Stage2!=None]
Std_Stage2=Std_Stage2[(Std_Stage2!=None)]

Std_Stage2=Std_Stage2**(1/2)


plt.plot(Std_Stage2[LCOE_Stage2>0], LCOE_Stage2[LCOE_Stage2>0], c='b',linestyle='-',linewidth=1, label="300MW Wind")



leg = plt.legend(loc='lower right')
leg.get_frame().set_alpha(1) # this will make the box totally transparent
leg.get_frame().set_edgecolor("k") # make the edges of the border white to match the background
leg.get_frame().set_linewidth(0.5)





#plt.title('Efficient Frontier: Different Portfolios')
plt.xlabel('\u03C3 (CF)')
plt.ylabel('LCOE [$/MWh]')
plt.yticks(np.arange(59, 150, step=10))
plt.ylim(top=150)
plt.ylim(bottom=50)
plt.grid(color='gray', linestyle='--', linewidth=0.5)
#plt.show()
plt.savefig('EfficientFrontier', dpi=700)
plt.close()  

