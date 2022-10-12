#!/usr/bin/env python


import matplotlib.pyplot as plt


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
V_ecoli = 8*(10**(-16))

# Constants from these papers:
# https://www.sciencedirect.com/science/article/pii/S0006349508001008#tblfn2
# https://www.frontiersin.org/articles/10.3389/fmicb.2020.00057/full
# https://www.sciencedirect.com/science/article/pii/S1631069105000685
# https://www.researchgate.net/figure/Degradation-rate-measurements-expressed-as-h-1-for-IPTG-ATc-and-HSL-in-all-the_tbl1_260022430
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108263

# In Ecoli, 0.002% of the total protein content is the repressor protein of he operon,
# average size of protein in Ecoli is 30 kDa and the average number of proteins in a cell is 2.5*10^6

k_Mc = 30
k_S1Myp = 0.5
OT = 1
k_sMR = 0.23
k_sR = 15
k_2R = 50
k_m2R = 10**(-3)

k_r = 960
k_mr = 2.4
k_dr1 = 3*(10**(-7))
k_mdr1 = 12

k_dr2 = 3*(10**(-7))
k_mdr2 = 4.8*(10**3)
k_s1MY = 0.5
k_sOMY = 0.01
k_sY = 30
k_p = 0.12

k_mp = 0.1
k_ft = 6*(10**4)
k_t = 0.92

#
k_leakyg = 0
#
k_leakyp = 0

l_MR = 0.462
l_MY = 0.462
l_R = 0.2
l_R2 = 0.2
l_Y = 0.2
l_Ylex = 0.2
l_I2R2 = 0.2
l_Iin = 0.001


t = np.arange(0, 100, 0.01)
x = np.arange(0, 1000000, 10000)
# x = [0, 10, 100, 1000, 2*10**3, 3*10**3, 4*10**3, 5*10**3, 6*10**3, 7*10**3, 8*10**3, 9*10**3, 10**4, 2*10**4, 3*10**4, 4*10**4, 5*10**4, 6*10**4, 7*10**4, 8*10**4, 9 *
#      10**4]  # , 10**5, 2*10**5, 3*10**5, 4*10**5, 5*10**5, 6*10**5, 7*10**5, 8*10**5, 9*10**5, 10**6, 2*10**6, 4*10**6, 6*10**6, 8*10**6, 10**7, 2*10**7, 4*10**7, 6*10**7, 8*10**7]  # , 10**8]
y = []

for i in x:
    print(i)

    def dD(D, tm):
        # Gene
        #Da = MR
        #Db = R
        #Dc =R2
        #Dd = O
        #De= I
        #Df = I2R2
        #Dg = MY
        #Dh = Y
        #Di = YIEX
        #OT= 50
        #Iex = 100

        Da, Db, Dc, Dd, De, Df, Dg, Dh, Di, Dj, Dk, Dl, Dm, Dn, Do, Dp, Dq, Dr, Ds = D
        d1 = (k_sMR - l_MR*Da)
        d2 = k_sR*Da - 2*k_2R*(Db)**2 + 2*k_m2R*Dc - l_R*Db
        d3 = k_2R*(Db)**2 - k_m2R*Dc - k_r*Dc*Dd + k_mr*Dd * \
            (50-Dd) - k_dr1*Dc*(De)**2 + k_mdr1*Df-l_R2
        d4 = -k_r*Dc*Dd + k_mr*(50-Dd) + k_dr2*(50-Dd)*(De)**2 - k_mdr2*Dd*Df
        d5 = -2*k_dr1*Dc*(De)**2 + 2*k_mdr1*Df - 2*k_dr2*(50-Dd)*(De)**2 + \
            2*k_mdr2*Dd*Df+k_ft*Di+k_t*(i-De) + 2*l_I2R2*Df + l_Ylex*Di
        d6 = k_dr1*Dc*(De)**2 - k_mdr1*Df+k_dr2*(50-Dd) * \
            (De)**2-k_mdr2*Dd*Df-l_I2R2*Df
        d7 = k_sOMY*(50-Dd) + k_s1MY*Dd - l_MY*Dg
        d8 = k_sY*Dg + (k_ft+k_mp)*Di - k_p*Dh*i-l_Y*Dh
        d9 = -(k_ft+k_mp)*Di + k_p*Dh*i - l_Ylex*Di

        # Plasmid
        #Dj = Mrp
        #Dk = Rp
        #Dl =R2p
        #Dn = Op
        #Do= Ip
        #Dp = I2R2p
        #Dq = Myp
        #Dr = Yp
        #Ds = YIEX

        d10 = (k_sMR - l_MR*Dj)
        d11 = k_sR*Dj - 2*k_2R*(Dk)**2 + 2*k_m2R*Dl - l_R*Dk
        d12 = k_2R*(Dk)**2 - k_m2R*Dl - k_r*Dl*Dm + k_mr*Dm * \
            (50-Dm) - k_dr1*Dl*(Dn)**2 + k_mdr1*Do-l_R2
        d13 = -k_r*Dl*Dm + k_mr*(50-Dm) + k_dr2*(50-Dm)*(Dn)**2 - k_mdr2*Dm*Do
        d14 = -2*k_dr1*Dl*(Dn)**2 + 2*k_mdr1*Do - 2*k_dr2*(50-Dm)*(Dn)**2 + \
            2*k_mdr2*Dm*Do+k_ft*Dr+k_t*(i-Dn) + 2*l_I2R2*Do + l_Ylex*Dr
        d15 = k_dr1*Dl*(Dn)**2 - k_mdr1*Do+k_dr2*(50-Dm) * \
            (Dn)**2-k_mdr2*Dm*Do-l_I2R2*Do
        d16 = k_sOMY*(50-Dm) + k_s1MY*Dm - l_MY*Dp
        d17 = k_sY*Dp + (k_ft+k_mp)*Dr - k_p*Dq*i-l_Y*Dq
        d18 = -(k_ft+k_mp)*Dr + k_p*Dq*i - l_Ylex*Dr
        # Ds is for C
        d19 = k_Mc*Dp

        return [d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19]
    D = odeint(dD, [20, 10, 0, 0, 0, 0, 0, 0, 20, 830,
               20, 10, 0, 0, 0, 0, 0, 0, 0], t)
    Da = D[:, 0]
    Db = D[:, 1]
    Dc = D[:, 2]
    Dd = D[:, 3]
    De = D[:, 4]
    Df = D[:, 5]
    Dg = D[:, 6]
    Dh = D[:, 7]
    Di = D[:, 8]
    Dj = D[:, 9]
    Dk = D[:, 10]
    Dl = D[:, 11]
    Dm = D[:, 12]
    Dn = D[:, 13]
    Do = D[:, 14]
    Dp = D[:, 15]
    Dq = D[:, 16]
    Dr = D[:, 17]
    Ds = D[:, 18]

    from scipy.signal import find_peaks
    peaks = find_peaks(Ds, height=1)
    y.append(Ds[9999])

x1 = [o/10**3 for o in x]
plt.plot(x1, y, linewidth=1, color='blue')
plt.xlabel('IPTG Concentration (uM)')
plt.ylabel('ClyA protein concentration (nM)')
print(y)
plt.show()
