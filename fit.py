#!/usr/bin/python3

import sys
import string
from scipy.optimize import leastsq
import numpy

filename=sys.argv[1]

def heaviside(x):
    return x * (x > 0)
#-------------------------------------------------------------#

def find_min_ene(ee):
    mine = 10000.0
    for i in range(len(ee)):
        ener=ee[i]
        if (ener < mine):
            mine = ener
            ix=i
    return ix

#::::::::::::#


def find_min_and_weights(vv,ee):
    mine = 10000.0
    minvol = 0.0
    for i in range(len(vv)):
        ener=ee[i]
        if (ener < mine):
            mine = ener
            minvol = vv[i]
    weights = []
    for i in range(len(vv)):
        tmp1 = (vv[i]/minvol)**1.5 - 1.0
        tmp2 = 1.0 - tmp1*tmp1
        weights.append(tmp2)
    return weights

#::::::::::::#

def parabola(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    a  = inpars[2]
    E = E0 + a*((vv/V0)**(2.0/3.0)-1.0)**2.0
    return E

def target0(inpars,vv,ee):
    err0 =  ee - parabola(inpars,vv)
    return err0

#............#

def target_bm(inpars,vv,ee,p):
    err =  p*(ee - birch_murnaghan(inpars,vv))
    tt=heaviside(-inpars[2]) 
    return err+tt*1000

def birch_murnaghan(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    eta = (V0/vv)**(2.0/3.0)
    tmp1 = (9.0/16.0)*B0*V0
    tmp2 = b1*(eta-1.0)**3.0
    tmp3 = (6.0-4.0*eta)*(eta-1.0)**2
    E = E0 + tmp1*(tmp2+tmp3)
    return E

#............#

def target_vinet(inpars,vv,ee,p):
    err =  p*(ee - vinet(inpars,vv))
    return err

def vinet(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    eta = (vv/V0)**(1.0/3.0)
    tmp1 = 2.0*B0*V0/(b1-1.0)**2.0
    tmp2 = 5.0 + 3.0*eta*(b1-1.0) - 3.0*b1
    tmp3 = -3.0*(b1-1)*(eta-1.0)/2.0
    E = E0 + tmp1*(2.0 - tmp2*numpy.exp(tmp3))
    return E

#............#

def target_murnaghan(inpars,vv,ee,p):
    err =  p*(ee - murnaghan(inpars,vv))
    tt=heaviside(-inpars[2])       
    return err+tt*1000

def murnaghan(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    eta = (V0/vv)**b1
    ww = vv * B0 
    tmp1 = ww/b1
    tmp2 = 1.0 + eta/(b1-1.0)
    tmp3 = B0*V0/(b1-1.0)
    E = E0 + tmp1*tmp2 - tmp3
    return E

def murnaghan_plot(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    eta = (V0/vv)**b1
    ww = B0*vv
    tmp1 = ww/b1
    tmp2 = 1.0 + eta/(b1-1.0)
    tmp3 = B0*V0/(b1-1.0)
    E = E0 + tmp1*tmp2 - tmp3
    return E

#............#

def target_pt(inpars,vv,ee,p):
    err =  p*(ee - poirier_tarantola(inpars,vv))
    return err

def poirier_tarantola(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    eta = numpy.log(V0/vv)
    tmp1 = B0*V0/2.0
    tmp2 = eta**2.0
    tmp3 = (b1-2.0)*eta**3.0
    E = E0 + tmp1*tmp2 + tmp1*tmp3/2.0
    return E

#............#

def target_birch6(inpars,vv,ee,p):
#    print inpars
    err =  p*(ee - birch6(inpars,vv))
#    print inpars[2],err
    tt=heaviside(-inpars[2])**2
    return err+tt*1000
#+tt

def birch6(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    a1 = inpars[4]
    a2 = inpars[5]
    a3 = inpars[6]
    eta = (V0/vv)**(2.0/3.0)
    tmp1 = 9.0*B0*V0/8.0
    tmp2 = (eta-1.0)**2.0
    tmp3 = 9.0*B0*V0*(b1-4.0)/16.0
    tmp4 = (eta-1.0)**3.0
    tmp5 = (eta-1.0)**4.0
    tmp6 = (eta-1.0)**5.0
    tmp7 = (eta-1.0)**6.0
    E = E0 + tmp1*tmp2 + tmp3*tmp4 + a1*tmp5 + a2*tmp6 + a3*tmp7
    return E

#............#

def target_birch5(inpars,vv,ee,p):
    err =  p*(ee - birch5(inpars,vv))
    return err

def birch5(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    a1 = inpars[4]
    a2 = inpars[5]
    eta = (V0/vv)**(2.0/3.0)
    tmp1 = 9.0*B0*V0/8.0
    tmp2 = (eta-1.0)**2.0
    tmp3 = 9.0*B0*V0*(b1-4.0)/16.0
    tmp4 = (eta-1.0)**3.0
    tmp5 = (eta-1.0)**4.0
    tmp6 = (eta-1.0)**5.0
    E = E0 + tmp1*tmp2 + tmp3*tmp4 + a1*tmp5 + a2*tmp6 
    return E

#............#

def target_birch4(inpars,vv,ee,p):
    err =  p*(ee - birch4(inpars,vv))
    return err

def birch4(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    a1 = inpars[4]
    eta = (V0/vv)**(2.0/3.0)
    tmp1 = 9.0*B0*V0/8.0
    tmp2 = (eta-1.0)**2.0
    tmp3 = 9.0*B0*V0*(b1-4.0)/16.0
    tmp4 = (eta-1.0)**3.0
    tmp5 = (eta-1.0)**4.0
    E = E0 + tmp1*tmp2 + tmp3*tmp4 + a1*tmp5 
    return E

#............#

def target_sj(inpars,vv,ee,p):
    err =  p*(ee - stabilized_jellium(inpars,vv))
    return err

def stabilized_jellium(inpars,vv):
    E0 = inpars[0]
    V0 = inpars[1]
    B0 = inpars[2]
    b1 = inpars[3]
    eta= ((V0/vv)**2)**0.5
    eta2 = eta**(2.0/3.0)
    eta3 = eta**(1.0/3.0)
    tmp1 = 9.0*B0*V0/2.0
    tmp2 = (b1-3.0)*eta
    tmp3 = (10.0-3.0*b1)*eta2
    tmp4 = -(11.0-3.0*b1)*eta3
    tmp5 = 4.0-b1
    E = E0 + tmp1*(tmp2+tmp3+tmp4+tmp5)
    return E

#::::::::::::::::::::#

def plot_and_error(func,inpars,ee,vv,p):
    error = 0.0
    plot = []
    ptot = 0.0
    for i in range(len(vv)):
        vvv = vv[i]
        eee = ee[i]
        aaa = func(inpars,vvv)
        error = error + p[i]*(eee-aaa)**2
        ptot = ptot + 1.0
        plot.append(aaa)
    error = error/ptot
    return (error,plot)

#-------------------------------------------------------------#



#read data
#print(filename)
data = numpy.loadtxt(filename)
vol = data[:,0]
ene = data[:,1]
# inital guess (use parabola)
reflat = float(sys.argv[2])
refbulk = float(sys.argv[3])
reflatyp = sys.argv[4]
#print(reflat,refbulk,reflatyp) 
if (reflatyp == "A1" or reflatyp == "A4" or reflatyp == "B1" or reflatyp == "B3"):
    refvol = 2.0*(reflat/2.0)**3.0
elif (reflatyp == "A2"):
    refvol = 4.0*(reflat/2.0)**3.0
b1 = 4.0
ix = find_min_ene(ene)
#refene = find_min_ene(ene,ix)
refene = ene[ix]
refvol = vol[ix]
print(refene,refvol)
print(refene,refvol,refbulk)
exit
print()
print(" +++++++++++++++++++++++++++++++++ SIMPLE FIT ++++++++++++++++++++++++++++++++++")
print()

wei = numpy.ones(len(vol))


print("Method                E_0 (Ha)    V_0 (A^3)    B_0 (GPa)     B' (Ha/A^6)   Error")
print("--------------------------------------------------------------------------------")

# fit with murnaghan
guess=(refene,refvol,refbulk/4359.7458849645852,b1)
outpars, ier = leastsq(target_murnaghan, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
#-saved guess------------
mineneguess=minene
minvolguess=minvol
bulkguess=outpars[2]
#--------------------
pe = plot_and_error(murnaghan_plot,outpars,ene,vol,wei)
murg_plot = pe[1]
print("Murnaghan         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))


# fit with birch_murnaghan
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1)
outpars, ier = leastsq(target_bm, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch_murnaghan,outpars,ene,vol,wei) 
bm_plot = pe[1]
print("Birch-Murnaghan   %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with Poirier-Tarantola
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1) 
outpars, ier = leastsq(target_pt, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(poirier_tarantola,outpars,ene,vol,wei) 
pt_plot = pe[1]
print("Poirier-Tarantola %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with vinet
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1) 
outpars, ier = leastsq(target_vinet, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(vinet,outpars,ene,vol,wei) 
vinet_plot = pe[1]
print("Vinet             %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with birch 4th-order
#guess=(refene,refvol,refbulk,b1,1.0)
guess=(mineneguess,minvolguess,bulkguess,b1,1.0) 
outpars, ier = leastsq(target_birch4, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch4,outpars,ene,vol,wei) 
b4_plot = pe[1]
print("Birch-4th         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with birch 5th-order
#guess=(refene,refvol,refbulk,b1,1.0,1.0)
#guess=(minene,minvol,bulk,b1,1.0,1.0)
guess=(mineneguess,minvolguess,bulkguess,b1,1.0,1.0) 
outpars, ier = leastsq(target_birch5, guess, args=(vol,ene,wei))
if (ier > 4):
    print("ERROR b5", ier)
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch5,outpars,ene,vol,wei) 
b5_plot = pe[1]
print("Birch-5th         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with birch 6th-order
#guess=(refene,refvol,refbulk,b1,1.0,1.0,1.0)
#guess=(minene,minvol,bulk,b1,1.0,1.0,1.0) 
guess=(mineneguess,minvolguess,bulkguess,b1,1.0,1.0,1.0)  
#print guess
outpars, ier = leastsq(target_birch6, guess, args=(vol,ene,wei))
if (ier > 4):
    print("ERROR b6", ier)
minene = outpars[0]
minene6 = minene
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch6,outpars,ene,vol,wei) 
b6_plot = pe[1]
print("Birch-6th         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with stabilized jellium
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1)  
outpars, ier = leastsq(target_sj, guess, args=(vol,ene,wei))
if (ier > 4): 
     print("ERROR sj", ier)
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(stabilized_jellium,outpars,ene,vol,wei) 
sj_plot = pe[1]
print("Stabilized-jell   %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))


print("Fitted curves saved in file fit.dat")
f = open("fit.dat","w")
f.write("#vol murnaghan birch-murnaghan poirier-tarantola vinet birch-4th birch-5th birch-6th stabilized-jell Ref\n")
for i in range(len(vol)):
    f.write("%10.3f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n" %(vol[i],murg_plot[i],bm_plot[i],pt_plot[i],vinet_plot[i],b4_plot[i],b5_plot[i],b6_plot[i],sj_plot[i],ene[i])) 
f.close()

f = open("fit0.dat","w")  
for i in range(len(vol)):   
  f.write("%10.3f %16.8f \n" %(vol[i],ene[i]-minene6))   
f.close()

print()
print(" +++++++++++++++++++++++++++++++ WEIGHTED FIT ++++++++++++++++++++++++++++++++++")
print()


wei = find_min_and_weights(vol,ene)
#print wei
   
print("Method                E_0 (Ha)    V_0 (A^3)    B_0 (GPa)     B' (Ha/A^6)   Error")
print("--------------------------------------------------------------------------------")

# fit with murnaghan
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1)  
outpars, ier = leastsq(target_murnaghan, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(murnaghan_plot,outpars,ene,vol,wei)
murg_plot = pe[1]
print("Murnaghan         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))


# fit with birch_murnaghan
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1)  
outpars, ier = leastsq(target_bm, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch_murnaghan,outpars,ene,vol,wei) 
bm_plot = pe[1]
print("Birch-Murnaghan   %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with Poirier-Tarantola
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1)  
outpars, ier = leastsq(target_pt, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(poirier_tarantola,outpars,ene,vol,wei) 
pt_plot = pe[1]
print("Poirier-Tarantola %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with vinet
#guess=(refene,refvol,refbulk,b1)
guess=(mineneguess,minvolguess,bulkguess,b1)  
outpars, ier = leastsq(target_vinet, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(vinet,outpars,ene,vol,wei) 
vinet_plot = pe[1]
print("Vinet             %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with birch 4th-order
#guess=(refene,refvol,refbulk,b1,1.0)
guess=(mineneguess,minvolguess,bulkguess,b1,1.0)  
outpars, ier = leastsq(target_birch4, guess, args=(vol,ene,wei))
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch4,outpars,ene,vol,wei) 
b4_plot = pe[1]
print("Birch-4th         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with birch 5th-order
#guess=(refene,refvol,refbulk,b1,1.0,1.0)
#guess=(minene,minvol,bulk,    b1,1.0,1.0)
guess=(mineneguess,minvolguess,bulkguess,b1,1.0,1.0)  
outpars, ier = leastsq(target_birch5, guess, args=(vol,ene,wei))
if (ier > 4):
    print("ERROR b5" , ier)
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch5,outpars,ene,vol,wei) 
b5_plot = pe[1]
print("Birch-5th         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))

# fit with birch 6th-order
#guess=(refene,refvol,refbulk,b1,1.0,1.0,1.0)
#guess=(minene,minvol,bulk,    b1,1.0,1.0,1.0)
guess=(mineneguess,minvolguess,bulkguess,b1,1.0,1.0,1.0)  
outpars, ier = leastsq(target_birch6, guess, args=(vol,ene,wei))
if (ier > 4):
    print("ERROR b6" , ier)
minene = outpars[0]
minene6 = minene
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(birch6,outpars,ene,vol,wei) 
b6_plot = pe[1]
print("Birch-6th         %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))


# fit with stabilized jellium
guess=(refene,refvol,refbulk,b1)
outpars, ier = leastsq(target_sj, guess, args=(vol,ene,wei))
if (ier > 4): 
     print("ERROR sj" , ier)
minene = outpars[0]
minvol = outpars[1]
bulk = outpars[2]*4359.7458849645852
pe = plot_and_error(stabilized_jellium,outpars,ene,vol,wei) 
sj_plot = pe[1]
print("Stabilized-jell   %17.7f  %10.3f    %8.2f     %8.2f      %8.2e" %(minene,minvol,bulk,outpars[3],pe[0]))



print("Fitted curves saved in file wfit.dat")
f = open("wfit.dat","w")
f.write("#vol murnaghan birch-murnaghan poirier-tarantola vinet birch-4th birch-5th birch-6th stabilized-jell Ref\n")
for i in range(len(vol)):
    f.write("%10.3f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n" %(vol[i],murg_plot[i],bm_plot[i],pt_plot[i],vinet_plot[i],b4_plot[i],b5_plot[i],b6_plot[i],sj_plot[i],ene[i])) 
f.close()


f = open("wfit0.dat","w")  
for i in range(len(vol)):    
 f.write("%10.3f %16.8f \n" %(vol[i],ene[i]-minene6))   
f.close()
