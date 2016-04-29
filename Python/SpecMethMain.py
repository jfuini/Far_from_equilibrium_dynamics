################# JOHN FUINI #######################
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import h5py

from scipy import empty, zeros, cos, sin, pi, optimize, exp
from numpy.polynomial.chebyshev import Chebyshev, chebval
from specmethmodules import *
from stringsformathematica import report
from matplotlib.backends.backend_pdf import PdfPages


#in order to interpret mathematica output correctly...
from numpy import power as Power
from numpy import e as E
from numpy import log as Log
from numpy import sqrt as Sqrt
from numpy import pi as Pi

############################## TO DO #####################################


############################## TO DO #####################################


setup = raw_input("Commander, which setup would you like to use?")
if setup in ('g','customfunc','custom'):
    N = 120 #240 for linear 60 for guassian
    N2 = 2*N/3
elif setup in ('l'):
    N = 60 #24 for linear 60 for guassian
    N2 = 2*N/3
else:
    print "ABORT ABORT ABORT!!!"
    sys.exit(0)    
    
    
x1 = 0
x2 = 1
grid = GLgrid(N,x1,x2)
coarsegrid = GLgrid(N2,x1,x2)
#N cardinals on the coarsegrid
d = dcard_M(N,x1,x2)
dd = d2card_M(N,x1,x2)

print "N:",N
print "N2:", N2
timestepchoice = 'rk4' # choose 'e' or 'rk4'
print "timestep choice:", timestepchoice



#### TIME STEPS
timesteps = 10000### <------------------------------ TIME STEPS
#### TIME STEPS
tstep = 0.001 # should be 0.001
print "Number of steps:", timesteps
horizontol = 10**(-3)




filt = "F" # use "T", "F", or "Multi"
print "Filter:", filt

#Lambda fixing option DO NOT USE VERY OLD
coordchange = False
backtracklamb = False



# initial guess for lamb
lamb = 0.0338098 #0.35#0.3284099042333

 


lambdafindchoice = raw_input("Should we loop to find the right lambda immediatly? (doesn't work with custom grid)")
if lambdafindchoice in ('yes','y'):
    lambdafindcountermax = 1000
    print "Ok, looping to find the right lambda to start with"
elif lambdafindchoice in ('no','n'):
    lambdafindcountermax = 0
    print "Fine. Not looping, just proceeding with your guess for lambda"
else:
    print "ruh roh"
    sys.exit(0)
if setup in ('custom'):
    found = True
    pass
else:
    maxwellchoice = raw_input("Should we turn on Maxwell fields?")
    if maxwellchoice in ('yes'):
        maxwelldetails = raw_input("Radial, Fxy, or both?")





########## LOOP INITIAL CONDITIONS TO GET LAMBDA RIGHT ON FIRST TRY ###############################


for lambdafindcounter in range(0,lambdafindcountermax):
    
    print "Attempt No.", lambdafindcounter, "to find the correct initial lamb"
    

    ########## LINEAR SETUP #####################
    # After we started subtracting logs from B, it is no longer correct to specify and B and then divide by u^3, so it is much easier to just specify a B_sub.

    if setup in ('l,lin, linear'):
        print 'B_sub = u, understood.'
        g = 1.0
        a2 = -1.0 
        amp =1.0
        uave = 0.0
        width = 0.0
        amp2 = 0.0
        uave2 = 0.0
        width2 = 0.0
        def initial_B(r):
            return amp*(1./r)**4
        def initial_B_sub(u):
            return initial_B(1./u + lamb)/Power(u,3) - (Power(fxy,2)*u*Log(1./u))/3. +  (4.*Power(fxy,2)*lamb*Power(u,2)*Log(1./u))/3. - (10.*Power(fxy,2)*Power(lamb,2)*Power(u,3)*Log(1./u))/3.
           
        
        #maxwellchoice = 'n'    
        if maxwellchoice in ('yes','y'):
            print "Roger, Maxwell fields are go."
            
            if maxwelldetails in ('r,radial'): 
                amp = 0.4
                f3max = (4./3*(-a2)**3)**(1./4)
                f3 = 0.0*f3max #6.9:, -15,  0.598342
                fxy = 0.0

            elif maxwelldetails in ('xy'):
                amp = 1.0
                f3 = 0.0
                fxy = .1
      
            else:
                sys.exit(0)
                
        elif maxwellchoice in ('no','n'):
            sys.exit(0)            
            print "Understood, negative on Maxwell fields."
            fxy = 0.0
            f3 = 0.0

    ############# GAUSSIAN SETUP #################################
                
    elif setup in ('gaussian, g, gauss'):
        print 'B is gaussian, wilco.'
        gausssetup = 'r'
        
        if gausssetup in ('u'):
            g = 1.0
            a2 = -1.0
            uave = 0.4
            width = 0.15 
            amp =  0.75
            def initial_B_sub(u): # FOR B WITHOUT LAMBDA - MUST START NUMERICS AT CORRECT LAMBDA# (only matters for initial B when rerunning with diff lambda -> NEED TO ADD IN lambda's SO THAT THIS FUNCTION SCALES ACCORDINGLY WITH RADIAL REPARAMETERIZATION
                #with lambda inside B
                #return (amp*u)/(Power(E,Power(u/(1. - lamb*u) - uave,2)/(2.*Power(width,2)))*Sqrt(2.*Pi)*(1. - lamb*u)*(1. + Power(u,4)/Power(1. - lamb*u,4))*width)
                #without lambda - initial B not rescaled when lambda scaled THIS IS THE WAY TO GO.
                return (amp*u)/(Power(E,Power(u - uave,2)/(2.*Power(width,2)))*Sqrt(2*Pi)*(1 + Power(u,4))*width)
        #    print "part 1", (amp*.3)/(Power(E,Power(.3/(1. - lamb*.3) - uave,2)/(2.*Power(width,2)))*Sqrt(2.*Pi)*(1. - lamb*.3)*(1. + Power(.3,4)/Power(1. - lamb*.3,4))*width)

        elif gausssetup in ('r'):
            g = 1.0  # DONT TOUCH
            a2 = -1.0# DONT TOUCH
            uave = 0.1# DONT TOUCH
            rave = 1./uave# DONT TOUCH
            width = 0.4 # DONT TOUCH
            amp =  1.0# DONT TOUCH
            ampLog = 1.0# DONT TOUCH
            uave2 = 0.6# DONT TOUCH
            rave2 = 1./uave2# DONT TOUCH
            width2 = 0.6# DONT TOUCH
            amp2 = 0.0# DONT TOUCH
                
            #def initial_B_sub(u):
                #return amp/Power(E,Power(rave - 1./u,2)/(2.*Power(width,2)))
            # The idea here is to define B, and not Bsub.  This is because B has good transformation properties, and we are trying to use the known lambda on the first time slice to make sure that the peak of B is set the same distance from the boundary in all cases.
            # defining B in lambda = 0 case.        
            def initial_B(r):
                #return amp/Power(E,Power(rave - r,2)/(2.*Power(width,2)))      
	         return amp/Power(E,Power(rave - r,2)/(2.*Power(width,2))) + ampLog*Log(r)*fxy**2/(3.*r**4) + amp2/Power(E,Power(rave2 - r,2)/(2.*Power(width2,2)))
           # def initial_B(u): 
               # return amp/Power(E,Power(rave - r,2)/(2.*Power(width,2)))            
#We then pull off Bsub from B(r + lambda), and we'll need to check at the end that the lambda we guess on the first time slice IS THE CORRECT LAMBDA that sets uh = 1
            def initial_B_sub(u):
                return initial_B(1./u + lamb)/Power(u,3) - (Power(fxy,2)*u*Log(1./u))/3. +  (4.*Power(fxy,2)*lamb*Power(u,2)*Log(1./u))/3. - (10.*Power(fxy,2)*Power(lamb,2)*Power(u,3)*Log(1./u))/3.
                

        else:
            print "no way broseidon, king of the brocean"
            sys.exit(0)
            

        
        if maxwellchoice in ('yes','y'):
            print "Roger, Maxwell fields are go."

            if maxwelldetails in ('r,radial'): #Maximum f3 = (4/3*a2^3)^(1/4)
                a2 = -1.0 #-3.577442 # -1.0
                amp = 0.005
                fxy = 0.0
                uave = 1.0#1.0/1.1 #0.2#100.0
                rave = 1./uave        
                width = 0.5
                f3max = (4./3*(-a2)**3)**(1./4)
                f3 = 0.8*f3max #6.9:, -15,  0.598342
                amp2 = 0.0#0.05
                width2 = 0.0
                uave2 = 1.0
                rave2 = 1./uave2
                
     
	
                  
    # NOTE for large FXY had to shift lambda by smaller amounts to make sure to not overshoot.
            elif maxwelldetails in ('xy'):
                a2 = -1.0#-1.675
                amp = 0.00005 #0.0005
                fxy = 1.5#1.901
                uave = 0.25#10./11.
                rave = 1./uave
                width = 0.5
                f3 = 0.0 ### DONT TOUCH
                amp2 = 0.00#0.05

                
            elif maxwelldetails in ('both','b'):
                fxy = 1.0
                f3 = 1.0
                
            ############### MULTI RUN CONTROL ##########################

            #This is an advanced setup and requires knowledge about good choices before hand.  The Control file contains all directions.
            #not sure how this works with the lambda loop.
            #elif maxwelldetails in ('multi'):
            #    print 'Multi-run Here we go.'
            #    g = 1.            
            #    amp = 0.000000005
            #    f3 = 0     
        
             #   a2choice = raw_input("a2?")
             #   a2 = float(a2choice)
             #   fxychoice = raw_input("fxy?")
             #   fxy = float(fxychoice)
                #Grid 1            
                #lamb = 0.018985 + 0.0587857*(fxy) - 0.0642652*(a2)# This needs to be setup serparetely depening on the paramaters of the multi run         
                #Grid 2
             #   lamb = -0.0458115 + 0.123288*(fxy) - 0.0606061*(a2)
                
                
            else:
                sys.exit(0)
        elif maxwellchoice in ('no','n'):
            print "Understood, negative on Maxwell fields."
            maxwelldetails = "NA"
            a2 = -1.0
            amp = 1.0 #0.182221
            fxy = 0.0
            f3 = 0.0
            
        else:
            print "Say again."
            sys.exit(0)
    
    elif setup in ('customfunc'):
        g = 1.0  # DONT TOUCH
        a2 = -1.0# DONT TOUCH
        uave = 0.1# DONT TOUCH
        rave = 1./uave# DONT TOUCH
        width = 0.4 # DONT TOUCH
        amp =  1.0# DONT TOUCH
        ampLog = 1.0# DONT TOUCH
        uave2 = 0.6# DONT TOUCH
        rave2 = 1./uave2# DONT TOUCH
        width2 = 0.6# DONT TOUCH
        amp2 = 0.0# DONT TOUCH
        f3 = 0.0
        fxy = 4.5
        def initial_B(r):
           	return (9.2757101)*(1./(r-0.789674))**4 + (-34.7932)*(1./(r-0.789674))**5 + (54.05604)*(1./(r-0.789674))**6 + (-39.394847)*(1./(r-0.789674))**7 + (11.0899369)*(1./(r-0.789674))**8
        def initial_B_sub(u):
            return initial_B(1./u + lamb)/Power(u,3) - (Power(fxy,2)*u*Log(1./u))/3. +  (4.*Power(fxy,2)*lamb*Power(u,2)*Log(1./u))/3. - (10.*Power(fxy,2)*Power(lamb,2)*Power(u,3)*Log(1./u))/3.
        

    ############# CUSTOM SETUP #################################
                
    
        
       
    else:
        print "Can't do that sir."
        sys.exit(0)


    if setup in ('custom'):
        pass
    else:
        B_sub = vectorize(initial_B_sub,N,x1,x2)
        B_sub[N] = 0 #### we can do this, because the limiting value is zero.  Blindly pluging in gives NAN.
        #print "Bsub is:", B_sub
        
        
    

    result = find_dots(grid,lamb,B_sub,g,a2,f3,fxy,N,x1,x2,filt = 'F', Nfilter = N2, horizonfind = False, dpass = d, ddpass = dd, lambdafix = True, sigmadottol = 10**(-11), backtracklamb = backtracklamb, horizontol = horizontol, pr = False, pl = False, lambdafindoption = True)
    templamb = result['lamb']
    print "difference", lamb - templamb
    if abs(lamb - templamb) < 10**(-10):
        print "Got it! Try lamb = ", lamb
        found = True
        break
    else:
        print "else!"
        found = False
        lamb = templamb



if found:
    print "Ok, here it is, lamb = ", lamb
else:
    print "Tried ", lambdafindcountermax, " times, and couldn't find lamb :("

if setup in ('custom'):
        maxwelldetails = "NA"
        print 'B_sub is custom made, double check the code - here we go.'
        a2 =  -0.41#-0.7 #-0.245333
        g = 1.            
        fxy = 1.56844 #4.2, 0.646632092249#3.4, 0.446936217096 #3.0 0.325508779913  #1.8, -0.078883 #1.2:  -0.300036199972 #0.4: -0.167026990718
        f3 = 0.0
        uave = 0.0
        amp = 0.0
        width = 0.
        amp2 = 0
        width2 = 0
        uave2 = 0
        amplog = 0
        B_sub = np.array([0.519833, 0.519763, 0.519552, 0.519201, 0.518709, 0.518078, 0.517308, \
0.516399, 0.515352, 0.514168, 0.512847, 0.511392, 0.509802, 0.508079, \
0.506225, 0.50424, 0.502127, 0.499887, 0.497521, 0.495031, 0.49242, \
0.489689, 0.48684, 0.483876, 0.480798, 0.477608, 0.47431, 0.470906, \
0.467397, 0.463787, 0.460078, 0.456273, 0.452375, 0.448386, 0.444309, \
0.440147, 0.435903, 0.43158, 0.427181, 0.422709, 0.418167, 0.413558, \
0.408886, 0.404153, 0.399363, 0.394519, 0.389623, 0.38468, 0.379693, \
0.374664, 0.369597, 0.364495, 0.359361, 0.354199, 0.349012, 0.343803, \
0.338575, 0.333331, 0.328074, 0.322808, 0.317536, 0.312259, 0.306983, \
0.301709, 0.29644, 0.29118, 0.285932, 0.280697, 0.275479, 0.27028, \
0.265104, 0.259952, 0.254827, 0.249732, 0.244669, 0.239641, 0.23465, \
0.229698, 0.224787, 0.219918, 0.215095, 0.21032, 0.205593, 0.200917, \
0.196293, 0.191724, 0.187211, 0.182755, 0.178358, 0.174022, 0.169746, \
0.165533, 0.161383, 0.157298, 0.153279, 0.149325, 0.14544, 0.141622, \
0.137873, 0.134194, 0.130584, 0.127045, 0.123576, 0.120177, 0.116849, \
0.113592, 0.110405, 0.107289, 0.104245, 0.10127, 0.098367, 0.0955337, \
0.0927703, 0.0900761, 0.0874505, 0.0848928, 0.0824025, 0.079979, \
0.0776215, 0.0753295, 0.0731022, 0.0709387, 0.0688381, 0.0667993, \
0.0648211, 0.0629023, 0.0610419, 0.0592388, 0.0574918, 0.0558, \
0.0541621, 0.0525772, 0.0510438, 0.0495605, 0.0481261, 0.0467392, \
0.0453982, 0.0441021, 0.0428496, 0.0416394, 0.0404703, 0.0393409, \
0.0382499, 0.0371956, 0.0361768, 0.0351919, 0.0342398, 0.0333192, \
0.0324292, 0.0315688, 0.0307367, 0.0299319, 0.0291528, 0.028398, \
0.0276661, 0.0269557, 0.0262659, 0.0255958, 0.0249446, 0.0243118, \
0.0236967, 0.0230985, 0.0225162, 0.0219487, 0.0213949, 0.0208537, \
0.0203243, 0.019806, 0.0192983, 0.018801, 0.0183135, 0.0178354, \
0.0173661, 0.0169051, 0.0164517, 0.0160053, 0.0155656, 0.0151323, \
0.0147055, 0.0142851, 0.013871, 0.0134631, 0.0130611, 0.0126645, \
0.0122728, 0.0118856, 0.0115026, 0.0111237, 0.010749, 0.0103789, \
0.0100139, 0.00965456, 0.00930121, 0.00895411, 0.00861324, \
0.00827833, 0.0079489, 0.00762437, 0.00730413, 0.00698765, \
0.00667459, 0.00636483, 0.00605849, 0.00575596, 0.00545783, \
0.00516479, 0.00487764, 0.00459717, 0.00432409, 0.00405902, \
0.00380241, 0.00355461, 0.00331578, 0.00308596, 0.0028651, \
0.00265307, 0.00244968, 0.00225473, 0.00206803, 0.0018894, \
0.00171869, 0.00155579, 0.00140065, 0.00125323, 0.00111356, \
0.000981676, 0.000857649, 0.000741566, 0.00063353, 0.000533648, \
0.000442029, 0.000358781, 0.000284005, 0.000217796, 0.00016024, \
0.000111411, 0.0000713726, 0.0000401776, 0.0000178664, 4.46804*10**(-6),
 0.])
        print B_sub
        B_sub[N] = 0.0

print "BEGINNING TIME STEPS.  initiallamb = ", lamb, ", f3 = ", f3, ", fxy = ", fxy, "."

#### REPORTS
numb4reports = 40000   #XXXX
reportb4 = timesteps/numb4reports
if reportb4 == 1:
    numb4reports = timesteps
if reportb4 == 0:
    reportb4 = 1
    numb4reports = timesteps+1
print "reportb4", reportb4, "numb4reports", numb4reports

numfuncreports = 100#want this many reports of entire radial function
reportfunc = timesteps/numfuncreports
if reportfunc == 0:
    reportfunc = 1
    numfuncreports = timesteps+1
print "reportfunct", reportfunc, "numfuncreports", numfuncreports





RK4_dvB_sub_list = zeros([4,N+1])
RK4_dvlamb_list = zeros(4)
reportb4_t_list = zeros(numb4reports)
b4_list = zeros(numb4reports)
lamb_list = zeros(numb4reports)
lambshift_list = zeros(numb4reports)
bdotbndry_list = zeros(numb4reports)
b4check_list = zeros(numb4reports)
horizon_list = zeros(numb4reports)
a2check_list = zeros(numb4reports)
Ahrzncheck_list = zeros(numb4reports)
sigmadotathor_list = zeros(numb4reports)
sigmadotsubscaleat0_list = zeros(numb4reports)
horizonalt_list = zeros(numb4reports)
dsigmasubscaleat0_list = zeros(numb4reports)
dAsubathor_list = zeros(numb4reports)
sigmasubscaleathor_list = zeros(numb4reports)
dsigmasubscaleathor_list = zeros(numb4reports)
ddsigmasubscaleathor_list = zeros(numb4reports)
Bsubathor_list = zeros(numb4reports)
dBsubathor_list = zeros(numb4reports)
ddBsubathor_list = zeros(numb4reports)

sigmaBCcheckD_list = zeros(numb4reports)
sigmaBCcheckN_list = zeros(numb4reports)
sigmadotBCcheckD_list = zeros(numb4reports)
sigmadotBCcheckN_list = zeros(numb4reports)
BdotBCcheckD_list = zeros(numb4reports)
BdotBCcheckN_list = zeros(numb4reports)
ABCcheckD_list = zeros(numb4reports)
ABCcheckN_list = zeros(numb4reports)

reportfunc_t_list = zeros(numfuncreports)
reportfunc_lamb_list = zeros(numfuncreports)
B_sub_list = zeros([numfuncreports, N+1])
Bdot_sub_list = zeros([numfuncreports, N+1])
dvB_sub_list = zeros([numfuncreports, N+1])
sigmadot_paul_list = zeros([numfuncreports, N+1])
sigma_sub_scale_list = zeros([numfuncreports, N+1])
A_sub_list = zeros([numfuncreports, N+1])
dsigma_sub_scale_list = zeros([numfuncreports, N+1])
Fvu_list = zeros([numfuncreports, N+1])


radialdest = '/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/SpecMethOutputRadialsetup' + setup + 'N' + str(N) + 'amp' + str(amp) + 'a2' + str(a2) + 'f3' + str(f3) + 'fxy' + str(fxy) + 'TS' + str(tstep) +'NS'+ str(timesteps) + 'F' + str(filt) + '.txt'
outputdest = '/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/SpecMethOutputsetup' + setup + 'N' + str(N) + 'amp' + str(amp)  + 'wid' + str(width) + 'uave' + str(uave) + 'a2' + str(a2) + 'f3' + str(f3) + 'fxy' + str(fxy) + 'TS' + str(tstep) +  'amp2' + str(amp2) + 'wid2' + str(width2) + 'uave2' + str(uave2) + 'NS'+ str(timesteps) + 'F' + str(filt) + '.h5'

#report(grid,comparedest,'w','grid')
if timestepchoice == 'rk4' and filt == 'T':
    print "Don't try to explicitally filter for rk4.  Rk4 has a filtering step built in."   
    sys.exit(0)

#### Python Plotting Stuff
finegrid = np.linspace(.5,1,100)
finegridrestricted = np.linspace(0,.1,100)
fig1 = plt.figure(figsize=(20,12))
plot1 = fig1.add_subplot(231)
plot2 = fig1.add_subplot(232)
plot3 = fig1.add_subplot(233)
plot4 = fig1.add_subplot(234)
plt.title('Sigmadot')
plot5 = fig1.add_subplot(235)
plot6 = fig1.add_subplot(236)
plot1.set_title('B_sub')
plot2.set_title('Sigmasubscale')
plot3.set_title('Fvu')
plot5.set_title('Sigmadot_sub')
plot6.set_title('A_sub')



if timestepchoice == 'rk4':
    for i in range(0,timesteps+1):
        if i%10 == 0 and i*10 >= 50:
            pr = False
        for k in range(0,4):
            if k == 0:
                print "Timestep/rkstep", i,"/",k,"   of",timesteps, "    Time:",tstep*i,"      a2 = ",a2, "      fxy = ", fxy
                #report([i,k],comparedest,'a', 'Time')
#                if i%100 == 0:
#                    coordchange = True
                result = find_dots(grid,lamb,B_sub,g,a2,f3,fxy,N,x1,x2,filt = 'F', Nfilter = N2, horizonfind = False, dpass = d, ddpass = dd, lambdafix = True, sigmadottol = 10**(-11), backtracklamb = backtracklamb, horizontol = horizontol, pr = False, destination = radialdest, pl = False, time = i)
#                if coordchange == True:
#                    coordchange = False
                RK4_dvlamb_list[k] = result['dvlamb']
                RK4_dvB_sub_list[k] = result['dvB_sub']
                
                
                #RK4_dvlamb_list[k] = result[1]
                #RK4_dvB_sub_list[k] = result[2]
                fixedlamb = result['lamb']
                print "sigmadot(1) = ", result['sigmadot at hor']
                if (i%reportb4 == 0) and (i/reportb4 < numb4reports):
                    reportb4_t_list[i/reportb4] = tstep*i
                    b4_list[i/reportb4] = result['b4']
                    lamb = result['lamb']
                    lamb_list[i/reportb4] = result['lamb']
                    lambshift_list[i/reportb4] = result['lambshift']
                    bdotbndry_list[i/reportb4] = result['Bdot_sub_scale at bndry']
                    b4check_list[i/reportb4] = result['b4check']
#                    horizon_list[i/reportb4] = result[6]
                    a2check_list[i/reportb4] = result['a2check']
                    Ahrzncheck_list[i/reportb4] = result['ABC']
#                    sigmadotvalcheck_list[i/reportb4] = result[9]
                    sigmadotsubscaleat0_list[i/reportb4] = result['sigmadot_sub_scale at bndry']
#                    horizonalt_list[i/reportb4] = result[15]
                    dsigmasubscaleat0_list[i/reportb4] = result['dsigma_sub_scale at bndry']
                    dAsubathor_list[i/reportb4] = result['dA_sub at hor']
                    sigmasubscaleathor_list[i/reportb4] = result['sigma_sub_scale at hor']
                    dsigmasubscaleathor_list[i/reportb4] = result['dsigma_sub_scale at hor']
                    dsigmasubscaleathor_list[i/reportb4] = result['ddsigma_sub_scale at hor']
                    Bsubathor_list[i/reportb4] = result['B_sub_scale at hor']
                    dBsubathor_list[i/reportb4] = result['dB_sub_scale at hor']
                    ddBsubathor_list[i/reportb4] = result['ddB_sub_scale at hor']
                    sigmadotathor_list[i/reportb4] = result['sigmadot at hor']
                    
                    ##### BC CONSISTENCY
                    sigmaBCcheckD_list[i/reportb4] = result['sigmaBCcheckD']
                    sigmaBCcheckN_list[i/reportb4] = result['sigmaBCcheckN']
                    sigmadotBCcheckD_list[i/reportb4] = result['sigmadotBCcheckD']
                    sigmadotBCcheckN_list[i/reportb4] = result['sigmadotBCcheckN']
                    BdotBCcheckD_list[i/reportb4] = result['BdotBCcheckD']
                    BdotBCcheckN_list[i/reportb4] = result['BdotBCcheckN']
                    ABCcheckD_list[i/reportb4] = result['ABCcheckD']
                    ABCcheckN_list[i/reportb4] = result['ABCcheckN']                     
#                    reportb4_t_list[i/reportb4] = tstep*i
#                    b4_list[i/reportb4] = result[3]
#                    lamb = result[0]
#                    lamb_list[i/reportb4] = result[0]
#                    bdotbndry_list[i/reportb4] = result[4]
#                    b4check_list[i/reportb4] = result[5]
#                    horizon_list[i/reportb4] = result[6]
#                    a2check_list[i/reportb4] = result[7]
#                    Ahrzncheck_list[i/reportb4] = result[8]
#                    sigmadotvalcheck_list[i/reportb4] = result[9]
#                    sigmadotsubscaleat0_list[i/reportb4] = result[11]
#                    horizonalt_list[i/reportb4] = result[15]
#                    dsigmasubscaleat0_list[i/reportb4] = result[17]
#                    dAsubathor_list[i/reportb4] = result[21]
                if (i%reportfunc == 0) and (i/reportfunc < numfuncreports):
                    reportfunc_t_list[i/reportfunc] = tstep*i
                    reportfunc_lamb_list[i/reportfunc] = result['lamb']
                    B_sub_list[i/reportfunc] = B_sub
                    dvB_sub_list[i/reportfunc] = result['dvB_sub'] 
                    Bdot_sub_list[i/reportfunc] = result['Bdot_sub_scale']
                    sigmadot_paul_list[i/reportfunc] = result['sigmadot_paul']
                    sigma_sub_scale_list[i/reportfunc] = result['sigma_sub_scale']
                    A_sub_list[i/reportfunc] = result['A_sub']
                    dsigma_sub_scale_list[i/reportfunc] = result['dsigma_sub_scale']
                    Fvu_list[i/reportfunc] = result['Fvu']
                      
#                    reportfunc_t_list[i/reportfunc] = tstep*i
#                    reportfunc_lamb_list[i/reportfunc] = result[0]
#                    B_sub_list[i/reportfunc] = B_sub
#                    dvB_sub_list[i/reportfunc] = result[2] 
#                    Bdot_sub_list[i/reportfunc] = result[12]
#                    sigmadot_sub_scale_list[i/reportfunc] = result[13]
#                    sigma_sub_scale_list[i/reportfunc] = result[14]
#                    A_sub_list[i/reportfunc] = result[16]
#                    dsigma_sub_scale_list[i/reportfunc] = result[18]
                    
                    ### Python Plotting

                    plt.figure(1)
                    plot1.plot(grid,B_sub_list[i/reportfunc],'-')
                    #plt.draw()
                    plot2.plot(grid,sigma_sub_scale_list[i/reportfunc],'-')
                    #plot3.plot(grid,result[19],'-')
                    #plot3.plot(finegridrestricted,cardfunct(B_sub,finegridrestricted,N,x1,x2),'-')
                    #plot4.plot(finegrid, result[20](finegrid),'o', finegrid, zeros([len(finegrid)]),'-')
                    plot5.plot(grid,sigmadot_paul_list[i/reportfunc],'-')
                    plot6.plot(grid,A_sub_list[i/reportfunc],'-')

                    
#                if result[10] == False:
#                    horizontol = horizontol*10

            elif k == 1:
                #report([i,k],comparedest,'a', 'Time')
#                print "Timestep/rkstep", i,"/",k , "----------------------------------"
                #report([lamb + 1./2.*tstep*RK4_dvlamb_list[0]],comparedest,'a','PYlamb')
                #report(B_sub + 1./2.*tstep*RK4_dvB_sub_list[0],comparedest,'a','PYBsubpaulAlt')
                
                result = find_dots(grid,lamb + 1./2.*tstep*RK4_dvlamb_list[0], B_sub + 1./2.*tstep*RK4_dvB_sub_list[0],g,a2,f3,fxy,N,x1,x2,filt = 'F', Nfilter = N2, horizonfind = False, dpass = d, ddpass = dd, lambdafix = False, destination=radialdest)
                RK4_dvlamb_list[1] = result['dvlamb']
                RK4_dvB_sub_list[1] = result['dvB_sub']                
                
                #RK4_dvlamb_list[1] = result[1]
                #RK4_dvB_sub_list[1] = result[2]
            elif k == 2:
                #report([i,k],comparedest,'a', 'Time')
#                print "Timestep/rkstep", i,"/",k , "----------------------------------"
                result = find_dots(grid,lamb + 1./2.*tstep*RK4_dvlamb_list[1], B_sub + 1./2.*tstep*RK4_dvB_sub_list[1],g,a2,f3,fxy,N,x1,x2,filt = 'F', Nfilter = N2, horizonfind = False, dpass = d, ddpass = dd, lambdafix = False, destination=radialdest)
                RK4_dvlamb_list[2] = result['dvlamb']
                RK4_dvB_sub_list[2] = result['dvB_sub']                
                
                #RK4_dvlamb_list[2] = result[1]
                #RK4_dvB_sub_list[2] = result[2]
            elif k == 3:
#                print "Timestep/rkstep", i,"/",k , "----------------------------------"
                #report([i,k],comparedest,'a', 'Time')
                result = find_dots(grid,lamb + tstep*RK4_dvlamb_list[2],B_sub + tstep*RK4_dvB_sub_list[2],g,a2,f3,fxy,N,x1,x2,filt = 'F', Nfilter = N2, horizonfind = False, dpass = d, ddpass = dd, lambdafix = False, destination=radialdest)
                RK4_dvlamb_list[3] = result['dvlamb']
                RK4_dvB_sub_list[3] = result['dvB_sub']                
                
#                RK4_dvlamb_list[3] = result[1]
#                RK4_dvB_sub_list[3] = result[2]
            
        dvB_sub_total = 1./6.*(RK4_dvB_sub_list[0] + 2.*RK4_dvB_sub_list[1] + 2.*RK4_dvB_sub_list[2] + RK4_dvB_sub_list[3])    
        #Filtering individually so don't need a total filter
        #print "Filtering..."
        #if (i%100 == 0):
        #    dvB_sub_total = filter_func(dvB_sub_total,N,N2,x1,x2) 
        B_sub = B_sub + tstep*dvB_sub_total
        #print "B_sub at bdnry =", B_sub[N]
        lamb = fixedlamb + 1./6.*tstep*(RK4_dvlamb_list[0] + 2.*RK4_dvlamb_list[1] + 2.*RK4_dvlamb_list[2] + RK4_dvlamb_list[3])





#### SIMPLE TIME STEP
#if timestepchoice == 'e':
#    flag = False
#    for i in range(0,timesteps+1):
#        print "Timestep:", i, "Time:",i*tstep
#        if filt == "Multi":
#            flag = "OFF"
#            if (i % 10 == 0):
#                filt = "T"
#                flag = "ON"
#        result = find_dots(grid,lamb,B_sub,g,a2,f3,fxy,N,x1,x2,filt = 'F', Nfilter = N2, horizonfind = True,dpass = d, ddpass = dd, fixlambdrift = True ,pr = False,destination=dest)
#        if flag == "ON":
#            filt = "Multi"
#            flag = "OFF"        
#        if (i%reportb4 == 0) and (i/reportb4 < numb4reports):
#            reportb4_t_list[i/reportb4] = tstep*i
#            b4_list[i/reportb4] = result[3]
#            lamb = result[0]
#            lamb_list[i/reportb4] = result[0]
#            bdotbndry_list[i/reportb4] = result[4]
#            b4check_list[i/reportb4] = result[5]
#            horizon_list[i/reportb4] = result[6]
#            a2check_list[i/reportb4] = result[7]
#            Ahrzncheck_list[i/reportb4] = result[8]
#            sigmadotvalcheck_list[i/reportb4] = result[9]
####  New observables havn't been added here.            
#        if (i%reportfunc == 0) and (i/reportfunc < numfuncreports):
#            print len(B_sub_list), len(lamb_list), i
#            reportfunc_t_list = tstep*i
#            B_sub_list[i/reportfunc] = B_sub
#            dvB_sub_list[i/reportfunc] = result[2] 
            
#        lamb = lamb + tstep*result[1]
#        B_sub = B_sub + tstep*result[2]






#### OUTPUT h5

output = h5py.File(outputdest,'w')
output.create_dataset('t_list', data = reportb4_t_list)
output.create_dataset('b4_list', data = b4_list)
output.create_dataset('lamb_list', data = lamb_list)
output.create_dataset('lambshift_list', data = lambshift_list)
output.create_dataset('dAsubathor_list', data = dAsubathor_list )
output.create_dataset('sigmasubscaleathor_list', data = sigmasubscaleathor_list )
output.create_dataset('t_list_for_radial_snapshots', data = reportfunc_t_list)
output.create_dataset('lamb_list_for_radial_snapshots', data = reportfunc_lamb_list)
output.create_dataset('B_sub_list', data = B_sub_list )
output.create_dataset('sigma_sub_scale_list', data = sigma_sub_scale_list)
output.create_dataset('sigmadot_paul_list', data = sigmadot_paul_list )
output.create_dataset('Fvu_list', data = Fvu_list )
output.create_dataset('Bdot_sub_list', data = Bdot_sub_list )
output.create_dataset('A_sub_list', data = A_sub_list )
output.create_dataset('dvB_sub_list', data = dvB_sub_list )
output.create_dataset('sigmaBCcheckD_list', data = sigmaBCcheckD_list )
output.create_dataset('sigmaBCcheckN_list', data = sigmaBCcheckN_list )
output.create_dataset('sigmadotBCcheckD_list', data = sigmadotBCcheckD_list )
output.create_dataset('sigmadotBCcheckN_list', data = sigmadotBCcheckN_list )
output.create_dataset('BdotBCcheckD_list', data = BdotBCcheckD_list )
output.create_dataset('BdotBCcheckN_list', data = BdotBCcheckN_list )
output.create_dataset('ABCcheckD_list', data = ABCcheckD_list )
output.create_dataset('ABCcheckN_list', data = ABCcheckN_list )
output.create_dataset('sigmadotathor_list', data = sigmadotathor_list )
output.create_dataset('a2', data = a2)
output.create_dataset('f3', data = f3)
output.create_dataset('fxy', data = fxy)
output.create_dataset('numfuncreports', data = numfuncreports)
output.close()
tempdeste = '/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/EnergyGrid'
tempdestf = '/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/FreeEnergyGrid'      
tempdestBTZ = '/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/BTZGrid'     

if maxwelldetails in ('multi'):
    ent = Pi*(1. + lamb_list[-1] + sigmasubscaleathor_list[-1])**3
    #print "Ent:", ent, "lastlamb", lamb_list[-1], "s", sigmasubscaleathor_list[-1]
    temperature = - (-2.*(1. + lamb_list[-1]) + dAsubathor_list[-1] + ((2.*Power(fxy,2)*(1. - 2.*lamb_list[-1] + 3.*Power(lamb_list[-1],2) - 4.*Power(lamb_list[-1],3)))/3.))/(4.*pi)
    free = -3./4*(a2) - temperature*ent
    report([fxy,-3./4*(a2), temperature], tempdeste,'a',' ')   
    report([fxy, free , temperature], tempdestf,'a',' ')
else:
    pass

#plt.figure(1)
#plt.savefig(pp, format='pdf')
#pp.close()     

#plt.figure(2)
#plt.plot(reportb4_t_list, exp(3.183*reportb4_t_list)*b4_list, 'o')
#plt.title('b4list modes')


#plt.figure(3)
#plt.plot(reportb4_t_list[len(reportb4_t_list)/3:], (lamb_list[len(reportb4_t_list)/3:]-0.000000000187217),'o')
#plt.title('lamblist')





#
#plt.figure(2)
#plt.plot(reportb4_t_list, b4_list, 'o')
#plt.title('b4list')
#
#
#
#plt.figure(4)
#plt.plot(reportb4_t_list, lamb_list,'o')
#plt.title('lamblist')
#
#plt.figure(5)
#plt.plot(reportb4_t_list, sigmaBCcheckD_list,'o')
#plt.title('sigmaBCcheckD')
#
#plt.figure(6)
#plt.plot(reportb4_t_list, sigmaBCcheckN_list,'o')
#plt.title('sigmaBCcheckN')
#
#plt.figure(7)
#plt.plot(reportb4_t_list, sigmadotBCcheckD_list,'o')
#plt.title('sigmadotBCcheckD')
#
#plt.figure(8)
#plt.plot(reportb4_t_list, sigmadotBCcheckN_list,'o')
#plt.title('sigmadotBCcheckN')
#
#plt.figure(9)
#plt.plot(reportb4_t_list, BdotBCcheckD_list,'o')
#plt.title('BdotBCcheckD')
#
#plt.figure(10)
#plt.plot(reportb4_t_list, BdotBCcheckN_list,'o')
#plt.title('BdotBCcheckN')
#
#plt.figure(11)
#plt.plot(reportb4_t_list, ABCcheckD_list,'o')
#plt.title('ABCcheckD')
#
#plt.figure(12)
#plt.plot(reportb4_t_list, ABCcheckN_list,'o')
#plt.title('ABCcheckN')
#plt.draw()
