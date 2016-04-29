############ JOHN FUINI ##########



import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

from scipy import empty, zeros, cos, pi, optimize, exp, absolute
from numpy.polynomial.chebyshev import Chebyshev, chebval, chebder
from stringsformathematica import report
from scipy.special import eval_chebyu, eval_chebyt

#in order to interpret mathematica output correctly...
from numpy import power as Power
from numpy import e as E
from numpy import log as Log
from numpy import sqrt as Sqrt
from numpy import pi as Pi


# MOVING THE CHEBYCHEV TOOLBOX TO WORK ON GENERAL DOMAIN (x1,x2)

def get_slope_intercept(x1,x2):
    m = 2./(x2-x1)
    b = -1. - m*x1
    return [m,b]



def GLgrid(N,x1,x2):
    grid = empty(N+1)    # include endpoints
    m = get_slope_intercept(x1,x2)[0]
    b = get_slope_intercept(x1,x2)[1]
    for i in range(0,N+1):
        grid[i] = 1/m*(cos(pi*i/N)-b)
    return grid

def coord_change_on_grid(lambdashift,grid):
    return grid/(1.0 + lambdashift*grid)

def p(j,N):
    if j == 0 or j == N or j == -N:
        return 2.0
    elif 0 < abs(j) and abs(j) < N:
        return 1.0
    else:
        print 'ERROR: bad j'
        sys.exit(0)

def c(j,N):
    if j == 0 or j == N or j == -N:
        return 2.0
    elif 0 < abs(j) and abs(j) < N:
        return 1.0
    else:
        print 'ERROR: bad j'
        sys.exit(0)

#  Scipy.special has a better way to get chebyshev's
#def T2(j,u,N,x1,x2):
#    m = get_slope_intercept(x1,x2)[0]
#    b = get_slope_intercept(x1,x2)[1]
#    if type(j) == int:
#        if j > N:
#            print "ERROR: j > N. This shouldn't be needed."
#            sys.exit(0) 
#        temp = zeros(N+1)
#        temp[j] = 1.
#        return chebval(m*u + b,temp)
#    else:
#        tempTarray = zeros(len(j))
#        for i in range(0,len(j)): 
#            if i > N:
#                print "ERROR: j > N. This shouldn't be needed."
#                sys.exit(0)
#            tempcoefarray = zeros(N+1)
#            tempcoefarray[j[i]] = 1.
#            tempTarray[i] = chebval(m*u+b,tempcoefarray)
#        return tempTarray

def T(j,u,N,x1,x2):
    m = get_slope_intercept(x1,x2)[0]
    b = get_slope_intercept(x1,x2)[1]
    return eval_chebyt(j,m*u+b)

def U(j,u,N,x1,x2):
    m = get_slope_intercept(x1,x2)[0]
    b = get_slope_intercept(x1,x2)[1]
    return eval_chebyu(j,m*u+b)

def dT(j,u,N,x1,x2):
    m = get_slope_intercept(x1,x2)[0]
    b = get_slope_intercept(x1,x2)[1]
    return m*j*U(j-1,u,N,x1,x2)

#def dT(j,u,N,x1,x2):
#    m = get_slope_intercept(x1,x2)[0]
#    b = get_slope_intercept(x1,x2)[1]
#    if type(j) == int:
#        if j > N:
#            print "ERROR: j > N. This shouldn't be needed."
#            sys.exit(0) 
#        temp = zeros(N+1)
#        temp[j] = 1.
#        temp2 = chebder(temp)
#        return m*chebval((m*u+b),temp2)
#    else:
#        tempTarray = zeros(len(j))
#        for i in range(0,len(j)): 
#            if i > N:
#                print "ERROR: j > N. This shouldn't be needed."
#                sys.exit(0)
#            tempcoefarray = zeros(N+1)
#            tempcoefarray[j[i]] = 1.
#            tempTarray[i] = m*chebval((m*u+b),chebder(tempcoefarray))
#        return tempTarray

#def T_array(x):
#    temp = empty(N+1)
#    i = 0
#    while i < N+1:
#        temp[i] = T(i,x)
#        i = i + 1
#    return temp

def cardone(j,u,N,x1,x2):
    m = get_slope_intercept(x1,x2)[0]
    b = get_slope_intercept(x1,x2)[1]
    grid = GLgrid(N,x1,x2)
    temp = zeros(N+1)
#    if grid[j] == u:
#        return 1.
    if (absolute((grid[j] - u)) < 1.0*10**(-13)):
        return 1.
    return (-1.)**(j+1)*(1.-(m*u+b)**2)/(c(j,N)*N**2*(m*u - m*grid[j]))*(1./m)*dT(N,u,N,x1,x2)
#    temppj = zeros(N+1)
#    tempcounterarray = zeros(N+1)
#    for i in range(0,N+1):
#        tempcounterarray[i] = i
#        temppj[i] = 1/p(i,N)
#    temp = np.multiply(np.multiply(T(tempcounterarray,grid[j],N),T(tempcounterarray,x,N)),temppj)
#    return (2.0/(N*p(j,N)))*np.sum(temp)
    
def card(j,u,N,x1,x2): # this guy will array like your cardinal function
    if type(j) == int and (type(u) == int or type(u) == float or type(u) == np.float64):
        return cardone(j,u,N,x1,x2)     
    elif type(j) != int and (type(u) != int and type(u) != float and type(u) != np.float64):
        print 'ERROR: cardinal called with 2 arrays'
        sys.exit(0)
    elif type(j) != int:
        temparray = zeros(len(j))
        for i in range(0,len(j)):
            temparray[i] = cardone(j[i],u,N,x1,x2)
        return temparray
    elif (type(u) != int and type(u) != float and type(u) != np.float64):
        temparray = zeros(len(u))
        for i in range(0,len(u)):
            temparray[i] = cardone(j,u[i],N,x1,x2)
        return temparray
    else:
        print 'ERROR: Not possible'
        sys.exit(0)


def cardarray(u,N,x1,x2): 
    if type(u) == np.ndarray and len(np.atleast_1d(u)) == 1:
        u = np.float64(u)
    if type(u) != int and type(u) != float and type(u) != np.float64:
        temparray = zeros([len(u),N+1])
        for i in range(0,len(u)):
            for j in range(0,N+1):
                temparray[i,j] = card(j,u[i],N,x1,x2)
    else:
        temparray = zeros([N+1])
        for i in range(0,N+1):
            temparray[i] = card(i,u,N,x1,x2)
    return temparray    


def cardfunct(v,u,N,x1,x2):
    #Sum = 0
    #for i in range(0,N+1):
    #    Sum = Sum + v[[i]]*card(i,u,N,x1,x2)
    #return Sum
    return np.dot(cardarray(u,N,x1,x2),v)    
        

def dcard_M(N,x1,x2):
    m = get_slope_intercept(x1,x2)[0]
    b = get_slope_intercept(x1,x2)[1]
    grid = GLgrid(N,x1,x2)
    tempMatrix = zeros([N+1,N+1])
    for i in range(0,N+1):
        for j in range(0,N+1):
            if i == 0 and j == 0: 
                tempMatrix[i,j] = (m)*(1. + 2.*N**2)/6.
            elif i == N and j == N:
                tempMatrix[i,j] = -(m)*(1. + 2.*N**2)/6.
            elif i == j:
                tempMatrix[i,j] = -(m)*(m*grid[j]+b)/(2.*(1. - (m*grid[j]+b)**2))
            else:
                tempMatrix[i,j] = (-1.)**(i+j)*(m)*p(i,N)/(p(j,N)*(m*grid[i]-m*grid[j]))
    return tempMatrix

def d2card_M(N,x1,x2):
    tempMatrix = np.dot(dcard_M(N,x1,x2),dcard_M(N,x1,x2))
    return tempMatrix

def vectorize(f,N,x1,x2):
    grid = GLgrid(N,x1,x2)
    temparray = zeros(N+1)
    for i in range(0,N+1):
        temparray[i] = f(grid[i])
    return temparray

def cards_on_grid(N,M,x1,x2): #Returns N cardinals on the M points.
    temparray = zeros([M+1,N+1])
    grid = GLgrid(M,x1,x2)
    for i in range(0,M+1):
        for j in range(0,N+1): 
            temparray[i,j] = card(j,grid[i],N,x1,x2)
    return temparray

def filter_func(v,N,N2,x1,x2): #N is number for finegrid, N2 is number for coarsegrid
    NtoN2 = cards_on_grid(N,N2,x1,x2)
    N2toN = cards_on_grid(N2,N,x1,x2)
    temp = np.dot(NtoN2,v)
    v = np.dot(N2toN,temp)
    return v
    

#def Lold(q2,q1,q0,N): #The prefactor functions will always be numbers times arrays (at grid points) of other functs, so the q's need to be arrays, or I could do them as general functions just evaluated at the points, but then I'd have to interpolate my solns. Lets keep them as arrays for now. heheh @ vectorize.
#    grid = GLgrid(N)
#    dC = dcard_M(N)
#    d2C = d2card_M(N)
#    tempMatrix = zeros([N+1,N+1])
#    for i in range(0,N+1):
#        for j in range(0,N+1):
#            tempMatrix[i,j] = 1.*q2[i]*d2C[i,j] + 1.*q1[i]*dC[i,j] + 1.*q0[i]*np.identity(N+1)[i,j]
#    return tempMatrix

def L(q2,q1,q0,N,x1,x2,dpass = 'None', ddpass = 'None'): 
    if (dpass == 'None') or (ddpass == 'None'):
        dC = dcard_M(N,x1,x2)
        d2C = d2card_M(N,x1,x2)
    else:
        dC = dpass
        d2C = ddpass
    tempMatrix = (q2*(d2C.T)).T + (q1*(dC.T)).T + (q0*np.identity(N+1)).T
    return tempMatrix

def spectral_solve(L,F,(BC1type,BC1pos,BC1val),(BC2type,BC2pos,BC2val),N,x1,x2,dpass = 'None', ddpass = 'None',cntrltype = None, cntrlBCloc = None, cntrlrowreplace = None, cntrlval = None ): #f should be vectorized, the positions of the BC should be given as a number for which the grid[BCpos] = actual position # Be careful with controlling the replacement.  If an automatic replacement conflicts with your direct control replacement you will not get the right answer, and this function will not complain.
    if (dpass == 'None') or (ddpass == 'None'):
        dC = dcard_M(N,x1,x2)
        d2C = d2card_M(N,x1,x2)
    else:
        dC = dpass
        d2C = ddpass
    grid = GLgrid(N,x1,x2)
    tempSoln = zeros(N+1)
    if cntrltype == 'D':
        for j in range(0,N+1):
            L[cntrlrowreplace,j] = 0    
        L[cntrlrowreplace, cntrlBCloc] = 1.
        F[cntrlrowreplace] = cntrlval
    elif cntrltype == 'N':
        for j in range(0,N+1):
            L[cntrlrowreplace,j] = dC[cntrlBCloc,j]
        F[BC1pos] = cntrlval
    elif cntrltype != None:
        print "ERROR: Please select N, D, or None for cntrltype"
        sys.exit(0)

    if BC1pos == cntrlrowreplace or BC2pos == cntrlrowreplace:
        print "ERROR: It seems that your automatic BC selections interfere with your direct control replacement."
        sys.exit(0)

    if (BC1type != 'D' and BC1type != 'N'and BC1type != 'na') or (BC2type != 'D' and BC2type != 'N' and BC2type != 'na'):
        print 'ERROR: Must chose Dirichel (D) or Nuemann (N) or (na) for both BCs'
        sys.exit(0)
    elif BC2type == 'na':
        if BC1type == 'D':
            for j in range(0,N+1):
                L[BC1pos,j] = 0
            L[BC1pos,BC1pos] = 1.
            F[BC1pos] = BC1val
        elif BC1type == 'N':
            for j in range(0,N+1):
                L[BC1pos,j] = dC[BC1pos,j]
            F[BC1pos] = BC1val 
        tempSoln = np.linalg.solve(L,F)
        return tempSoln
    elif BC1type == 'N' and BC2type == 'D':
        print "ERROR: Don't do that.  Make BC1 Dirchlet, and BC2 Nuemann." 
        sys.exit(0)
    elif BC1type == 'D' and BC2type == 'D':
        if BC1pos == BC2pos:
            print "ERROR: Can't put two Dirchet conditions at the same point smart guy..."
            sys.exit(0)
        # So lets remove the row that is where the BC is.
        for j in range(0,N+1):
            L[BC1pos,j] = 0
            L[BC2pos,j] = 0
        L[BC1pos,BC1pos] = 1.
        L[BC2pos,BC2pos] = 1. 
        F[BC1pos] = BC1val
        F[BC2pos] = BC2val
        tempSoln = np.linalg.solve(L,F)
        return tempSoln
    elif BC1type == 'D' and BC2type == 'N':
        BC2remove = BC2pos
        if BC1pos == BC2pos:
            if BC1pos == N:
                BC2remove = BC1pos - 1
            else:
                BC2remove = BC1pos + 1
        for j in range(0,N+1):
            L[BC1pos,j] = 0
            L[BC2remove,j] = dC[BC2pos,j] 
        L[BC1pos,BC1pos] = 1.
        F[BC1pos] = BC1val
        F[BC2remove] = BC2val   #change this all out for a new set of removal and BCloc.
        tempSoln = np.linalg.solve(L,F)
        return tempSoln
    elif BC1type == 'N' and BC2type == 'N':
        if BC1pos == BC2pos:
            print "ERROR: Can't put two Nuemann conditions at the same point captain..." 
        for j in range(0,N+1):
            L[BC1pos,j] = dC[BC1pos,j]
            L[BC2pos,j] = dC[BC2pos,j]
        F[BC1pos] = BC1val
        F[BC2pos] = BC2val
        tempSoln = np.linalg.solve(L,F)
    elif BC1type == 'na' and BC2type == 'na':
        tempsoln = np.linalg.solve(L,F)


################################################END BEFORE PHYSICS SPECIFIC PROBLEM################################################

#### Find the Dots for future time stepping

def find_dots(grid,lamb, B_sub, g, a2, f3, fxy, N, x1, x2,filt = 'None',Nfilter = 'None', horizonfind = False,dpass = 'None', ddpass = 'None', lambdafix = False, sigmadottol = "None", backtracklamb = False, horizontol = 'None', pr=False, destination='Null',pl=False , time = '0', lambdafindoption = False):
#    if fixlambdrift == True and coordchange == True:
#        print "Do not try to fix lambda independently AND do the coord shift (the former is already done in the latter)"
#        sys.exit(0)
#    if fixlambdrift == True and backtracklamb == True:
#        print "Do not try to fix lambda independently AND do the backtracklamb fix"
#        sys.exit(0)
#    if backtracklamb == True and coordchange == True:
#        print "Do not try to do the coord shift AND do the backtracklamb fix"
#        sys.exit(0)
    res = -1
    #report([lamb],'/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/ComparetoPaul.txt','a','PYlamb')

    finegrid = np.linspace(.01,1.0,1000)
    foundhorizon = True

    if (dpass == 'None') or (ddpass == 'None'):
        dC = dcard_M(N,x1,x2)
        d2C = d2card_M(N,x1,x2)
    else:
        dC = dpass
        d2C = ddpass

    if backtracklamb == True:
        backtrackcounter = 1
    else: 
        backtrackcounter = 0

    comparedest = '/phys/users/fuini/Research/Dropbox/Yaffe/Python/Output/ComparetoPaul.txt'
    for i in range(0,backtrackcounter+1):
        
    #B_sub1
#        B = grid**3*B_sub # not true anymore
       # print "B_sub before dB is", B_sub
        dB_sub = np.dot(dC,B_sub)
        b4 = np.dot(dC,B_sub)[N]
       # print "dB_sub IS", dB_sub
       # print "b4= ", b4
        
        B = B_sub*Power(grid,3) + ((Power(fxy,2)*Power(grid,4))/(3.*Power(g,2)) - (4.*Power(fxy,2)*Power(grid,5)*lamb)/(3.*Power(g,2)) + (10.*Power(fxy,2)*Power(grid,6)*Power(lamb,2))/(3.*Power(g,2)))*Log(1./grid)
        B[N] = 0 #know this analytically, may as well enforce it.        
        if lambdafix == True:
            print "b4: ",b4
#        print "lambda:" ,lamb

        #report(B_sub, comparedest, 'a', 'PYBpaul')
 

    #### SIGMA
        
    #### new lambda fix
        bub = 0
        lambbefore = lamb
        lambshift = 0
        while True:
            bub = bub + 1
            
            vq2 = vectorize(lambda u: 1.,N,x1,x2)
            vq1 = 12./grid        
            vq0fxy0 = (60. + Power(grid,6)*Power(3.*B_sub + dB_sub*grid,2))/(2.*Power(grid,2))
            vq0fxy2 = (Power(grid,5)*(3.*B_sub + dB_sub*grid)*(-1. + 2.*grid*lamb*(2. - 5.*grid*lamb) + 4.*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid)))/(3.*Power(g,2))
            vq0fxy4 = (Power(grid,6)*Power(-1. + 2.*grid*lamb*(2. - 5.*grid*lamb) + 4.*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid),2))/(18.*Power(g,4))
            vq0 = vq0fxy0 + vq0fxy2*fxy**2 + vq0fxy4*fxy**4
            vF0 = -(Power(3.*B_sub + dB_sub*grid,2)*(1. + grid*lamb))/(2.*Power(grid,2))
            vF2 = -((3.*B_sub + dB_sub*grid)*(1. + grid*lamb)*(-1. + 2.*grid*lamb*(2. - 5.*grid*lamb) + 4.*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid)))/(3.*Power(g,2)*grid)
            vF4 = -((1. + grid*lamb)*Power(-1. + 2.*grid*lamb*(2. - 5.*grid*lamb) + 4.*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid),2))/(18.*Power(g,4))
            vF = vF0 + vF2*fxy**2 + vF4*fxy**4
         
         



            LM = L(vq2,vq1,vq0,N,x1,x2,dpass = dC, ddpass = d2C)

            sigma_sub_scale = spectral_solve(LM,vF,('D',N,0),('N',N,0),N,x1,x2, dpass = dC, ddpass = d2C)
            dsigma_sub_scale = np.dot(dC,sigma_sub_scale)
            sigma = 1./grid + lamb + (grid**5)*sigma_sub_scale
            
            sigma_paul = lamb + grid**5*sigma_sub_scale
            #report(sigma_paul, comparedest, 'a', 'PYSigmapaul')

        #### FVU

        # from the radial derivative equation for Fvr, we can easily solve for Fvr as a function of s

        #Fvr = f3*(sigma[N]**3)*sigma[N]**(-3) Can't figure out how to use this yet

#            vq2 = vectorize(lambda u: 0,N,x1,x2)
#            vq1 = vectorize(lambda u: 1.,N,x1,x2)
#            vq0 =(-1. + 2.*lamb*grid + 17.*sigma_sub_scale*(grid**6) + 3.*dsigma_sub_scale*(grid**7))/(grid*(1. + lamb*grid + sigma_sub_scale*(grid**6)))
#            vF = vectorize(lambda u: 0,N,x1,x2)

#            LM = L(vq2,vq1,vq0,N,x1,x2,dpass = dC, ddpass = d2C)

#            Fvu = spectral_solve(LM,vF,('D',N,0),('na',N,-f3),N,x1,x2,dpass = dC, ddpass = d2C)
            

            # Don't need to solve spectrally.
            Fvu = -f3/(grid**2*sigma**3)
            Fvu[N] = 0 



        #### SIGMA DOT


    ######## If we want to exactly match pauls choice of scalings then... ######################

            vq2 = vectorize(lambda u: 0,N,x1,x2)
            vq1 = vectorize(lambda u: 1,N,x1,x2)
            vq0 = (-1. + 2.*dsigma_sub_scale*Power(grid,7) + grid*lamb + 11.*Power(grid,6)*sigma_sub_scale)/(grid*(1. + grid*lamb + Power(grid,6)*sigma_sub_scale))            
            vF0 = (2.*(1. + grid*lamb + Power(grid,6)*sigma_sub_scale)*(Power(Fvu,2) - 6.*Power(g,2)*Power(grid,2)*(dsigma_sub_scale*grid + 5.*sigma_sub_scale)))/(3.*Power(g,2))            
            #vF2 = (2.*(1./(Power(E,2.*B_sub*Power(grid,3))*Power(1./grid,(2.*Power(fxy,2)*Power(grid,4)*(1. + 2.*grid*lamb*(-2. + 5.*grid*lamb)))/(3.*Power(g,2)))) + (-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb)))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3) - 2.*Power(grid,4)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2)*(10.*Power(lamb,4) + dsigma_sub_scale*Power(grid,3)*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb))) + Power(grid,2)*(-6. + grid*lamb*(13. + 3.*grid*lamb*(-7. + 10.*grid*lamb)))*sigma_sub_scale)*Log(1./grid)))/(3.*Power(g,2)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3))
# this vF2 sends 1/u to u to negative power            
            vF2 = (2.*(1./(Power(E,2.*B_sub*Power(grid,3))*Power(grid,-(2.*Power(fxy,2)*Power(grid,4)*(1. + 2.*grid*lamb*(-2. + 5.*grid*lamb)))/(3.*Power(g,2)))) + (-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb)))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3) - 2.*Power(grid,4)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2)*(10.*Power(lamb,4) + dsigma_sub_scale*Power(grid,3)*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb))) + Power(grid,2)*(-6. + grid*lamb*(13. + 3.*grid*lamb*(-7. + 10.*grid*lamb)))*sigma_sub_scale)*Log(1./grid)))/(3.*Power(g,2)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3))
            vF = vF0 + vF2*fxy**2
           
            LM = L(vq2,vq1,vq0,N,x1,x2,dpass = dC, ddpass = d2C)
            #sigmadot_paul = spectral_solve(LM,vF,('N', N, a2),('na',N,'None'),N,x1,x2,dpass = dC, ddpass = d2C)
            sigmadot_paul = spectral_solve(LM,vF,('D',N,'0'),('N', N, a2),N,x1,x2,dpass = dC, ddpass = d2C)            
            dsigmadot_paul = np.dot(dC, sigmadot_paul)

            #### back to my variables
            sigmadot_sub_scale = 1./2.*(sigmadot_paul - ((-((grid**3)*sigma_sub_scale*(2. + 2.*grid*lamb + (grid**6)*sigma_sub_scale)))))
#           # report(sigmadot_paul, comparedest, 'a', 'PYsigmadot_paul')
            #temp = 2.*sigmadot_sub_scale + (-((grid**3)*sigma_sub_scale*(2. + 2.*grid*lamb + (grid**6)*sigma_sub_scale)))

 
    ############################################################################################


           # sigmadot =   1./(2.*grid**2) + lamb/grid + Power(lamb,2)/2. + (-((fxy**2)*grid**2)/(3.*(g**2)) + (2.*(fxy**2)*grid**3*lamb)/(3.*(g**2)))*Log(1./grid) + grid*sigmadot_sub_scale

            def sigmadotfunct(x): 
                return (cardfunct(sigmadot_paul,x,N,x1,x2)*x)/2. + Power(lamb + 1./x + cardfunct(sigma_sub_scale,x,N,x1,x2)*Power(x,5),2)/2. + (-(Power(fxy,2)*Power(x,2))/(3.*Power(g,2)) + (2.*Power(fxy,2)*lamb*Power(x,3))/(3.*Power(g,2)) - (Power(fxy,2)*Power(lamb,2)*Power(x,4))/Power(g,2) + (4.*Power(fxy,2)*Power(lamb,3)*Power(x,5))/(3.*Power(g,2)))*Log(1./x)
                
            if bub in range(0):
                finegrid = np.linspace(.3,1.0,100)
                tempmax = max(abs(sigmadotfunct(finegrid[50:100])))
                tempmin = min(abs(sigmadotfunct(finegrid[50:100])))
                fig1 = plt.figure(figsize=(20,12))
                plot1 = fig1.add_subplot(231)
                plot1.plot(finegrid[50:100],sigmadotfunct(finegrid[50:100]),'o',finegrid[50:100], zeros(len(finegrid[50:100])),'-')
                plt.axis([.5, 1.01,-abs(tempmax + .2*tempmax), abs(tempmax + .2*tempmax)])
                plt.show() 
    #### FINDING HORIZON
             
            if horizonfind == True:



#                 gridforhorizon = np.linspace(.99, 1.01, 50)
#                 maxendzeroSigmadot = 1.001  
               # Finding first zero of Sigma
#                 for i in range(len(gridforhorizon)):
#                     if (sigmadotfunct(gridforhorizon[i])) < 0:
#                        maxendzeroSigmadot = gridforhorizon[i]
    #                    print 'found it! pos -> neg!' ; 
#                        break


    #             try:
    #                 horizon = sp.optimize.broyden1(sigmadotfunct, 0.9, f_tol=horizontol) 
    #                 horizon = np.float64(horizon)  
    #                 print "HORIZON:", horizon
    #                 print "TOLERANCE:", horizontol

    #             except Exception:
    #                 print "Horizon Not found!"
    #                 foundhorizon = False
    #                 horizon = 1.
    #                 print "Sigmadot at 1 is:", sigmadotfunct(1)                
    #                 report(sigmadotfunct(finegrid),destination,'a','ErrorSigmadotfunct')
    #                 pass

                 try:
#                     horizonalt = optimize.brentq(sigmadotfunct, .99, maxendzeroSigmadot) 
                     print "HORIZON alt", horizonalt
                     horizon = horizonalt
                 except Exception:
                     horizonalt = -1
                     print "HORIZON alt", horizonalt
                     horizon = 1.0
                     pass
            else:
                horizon = None 
                horizonalt = None
               

#            if backtracklamb == True and i == 0:
#                print "***Back track lambda fix***"
#                lambshift = (1. - horizon)/(horizon*1.)
#                lamb = lamb + lambshift
#                print "lambshift = ", lambshift
#                B_sub = cardfunct(B_sub,coord_change_on_grid(lambshift,grid),N,x1,x2)
            
            #### Fixing lambda
    
            
            if lambdafix != True:
                break
 
            #res1 = sigmadotfunct(1)
            res = sigmadot_paul[0]/2. + (1. + lamb + sigma_sub_scale[0])**2/2.
            # withouth MORE LOGS Lp2 = -((Power(fxy,2)*(1. - 2.*lamb))/(3.*Power(g,2)) + (3.*(dsigmadot_paul[0] + sigmadot_paul[0]) + 6.*(1. + lamb + sigma_sub_scale[0])*(-1. + dsigma_sub_scale[0] + 5.*sigma_sub_scale[0]))/6.)
            #lambshift1 = (res)/((np.dot(dC,sigmadot_paul)[0]*grid[0] + sigmadot_paul[0] + (2.*(1. + grid[0]*lamb + Power(grid[0],6)*sigma_sub_scale[0])*(-1. + np.dot(dC, sigma_sub_scale)[0]*Power(grid[0],7) + 5.*Power(grid[0],6)*sigma_sub_scale[0]))/Power(grid[0],3))/2.)
            Lp = (Power(fxy,2)*(-2. + 2.*lamb*(2. + lamb*(-3. + 4.*lamb))) - 3.*Power(g,2)*(-2. + dsigmadot_paul[0] - 2.*lamb + sigmadot_paul[0] + 8.*sigma_sub_scale[0] + 2.*(5.*sigma_sub_scale[0]*(lamb + sigma_sub_scale[0]) + dsigma_sub_scale[0]*(1. + lamb + sigma_sub_scale[0]))))/(6.*Power(g,2))
            #print "RES:", res,"(-L):", -Lp

#### METHOD PAUL
            if abs(res) < sigmadottol:
                print "fixed lamb!!", lambbefore,"->", lamb, "   Shift:", lamb - lambbefore
                print "residual:", res
                if time == 0:
                    if lambdafindoption == False:
                        if abs(lamb - lambbefore) > 10**(-10):
                            print "Sir, the initial lambda wasn't correct.  I advise you get a good starting lambda from the Lambda Find program"
                            #sys.exit(0)
                            # RECENT NOTE:  This should be here, unless using custom, then I want it off. :(
                break
            lambshift = -(1./3)*(res)/Lp # this has a slowing factor out front, to make convergence better.  Should be investigated as to why it's needed.



            
#            print "fixing lamb...", lamb, "->", lamb + lambshift, "res before:", res #"dlamb = ", lambshift,
            lamb = lamb + lambshift
            #lamb = lamb + lambshift
            
   
            
    #### Fix Horizon drift
    #    if coordchange == True:
    #        print "Shifting coords..."
    #        lambshift = (1. - horizon)/(horizon*1.)
    #        lamb = lamb + lambshift
    #        print "lambshift = ", lambshift
    #        # Shift known functions to correspond correctly to horizon at 1.
    #        B_sub = cardfunct(B_sub,coord_change_on_grid(lambshift,grid),N,x1,x2)
    #        sigma_sub_scale = cardfunct(sigma_sub_scale,coord_change_on_grid(lambshift,grid),N,x1,x2)
    #        Fvu = cardfunct(Fvu,coord_change_on_grid(lambshift,grid),N,x1,x2)
    #        sigmadot_sub_scale = cardfunct(sigmadot_sub_scale,coord_change_on_grid(lambshift,grid),N,x1,x2)


    #        # derived functions reinitialized correctly
    #        B = grid**3*B_sub
    #        dB = np.dot(dC,B) 
    #        b4 = np.dot(dC,B_sub)[N]
    #        dsigma_sub_scale = np.dot(dC,sigma_sub_scale)
    #        sigma = 1./grid + lamb + grid**5*sigma_sub_scale
    #        sigmadot = 1./(2.*grid**2) + lamb/grid + Power(lamb,2)/2. + (-((fxy**2)*grid**2)/(3.*(g**2)) + (2.*(fxy**2)*grid**3*lamb)/(3.*(g**2)))*Log(1./grid) + grid*sigmadot_sub_scale
        

               

    
#### BDOT Now using Bdot_sub_scale (has logs removed).
  

### without MORE LOGS and with out subtracted logs  #vq2 = vectorize(lambda x: 0,N,x1,x2)
                                                    #vq1 = vectorize(lambda x: 1,N,x1,x2)
                                                    #vq0 = (1. + 3.*dsigma_sub_scale*(grid**7) + 4.*grid*lamb + 19.*(grid**6)*sigma_sub_scale)/(2.*grid + 2.*(grid**2)*lamb + 2.*(grid**7)*sigma_sub_scale)
                                                    #vF = (-3.*(3.*B_sub + dB_sub*grid)*((grid**3)*sigmadot_paul + Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,2)))/(4.*grid*(1 + grid*lamb + (grid**6)*sigma_sub_scale)) + (Power(fxy,2)*(4./Power(E,2.*B_sub*(grid**3)) - 3.*(grid**3)*(3.*B_sub + dB_sub*grid)*(-1. + 2.*grid*lamb)*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,3)*Log(1./grid)))/(6.*Power(g,2)*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,4))

    vq2 = vectorize(lambda x: 0,N,x1,x2)
    vq1 = vectorize(lambda x: 1,N,x1,x2)
    vq0 = (1. + 3.*dsigma_sub_scale*Power(grid,7) + 4.*grid*lamb + 19.*Power(grid,6)*sigma_sub_scale)/(2.*grid + 2.*Power(grid,2)*lamb + 2.*Power(grid,7)*sigma_sub_scale)    
    vF0 = (-3.*(3.*B_sub + dB_sub*grid)*(1. + 2.*grid*lamb + Power(grid,2)*Power(lamb,2) + Power(grid,3)*sigmadot_paul + 2.*Power(grid,6)*sigma_sub_scale + 2.*Power(grid,7)*lamb*sigma_sub_scale + Power(grid,12)*Power(sigma_sub_scale,2)))/(4.*grid*(1. + grid*lamb + Power(grid,6)*sigma_sub_scale))   
    #vF2 = -(-8./(Power(E,2.*B_sub*Power(grid,3))*Power(1./grid,(2.*Power(fxy,2)*Power(grid,4)*(1. - 4.*grid*lamb + 10.*Power(grid,2)*Power(lamb,2)))/(3.*Power(g,2)))) - Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3)*(-5. + 10.*grid*lamb - 15.*Power(grid,2)*Power(lamb,2) + 30.*Power(grid,5)*Power(lamb,2)*sigmadot_paul + Power(grid,3)*(80.*Power(lamb,3) + 3.*sigmadot_paul) + 2.*Power(grid,4)*(55.*Power(lamb,4) - 6.*lamb*sigmadot_paul) - 2.*Power(grid,6)*sigma_sub_scale + 6.*Power(grid,7)*lamb*sigma_sub_scale - 12.*Power(grid,8)*Power(lamb,2)*sigma_sub_scale + 140.*Power(grid,9)*Power(lamb,3)*sigma_sub_scale + 3.*Power(grid,12)*Power(sigma_sub_scale,2) - 12.*Power(grid,13)*lamb*Power(sigma_sub_scale,2) + 30.*Power(grid,14)*Power(lamb,2)*Power(sigma_sub_scale,2)) + 6.*Power(grid,3)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3)*(3.*B_sub*(-1. + 2.*grid*lamb - 3.*Power(grid,2)*Power(lamb,2) + 4.*Power(grid,3)*Power(lamb,3)) + dB_sub*grid*(-1. + 2.*grid*lamb - 3.*Power(grid,2)*Power(lamb,2) + 4.*Power(grid,3)*Power(lamb,3)) + 2.*(55.*grid*Power(lamb,4) + dsigma_sub_scale*Power(grid,4)*(-1. + 3.*grid*lamb - 6.*Power(grid,2)*Power(lamb,2) + 10.*Power(grid,3)*Power(lamb,3)) + sigmadot_paul - 5.*Power(grid,3)*sigma_sub_scale + Power(grid,9)*Power(sigma_sub_scale,2) + 5.*Power(lamb,3)*(7. + 24.*Power(grid,6)*sigma_sub_scale) - 5.*grid*lamb*(sigmadot_paul + Power(grid,3)*sigma_sub_scale*(-3. + Power(grid,6)*sigma_sub_scale)) + 15.*Power(grid,2)*Power(lamb,2)*(sigmadot_paul + Power(grid,3)*sigma_sub_scale*(-2. + Power(grid,6)*sigma_sub_scale))))*Log(1./grid))/(12.*Power(g,2)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,4))
 ### this vf2 sends 1/u to u to negative power  
    vF2 = -(-8./(Power(E,2.*B_sub*Power(grid,3))*Power(grid,-(2.*Power(fxy,2)*Power(grid,4)*(1. - 4.*grid*lamb + 10.*Power(grid,2)*Power(lamb,2)))/(3.*Power(g,2)))) - Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3)*(-5. + 10.*grid*lamb - 15.*Power(grid,2)*Power(lamb,2) + 30.*Power(grid,5)*Power(lamb,2)*sigmadot_paul + Power(grid,3)*(80.*Power(lamb,3) + 3.*sigmadot_paul) + 2.*Power(grid,4)*(55.*Power(lamb,4) - 6.*lamb*sigmadot_paul) - 2.*Power(grid,6)*sigma_sub_scale + 6.*Power(grid,7)*lamb*sigma_sub_scale - 12.*Power(grid,8)*Power(lamb,2)*sigma_sub_scale + 140.*Power(grid,9)*Power(lamb,3)*sigma_sub_scale + 3.*Power(grid,12)*Power(sigma_sub_scale,2) - 12.*Power(grid,13)*lamb*Power(sigma_sub_scale,2) + 30.*Power(grid,14)*Power(lamb,2)*Power(sigma_sub_scale,2)) + 6.*Power(grid,3)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,3)*(3.*B_sub*(-1. + 2.*grid*lamb - 3.*Power(grid,2)*Power(lamb,2) + 4.*Power(grid,3)*Power(lamb,3)) + dB_sub*grid*(-1. + 2.*grid*lamb - 3.*Power(grid,2)*Power(lamb,2) + 4.*Power(grid,3)*Power(lamb,3)) + 2.*(55.*grid*Power(lamb,4) + dsigma_sub_scale*Power(grid,4)*(-1. + 3.*grid*lamb - 6.*Power(grid,2)*Power(lamb,2) + 10.*Power(grid,3)*Power(lamb,3)) + sigmadot_paul - 5.*Power(grid,3)*sigma_sub_scale + Power(grid,9)*Power(sigma_sub_scale,2) + 5.*Power(lamb,3)*(7. + 24.*Power(grid,6)*sigma_sub_scale) - 5.*grid*lamb*(sigmadot_paul + Power(grid,3)*sigma_sub_scale*(-3. + Power(grid,6)*sigma_sub_scale)) + 15.*Power(grid,2)*Power(lamb,2)*(sigmadot_paul + Power(grid,3)*sigma_sub_scale*(-2. + Power(grid,6)*sigma_sub_scale))))*Log(1./grid))/(12.*Power(g,2)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,4))    
    vF4 = -(Power(grid,4)*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb)))*Log(1./grid)*(-1. + 2.*grid*lamb*(2. - 5.*grid*lamb) + 4.*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid)))/(6.*Power(g,4)*(1. + grid*lamb + Power(grid,6)*sigma_sub_scale)) 
    vF = vF0 + vF2*fxy**2 + vF4*fxy**4
   

    LM = L(vq2,vq1,vq0,N,x1,x2,dpass = dC, ddpass = d2C)
    Bdot_sub_scale = spectral_solve(LM,vF,('D',N,0),('N',N,(-2.*b4 + fxy**2/(6*g**2))),N,x1,x2,dpass = dC, ddpass = d2C)#(-2.*b4)
    dBdot_sub_scale = np.dot(dC,Bdot_sub_scale)




    Bdot = Bdot_sub_scale*(grid**2) + ((-2.*Power(fxy,2)*(grid**3))/(3.*Power(g,2)) + (2.*Power(fxy,2)*(grid**4)*lamb)/Power(g,2) - (4.*Power(fxy,2)*(grid**5)*Power(lamb,2))/Power(g,2) + (20.*Power(fxy,2)*(grid**6)*Power(lamb,3))/(3.*Power(g,2)))*Log(1./grid)

  



#### A

#### without MORE LOGS    #vq2 = vectorize(lambda x: 1.,N,x1,x2)
                            #vq1 = 2./grid
                            #vq0 = vectorize(lambda x: 0.,N,x1,x2)
                           #vF = (9.*Bdot_su(-27*B_sub*Power(E,2*B_sub*(grid**3))*Power(g,2)*Power(1 + grid*lamb + (grid**6)*sigma_sub_scale,3)*(1 + 2*grid*lamb + (grid**2)*Power(lamb,2) + (grid**3)*sigmadot_paul + 2*(grid**6)*sigma_sub_scale + 2*(grid**7)*lamb*sigma_sub_scale + Power(grid,12)*Power(sigma_sub_scale,2)) - grid*(9*dB_sub*Power(E,2*B_sub*(grid**3))*Power(g,2)*Power(1 + grid*lamb + (grid**6)*sigma_sub_scale,3)*(1 + 2*grid*lamb + (grid**2)*Power(lamb,2) + (grid**3)*sigmadot_paul + 2*(grid**6)*sigma_sub_scale + 2*(grid**7)*lamb*sigma_sub_scale + Power(grid,12)*Power(sigma_sub_scale,2)) - 8*Power(fxy,2)*(1 + Power(E,2*B_sub*(grid**3))*(-1 + 3*grid*lamb - 6*(grid**2)*Power(lamb,2) + 10*(grid**3)*Power(lamb,3))*Power(1 + grid*lamb + (grid**6)*sigma_sub_scale,4))) - 6*Power(E,2*B_sub*(grid**3))*Power(fxy,2)*grid*Power(1 + grid*lamb + (grid**6)*sigma_sub_scale,3)*(3*B_sub*(grid**3)*(-1 + 2*grid*lamb - 3*(grid**2)*Power(lamb,2) + 4*(grid**3)*Power(lamb,3)) + dB_sub*(grid**4)*(-1 + 2*grid*lamb - 3*(grid**2)*Power(lamb,2) + 4*(grid**3)*Power(lamb,3)) + 2*(-1 + 3*grid*lamb - 6*(grid**2)*Power(lamb,2) + 10*(grid**3)*Power(lamb,3) + 40*(grid**4)*Power(lamb,4) + dsigma_sub_scale*(grid**7)*(-1 + 3*grid*lamb - 6*(grid**2)*Power(lamb,2) + 10*(grid**3)*Power(lamb,3)) - 7*(grid**6)*sigma_sub_scale + 23*(grid**7)*lamb*sigma_sub_scale - 50*Power(grid,8)*Power(lamb,2)*sigma_sub_scale + 90*Power(grid,9)*Power(lamb,3)*sigma_sub_scale))*Log(1/grid))/(12.*Power(E,2*B_sub*(grid**3))*Power(g,2)*grid*Power(1 + grid*lamb + (grid**6)*sigma_sub_scale,4))b*Power(g,2)*(grid**3)*(3.*B_sub + dB_sub*grid)*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,2) - 2.*(-7.*Power(Fvu,2)*grid*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,2) + 9.*Power(g,2)*((grid**3)*(dsigma_sub_scale*grid + 5.*sigma_sub_scale)*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,2) + sigmadot_paul*(-1. + dsigma_sub_scale*(grid**7) + 5.*(grid**6)*sigma_sub_scale))))/(3.*Power(g,2)*grid*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,2)) + (2.*Power(fxy,2)*(5. - Power(E,2.*B_sub*(grid**3))*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,2)*(5. + 2.*(grid**3)*Power(lamb,3)*(-7. + 12.*Log(1./grid)) + (grid**6)*(sigma_sub_scale*(10. + (grid**6)*sigma_sub_scale*(5. - 6.*Log(1/grid)) - 42.*Log(1./grid)) - 6.*dsigma_sub_scale*grid*Log(1./grid)) + (grid**2)*Power(lamb,2)*(-23. + 42.*Log(1/grid) + 4.*(grid**6)*sigma_sub_scale*(-7. + 12.*Log(1./grid))) + 2.*grid*lamb*(-2. + (grid**6)*(6.*dsigma_sub_scale*grid*Log(1./grid) + sigma_sub_scale*(-9. + 48.*Log(1./grid) + (grid**6)*sigma_sub_scale*(-7. + 12.*Log(1./grid))))))))/(3.*Power(E,2*B_sub*(grid**3))*Power(g,2)*Power(1. + grid*lamb + (grid**6)*sigma_sub_scale,4))
    
    vq2 = vectorize(lambda x: 1.,N,x1,x2)
    vq1 = 2./grid
    vq0 = vectorize(lambda x: 0.,N,x1,x2)
    vF0 = (9.*Power(g,2)*(2.*sigmadot_paul + Power(grid,3)*(Bdot_sub_scale*(3.*B_sub + dB_sub*grid)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) - 2.*(dsigma_sub_scale*grid + 5.*sigma_sub_scale)*(1. + 2.*grid*lamb + Power(grid,2)*Power(lamb,2) + Power(grid,3)*sigmadot_paul + 2.*Power(grid,6)*(1. + grid*lamb)*sigma_sub_scale + Power(grid,12)*Power(sigma_sub_scale,2)))) + 14.*grid*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2)*Power(Fvu,2))/(3.*Power(g,2)*grid*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2))
    #vF2 = -(-10./(Power(E,2.*B_sub*Power(grid,3))*Power(1./grid,(2.*Power(fxy,2)*Power(grid,4)*(1. - 4.*grid*lamb + 10.*Power(grid,2)*Power(lamb,2)))/(3.*Power(g,2)))) + (10. - 28.*grid*lamb + 54.*Power(grid,2)*Power(lamb,2) - 88.*Power(grid,3)*Power(lamb,3) + 3.*Bdot_sub_scale*Power(grid,3)*(1. - 4.*grid*lamb + 10.*Power(grid,2)*Power(lamb,2)))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,4) - 6.*Power(grid,3)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2)*(2.*Bdot_sub_scale*(1. - 5.*grid*lamb + 15.*Power(grid,2)*Power(lamb,2))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) + 3.*B_sub*(-1. + 3.*grid*lamb - 6.*Power(grid,2)*Power(lamb,2) + 10.*Power(grid,3)*Power(lamb,3))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) + grid*(dB_sub*(-1. + 3.*grid*lamb - 6.*Power(grid,2)*Power(lamb,2) + 10.*Power(grid,3)*Power(lamb,3))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) - 2.*(20.*grid*Power(lamb,5) + dsigma_sub_scale*Power(grid,3)*(-1. + 2.*grid*lamb - 3.*Power(grid,2)*Power(lamb,2) + 4.*Power(grid,3)*Power(lamb,3)) + 20.*Power(grid,5)*Power(lamb,3)*sigma_sub_scale*(2. + Power(grid,6)*sigma_sub_scale) + 4.*Power(grid,3)*lamb*sigma_sub_scale*(4. + Power(grid,6)*sigma_sub_scale) - Power(grid,2)*sigma_sub_scale*(7. + Power(grid,6)*sigma_sub_scale) + 10.*Power(lamb,4)*(3. + 4.*Power(grid,6)*sigma_sub_scale) - Power(grid,4)*Power(lamb,2)*sigma_sub_scale*(27. + 10.*Power(grid,6)*sigma_sub_scale))))*Log(1./grid))/(3.*Power(g,2)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,4))  
   # this vf2 sends 1/u to u to negative power 
    vF2 = -(-10./(Power(E,2.*B_sub*Power(grid,3))*Power(grid,-(2.*Power(fxy,2)*Power(grid,4)*(1. - 4.*grid*lamb + 10.*Power(grid,2)*Power(lamb,2)))/(3.*Power(g,2)))) + (10. - 28.*grid*lamb + 54.*Power(grid,2)*Power(lamb,2) - 88.*Power(grid,3)*Power(lamb,3) + 3.*Bdot_sub_scale*Power(grid,3)*(1. - 4.*grid*lamb + 10.*Power(grid,2)*Power(lamb,2)))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,4) - 6.*Power(grid,3)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2)*(2.*Bdot_sub_scale*(1. - 5.*grid*lamb + 15.*Power(grid,2)*Power(lamb,2))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) + 3.*B_sub*(-1. + 3.*grid*lamb - 6.*Power(grid,2)*Power(lamb,2) + 10.*Power(grid,3)*Power(lamb,3))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) + grid*(dB_sub*(-1. + 3.*grid*lamb - 6.*Power(grid,2)*Power(lamb,2) + 10.*Power(grid,3)*Power(lamb,3))*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,2) - 2.*(20.*grid*Power(lamb,5) + dsigma_sub_scale*Power(grid,3)*(-1. + 2.*grid*lamb - 3.*Power(grid,2)*Power(lamb,2) + 4.*Power(grid,3)*Power(lamb,3)) + 20.*Power(grid,5)*Power(lamb,3)*sigma_sub_scale*(2. + Power(grid,6)*sigma_sub_scale) + 4.*Power(grid,3)*lamb*sigma_sub_scale*(4. + Power(grid,6)*sigma_sub_scale) - Power(grid,2)*sigma_sub_scale*(7. + Power(grid,6)*sigma_sub_scale) + 10.*Power(lamb,4)*(3. + 4.*Power(grid,6)*sigma_sub_scale) - Power(grid,4)*Power(lamb,2)*sigma_sub_scale*(27. + 10.*Power(grid,6)*sigma_sub_scale))))*Log(1./grid))/(3.*Power(g,2)*Power(1. + grid*lamb + Power(grid,6)*sigma_sub_scale,4))      
    vF4 = (2.*Power(grid,4)*(-1. + grid*lamb*(3. + 2.*grid*lamb*(-3. + 5.*grid*lamb)))*Log(1./grid)*(-1. + 2.*grid*lamb*(2. - 5.*grid*lamb) + 4.*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid)))/(3.*Power(g,4))
    vF = vF0 + vF2*fxy**2 + vF4*fxy**4

    LM = L(vq2,vq1,vq0,N,x1,x2,dpass = dC, ddpass = d2C)

    correct_BC = ((Bdot[0]**2)*sigma[0])/((fxy**2)/(3.*Power(E,2.*B[0])*(g**2)*Power(sigma[0],3)) - 2.*sigma[0] + (Power(Fvu[0],2)*sigma[0]*Power(grid[0],4))/(3.*(g**2)))



#    correct_BC = -0.031294977000893336

#### A - method 1

    #A_sub_temp = spectral_solve(LM,vF,('D',N,(2)),('na',N,0),N,x1,x2) # any value can go in there, and regularity is already forced by spectral method

    #Now we impose the boundary condition at the horizon.  We can add a constant to fix the value at the boundary, because a constant is a homogenous solution to the equation.

    #A_temp = lamb**2 + grid**(-2) + (2.*lamb)/grid + ((-2.*fxy*2*grid**2)/(3.*(g**2)) + (4.*(fxy**2)*lamb*(grid**3))/(3.*(g**2)))*Log(1./grid) + A_sub_temp

    #print "Correct BC:", correct_BC

    #A = A_temp - A_temp[0] + correct_BC  #adding a constant function, by adding the constant to each coefficient.

    #A_sub = A_sub_temp - A_temp[0] + correct_BC

    #print A[0]


#### Method 2

  # without MORE LOGS  correct_BC_sub = correct_BC - (Power(lamb,2) + Power(grid[0],-2) + (2.*lamb)/grid[0] + ((-2.*(fxy**2)*Power(grid[0],2))/(3.*(g**2)) + (4.*(fxy**2)*lamb*Power(grid[0],3))/(3.*(g**2)))*Log(1./grid[0]))

    correct_BC_sub = correct_BC + (-3.*Power(g,2) - 6.*Power(g,2)*grid[0]*lamb - 3.*Power(g,2)*Power(grid[0],2)*Power(lamb,2) + 2.*Power(fxy,2)*Power(grid[0],4)*Log(1./grid[0]) - 4.*Power(fxy,2)*Power(grid[0],5)*lamb*Log(1./grid[0]) + 6.*Power(fxy,2)*Power(grid[0],6)*Power(lamb,2)*Log(1./grid[0]) - 8.*Power(fxy,2)*Power(grid[0],7)*Power(lamb,3)*Log(1./grid[0]))/(3.*Power(g,2)*Power(grid[0],2))
    #correct_BC_sub = correct_BC - Power(1 + lamb,2) #simpler form
    A_sub = spectral_solve(LM,vF,('D',0,correct_BC_sub),('N',N,0),N,x1,x2,dpass = dC, ddpass = d2C) # This should have the correct value immediatly
    dA_sub = np.dot(dC, A_sub)

#Old Way A_sub = spectral_solve(LM,vF,('na',None,None),('na',None,None),N,x1,x2,dpass = dC, ddpass = d2C,cntrltype = 'D', cntrlBCloc = 0, cntrlrowreplace = N, cntrlval = correct_BC_sub) # This should have the correct value immediatly

    
    A = Power(grid,-2) + (2.*lamb)/grid + Power(lamb,2) + A_sub + ((-2.*Power(fxy,2)*(grid**2))/(3.*Power(g,2)) + (4.*Power(fxy,2)*(grid**3)*lamb)/(3.*Power(g,2)) -  (2.*Power(fxy,2)*(grid**4)*Power(lamb,2))/Power(g,2) + (8.*Power(fxy,2)*(grid**5)*Power(lamb,3))/(3.*Power(g,2)))*Log(1./grid)
    
  

    dvlamb = -A_sub[N]/2. 


    dvB_sub = (6.*Bdot_sub_scale*Power(g,2) + grid*(4.*Power(fxy,2)*(-1. + grid*(3.*lamb + 2.*grid*(dvlamb - 5.*dvlamb*grid*lamb + Power(lamb,2)*(-3. + 5.*grid*lamb))))*Log(1./grid) + ((9.*B_sub*Power(g,2) - grid*(-3.*dB_sub*Power(g,2) + Power(fxy,2)*(1. + 2.*grid*lamb*(-2. + 5.*grid*lamb))) + 4.*Power(fxy,2)*grid*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid))*(3.*Power(g,2)*(A_sub*Power(grid,2) + Power(1. + grid*lamb,2)) + 2.*Power(fxy,2)*Power(grid,4)*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb)))*Log(1./grid)))/(3.*Power(g,2)*grid)))/(6.*Power(g,2)*grid)

#    dvB_sub = (18.*Bdot_sub_scale + (24.*dvlamb*Power(fxy,2)*Power(grid,3)*(1. - 5.*grid*lamb)*Log(1./grid))/Power(g,2) + (12.*Power(fxy,2)*grid*(-1. + grid*lamb*(3. + 2.*grid*lamb*(-3. + 5.*grid*lamb)))*Log(1./grid))/Power(g,2) + ((9.*B_sub*Power(g,2) - grid*(-3.*dB_sub*Power(g,2) + Power(fxy,2)*(1. + 2.*grid*lamb*(-2. + 5.*grid*lamb))) + 4.*Power(fxy,2)*grid*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*Log(1./grid))*(3.*Power(g,2)*(A_sub*Power(grid,2) + Power(1. + grid*lamb,2)) + 2.*Power(fxy,2)*Power(grid,4)*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb)))*Log(1./grid)))/Power(g,4))/(18.*grid)
    
#    dvB_sub0logs = (6.*Bdot_sub_scale*Power(g,2) + (A_sub*Power(grid,2) + Power(1. + grid*lamb,2))*(9.*B_sub*Power(g,2) - grid*(-3.*dB_sub*Power(g,2) + Power(fxy,2)*(1. + 2.*grid*lamb*(-2. + 5.*grid*lamb)))))/(6.*Power(g,2)*grid)
#    dvB_sub1logs = (Power(fxy,2)*Power(grid,2)*(12.*dvlamb*Power(g,2)*(1. - 5.*grid*lamb) + 6.*A_sub*Power(g,2)*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb)) + grid*(9.*B_sub*Power(g,2)*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb))) - Power(fxy,2)*grid*(1. + 2.*grid*lamb*(-2. + 5.*grid*lamb))*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb))) + 3.*Power(g,2)*(10.*Power(lamb,3)*(7. + 3.*grid*lamb) + dB_sub*grid*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb)))))))/(9.*Power(g,4))
#    dvB_sub2logs = (4.*Power(fxy,4)*Power(grid,4)*(1. + 5.*grid*lamb*(-1. + 3.*grid*lamb))*(-1. + grid*lamb*(2. + grid*lamb*(-3. + 4.*grid*lamb))))/(9.*Power(g,4))
#    dvB_sub = dvB_sub0logs + Log(1./grid)*dvB_sub1logs + (Log(1./grid)**2)*dvB_sub2logs

    dvB_sub[N] = 0 #can set this to zero.  Just expand out near boundary to see this.
    #report(dvB_sub, comparedest, 'a', 'PYdvBpaul')
    #report((np.dot(dC,B_sub)), comparedest, 'a', 'PYdBpaul')
    #report([dvlamb], comparedest, 'a', 'PYdvlamb')



    

####  FILTERING/ALIASING PROBLEM

    #define collocation points rougly 2/3 of original. #DONT DO FILTERING HERE ANYMORE
#    if filt == "T": 
#        print "filtering!"
#        dvB_sub = filter_func(dvB_sub,N,Nfilter,x1,x2)        


####  FIXING LAMBDA DRIFT
#    if fixlambdrift == True:
#        print "fixing lamb"
#        lambold = lamb
#        lamb = lamb + (1. - horizon)/(horizon*1.) 

#### REPORT
#    if pr:
#        report(B_sub,destination,'a','PYB_subCoefs' + str(time))
#        report(sigma_sub_scale,destination,'a','PYsigma_sub_scaleCoefs' + str(time))
 #       report(dsigma_sub_scale,destination,'a','dSigma_sub_scale')
 #       report(Fvu,destination,'a','Fvu')
#        report(sigmadot_sub_scale,destination,'a','PYSigmadotsubscaleCoefs' + str(time))
#        report(Bdot_sub, destination,'a','PYBdotsubCoefs' + str(time))
#        report(A_sub,destination,'a','PYA_subCoefs' + str(time))
#        report(dvB_sub,destination,'a','PYdvB_subCoefs' + str(time)) 
#        report([dvlamb],destination,'a','dvlamb' + str(time)) 
#        report([dvlamb],destination,'a','b4' + str(time)) 
 

#### PLOT
    if pl:
        finegrid = np.linspace(.01,1.01,100)
        fig1 = plt.figure(figsize=(20,12))
        plot1 = fig1.add_subplot(231)
        plot1.plot(finegrid,cardfunct(Bdot_sub_scale,finegrid,N,x1,x2),'o')
        plt.xlabel(sigmadotfunct(1)) 
#        plt.axis([.9, 1.01, -.05, 1])
        plt.title('Bdot')

        plot2 = fig1.add_subplot(232)
        plot2.plot(finegrid,cardfunct(Fvu,finegrid,N,x1,x2))
        plt.axis([0, 1.01, None, None])
        plt.title('Fvu')

#        plot3 = fig1.add_subplot(233)
#        plot3.plot(finegrid,cardfunct(sigmadot_sub_scale,finegrid,N,x1,x2))
#        plt.title('sigmadot_sub_scale')

#        plot4 = fig1.add_subplot(234)
#        plot4.plot(finegrid,finegrid**2*cardfunct(Bdot_sub,finegrid,N,x1,x2))
#        plt.title('Bdot_sub')        

#        plot5 = fig1.add_subplot(235)
#        plot5.plot(finegrid,cardfunct(A_sub,finegrid,N,x1,x2))
#        plt.title('A_sub')
        #fig2 = plt.figure(figsize=(20,12))
        #myplot = fig2.add_subplot(231)
        #myplot.plot([1,2,3],[2,4,6],'o')
        #plt.title('My Plot')

#        plt.draw()
#        plt.show()
#        plt.close()
    

    b4check = -1./2.*np.dot(dC,Bdot_sub_scale)[N]
    a2check = 2.*np.dot(dC,sigmadot_sub_scale)[N]
    templist1 =  [lamb, dvlamb, dvB_sub, b4, Bdot_sub_scale[N], b4check, horizon, a2check, A[0]-correct_BC, -1, foundhorizon, sigmadot_sub_scale[N], Bdot_sub_scale, sigmadot_sub_scale, sigma_sub_scale, horizonalt, A_sub, dsigma_sub_scale[N],dsigma_sub_scale, Fvu, -1, dA_sub[0]]
#    print "is this zero?"
#    val = b4 == templist1[3]
#    print('%.30f' % val)
#    print np.dtype(val)
#    print np.dtype(templist1[3])
 
    #return (lamb, dvlamb, dvB_sub, b4, Bdot_sub_scale[N], b4check, horizon, a2check, A[0]-correct_BC, -1, foundhorizon, sigmadot_sub_scale[N], Bdot_sub_scale, sigmadot_sub_scale, sigma_sub_scale, horizonalt, A_sub, dsigma_sub_scale[N],dsigma_sub_scale, Fvu, -1, dA_sub[0])
    return {'lamb':lamb, 'dvlamb': dvlamb, 'lambshift': lamb-lambbefore, 'dvB_sub':dvB_sub,'b4': b4, 'sigmaBCcheckD': sigma_sub_scale[N], 'sigmaBCcheckN': dsigma_sub_scale[N], 'sigmadotBCcheckD': sigmadot_paul[N],'sigmadotBCcheckN': dsigmadot_paul[N]-a2,'BdotBCcheckD': Bdot_sub_scale[N], 'BdotBCcheckN': dBdot_sub_scale[N]-(-2.*b4 + fxy**2/(6*g**2)),'ABCcheckD': A_sub[0] - correct_BC_sub, 'ABCcheckN': np.dot(dC,A_sub)[N], 'Bdot_sub_scale at bndry':Bdot_sub_scale[N], 'b4check':b4check, 'hor': horizon, 'a2check':a2check,'ABC': A[0]-correct_BC,'none': -1,'foundhor': foundhorizon, 'sigmadot_sub_scale at bndry': sigmadot_sub_scale[N], 'Bdot_sub_scale':Bdot_sub_scale, 'sigmadot_paul': sigmadot_paul, 'sigma_sub_scale': sigma_sub_scale, 'horalt':horizonalt, 'A_sub':A_sub, 'dsigma_sub_scale at bndry':dsigma_sub_scale[N], 'dsigma_sub_scale': dsigma_sub_scale, 'Fvu':Fvu, 'sigmadot at hor': res, 'dA_sub at hor': dA_sub[0], 'sigma_sub_scale at hor': sigma_sub_scale[0], 'dsigma_sub_scale at hor': dsigma_sub_scale[0], 'ddsigma_sub_scale at hor': (np.dot(dC,dsigma_sub_scale)[0]), 'B_sub_scale at hor': B_sub[0], 'dB_sub_scale at hor': dB_sub[0], 'ddB_sub_scale at hor': (np.dot(dC,dB_sub)[0])} 

