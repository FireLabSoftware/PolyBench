from math import log
def quadratic(aq,bq,cq):
    ''' solutions to a quadratic equation with a*x^2+b*x+c==0; aq, bq, and cq are a,b, and c
    returns (s1,s2) where s1 and s2 are the different solutions; returns "NoRealSolution" if there is no
    real solution'''
    dia1=bq**2-4*aq*cq
    if dia1<0:
        return ("No_real_solution",(aq,bq,cq))
    if aq==0:
        if bq!=0:
            return (-cq/float(bq),-cq/float(bq))
        else:
            return ("No_real_solution",(aq,bq,cq))
    s1=(-bq+(dia1)**0.5)/(2*aq)
    s2=(-bq-(dia1)**0.5)/(2*aq)
    return (s1,s2)
'''find a solution s for f(s)=goal so that the margin of error abs(log(f(s)/goal,2)<logerror and abs(f(s)-goal)<error
    f(unction), g(oal),s(tartInput,e(ndInput),t(olerance),l(ogTolerance),m(inInput),h(ighestInput),c(cycles)         The maximum number of test cycles
    '''
def solve(**d):
    '''f(unction), g(oal),s(tartInput),e(ndInput),t(olerance),l(og10Tolerance),m(inInput),h(ighestInput),c(cycles)
        Find a solution s for f(s)=goal so that the margin of error abs(log(f(s)/goal,2)<logerror and abs(f(s)-goal)<error
        This will only work well if the function is monotonic for the range of inputs between s and e and if there is indeed a solution'''
    ## set defauls for a few parameters
    c=100.0
    g=0.0
    
    
    for n in d:
        n1=n[0].lower()
        d[n1]=d[n]
        if n1=='f':f=d[n]  ## the function to solve
        if n1=='g':g=float(d[n])  ## the intended value of the function
        if n1=='s':s=d[n]  ## first (lower) value to test
        if n1=='e':e=d[n]  ## first (upper) value to test
        if n1=='t':t=d[n]  ## tolerance... the maximum abs difference between f(x) and the intended value g
        if n1=='l':l=d[n]  ## log tolerance ... the maximum abs log10 value of f(x)/g
        if n1=='m':m=d[n]  ## the minimum value that can be fed into the function f
        if n1=='h':h=d[n]  ## the maximum value that can be fed into the function f
        if n1=='c':c=d[n]  ## the maximum number of cycles
    ## set default guesses if the user doesn't provide s and e
    if not('s' in d):
        if ('m' in d):
            d['s']=d['m']
        elif 'e' in d:
            d['s']=d['e']-1.0
        else:
            d['s']=0.0
        s=d['s']
    if not('e' in d):
        if ('h' in d):
            d['e']=d['h']
        else:
            d['e']=d['s']+1.0
        e=d['e']
    ## set default guesses if the user doesn't provide t or l
    if not('t' in d) and not('l' in d):
        d['t']=0.01
        t=d['t']
        if g>0.0:
            d['l']=0.01
            l=d['l']
        
    if g<=0.0 and 'l' in d:
        print("Trouble-- you are seeking a zero value with a log tolerance")
    a=False  ## are we at the goal?
    i1=s     ## current lower bound test value
    i2=e     ## current upper bound test value
    r1=f(i1) ## current lower bound result
    r2=f(i2) ## current upper bound result
    oldr='Aardvark' ## the most recent test value, to test for convergence
    cycle=0
    while not(a):
        if r2<=r1:
            ## first sort the two test values i1, i2 so r2 is the one with the higher function value
            ## if the values are equal, switch the order (since at that point we are at a flat part
            ## of the curve and need to be checking both sides to see where things get better
            r1,r2,i1,i2=r2,r1,i2,i1
        if r1<g and r2>g:
            i3=(i1+i2)/2.0
            r3=f(i3)
            if r3<g:
                r1=r3
                i1=i3
            else:
                r2=r3
                i2=i3
        elif r2<g:   ## r1<r2<g
            i3=i2+2*(i2-i1)
            if 'm' in d: i3=max(m,i3)
            if 'h' in d: i3=min(h,i3)
            r3=f(i3)
            r1,r2,i1,i2=r2,r3,i2,i3
        else:        ## g<r1<r2
            i3=i1-2*(i2-i1)
            if 'm' in d: i3=max(m,i3)
            if 'h' in d: i3=min(h,i3)
            r3=f(i3)
            r1,r2,i1,i2=r3,r1,i3,i1
        a=True
        if ('l' in d) and ( g==0.0 or r3/g<=0 or abs(log(r3/g,10))>l ):
            a=False
        if ('t' in d) and (abs(r3-g)>t):
            a=False
        cycle+=1
        if cycle>c:
            print('Cycle count hit '+str(c)+' , giving current value')
            a=True
        oldr=r3
    return i3
class BRspecies():
    '''Information for incidence of a single species in two different datasets'''
    def __init__(self,P1,P2,T1,T2,ts1=2,ts2=2,r=1):
        '''initialize the object;
            P1 and P2 are the incidences of the the species in datasets 1 and 2
            T1 and T1 are total incidences in 1 and 2
            ts1 and ts2 are estimates of total species in 1 and 2
            r1 is a regularization constant (involves bumping up each sample to avoid division by zero'''
        if (P1<0) or (P2<0):
            print(str(("Warning, data point with success count(s)<0. Setting success count to zero.  P1,P2,T1,T2=",P1,P2,T1,T2)))
            P1=max(0,P1)
            P2=max(0,P2)
        if (P1>T1) or (P2>T2):
            print(str(("Warning, data point with more successes than trials. Will set success count to total trials.  P1,P2,T1,T2=",P1,P2,T1,T2)))
            P1=min(P1,T1)
            P2=min(P2,T2)
        self.P1=P1
        self.P2=P2
        self.T1=T1
        self.T2=T2
        self.P1c=P1+r ## these are regularizations that avoid any problem when dividing by zero.  
        self.P2c=P2+r ## with reasonably large datasets, results should be reasonably independent of regularization values
        self.T1c=T1+ts1*r  ##r==0 avoids any regularization
        self.T2c=T2+ts2*r
        self.milp1=float(self.P1c)/self.T1c
        self.milp2=float(self.P2c)/self.T2c
    def logPprob(self,p1,p2):
        if p1<=0.0 or p1>=1.0 or p2<=0.0 or p2>=1.0:
            return -1e99
        '''gives the log2 posterior probability of the observed results for individual dataset probabilities p1 and p2'''
        pp1=self.P1*log(p1,2)+(self.T1c-self.P1c)*log(1-p1,2)
        pp2=self.P2*log(p2,2)+(self.T2c-self.P2c)*log(1-p2,2)
        return pp1+pp2
    def logPi(self):
        '''provides a maximum likelihood log p value for the event with p1 and p2 being essentially the observed incidences'''
        return self.logPprob(self.milp1,self.milp2)
    def ratap1p2(self,a):
        '''Mnemonic: Ratio a, p1 p2.
        Provides a maximum likelihood value for p1 with the assumption that p2=a*p1
        Output is a 3-tuple:
            p1,p2 (a*p1) and log-probability
        Here's a brief summary of the math
        For observations P1,P2,T1,T2 (un-regularized) and probabilities p1,p2=p1,ap1
        Pprob=p1^P1 * (1-p1)^(T1-P1) * (a*p1)^P2 * (1-a*p1)^(T2-P2)
        To find a maximum likelihood p1 under these conditions, we set dPprob/dp1==0
        Doing the term-by-term differential gives
        dPprob/dp1 = Pprob * ( P1/p1 - (T1-P1)/(1-p1) + P2/p1 - a*(T2-P2)/(1-a*p1)
        Expressing this as an equation in p1, we get
        (P1+P2)*(1-p1)*(1-a*p1) = (T1-P1)*p1*(1-a*p1)+a*(T2-P2)*p1*(1-p1)
        a*(T1+T2)*p1^2 - (P1+a*P1+P2+a*P2+(T1-P1)+a*(T2-P2))*p1 + (P1+P2)
        a*(T1+T2)*p1^2 - (T1+a*P1+P2+a*T2)*p1 + (P1+P2)
        The zero points of this can be solved as a quadratic formula'''
        aq1=a*(self.T1c+self.T2c)
        bq1=-(self.T1c+a*self.P1c+self.P2c+a*self.T2c)
        cq1=self.P1c+self.P2c
        p1L,p1H=quadratic(aq1,bq1,cq1)
        p2L=p1L*a
        p2H=p1H*a
        logPaL=self.logPprob(p1L,p2L)
        logPaH=self.logPprob(p1H,p2H)
        if logPaL>logPaH:
            return (p1L,p2L,logPaL)
        else:
            return (p1H,p2H,logPaH)
    
    def logPa(self,a):
        '''this gives log of conditional probability given
            observations P1/T1, P2/T2 and
            maximum likelihood values obtained from ratap1p2 with ratio a'''
        return self.ratap1p2(a)[2]
    def logDifa(self,a):
        '''this gives a log probability ratio between
            i) a model in which maximum liklihood values are chosen for p1 and p2
            ii) a model in which p2 is constrained to be equal to a*p1
        The result is in log2 scale and is most positive if the independent model is much better
        The independent model will always be at least as good, so A.logDifa() is always at minimum 0.
        The lowest value for logDifa is when p2=a*p1, in which case the models are essentially identical and logDifa is zero'''
        logPa1=self.logPa(a)
        logPi1=self.logPi()
        return logPa1-logPi1
    def siga(self,a0):
        '''A.siga(a) Gives a significance value for any given set of observations (P1/T1 and P2/T2 for two samples) being
        at least a-fold different.  Essentially this is the log base 2 of the ratio between the conditional probability of the
        observations given a model in which optimal incidences are chosen for the two samples with no linkage and a model in which
        the instances are constrained to differ by a factor of a in either direction.
        a should be a value of >1 for this check, what is returned is a >=0 number which is a strongly positive value if the experimental
        result is difficult to account for without an anomoly of >a-fold in the underlying instances in the two samples
        '''
        a1=float(a0)
        if a0<=1.0:
            return 0.0 
        if self.milp2>a1*self.milp1:
            return self.logDifa(a1)
        elif self.milp2<self.milp1/a1:
            return self.logDifa(1/a1)
        else:
            return 0.0
    def maxd(self,u0,minLogDif=0.0001):
        '''A.maxd(u) gives the least significant (closest to zero) log2-fold \difference
            between incidence in samples P1 and P2 for which the
            simple probability of having a lesser fold difference is less than u
            minLogDif is the granularity of the reported analysis, this is the smallest log-2-fold difference that one might
            consider significant'''
        u1=log(u0,2)
        at0=0.0
        at1=log(self.milp2/self.milp1,2)
        
        if at1>0:
            atdel=at1-at0
            while atdel>minLogDif:
                atmed=(at1+at0)/2.0
                if self.logDifa(2.0**atmed)<u1:
                    at0=atmed
                else:
                    at1=atmed
                atdel=abs(at1-at0)
            return at0
        elif at1<0:
            at2=-at1
            atdel=at2-at0
            while atdel>minLogDif:
                atmed=(at2+at0)/2.0
                if self.logDifa(0.5**atmed)<u1:
                    at0=atmed
                else:
                    at2=atmed
                atdel=abs(at2-at0)
            return at0
        else:
            return 0.0
            
           
    def valNewP2(self,a,P2new):
        '''starting with the given set of values, gives the log2Difa with a different value of P2'''
        test=BRspecies(self.P1,P2new,self.T1,self.T2)
        return test.logDifa(a)
    def maxP2givena(self,a,u):
        '''starting with a [ratio between P2 and P1] and u (minimum probability), provide a maximum value of P2'''
        u1=log(u,2)
        x=solve(f=lambda x1:self.valNewP2(a,x1),s=a*self.T2*self.P1/(1.0*self.T2),m=0.01,e=10+a*self.T2*self.P1/(1.0*self.T2),t=0.0001,g=u1)
        return x
    
    def minP2givena(self,a,u):
        '''starting with a [ratio between P2 and P1] and u (minimum probability), provide a maximum value of P2'''
        u1=log(u,2)
        x=solve(f=lambda x1:self.valNewP2(a,x1),m=0.01,e=a*self.T2*self.P1/(1.0*self.T2),s=-10+a*self.T2*self.P1/(1.0*self.T2),t=0.0001,g=u1)
        return x
                
def ConservativeEnrichment(PExpt,PRef,TExpt,TRef,FDR=0.05,Sp=2):
    '''This function returns the most conservative enrichment of a species P in an experimental relative to a reference dataset,
        i.e., the least enrichment that is compatible with the parameters and data provided.
        (PExpt,TExpt) and (PRef,TRef) are the instances of a specific species [P] and total instances [T] in experiment and reference respectively
        FDR is the desired False Discovery rate (species that might be mis-identified per experiment)
        Sp is a rough estimate of the total number of observed species-- Sp is for regularization and need not be precise
        The value returned is the most conservative value of the enrichment in the experimental sample, or 1 if there is no
        evidence for such an enrichment under the statistical conditions of the query'''
    D1=BRspecies(PExpt,PRef,TExpt,TRef,Sp,Sp,0.5)
    if D1.milp2/D1.milp1<=1.0:
        return 1.0
    s1=D1.maxd(FDR)
    if s1>0:
        return 0.5**s1
    else:
        return 1.0
def ConservativeDeEnrichment(PExpt,PRef,TExpt,TRef,FDR=0.05,Sp=2):
    '''This function returns the most conservative de-enrichment of a species P in an experimental relative to a reference dataset,
        i.e., the least de-enrichment that is compatible with the parameters and data provided.
        (PExpt,TExpt) and (PRef,TRef) are the instances of a specific species [P] and total instances [T] in experiment and reference respectively
        FDR is the desired False Discovery rate (species that might be mis-identified per experiment)
        Sp is a rough estimate of the total number of observed species-- Sp is for regularization and need not be precise
        The value returned is the most conservative value of the deenrichment in the experimental sample, or 1 if there is no
        evidence for such a de-enrichment under the statistical conditions of the query'''
    D1=BRspecies(PRef,PExpt,TRef,TExpt,Sp,Sp,0.5)
    if D1.milp1/D1.milp2>=1.0:
        return 1.0
    s1=D1.maxd(FDR)
    if s1>0:
        return 2.0**s1
    else:
        return 1.0
def ConservativeFoldDifference(PExpt,PRef,TExpt,TRef,FDR=0.05,Sp=2):
    '''This function returns the most conservative differential evaluation of a species P in an experimental relative to a reference dataset,
        i.e., the least difference in incidence that is compatible with the parameters and data provided.
        (PExpt,TExpt) and (PRef,TRef) are the instances of a specific species [P] and total instances [T] in experiment and reference respectively
        FDR is the desired False Discovery rate (species that might be mis-identified per experiment)
        Sp is a rough estimate of the total number of observed species-- Sp is for regularization and need not be precise
        The value returned is the most conservative value of the difference in the experimental sample, or 1 if there is no
        evidence for such a de-enrichment under the statistical conditions of the query'''
    v1=ConservativeEnrichment(PExpt,PRef,TExpt,TRef,FDR=FDR,Sp=Sp)
    v2=ConservativeDeEnrichment(PExpt,PRef,TExpt,TRef,FDR=FDR,Sp=Sp)
    if v1!=1.0:
        return v1
    elif v2!=1.0:
        return v2
    else:
        return 1.0
    
def maxD(PExpt,PRef,TExpt,TRef,FDR=0.05,Sp=2):
    '''This function returns the most conservative differential evaluation of a species P in an experimental relative to a reference dataset,
        i.e., the least difference in incidence that is compatible with the parameters and data provided.
        (PExpt,TExpt) and (PRef,TRef) are the instances of a specific species [P] and total instances [T] in experiment and reference respectively
        FDR is the desired False Discovery rate (species that might be mis-identified per experiment)
        Sp is a rough estimate of the total number of observed species-- Sp is for regularization and need not be precise
        The value returned is the most conservative value of the difference in the experimental sample, or 1 if there is no
        evidence for such a de-enrichment under the statistical conditions of the query'''
    D1=BRspecies(PExpt,PRef,TExpt,TRef,Sp,Sp,0.5)
    return D1.maxd(FDR)
    
    
    
    
         
            
        
                
                
        
        
    
    
