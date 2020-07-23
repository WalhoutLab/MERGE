import numpy as np
import pandas as pd
import re

from pylab import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


#########################################
def normal(x,mu,sigma,N):
    '''Returns a normal distribution scaled with <N> as the population size.'''
    deltax=x[2]-x[1];
    coef=deltax*N*(1./(sqrt(2.*pi*sigma**2)));
    return coef*exp(-(x-mu)**2/2/sigma**2)

def superimpose_bimodal(x,mu1,sigma1,mu2,sigma2,Np):
    '''Superimposition function for bi-modal curve fit. Returns a bimodal distribution with two subpopulations,
    one defined by <mu1> and <sigma1> and scaled with global Ndetected-<Np> serving 
    as the subpopulation size, and the other defined by 
    <mu2> and <sigma2> and scaled with <Np>. <Np> is the subpopulation size for curve 2.'''    
    return normal(x,mu1,sigma1,Ndetected-Np)+normal(x,mu2,sigma2,Np);
    
def bimodal(datao,expectedStats=[0,1.5,5.,2.5],Nbin=100,figsize=(12,6),xlims=None,showPlot=True):
    ''' 
    Fits the superimposition of two Gaussian curves to a histogram of data in <datao> with the number of bins indicated by <Nbin>.
    
    <expectedStats> is an estimate of the mean (mu) and standard deviation (sigma) of the two curves to be fitted (as [mu1,sigma1,mu2,sigma2]). 
    If <showPlot> is True, the histogram and fit are plotted. The cfunction should be run iteratively to change <expectedStats> in case of misfit.
    <figsize> is the width and height of figure to be plotted in the respective order.
    <xlims> is the range of the values (bins) used in the histogram. If not provided (default), this range is automatically calculated.
    
    Returns a dictionary of fitted parameters and statistics:
        -mu and sigma are as defined above
        -N indicates subpopulation size for the corresponding curves
        -n is number of bins
        -RSS is Residual Sum of Squares
        -Rsquared is R2 of the fit
        
    USAGE
    
        Stats=bimodal(data,expectedStats=[0,1.5,5.,2.5],Nbin=100,figsize=(12,6),xlims=None,showPlot=True);
    '''
    global Ndetected
    #Static Input
    sumcolor='black';
    bgcolor='black';
    poscolor='black';
    bgline=':';
    posline='--';

    #Data processing and initial plotting
    V=np.array(datao); 
    V=V[np.nonzero(V)];
    data=np.log2(V);
    if xlims is None:
        txlim=(-(int(abs(min(data))/5.)+1)*5.,(int(max(data)/5.)+1)*5.);
    else:
        txlim=(xlims[0],xlims[1]);
    figure(figsize=figsize);
    y,x,_=hist(data,Nbin,range=txlim,alpha=.3,label='data');
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    #Fitting
    Ndetected=len(data);
    expected=tuple(expectedStats+[Ndetected/2.]);
    params,_=curve_fit(superimpose_bimodal,x,y,expected);
    yhat=superimpose_bimodal(x,*params);
    resid=y-yhat;
    ss_res=np.sum(resid**2);
    rss=ss_res/len(y);
    ss_tot=np.sum((y-np.mean(y))**2);
    Rsquared=1.-(ss_res/ss_tot);
    Dparams={'mu1':params[0],'sigma1':params[1],'mu2':params[2],'sigma2':params[3],\
             'N1':Ndetected-params[4],'N2':params[4],'Rsquared':Rsquared,'RSS':rss,'n':len(y)};            
    #Plot the fits
    plot(x,superimpose_bimodal(x,*params),color=sumcolor,lw=3,label='superimposed');
    plot(x,normal(x,params[0],params[1],Ndetected-params[-1]),color=bgcolor,lw=3,linestyle=bgline,label='curve 1');
    plot(x,normal(x,*params[2:]),color=poscolor,lw=3,linestyle=posline,label='curve 2');
    xlim(txlim[0],txlim[1]);
    legend();
    if showPlot:
        show();
    else:
        close();
    return Dparams
    
    
def superimpose_trimodal(x,mu1,sigma1,mu2,sigma2,mu3,sigma3,Nn,Np):
    '''Superimposition function modified for tri-modal curve fit. See superimpose_bimodal(). <Nn> is the subpopulation size for curve 1.'''    
    return normal(x,mu1,sigma1,Nn)+normal(x,mu2,sigma2,Ndetected-(Np+Nn))+normal(x,mu3,sigma3,Np);

def trimodal(datao,expectedStats=[-5,2,1,2,5.,1],Nbin=100,figsize=(12,6),xlims=None,showPlot=True):
    ''' 
    Fits the superimposition of three Gaussian curves to a histogram of data in <datao> with the number of bins indicated by <Nbin>.
    
    <expectedStats> is an estimate of the mean (mu) and standard deviation (sigma) of the three curves to be fitted (as [mu1,sigma1,mu2,sigma2,mu3,sigma3]). 
    If <showPlot> is True, the histogram and fit are plotted. The cfunction should be run iteratively to change <expectedStats> in case of misfit.
    <figsize> is the width and height of figure to be plotted in the respective order.
    <xlims> is the range of the values (bins) used in the histogram. If not provided (default), this range is automatically calculated.
    
    Returns a dictionary of fitted parameters and statistics:
        -mu and sigma are as defined above
        -N indicates subpopulation size for the corresponding curves
        -n is number of bins
        -RSS is Residual Sum of Squares
        -Rsquared is R2 of the fit
        
    USAGE
    
        Stats=trimodal(data,expectedStats=[-5,2,1,2,5.,1],Nbin=100,figsize=(12,6),xlims=None,showPlot=True);
    '''
    global Ndetected
    #Static Input
    sumcolor='black';
    bgcolor='red';
    poscolor='green';
    midcolor='grey';
    bgline=':';
    posline='-';
    midline='--';
    
    #Data processing and initial plotting
    V=np.array(datao); 
    V=V[np.nonzero(V)];
    data=np.log2(V);
    if xlims is None:
        txlim=(-(int(abs(min(data))/5.)+1)*5.,(int(max(data)/5.)+1)*5.);
    else:
        txlim=(xlims[0],xlims[1]);
    figure(figsize=figsize);
    y,x,_=hist(data,Nbin,range=txlim,alpha=.3,label='data');
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    #Fitting
    Ndetected=len(data);
    expected=tuple(expectedStats+[Ndetected/4.,Ndetected/3.]);
    params,_=curve_fit(superimpose_trimodal,x,y,expected);
    print params
    yhat=superimpose_trimodal(x,*params);
    resid=y-yhat;
    ss_res=np.sum(resid**2);
    rss=ss_res/len(y);
    ss_tot=np.sum((y-np.mean(y))**2);
    Rsquared=1.-(ss_res/ss_tot);
    Dparams={'mu1':params[0],'sigma1':params[1],'mu2':params[2],'sigma2':params[3],\
    'mu3':params[4],'sigma3':params[5],'N1':params[6],'N3':params[7],\
    'N2':Ndetected-(params[6]+params[7]),'Rsquared':Rsquared,'RSS':rss,'n':len(y)};            
    #Plot the fits
    plot(x,superimpose_trimodal(x,*params),color=sumcolor,lw=3,label='superimposed');
    plot(x,normal(x,params[0],params[1],params[6]),color=bgcolor,lw=3,linestyle=bgline,label='curve 1');
    plot(x,normal(x,params[2],params[3],Ndetected-(params[6]+params[7])),color=midcolor,lw=3,linestyle=midline,label='curve 2');
    plot(x,normal(x,params[4],params[5],params[7]),color=poscolor,lw=3,linestyle=posline,label='curve 3');
    xlim(txlim[0],txlim[1]);
    legend();
    if showPlot:
        show();
    else:
        close();
    return Dparams
    
    
def subTable(T,geneList,map=None):
    '''
    Extracts a subset of rows, given in <geneList>, from table <T>. 
    
    If the names in geneList have a different nomenclature than the rows of the table, than a mapper (<map>) needs to be provided to go from <geneList> to rows of <T>.

    USAGE
    
        Tsub=subTable(T,geneList,map=None)
        
    '''
    LL=[[],[]];
    if map:
        for i in T.index:
            if map.has_key(i):
                if map[i] in geneList:
                    LL[0].append(i);
                    LL[1].append(map[i]);
    else:
        LL=[geneList,geneList];

    Tsub=pd.DataFrame(index=LL[1],columns=T.columns);
    for i in range(len(Tsub)):
        for col in Tsub.columns:
            Tsub[col][LL[1][i]]=T[col][LL[0][i]];
    return Tsub
        
 
def categorize_absCutoff(T,cutoffs,excludeCols=[]):
    '''
    Categorizes genes (rows) of expression table <T> according to <cutoffs>, list-like object with [rare cutoff,low cutoff,high cutoff].
    
    Genes with expression levels less than rare cutoff are categorized as Rare.
    Other genes with expression levels less than low cutoff are categorized as Low.
    Genes with expression levels greater than high cutoff are categorized as High.
    All other genes are categorized as Moderate.
    
    <excludeCols> list indicates which columns of <T> should not be categorized.
    
    USAGE
    
        Tcat=categorize_absCutoff(T,cutoffs,excludeCols=[])
    '''
    Lcols=[];
    for col in T.columns:
        if not col in excludeCols:
            Lcols.append(col);
    Tcat=pd.DataFrame(index=T.index,columns=Lcols);
    for i in Tcat.index:
        for col in Tcat.columns:
            val=T[col][i];
            if val<=0:
                Tcat[col][i]='Rare';
            elif val<cutoffs[0]:
                Tcat[col][i]='Rare';
            elif val<cutoffs[1]:
                Tcat[col][i]='Low';
            elif val>cutoffs[2]:
                 Tcat[col][i]='High';
            else:
                Tcat[col][i]='Moderate';       
    return Tcat
    
def stackedCat(Tcat):
    '''Plots a stacked bar graph of categories (Rare, Low, Moderate, High) for each column of <Tcat>.'''
    h,m,l,r=[],[],[],[];
    N=len(Tcat.columns);
    ind = np.arange(N); 
    width = 0.8;
    
    for col in Tcat.columns:
        s=pd.value_counts(Tcat[col]);
        h.append(s['High']);
        m.append(s['Moderate']);
        l.append(s['Low']);
        r.append(s['Rare']);

    pr = plt.bar(ind, r, color='red',width=width);
    pl = plt.bar(ind, l, color='orange',bottom=r,width=width);
    pm = plt.bar(ind, m, color='grey',bottom=[r[i]+l[i] for i in range(N)],width=width);
    ph = plt.bar(ind, h, color='green',bottom=[r[i]+l[i]+m[i] for i in range(N)],width=width);

    plt.ylabel('Number of genes');
    plt.xticks([i+0.5 for i in ind], Tcat.columns);
    plt.xticks(rotation=90, fontsize=8);
    plt.legend((pr[0],pl[0],pm[0],ph[0]), ('Rare', 'Low','Moderate','High'));

    plt.show()   

def relativeExp(Texpo,Tcato,tao_rel,fc_mid=1.5,fc_end=4.,fcHigh=None,fcLow=None):
    '''
    Recategorizes some Moderate genes in <Tcato>, a table of categorized genes (rows) in tissues or conditions (columns) as Low or High based on relative expression levels.
    
    Expression levels are provided as <Texpo>, a table that matches <Tcato> in rows and columns, but has expression levels instead of categories.
    
    A heuristic method is used such that the expression profile of each gene (row) is first obtained by sorting the row in <Texpo> from low to high expression.
    Then fold changes from one tissue or condition to the next is monitored and significant jumps are tracked.
    For an increase to be considered significant, the fold change (FC) should be greater than a threshold and the larger number should be greater than <tao_rel>.
    If the increase is  greater than an FC threshold (FC for low), then the lower value and all below are labeled Low.
    If the increase is  greater than an FC threshold (FC for high), the higher value and all above are labeled High.
    Variable (High and Low) labeling of the same value in different steps result in no categorization for that value.

    FC for low and FC for high values for every increment (for N columns, there are N-1 increments) can be provided as lists (<fcLow> and <fcHigh>).
    If low vs high thresholds are to be the same and the same threshold is going to be used for all middle increments and another threshold for terminal increments, then <fc_mid> and <fc_end> thresholds can be provided instead, as single thresholds for the respective increment sets.  
    If the same threshold is to be used throughout, then just enter this value for both <fc_mid> and <fc_low>.
    If <fcHigh> or <fcLow> list is provided, <fc_mid> and <fc_end> will be automatically null for the thresholds determined by the list.
    
    In the end, a Moderate gene in a tissue can be labeled as Low only if its expression level is less than <tao_rel> and labeled as High only if its expression level is higher than <tao_rel>. 
    
    USAGE
    
        Tcat_final=relativeExp(Texp,Tcat,tao_rel,fc_mid=1.5,fc_end=4.,fcHigh=None,fcLow=None)
 
    '''
    Tcatf=Tcato.copy();
    Texp=Texpo.drop('ave',axis=1);

    n=len(Texp.columns);
    if not fcHigh:
        fcHigh=[fc_end]+(n-2)*[fc_mid];
    elif not len(fcHigh)==n-1:
        raise Exception('fcHigh vector must of size n-1, where n is the number of conditions');
    if not fcLow:
        fcLow=(n-2)*[fc_mid]+[fc_end];
    elif not len(fcLow)==n-1:
        raise Exception('fcLow vector must of size n-1, where n is the number of conditions');
    for gene in Texp.index:
        row=Texp.loc[gene,Texp.columns].copy(); ###
        row=row.sort_values(); ###
        S=[set() for i in range(n)];
        fc=[0];
        for k in range(1,n):
            if row[k-1]:
                fc.append(row[k]/row[k-1]);
            else:
                fc.append(inf);       
        for k in range(1,n):
            if row[k]>tao_rel:
                #check relatively high
                if fc[k]>fcHigh[k-1]:
                    for i in range(0,k):
                        S[i].add('not high');
                    for i in range(k,n):
                        S[i].add('high');   
                #check relatively low
                if fc[k]>fcLow[k-1]:
                    for i in range(0,k):
                        S[i].add('low');
                    for i in range(k,n):
                        S[i].add('not low');                   
        for i in range(n):
            cond=row.index[i];
            if Tcatf[cond][gene]=='Moderate':
                if row[i]>tao_rel:
                    if S[i].issubset(['high','not low']) and 'high' in S[i]:
                        Tcatf [cond][gene]='High';
                elif S[i].issubset(['low','not high']) and 'low' in S[i]:
                    Tcatf [cond][gene]='Low';    
    return Tcatf

def plotCatExp(gene,Texp,Tcato,Tcatf):
    '''
    Plots a bar chart that shows ascending tissue (or condition) profile of <gene> based on expression level table <Texp> with genes as rows and tissues (coditions) as columns.
    Categorization tables from categorize_absCutoff() function (<Tcato>) and relativeExp() function (<Tcatf>) are to be provided.
    
    Bars are colored according to category ({'High':'green','Moderate':'gray','Low':'orange','Rare':'red'}).
    Categories that differ between <Tcato> and <Tcatf> are indicated by hatches.
    
    USAGE
        
        plotCatExp(gene,Texp,Tcato,Tcatf)
    '''
    width=0.5;
    dcolor={'High':'green','Moderate':'gray','Low':'orange','Rare':'red'};

    x=Texp.loc[gene,Texp.columns];
    x=x.sort_values();  
    ind=np.arange(len(x));      
    rects = plt.bar(ind, x.values, width=width);
    plt.xticks([i+width/2 for i in ind],x.index);
    plt.xticks(rotation=90, fontsize=8);
    for i in ind:
        cond=x.index[i];
        rects[i].set_facecolor(dcolor[Tcatf[cond][gene]]);
        if not Tcatf[cond][gene]==Tcato[cond][gene]:
            rects[i].set_hatch('/');
    plt.show(); 

  
 
