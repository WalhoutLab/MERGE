###############################################################################
#########           EFFICIENCY         ########################################
###############################################################################


def getEfficiencyModel_dual(TmodelFile='Tmultiple_2',cpd=None,consume=False,\
                       stripCoef=('DGR0007_L',6.0),overwrite_sideR=None,overwrite_sideR_ind=None,
                       overwrite_storageR=None,addReactants=[],addProducts=[],bacIntake=1.):                      
    '''Creates, constrains, and returns a model for subsequent use in next generation
    flux or production potential predictions. If metabolite <cpd> is provided, adds a
    sink reaction (SINK) that drains (<consume>=False) or introduces (<consume>=True)
    this metabolite. Model table and basic constraints defining nutritional
    state are in .pkl files starting with <TmodelFile> and <DconFile>, respectively.
    
    <stripCoef> tuple is added to cut the dependence on energy because of
    bacterial digestion. This tuple has the name of the atp-dependent reaction that
    degrades bacterial whole mass into macromolecules as its first element and the
    coefficient of atp therein as the second. This coefficient is replaced by
    0.0 in the resulting model. If this tuple is None, then no action will be
    taken to modify the digestion energy. 
    
    If <cpd> exists and production or consumption of the compound requires other metabolites
    then these should be defined as a list in <addReactants> and/or <addProducts>.
    
    IMPORTANT:Model is converted to a single direction model. 
    
    Note that, because M (model) is made single dir, gurobi object under M sees the
    SINK reaction as SINKf (<consume>=False).   
    
        Created     :   10/2018
        Modified    :   10/2019
    '''
    #STATIC
    Db={'DMN0033': [0.0, 1.0],\
     'EX00001': [-1000.0, 1000.0],\
     'EX00007': [-1000.0, 1000.0],\
     'EX00009': [-1000.0, 1000.0],\
     'EX00080': [-1000.0, 1000.0],\
     'EXC0050': [-bacIntake, 0.0],\
     'RM00112': [-1000.0, 1000.0],\
     'RMC0005': [0.0, 1000.0]};
     
    #Preps
    if overwrite_sideR is None:
        sideRf=sideR;
    else:
        sideRf=overwrite_sideR;
    if overwrite_sideR_ind is None:
        sideR_indf=sideR_ind;
    else:
        sideR_indf=overwrite_sideR_ind;
    if overwrite_storageR is None:
        storageRf=storageR;
    else:
        storageRf=overwrite_storageR;      

    addReactants_=[x for x in addReactants];
    addProducts_=[x for x in addProducts];

    Tmodel=op(TmodelFile);
    if stripCoef:
        for col in ['BIGG','KEGG']:
            Tmodel[stripCoef[0],col]=re.sub('\('+str(stripCoef[1])+'\)','(0.0)',Tmodel[stripCoef[0],col]);
    M=myFlux.mnmGurobi(Tmodel)
    ##Modify model with Sink
    if not cpd is None:            
        if consume:
            addProducts_.append(cpd);
            rxn=' + '.join(addReactants_)+' --> '+' + '.join(addProducts_);
        else:
            addReactants_.append(cpd);
            rxn=' + '.join(addReactants_)+' --> '+' + '.join(addProducts_);
        if rxn[0]==' ':
            rxn=rxn[1:];
        if rxn[-1]==' ':
            rxn=rxn[:-1];
            
        ID='SINK';
        M.maddRxn({'ID':ID,'BIGG':rxn});            
        M.mupdate();
    ##Add given constraints
    Dcon={};
    for rid in Db:
        if M.Drid2row.has_key(rid+'_L'):
            Dcon[rid+'_L']=Db[rid];
        elif M.Drid2row.has_key(rid+'_E'):
            Dcon[rid+'_E']=Db[rid];
        elif M.Drid2row.has_key(rid+'_I') and M.Drid2row.has_key(rid+'_X'):
            Dcon[rid+'_I']=Db[rid];
            Dcon[rid+'_X']=Db[rid];
        else:
            raise Exception('The following reaction belongs to a category that is not addressed by this code: '+rid);
    
    M.setConstraints(Dcon);
    M.setgurobi(isSingleDir=True);
    ##Constraints on side metabolites
    Dside={};
    M.faddD2Dcon(Dside,'sideTotal',sign='<=');
    Dside['sideTotal']['RHSvar']=['EXC0050_Lr'];
    Dside['sideTotal']['RHScoef']=[sideRf*bacMW/dividedBy];
    Dside['sideTotal']['LHSvar']=['EXC9998_Ef'];
    Dside['sideTotal']['LHScoef']=[1.];
    for rid in M.Drid2row:
        mw=0.
        for i in range(len(M.Drbigg[rid].reactants[0])):
            if M.Drbigg[rid].reactants[0][i]=='sideMet':
                if not rid=='EXC9998_E':
                    mw=dividedBy*M.Drbigg[rid].reactants[1][i];
        if mw:
            cname='side_'+rid;
            M.faddD2Dcon(Dside,cname,sign='<=');
            Dside[cname]['RHSvar']=['EXC0050_Lr'];
            Dside[cname]['RHScoef']=[sideR_indf*bacMW];
            Dside[cname]['LHSvar']=[rid+'r'];
            Dside[cname]['LHScoef']=[mw];
           
    M.g.maddAnyCon(Dside);
    ##Constraints on storage metabolites
    Dstorage={};
    M.faddD2Dcon(Dstorage,'storageTotal',sign='<=');
    Dstorage['storageTotal']['RHSvar']=['EXC0050_Lr'];
    Dstorage['storageTotal']['RHScoef']=[storageRf*bacMW/dividedBy];
    Dstorage['storageTotal']['LHSvar']=['EXC9999_Ef'];
    Dstorage['storageTotal']['LHScoef']=[1.];
    M.g.maddAnyCon(Dstorage);   
    return M


class fluxEfficiency_dual():
    para={'norder':1.,'dorder':1.5,'allowedFlux':1.0,'nongenePen':1.0,'Stat':'max',\
    'minExp':1.,'maxPen':None,'maxDistance':25};
    dtoggle={'f':'r','r':'f'};
    provider='Intestine';
    def __init__(self,M,Texp,**kwargs):
        '''Class for calculation of flux potentials and efficiencies. 
        
        INPUT: 
        
        <M> is a model from getEfficiencyModel_dual() 
        
        <Texp> is a table of genes (rows) and conditions (tissues) (columns) in 
        a non-dual format. 
        
        Arguments in <kwargs> can fill in self.para dictionary.
        
        '''
        self.conds=Texp.cols;
        if self.conds.count('supercond'):
            raise Exception('You cannot use "supercond" as a name in your experiments. This word is spared for "super condition".')
        self.Texp=Texp;
        self.M=M;
        
        self.non_provider=list(set(self.conds).difference(self.provider));
        self.Dcon_ori=M.getAllConstraints(False).copy();
        
        self.setParams(**kwargs);
        
        self.Archive={};
        self.Dbugs={};
        self.DTflux={};
        self.DDpotential={};
        self.DDefficiency={};
        self.DDconstraints=self._makeEmptyNestedDict(True);
        self.Dlogtables={};
        

    def setParams(self,**args):
       '''Changes parameters in the self.para dictionary.'''
       for key, value in args.items():
            if self.para.has_key(key):
                self.para[key]=value;
            else:
                print 'Argument '+key+' is not understood. Available args are:\n'+str(self.para.keys());
                
    def setPenalties(self,Dspc_pen={},DDspc_pen={},DDspc_minPen={}):
        '''Sets penalties by first determining default penalties from
        convertExp2Penalty() function, the insering special penalties in the
        dictionaries in the input:
        
            <Dspc_pen>: rxn -> penalty for all conditions.
            <DDspc_pen>: cond -> rxn -> penalty (nested dictionary).
            <DDspc_minPen>: cond -> rxn -> penalty (nested dictionary). Penalty
                            will be minimum of the value from this dictionary and
                            that in the current table.
        
        Finally, penalty table is converted to a unidirectional-reaction table.
        '''
        ##get default penalties
        Tpen=self.convertExp2Penalty();
        self.Archive['Tpeno']=Tpen.mCopy();
        if self.para['maxPen']:
            for row in Tpen.rows:
                for col in Tpen.cols:
                    xo=Tpen[row,col];
                    if not np.isnan(xo):
                        Tpen[row,col]=min(self.para['maxPen'],xo);

        ##enter super condition
        Tpen.minsertCol(Tpen.rowno*[1.],Tpen.cols[-1],colName='supercond');
        ##non-gene reactions
        Snongene=self.getNongeneRxns(self.M);
        for rid in Snongene:
            for col in Tpen.cols:
                xo=Tpen[rid,col];
                if not np.isnan(xo):
                    Tpen[rid,col]=self.para['nongenePen'];
                 
        ##enter special penalties
        #special reactions for all
        for rid in Dspc_pen:
             for col in Tpen.cols:
                Tpen[rid,col]=Dspc_pen[rid];           
        #special reactions for special conditions
        for cond in DDspc_pen:
            for rid in DDspc_pen[cond]:
                 for col in Tpen.cols:
                    Tpen[rid,col]=DDspc_pen[cond][rid];                  
        #special reactions for special conditions - minimum penalty rules
        for cond in DDspc_minPen:
            for rid in DDspc_minPen[cond]:
                 for col in Tpen.cols:
                    Tpen[rid,col]=max(Tpen[rid,col],DDspc_minPen[cond][rid]);                
        self.Tpen_ori=Tpen.mCopy();
        
        ##convert to single-dir model penalties
        N=len(self.M.LLsolver[0]);
        LL=myList.fmake2Dlist([Tpen.colno,N],1);
        Tpenf=myTable.cTable(LL);
        Tpenf.setrowmap(self.M.LLsolver[0]);
        Tpenf.setcolmap(Tpen.cols);
        for row in Tpenf.rows:
            for col in Tpenf.cols:
                rxnroot=row[:-1];
                Tpenf[row,col]=Tpen[rxnroot,col];
                                
        self.Tpen=Tpenf.mCopy();
        self.Tcoef=Tpenf.mCopy(); #for use during FBA, this is where distance would be incorporated
        

    def incorporateDistance(self,Ddist,quiet=False):
        '''Multiplies penalties by 1/(d+1)^self.para['dorder'], where d for a particular rxn is from <Ddist>. 
        If there is a self.para['maxDistance'],then d-values greater than this number will be converted to 
        this number. If <Ddist> has inf AND self.para['maxDistance'] is None AND self.para['dorder'] is 0,
        then there is likely gonna be a problem, so distance order of zero must be used together
        with a maximum distance set.'''
        self.Tcoef=self.Tpen.mCopy(); #in case some other distance was previously incorporated
        self.currentDistances={};
        ## determine upper limit
        if self.para['maxDistance'] is None:
            Ldist=[];
            for x in Ddist:
                if not Ddist[x]==inf:
                    Ldist.append(Ddist[x]);
            maxDist=np.max(Ldist)+1;
            if not quiet:
                print 'Since no upper boundary was given for distance, maximum of current integers plus 1 is used as limit, which is '+str(maxDist);
        else:
            maxDist=self.para['maxDistance'];
        ##incorporate distances into coefficients table
        for rxn in self.Tcoef.rows:
            if Ddist.has_key(rxn):
                d=Ddist[rxn];
                if d==inf:
                    d=maxDist;
            else:
                d=maxDist;
            for col in self.Tcoef.cols:
                self.Tcoef[rxn,col]=self.Tcoef[rxn,col]/((d+1)**self.para['dorder']);
            self.currentDistances[rxn]=d;
            
    def incorporateNetwork(self,Tnetworko,threshold,merge2one=False):
        '''Receives a reaction to condition table (<Tnetworko>, with unidirectional reactions) that must
        have information about confidence levels in reaction fluxes and for each
        condition, sets the constraints of reactions less than <threshold>
        to (0,0) using self.DDconstraints. A self.Tnetwork table will be
        generated, with 0 in constrained and 1 in non-constrained reactions for each
        condition.
        
        If <merge2one>, only reactions with confidence less than <threshold> in all
        conditions will be constrained.
        '''
        LL=myList.fmake2Dlist((len(self.conds),len(self.Dcon_ori)),1);
        T=myTable.cTable.fmakeTable(LL,self.Dcon_ori.keys(),self.conds);
        for row in Tnetworko.rows:
            for col in T.cols:
                if Tnetworko[row,col]>=threshold:
                    T[row,col]=1;
                else:
                    T[row,col]=0;
        if merge2one:
            for row in T.rows:
                maxval=np.max([T[row,col] for col in T.cols]);
                for col in T.cols:
                    T[row,col]=maxval;
        self.Tnetwork=T;
        for row in self.Tnetwork.rows:
            for col in self.Tnetwork.cols:
                if not self.Tnetwork[row,col]:
                    self.DDconstraints[col][row]=[0.,0.];
                
               
    def convertExp2Penalty(self):
        '''Expression table of genes is converted to a penalty table of 
        reactions. The penalty score of a reaction for a condition is
        a function of relative expression level of that gene in that
        condition in GPR blocks. Relative value is determined based on
        max, mean or medan expression (in a block, then minimum value is 
        used across blocks) depending on self.para['Stat']. Penalty score
        is 1/(rel.value)^self.para['norder'].
        
        Penalty table is adjusted to dual model (_I and _L reactions have intestine
        penalty for all tissues and _X reactions have 0 for intestine). 
        '''
        Texpmod=self.Texp.mCopy();
        for row in Texpmod.rows:
            for col in Texpmod.cols:
                Texpmod[row,col]=Texpmod[row,col]+self.para['minExp'];
        #GPR blocks
        Drid2blocks=self.getRxnBlocks(self.M);         
        #Block sums of expression
        LL=myList.fmake2Dlist([Texpmod.colno,self.M.T.rowno],[]);
        TLrxnsums=myTable.cTable(LL);
        TLrxnsums.setcolmap(Texpmod.cols);
        TLrxnsums.setrowmap(self.M.T.getcol('ID'));
        for rid in TLrxnsums.rows:
            for cond in TLrxnsums.cols:
                Lsums=[];
                for Lblock in Drid2blocks[rid]:
                    Lsums.append(self._getBlockSums(Lblock,Texpmod,cond));
                TLrxnsums[rid,cond]=Lsums;
        #Sums -> reaction penalties
        Tpen=TLrxnsums.mCopy();
        for rid in TLrxnsums.rows:
            LL=[];
            for cond in Tpen.cols:
                LL.append(TLrxnsums[rid,cond]);
                A=np.array(LL);
            if A.any():
                if self.para['Stat']=='mean':
                    Aplus=np.concatenate((A,np.reshape(np.mean(A,0),(1,A.shape[1]))),0);
                elif self.para['Stat']=='median':
                    Aplus=np.concatenate((A,np.reshape(np.median(A,0),(1,A.shape[1]))),0);
                elif self.para['Stat']=='max':
                    Aplus=np.concatenate((A,np.reshape(np.max(A,0),(1,A.shape[1]))),0);
                else:
                    raise Exception('Stat parameter is not understood. Use mean, median, or max');
                for j in range(Aplus.shape[1]):
                    for i in range(Tpen.colno):
                        Aplus[i,j]=Aplus[i,j]/Aplus[-1,j];            
                for i in range(Tpen.colno):
                    x=np.min(Aplus[i,:]);
                    Tpen[rid,Tpen.cols[i]]=1./(x**self.para['norder']);  
            else:
                for col in Tpen.cols:
                    Tpen[rid,col]=0;  
        #Dual model correction
        for rid in Tpen.rows:
            if rid[-1]=='I' or rid[-1]=='L':
                for cond in self.non_provider:
                    Tpen[rid,cond]=Tpen[rid,self.provider];
            elif rid[-1]=='X':
                Tpen[rid,self.provider]=0;
        return Tpen    

        
    def calculatePotential(self,rxnid,cond,Ddistance,Dcon_temp_,tabulateFluxes=True,takeLog=True,quiet=False):
        '''
        Calculates flux potential for <rxnid> reaction in the model for condition <cond>.
        
        If <Ddistance> (rxn -> distance from <rxnid>) is provided, then distances are 
        incorporated to modify penalties. If <tabulateFluxes>, then self.DTflux dictionary
        will map <rxnid> to a flux table of conditions. If <takeLog>, self.Dlogtables 
        will map condition to a detailed data of flux calculation (penalties, distances etc.).
        
        Constraints in <Dcon_temp> (typical constraints dictionary) are used but default is restored in
        in the pertaining reactions in the end.
        '''
        ##Preps
        Dcon_temp=deepcopy(Dcon_temp_);
        reverseid=rxnid[:-1]+self.dtoggle[rxnid[-1]];
        Dcon_restore={};
        if self.Dcon_ori.has_key(reverseid) and not Dcon_temp.has_key(reverseid):
            Dcon_temp[reverseid]=[0.,0.];
        for x in self.DDconstraints[cond]:
            if not Dcon_temp.has_key(x):
                Dcon_temp[x]=self.DDconstraints[cond][x];
        for x in Dcon_temp:
            Dcon_restore[x]=self.Dcon_ori[x];
        ##Incorporate distance
        if not Ddistance is None:
            self.incorporateDistance(Ddistance,quiet=quiet);
        ##Set temporary constraints
        self.M.g.meditCons(Dcon_temp);
        ##Set penalty constraints
        LL=[self.Tcoef.rows,self.Tcoef.getcol(cond)];
        self.M.constrain('cWeighedSum',LL,self.para['allowedFlux'],'<=');
        ##Solve
        self.M.fba([[rxnid],[1.]]);
        if self.DDpotential.has_key(rxnid):
            self.DDpotential[rxnid][cond]=self.M.g.optVal;
        else:
            self.DDpotential[rxnid]={cond: self.M.g.optVal};
        ##Write down flux distribution
        if tabulateFluxes:
            self.M.setDflux();
            #Flux table
            if not self.DTflux.has_key(rxnid):
                self.DTflux[rxnid]=myTable.cTable.fmakeTable(\
                myList.fmake2Dlist([self.Tcoef.colno,self.M.T.rowno],0),\
                self.M.T.getcol('ID'),self.Tcoef.cols);
            for row in self.DTflux[rxnid].rows:
                self.DTflux[rxnid][row,cond]=self.M.T[self.M.Drid2row[row],'FLUX'];
        ##Take Log
        if takeLog:
            self.Dlogtables[cond]=myTable.cTable.fmakeTable(myList.fmake2Dlist([6,self.Tcoef.rowno],0),\
            self.Tcoef.rows,['pen','dist','coef','flux','contr','bigg']);   
            for row in self.Dlogtables[cond].rows:
                self.Dlogtables[cond][row,'pen']=self.Tpen[row,cond];
                if Ddistance:
                    self.Dlogtables[cond][row,'dist']=self.currentDistances[row];
                self.Dlogtables[cond][row,'coef']=self.Tcoef[row,cond];
                self.Dlogtables[cond][row,'flux']=self.M.g.Dflux[row];
                self.Dlogtables[cond][row,'contr']=self.Dlogtables[cond][row,'coef']*self.Dlogtables[cond][row,'flux'];
                self.Dlogtables[cond][row,'bigg']=self.M.T[self.M.Drid2row[row[:-1]],'BIGG'];
        ##Restore temporary constraints
        self.M.g.M.remove(self.M.g.M.getConstrByName('cWeighedSum')); 
        self.M.g.meditCons(Dcon_restore);
        self.M.g.M.update();
        
        if not quiet:
            print 'Flux potential for the objective reaction in condition '+cond+' is '+str(self.DDpotential[rxnid][cond]);
            
        return self.DDpotential[rxnid][cond]

        
    def calculateEfficiency(self,quiet=True):
        '''Calculates potential_cond/potential_supercond for each reaction self.DDpotential. Ifsupercond potential has not 
        been calculated for a particular reaction, than this reaction is skipped.'''
        for x in self.DDpotential:
            if not self.DDpotential[x].has_key('supercond'):
                print 'Skipping '+x+' data since supercond calculations are not available';
            elif self.DDpotential[x]['supercond']==0.:
                print 'WARNING: supercondition gives no yield for '+x+'.';
                self.DDefficiency[x]=None;
            else:
                self.DDefficiency[x]={};
                for cond in self.DDpotential[x]:
                    if not cond=='supercond':
                        self.DDefficiency[x][cond]=self.DDpotential[x][cond]/self.DDpotential[x]['supercond'];
        if not quiet:
            for x in self.DDefficiency:
                for cond in self.DDefficiency[x]:
                    print x+'\t'+cond+'\t'+str(self.DDefficiency[x][cond])
        
    def getFluxModel(self,rxnid,cond):
        M=myFlux.mnmGurobi(self.M.T);
        Dflux={};
        for row in self.DTflux[rxnid].rows:
            Dflux[row]=self.DTflux[rxnid][row,cond];
        M.setFlux(Dflux);
        M.setColsOff(['TAG','KEGGID','EC','KEGG','PATHWAY']);
        return M
    
                       
    def _makeEmptyNestedDict(self,addSuper=False):
        DD={};
        for cond in self.conds:
            DD[cond]={};
        if addSuper:
            DD['supercond']={};
        return DD
    
    @staticmethod
    def getNongeneRxns(M):
        Dg2r=M.getDgene2rxn(False);
        S=set();
        for g in Dg2r:
            S=S.union(Dg2r[g]);
        return set(M.T.getcol('ID')).difference(S);
    
    @staticmethod
    def getRxnBlocks(M):
        '''Creates a dictionary from reaction IDs in model M to double lists of GPR blocks
        where different blocks (lists inside the double list) are gene subsets connected with AND,
        and within each block, subsets are connected with ORs, but within this secondary sets,
        if a set has multiple genes, they are connected by ANDs. 
        
        Directly adopted from an earlier script with a step added to correct some exceptions in the end, 
        which was also used elsewhere.

            Created     :   07/03/18
            Modified    :   09/19/19
        '''
        Drid2blocks={}; #rid --> LL with LL is a 2D list of blocks defining the reaction GPR
        DDDrpg=M.getResolvedGeneGroups_adv(); #nested dictionary of reaction to gene associations
        for rid in M.DIrxn:
            DD=DDDrpg[rid];
            Sneg=set();
            Spos=set(); #these sets are made just to confirm assumptions in the code
            LL=[];
            for first in DD:
                for second in DD[first]:
                    L=DD[first][second];
                    if L:
                        LL.append(L);
                        Sblockgenes=set();
                        for S in L:
                            Sblockgenes=Sblockgenes.union(S);
                        if first<0 or second<0:
                            for S in L:
                                Sneg=Sneg.union(S);
                        elif first>0 and second>0:
                            Spos=Spos.union(Sblockgenes);
            if Sneg and Spos:
                if DD[-1]:
                    raise Exception('Reaction '+rid+' seems to have both negative and positive RPG cells filled in the first level, which is against a key assumption in this code');
            elif DD.has_key(-2):
                if Spos:
                    raise Exception('Reaction '+rid+' seems to have both positive RPG cells and a -2 key at the top level, which is against a key assumption in this code');
                else:
                    LL=[[set([g]) for g in Sneg]];
            Drid2blocks[rid]=LL;
        #correction for multiple & connected genes within a single block with no other sets      
        for rid in Drid2blocks:
            LL=Drid2blocks[rid];
            LLmod=[];
            for L in LL:
                if len(L)==1 and len(L[0])>1:
                    for g in L[0]:
                        LLmod.append([set([g])]);
                else:
                    LLmod.append(L);
            Drid2blocks[rid]=myList.fcopy2D(LLmod);  
        return Drid2blocks
        
    @staticmethod
    def _getBlockSums(Lblock,Texp,cond):
        '''UNDER CONSTRUCTION'''
        s=0.;
        for S in Lblock:
            minS=inf;
            for g in S:
                minS=min(minS,Texp[g[:-2],cond]);
            s=s+minS;
        return s

def getExpressionTable(orifile='Data/Caoetal/Ttissue'):
    '''Converts original tissue expression table to a model-specific table from genes to tissues.
    Returns this table.'''
    Tori=op(orifile);
    Dw2n=op('Dw2n_currentModel');
    L=list(set(Tori.rows).intersection(Dw2n.keys()));
    Texp=Tori.getTableByLists(Dcao_tissue.values(),L);
    Texp.setrowmap([Dw2n[w] for w in Texp.rows]);
    return Texp

def getConfidencetables(saveFolder='Data/',saveSfx='',sensDictAt='Data/Dsensitivity_tables',\
    tolZero=1E-6,tolDeltaNfit=0,tolDeltaLow=1E-5):
    '''Uses information in dictionary from <sensDictAt> to generate and save and return two
    confidence tables: one for dual model format (long, to be used as dual model network), the other 
    for regular model format (to be used as confidence heatmap). 
    
    Nonzero flux is determined based on <tolZero>, highest confidence level is determined based
    on <tolDeltaNfit>, and medium level based on <tolDeltaLow>.'''
    #STATIC
    provider='Intestine';
    
    M=getEfficiencyModel_dual();   #to get a single-dir model 
    Dsensori=op(sensDictAt);
    Snon_provider=set(Dcao_tissue.values()).difference([provider]);
    Dcomp2ts={'I':[provider],'L':[provider],'X':Snon_provider,'E':Dcao_tissue.values()};
    #Modification of Dsens
    Dsens={};
    Dtoggle={'f':1.,'r':-1.};
    for ts in Dsensori:
        T=Dsensori[ts].mCopy();
        for row in T.rows:
            vfull=T[row,'integrated_full'];
            vpartial=T[row,'integrated_partial'];
            if not vfull=='na':
                for pfx in ['f','r']:
                    x=T[row,pfx+'_dNfalse'];
                    if x=='na':
                        vfullplus=Dtoggle[pfx]*vfull;
                        vpartialplus=Dtoggle[pfx]*vpartial;
                        if vpartialplus>tolZero or vfullplus>tolZero:
                            for col in [pfx+'_'+para for para in ['dNfalse','dlowFlux','dminFlux']]:
                                T[row,col]=0.;
                        else:
                            for col in [pfx+'_'+para for para in ['dNfalse','dlowFlux','dminFlux']]:
                                T[row,col]=9999.;   
                pfx='z'; #the same modification repeated for zeros    
                x=T[row,pfx+'_dNfalse'];
                if x=='na':
                    vfullplus=abs(vfull);
                    vpartialplus=abs(vpartial);
                    if vpartialplus<tolZero or vfullplus<tolZero:
                        for col in [pfx+'_'+para for para in ['dNfalse','dlowFlux','dminFlux']]:
                            T[row,col]=0.;
                    else:
                        for col in [pfx+'_'+para for para in ['dNfalse','dlowFlux','dminFlux']]:
                            T[row,col]=9999.; 
        Dsens[ts]=T.mCopy();
    #Confidence table (long);
    Drxnsets={};
    for ts in Dsens:
        Drxnsets[ts]=set(Dsens[ts].rows);
    Lridd=M.LLsolver[0];
    Lts=Dsens.keys();
    Lridd.sort();
    Lts.sort();
    Sridd=set(Lridd);
    n=len(Lridd);
    m=len(Lts);
    LL=myList.fmake2Dlist([m,n],0);
    T=myTable.cTable(LL);
    T.setrowmap(Lridd);
    T.setcolmap(Lts);
    Sshort_ridd=set();
    for ridd in T.rows:
        rid=ridd[:-3];
        thisdir=ridd[-1];
        Sshort_ridd.add(rid+thisdir); #for later use to make a short table.
        comp=ridd[-2];
        if thisdir=='f':
            multiplier=1.;
            otherdir='r'; #unused if this does not exist
        else:
            multiplier=-1.;
            otherdir='f'; #unused if this does not exist
        if not Sridd.issuperset([ridd[:-1]+otherdir]):
            otherdir=None;
        for ts in Dcomp2ts[comp]:
            if Drxnsets[ts].issuperset([rid]):
                vofd=Dsens[ts][rid,'integrated_full'];
                if not vofd=='na':
                    vofd_plus=multiplier*vofd;  
                    if vofd_plus<tolZero:
                        if Dsens[ts][rid,thisdir+'_dNfalse']>tolDeltaNfit:
                            T[ridd,ts]=-3;
                        elif Dsens[ts][rid,thisdir+'_dlowFlux']>tolDeltaLow:
                            T[ridd,ts]=-2;
                        else:
                            T[ridd,ts]=-1;
                    elif otherdir:
                        if min(Dsens[ts][rid,'z_dNfalse'],Dsens[ts][rid,otherdir+'_dNfalse'])>tolDeltaNfit:
                            T[ridd,ts]=3;
                        elif min(Dsens[ts][rid,'z_dlowFlux'],Dsens[ts][rid,otherdir+'_dlowFlux'])>tolDeltaLow:
                            T[ridd,ts]=2;
                        else:
                            T[ridd,ts]=1;  
                    else:
                        if Dsens[ts][rid,'z_dNfalse']>tolDeltaNfit:
                            T[ridd,ts]=3;
                        elif Dsens[ts][rid,'z_dlowFlux']>tolDeltaLow:
                            T[ridd,ts]=2;
                        else:
                            T[ridd,ts]=1;  
    #Correction for intestinal reactions
    for row in T.rows:
        if row[-2]=='I' or row[-2]=='L':
            for col in Snon_provider:
                T[row,col]=T[row,provider];
    #Making the short confidence table
    Lrid=list(Sshort_ridd);
    Lrid.sort();
    Tshort=T.fmakeTable(myList.fmake2Dlist([T.colno,len(Lrid)],0),Lrid,T.cols);
    for row in T.rows:
        ridd=row[:-3]+row[-1];
        comp=row[-2];
        if comp=='I':
            Tshort[ridd,provider]=T[row,provider];
        elif comp=='X':
            for col in Snon_provider:
                Tshort[ridd,col]=T[row,col]; #test        
        elif comp=='E':
            for col in T.cols:
                Tshort[ridd,col]=T[row,col]; #test
        elif not ridd[0]=='T': #This may not be important as L reactions are not part of the common network
            Tshort[ridd,provider]=T[row,provider];
    #Converting zeros to -3 for intestine for X reactions 
    for row in T.rows:
        if row[-2]=='X':
            T[row,provider]=-3;
    #Acknowledging that sensitivity analysis is lacking for L transports
    for row in T.rows:
        if row[0]=='T' and row[-2]=='L':
            for col in T.cols:
                T[row,col]=0;
    
    sv(T,saveFolder+'Tconf_dual'+saveSfx);
    sv(Tshort,saveFolder+'Tconf'+saveSfx);
    return T,Tshort


def getDistanceDict(name,X,Di,Dj):
    '''Picks a row (at <name>) from <X> matrix and converts to a dictionary
    using row and column indices in <Di> and <Dj>. I no rows are present for name,
    then an empty dictionary is returned.'''
    D={};
    if Di.has_key(name):
        i=Di[name];
        for x in Dj:
            D[x]=X[i,Dj[x]];
    return D
    
    
def metabolite2terminals_dualModel(M,comp):
    '''Derives, for every metabolite in M, input, output and transport reactions categorized
    as "uptake","exchange","demand", and "transport" in a dictinary of category -> list.
    The model must have TAGs in modified model format (u for uptake etc.).'''
    DDcpd2io={};
    if comp=='I' or comp=='L':
        Scomp=set('IL');
    else:
        Scomp=set(comp);
    for cpd in M.Dbigg2cel:
        cpdRoot=cpd[:-2];
        DDcpd2io[cpdRoot]={'uptake':[],'exchange':[],'demand':[],'transport':[],'sink':[]}; #overwriting is allowed
    for rid in M.Drbigg:
        if Scomp.issuperset(rid[-1]):
            Tag=M.T[M.Drid2row[rid],'TAG'];
            if Tag=='u':
                if rid[2]=='K':
                    for spc in M.Drbigg[rid].reactants[0]:
                        DDcpd2io[spc[:-2]]['sink'].append(rid);
                else:
                    for spc in M.Drbigg[rid].reactants[0]:
                        DDcpd2io[spc[:-2]]['uptake'].append(rid);
            elif Tag=='e' or Tag=='ae':
                for spc in M.Drbigg[rid].reactants[0]:
                    DDcpd2io[spc[:-2]]['exchange'].append(rid);   
            elif Tag=='s' or Tag=='d':
                for spc in M.Drbigg[rid].reactants[0]:
                    DDcpd2io[spc[:-2]]['demand'].append(rid); 
            elif rid[:3]=='TCE' or rid[:3]=='TEC':            
                for spc in M.Drbigg[rid].reactants[0]:
                    DDcpd2io[spc[:-2]]['transport'].append(rid);  
    return DDcpd2io 


def efficiency(ID,dorder,Lgiven=None,sensTable='Data/Tconf_dual',**kwargs):   
    '''Given a dual-model reaction (with direction in the end as "f" or "r") or a dual-model
    metabolite (again with "f" or "r" in the end, representing production and consumption,
    respectively) as <ID>, and using distance of <dorder>, determines flux efficiencies
    in all relevant tissues, or thiose in <Lgiven> if this is provided, for this reaction
    or metabolite. To do this, an instance of of fluxEfficiency_dual is generated, efficiencies are
    calculated, and the instance is returned. Tissue networks are informed by the table at <sensTable>.
    Make this parameter None, to avoid defining tissue networks. Gene expression levels are obtained
    using default input for getExpressionTable().
    
    For other arguments, see "Dpara" in DYNAMIC.
    '''
    ###DYNAMIC
    Dpara={'uptakePen':1.,'exchangePen':0.,'bacUptakePen':0.,'cpdSinkPen':1.,'norder':1,'uptakeDist':0,\
    'mergenetworks2one':False,'maxDistance':31,'X2Idist':0,'allowedFlux':1.,'nongenePen':1.,'maxBacIntake':1000.};
    ###STATIC
    provider='Intestine';
    rxnDistanceFile='/project/umw_marian_walhout/data/yilmazl/Tissue/distance/distance_rxn2rxn';
    cpdProductionDistanceFile='/project/umw_marian_walhout/data/yilmazl/Tissue/distance/distance_metabolite_production';
    cpdConsumptionDistanceFile='/project/umw_marian_walhout/data/yilmazl/Tissue/distance/distance_metabolite_consumption';
    bacterialUptakeRxn='EXC0050_Lr';
    byproducts_noncomp=['co2','amp','nadp','nadph','ppi','o2','nadh','nad','pi','adp','coa', 'atp', 'h2o', 'h', 'gtp', 'gdp','etfrd','etfox']+['na1','oh1','cl','i'];
    dgrPen=0.; #uniform penalty for degradation reactions
    ###PREPS
    #identify type of entry and rearrange the input for algorithm
    if ID[-2]==']':
        cpdORrxn='cpd';
        consume={'f':False,'r':True}[ID[-1]];
        ID=ID[:-1];
    else:
        consume=None;
        cpdORrxn='rxn';
    #compartment
    if cpdORrxn=='rxn':
        comp=ID[-2];
    else:
        comp=ID[-4];     
    #Conditions to be evaluated
    if Lgiven:
        conds=Lgiven;
    else:
        if comp=='I' or comp=='L':
            conds=[provider];
        elif comp=='X':
            conds=list(set(Dcao_tissue.values()).difference([provider]));
        else:
            conds=Dcao_tissue.values();

    #Modifications
    if type(dorder)==str:
        dorder=float(dorder);
    if type(consume)==str:
        consume=eval(consume);
    for key, value in kwargs.items():
        if Dpara.has_key(key):
            Dpara[key]=value;
        else:
            raise Exception('Argument '+key+' is not understood. Available args are:\n'+str(Dpara.keys()));    

    uptakePen=Dpara['uptakePen'];
    exchangePen=Dpara['exchangePen'];
    bacUptakePen=Dpara['bacUptakePen'];
    cpdSinkPen=Dpara['cpdSinkPen'];
    uptakeDist=Dpara['uptakeDist'];
    mergenetworks2one=Dpara['mergenetworks2one'];
    maxDistance=Dpara['maxDistance'];
    norder=Dpara['norder'];
    X2Idist=Dpara['X2Idist'];
    allowedFlux=Dpara['allowedFlux'];
    nongenePen=Dpara['nongenePen'];
    maxBacIntake=Dpara['maxBacIntake'];

    if cpdORrxn=='cpd':
        cpd=ID;
        rxnid='SINKf';
    elif cpdORrxn=='rxn':
        cpd=None;
        rxnid=ID;
    else:
        raise Exception('cpdORrxn input not understood. It should be either "cpd" or "rxn".');
        
    ##Distance dictionary
    if cpd is None:
        X,Di=op(rxnDistanceFile);
        Ddist=getDistanceDict(rxnid,X,Di,Di);
        Dj=Di;
    elif not consume:
        X,Di,Dj=op(cpdProductionDistanceFile);
        Ddist=getDistanceDict(cpd,X,Di,Dj);
        Ddist[rxnid]=0;
    else:
        X,Di,Dj=op(cpdConsumptionDistanceFile);
        Ddist=getDistanceDict(cpd,X,Di,Dj);
        Ddist[rxnid]=0; 
    #special distances
    if comp=='X':
        for k in Dj:
            if k[-2]=='I':
                Ddist[k]=X2Idist;
    for x in Dj:
        if x[:2]=='UP':
            Ddist[x]=uptakeDist;
    Ddist[bacterialUptakeRxn]=uptakeDist;

    ##Metabolite-run-specific
    if cpd is None:
        M=getEfficiencyModel_dual(bacIntake=maxBacIntake);
        print ID;
    elif cpd[-8:-5]=='coa' and not consume:
        M=getEfficiencyModel_dual(cpd=cpd,consume=consume,addReactants=['h2o'+cpd[-5:]],addProducts=['coa'+cpd[-5:],'h'+cpd[-5:]],bacIntake=maxBacIntake);
    else:
        M=getEfficiencyModel_dual(cpd=cpd,consume=consume,bacIntake=maxBacIntake);
    
    ##Sensitivity Table
    if sensTable:
        Tsens=op(sensTable);

    ##Experimental data
    Texp=getExpressionTable();

    ##Special penalties
    Dspc_pen={};
    #exchange, uptake, and degradation
    for rid in M.T.getcol('ID'):
        if rid[:2]=='UP':
            Dspc_pen[rid]=uptakePen;
        elif rid[:2]=='EX' and not rid==bacterialUptakeRxn[:-1]:
            Dspc_pen[rid]=exchangePen;
        elif rid[:3]=='DGR':
            Dspc_pen[rid]=dgrPen;
    Dspc_pen[bacterialUptakeRxn[:-1]]=bacUptakePen;
    if not cpd is None:
        Dspc_pen['SINK']=cpdSinkPen;
    #nullify intestinal transport for X
    DDspc_pen={};
    for cond in conds:
        DDspc_pen[cond]={};
    if not conds==[provider]:
        for rid in M.T.getcol('ID'):
            if rid[-1]=='I' or rid[-1]=='L':
                if rid[:3]=='TCE' or rid[:3]=='TEC':
                    for cond in set(conds).difference([provider]):
                        DDspc_pen[cond][rid]=1E-6;

    ##Special constraints
    rxnroot=rxnid[:-1];
    Tag=M.T[M.Drid2row[rxnroot],'TAG'];
    if Tag=='e' or Tag=='ae' or Tag=='s' or Tag=='d' or Tag=='u' or rxnid[:3]=='TCE' or rxnid[:3]=='TEC' or rxnroot=='SINK':
        Lspc=M.Drbigg[rxnroot].reactants[0]+M.Drbigg[rxnroot].products[0];
    else:
        Lspc=[];
    Sspc=set([spc for spc in set([s[:-2] for s in Lspc]).difference(byproducts_noncomp)]);
    if Sspc:
        Dspcon={};
        Dcpd2io=metabolite2terminals_dualModel(M,comp); 
        for spc in Sspc:
            if rxnid[:3]=='TCE' or rxnid[:3]=='TEC':
                Snull=set(Dcpd2io[spc]['sink']+Dcpd2io[spc]['demand']+Dcpd2io[spc]['transport']).difference([rxnroot]);
            elif Tag=='d' or Tag=='s':
                Snull=set(Dcpd2io[spc]['sink']+Dcpd2io[spc]['uptake']+Dcpd2io[spc]['transport']).difference([rxnroot]);
            elif Tag=='u':
                Snull=set(Dcpd2io[spc]['demand']+Dcpd2io[spc]['transport']).difference([rxnroot]);
            elif Tag=='e' or Tag=='ae':
                Snull=set(Dcpd2io[spc]['sink']+Dcpd2io[spc]['uptake']).difference([rxnroot]);
            elif rxnroot=='SINK':
                Snull=set(Dcpd2io[spc]['sink']+Dcpd2io[spc]['uptake']+Dcpd2io[spc]['demand']+Dcpd2io[spc]['exchange']+Dcpd2io[spc]['transport']).difference([rxnroot]);
            else:
                raise Exception('Categorization of related terminal reactions failed. Check the code!')
            for rxn in Snull:
                if M.LLsolver[0].count(rxn+'f'):
                    Dspcon[rxn+'f']=[0.,0.];
                if M.LLsolver[0].count(rxn+'r'):
                    Dspcon[rxn+'r']=[0.,0.];
    else:
        Dspcon={};
    
    ##Creation of effciency instance
    eff=fluxEfficiency_dual(M,Texp,norder=norder,dorder=dorder,maxDistance=maxDistance,allowedFlux=allowedFlux,nongenePen=nongenePen); #np.max(X[X!=inf])+1    
    eff.setPenalties(Dspc_pen=Dspc_pen,DDspc_pen=DDspc_pen); 
    if sensTable:
        eff.incorporateNetwork(Tsens,-1.1,mergenetworks2one);
    
    ##Correction for supercond in the case of sole provider
    if conds==[provider]:
        for ridd in eff.DDconstraints[provider]:
            if ridd[-2]=='X':
                eff.DDconstraints['supercond'][ridd]=eff.DDconstraints[provider][ridd];
    
    ##Caclculation
    for cond in conds+['supercond']: 
        print cond
        print Dspcon
        if cpdORrxn=='cpd':
            _=eff.calculatePotential('SINKf',cond,Ddist,Dspcon);
        else:
            _=eff.calculatePotential(ID,cond,Ddist,Dspcon);
    eff.calculateEfficiency(True);
    return eff
    

def generalEfficiency(ID,dorder,startk=1,saveToFolder=None,saveLogToSubfolder=None,sensTable='Data/Tconf_dual',**kwargs):
    '''Uses efficiency() to establish and save (if not <saveToFolder> is None) a 
    tissue -> flux efficiency dict ("DDeff") and a tissue -> logtable dict ("DDlog").
    
    ID can be reaction ID or metabolite  ID as in modmodel, singledir format, with "f" and "r"
    also used in metabolite names as in efficiency() (to indicate production and consumption, respectively)
   (e.g., RM04432f, collg[c]f). This will be converted to intestine or other tissue element
    using _I or _X, respectively, except that, if ID denotes a TCE type reaction, then intestinal ID
    is labeled by _L.
    
    dorder is a list or number. For every number given, the efficiency calculation will be repeated.
    Every step i in this loop produces two tuples ("DDeff",dorder_i) and ("DDlog",dorder_i). Tuples matching different
    dorders are returned as two separate lists, "Leff" and "Llog".
    
    Given valid <saveToFolder>, which should be the path to a folder, tuples in the above mentioned "Leff" list 
    are gonna be saved after mod<ID>_i, where modID has [ and ] replaced by __ if ID is a metabolite ID
    and i is order of dorder used in dorder list. Given valid <saveToSubfolder>, which should be a subfolder
    under <SaveToFolder>, the same is done for "Llog".
    
    If <startk> is different from 1, the _i in saved filenames will start from <startk>.
    '''
    #STATIC
    modelFile='Tmultiple_2'; #model needed for TCE-type transport to determine the side (lumen or other tissue) of transport
    Stransby=set(['atp','h2o','na1','oh1','cl','i','h']); #also for deciphering intestine transport
    ##
    startk=int(startk);
    Dir=ID[-1];
    if not Dir=="f" and not Dir=="r":
        raise Exception('Invalid ID. The name should have "f" or "r" in the end.')
    if not type(dorder)==list:
        Ldorder=[float(dorder)];
    else:
        Ldorder=[float(d) for d in dorder];
        
    if ID[-2]==']':
        iscpd=True;
        IDi=ID[:-4]+'_I'+ID[-4:];
        IDx=ID[:-4]+'_X'+ID[-4:];
    else:
        iscpd=False;
        IDi=ID[:-1]+'_I'+ID[-1];
        IDx=ID[:-1]+'_X'+ID[-1];
        
    ##direction of intestine transport
    if ID[:3]=='TCE' or ID[:3]=='TEC':
        M=loadModel(modelFile);
        rid=IDi[:-1];
        rid_L=rid[:-1]+'L';
        if M.Drbigg.has_key(rid_L): #means lumen counterpart exists
            Sby=set();
            for sfx in ('_I','_E'):
                for cpd in Stransby:
                    Sby.add(cpd+sfx);
            r=M.Drbigg[rid];
            r.mlocalize();
            Sbyloc=set();
            Sloc=set();
            for i in range(len(r.reactants[0])):
                spc=r.reactants[0][i];
                loc=r.rlocs[i];
                if Sby.issuperset([spc]):
                    Sbyloc.add(loc);
                else:
                    Sloc.add(loc);
            if len(Sloc)==1:
                loc=list(Sloc)[0];
            elif not Sloc and len(Sbyloc)==1:
                loc=list(Sbyloc)[0];
            else:
                loc='na';
            if (loc=='e' and Dir=='f') or (loc=='c' and Dir=='r'):
                IDi=rid_L+Dir;
   
    Leff=[];
    Llog=[];
    for dorder in Ldorder:
        teff=({},dorder);
        tlog=({},dorder);
        effi=efficiency(IDi,dorder,sensTable=sensTable,**kwargs);
        effx=efficiency(IDx,dorder,sensTable=sensTable,**kwargs);
        if not iscpd:
            if effi.DDefficiency[IDi]:
                teff[0].update(effi.DDefficiency[IDi]);
            if effx.DDefficiency[IDx]:
                teff[0].update(effx.DDefficiency[IDx]);
        else:
            if effi.DDefficiency['SINKf']:
                teff[0].update(effi.DDefficiency['SINKf']);
            if effx.DDefficiency['SINKf']:
                teff[0].update(effx.DDefficiency['SINKf']);        
        tlog[0].update(effi.Dlogtables);
        tlog[0].update(effx.Dlogtables);
        tlog[0]['supercond_I']=effi.Dlogtables['supercond'];
        Leff.append(teff);
        Llog.append(tlog);
        
    if saveToFolder and not saveLogToSubfolder:
        fname=re.sub('[\[\]]','__',ID);
        for i in range(len(Leff)):
            sv(Leff[i],saveToFolder+fname+'_'+str(i+startk));
    elif saveToFolder and saveLogToSubfolder:
        fname=re.sub('[\[\]]','__',ID);
        for i in range(len(Leff)):
            sv(Leff[i],saveToFolder+fname+'_'+str(i+startk));
            sv(Llog[i],saveToFolder+saveLogToSubfolder+fname+'_'+str(i+startk));
    elif saveLogToSubfolder:
        print 'WARNINGS: Logs are not saved because saveLogToSubfolder is provided without saveToFolder. Read the code and use the correct guidance for saving data';
        
    return Leff,Llog
 