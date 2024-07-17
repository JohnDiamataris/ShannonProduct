from __future__ import print_function
from sympy import Matrix,symbols,Function,Symbol,cse,expand
from scipy.linalg import orth
from sympy.solvers import solve,solve_poly_system,nsolve
from sympy.core.symbol import Dummy
from sympy.parsing.sympy_parser import parse_expr
from sympy.abc import x, y, z
import numpy as np

#sympifyexpr = True
class phen_to_gen() :
# From Phenotypes to Genotypes based on Shannon algebra 
    def __init__(self,alleletypes):

        self.solved=False
        self.used_non_eq=[]
        self.symetryeq=[]        
        self.genotypes=[]
        self.setupgen(alleletypes)
        self.set_gent()
#        Cret=self.rec_gent()
        self.genotypes_dic={}
        for igen in range(len(self.genotypes)):
#            strloc=str(self.genotypes[igen][0])
#            for ilock in self.genotypes[igen][1:] :
#                strloc += ";"+str(ilock)
            strloc= self.gen_to_string(self.genotypes[igen])
            self.genotypes_dic[strloc]=igen

    def checksympl(self,sympifyexpr):
        self.sympifyexpr = sympifyexpr

    def initialize(self,A,B):
        self.A=A
        self.L=[[0 for x in range(len(self.genotypes_dic))] for y in range(len(self.aplotypes_dic))] 
        for iapl in self.aplotypes_dic :
                indaplo=self.aplotypes_dic[iapl]
                for igen in range(len(self.genotypes)):
                    strloc= self.gen_to_string(self.genotypes[igen])
                    intexgen=self.genotypes_dic[strloc]
                    coppys=strloc.count(iapl)
                    self.L[indaplo][intexgen]=coppys/(2.)        
#                    self.L[indaplo][intexgen]=coppys/(2.*len(self.aplotypes))         
#  check
        sumtest=[]
        for irow in self.L:
            sumtest.append(0.)
        int_row=-1    
        for jrow in self.L:
            int_row=int_row+1
            for icol in range(len(jrow)):
                sumtest[int_row]=sumtest[int_row]+jrow[icol]

        self.L=Matrix(self.L)
#  the left and the rigth vector space of L is an messure of the independent geotypes and applotypes        
        self.B=B
        self.ncols=A.cols
        self.nrows=A.rows
        self.Lncols=self.L.cols
        if (self.ncols !=self.Lncols) :
            print ("oupts ")
            exit()

    def aplotypes_dicD(self):
        return self.aplotypes_dic

    def gen_to_aplo_List(self):
        return self.L.tolist()

    def gen_to_string(self,igenotype):
        if  type(igenotype)==list: 
            strloc=str(igenotype[0])
            for ilock in igenotype[1:] :
                strloc += ";"+str(ilock)
        else :
            strloc=igenotype
        return strloc

    def phen_gent(self):
        
        pass
    def genotypesL(self):
        return self.genotypes

    def genotypesM(self):
        return Matrix(self.genotypes)

    def genotypesD(self):
        return self.genotypes_dic

    def solve(self) :       
        self.ncols=self.A.cols
        self.nrows=self.A.rows      
        self.Bnrows=self.B.rows
        if (self.nrows !=self.Bnrows) :
            print ("oupts ")
            exit()

        Xw=self.A.pinv_solve(self.B)
        X=self.A.pinv_solve(self.B, arbitrary_matrix=Matrix.zeros(self.ncols,1))
        self.solved=True
        return (Xw,X)

    def write_to_exfile(self,exe_filepath,w,wl='wlist'):
        exe_file = open(exe_filepath,'w+')
#        sttr = "    def loadw():\n"
        sttr = "from sympy import symbols\n"
        sttr +="global "+wl+"\n"
        sttr +=wl+"=[]\n"            
        for iw in w:
            sttr +="global "+str(iw)+"\n"
            sttr +=str(iw)+"= symbols('"+str(iw)+"')"+"\n"
            sttr +=wl+".append("+str(iw)+")\n"
        exe_file.write(sttr)
        exe_file.close()
        execfile(exe_filepath,globals())

    def allpos(self,list):
        for i in list:
            if i.is_Number :
                if i < 0. :
                    return False
        return True

# example
    def apop(self,B,SelfCons=True) :
#    ''' Analitical allele frequency estimation based on phynotipic frequency '''         
        self.B = B
#   Lamda ij_kl ... genotype        
        Xw,X=self.solve()
        lw=self.L*Xw
        l=self.L*X
#    lamda i.... (allele frequensy in the population)
        n=lw.rows
        one=Matrix.ones(1,n)
        w = symbols('w:{0}_:{1}'.format(self.nrows, self.ncols), cls=Dummy)
        self.write_to_exfile("temp_w.py",w)
#        f=[]
#        for iS in self.Sim :
#            i_,j_=iS
#            f.append(parse_expr(str(Xw[i_]-Xw[j_])))
#        for i_ in range(len(lw)) :
#            for j_ in range(len(lw)) :
#                f.append(parse_expr(str(lw[i_]*lw[j_]-Xw[i_*(n-1)+j_])))                    
        if SelfCons :
            f,fall=self.set_non_lin_eqSelfCons(Xw,lw)
        else :             
            f,fall=self.set_non_lin_eq(Xw,lw)
        fs=tuple(f)
        Sep_solutions=[]
        solutions=solve(fs)
        self.nfreedom=self.A.cols-self.A.rows
        if len(solutions)==0 :
#            for ifs in fs :
#                Sep_solutions.append(solve(ifs))
#            solutions=Sep_solutions[0]
            for i in range(len(fs)-1) :
                Sep_solutions.append(solve(fs[i:i+self.nfreedom]))
            solutions=Sep_solutions[0]
        global wlist
        
        S_lamdas=[]
        S_lamdasAN=[]
        S_lamdasw=[]
        for isol in solutions :
            lamda=[]
            lamdaw=[]
            sub=[]
            for ilist in wlist :
                if ilist in isol.keys():
                    iso=isol[ilist]
                    sub.append((ilist,iso))
            for iXw in Xw:
                    try :
                        lamda.append(parse_expr(str(iXw)).subs(sub))
                        lamdaw.append(parse_expr(str(iXw)))
                    except :
                        pass
            if (sum(lamda)-1.)**2<0.000000001 and self.allpos(lamda) :             
                S_lamdas.append(self.L*Matrix(lamda))
                S_lamdasw.append(self.L*Matrix(lamdaw))            
            else :
                S_lamdasAN.append(self.L*Matrix(lamda))
        self.solved=True        
        return S_lamdas,S_lamdasw,S_lamdasAN

    ''' setup the set of non linear equations nased on Xw,lw , self.L'''
    def set_non_lin_eq(self,Xw,lw):        
        self.used_non_eq=[]
        min_ttf=[]
        ttf=[]
        for ix in range(len(Xw.tolist())):
            strloc= self.gen_to_string(self.genotypes[ix])
            refalel=0
            strtext=""
            strtextEQ=""
            use_firstAllel=False
            for i_ in range(len(lw)) :
                if self.L.tolist()[i_][ix] > 0 :
                    if self.L.tolist()[i_][ix] > 0.6 :
# select only the first row of idipendend non linear equations based on 0 allel 
                        if i_==refalel :
                            xy=self.genotypes[ix]
                            ifm=[]
                            if type(xy)==list:
                                for ixy in xy :
                                    ifm.append(ixy.split(","))
                            else :
                                ifm.append(xy.split(","))
                            test = False
                            for itest in range(len(ifm)):
                                test = self.aplotypes_dic[ifm[itest][0]]==refalel
                                if test :
                                    break 
                            if test :
                                use_firstAllel=True
                        if strtext=="":
                            strtext=strtext+"("+str(lw[i_])+")*("+str(lw[i_])+")"
                            strtextEQ=strtextEQ+"("+self.aplotypeList[i_]+")^2"
                        else :    
                            strtext=strtext+"*("+str(lw[i_])+")*("+str(lw[i_])+")"
                            strtextEQ=strtextEQ+"*("+self.aplotypeList[i_]+")^2"
                    else :
                        if i_==refalel :
                            xy=self.genotypes[ix]
                            ifm=[]
                            if type(xy)==list:
                                for ixy in xy :
                                    ifm.append(ixy.split(","))
                            else :
                                ifm.append(xy.split(","))
                            test = False
                            for itest in range(len(ifm)):
                                test = self.aplotypes_dic[ifm[itest][0]]==refalel
                                if test :
                                    break 
                            if test :
                                use_firstAllel=True
                        if strtext=="":
                            strtext=strtext+"("+str(lw[i_])+")"
                            strtextEQ=strtextEQ+"("+self.aplotypeList[i_]+")"                            
                        else :    
                            strtext=strtext+"*("+str(lw[i_])+")"
                            strtextEQ=strtextEQ+"*("+self.aplotypeList[i_]+")"                            
            strtext=strtext+"+("+str(-Xw[ix])+")"
            strtextEQ=strtextEQ+"-["+self.gen_to_string(self.genotypes[ix])+"]"                            

            if use_firstAllel :
               min_ttf.append(parse_expr(strtext))
               self.used_non_eq.append(strtextEQ)            
            ttf.append(parse_expr(strtext))
        return min_ttf,ttf   

    ''' Implement Estimation of Genotype frequency via shannon product,n_Xw=len(Xw.tolist())'''
    def shannon_product(self,n_Xw,lw):
        pass

    ''' Estimation of Genotype frequency via genetic equlibrium of shannon product,n_Xw=len(Xw.tolist()) '''
    def est_geno_from_alle(self,n_Xw,lw,as_equl=True):        
        self.aplot_genot_eqD={}
        Xw_est=[]
#        for iw in range(len(Xw.tolist())):
        for iw in range(n_Xw):
            Xw_est.append(0.)
        if as_equl :
            for ix in range(n_Xw):
                strtext=""
                strtextEQ=""
                for i_ in range(len(lw)) :
                    if self.L.tolist()[i_][ix] > 0 :
                        if self.L.tolist()[i_][ix] > 0.6 :
    # select only the first row of idipendend non linear equations based on 0 allel 
                            if strtext=="":
                                strtext=strtext+"("+str(lw[i_])+")*("+str(lw[i_])+")"
                                syb_text=((lw[i_])*(lw[i_]))
                                strtextEQ=strtextEQ+"("+self.aplotypeList[i_]+")^2"
                            else :    
                                strtext=strtext+"*("+str(lw[i_])+")*("+str(lw[i_])+")"
                                syb_text=syb_text*((lw[i_])*(lw[i_]))
                                strtextEQ=strtextEQ+"*("+self.aplotypeList[i_]+")^2"
                        else :
                            if strtext=="":
                                strtext=strtext+"("+str(lw[i_])+")"
                                syb_text=(lw[i_])
                                strtextEQ=strtextEQ+"("+self.aplotypeList[i_]+")"                            
                            else :    
                                strtext=strtext+"*("+str(lw[i_])+")"
                                syb_text=syb_text*(lw[i_])
                                strtextEQ=strtextEQ+"*("+self.aplotypeList[i_]+")"                            
    # estimation of genotype frequency based on genetic equlibrium
                Xw_est[ix]=syb_text
    #            strtext=strtext+"+("+str(-Xw[ix])+")"
                strgenotype=self.gen_to_string(self.genotypes[ix])
                self.aplot_genot_eqD[strgenotype] = [syb_text,strtext,strtextEQ]
        else :
            # implement shannon product
            self.shannon_product(n_Xw,lw)
# implement
            exit()
# estimation of phenotypes based on genetic equlibrium
        Xw_est=Matrix(Xw_est)
        return Xw_est

#  solve for one variable so you can substitute  
    def subvarequat(self,xvar,dic,Ndic={}) :
        maxv = -1
        for ivar in dic.keys() :
            co = abs(dic[ivar])
            if co > maxv and ivar not in Ndic.keys():
                maxv = co
                smaxv =  dic[ivar]
                refv = ivar

        for ivar in dic.keys() :
            if abs(dic[ivar])> 1e-14 :
                dic[ivar] = -dic[ivar]/smaxv
            else :
                dic[ivar] = 0
        sc =xvar/smaxv
        for ivar in dic.keys() :
            if ivar != refv :
                sc = sc + dic[ivar]*ivar
        Ndic[refv]=parse_expr(str(sc))
        return Ndic 



    def subvareigenv(self,xvar,dic,Ndic={}) :
        maxv = -1
        for ivar in dic.keys() :
            co = abs(dic[ivar])
            if co > maxv and ivar not in Ndic.keys():
                maxv = co
                smaxv =  dic[ivar]
                refv = ivar

        for ivar in dic.keys() :
            if abs(dic[ivar])> 1e-14 :
                dic[ivar] = -dic[ivar]/smaxv
            else :
                dic[ivar] = 0
        sc =xvar/smaxv
        for ivar in dic.keys() :
            if ivar != refv :
                sc = sc + dic[ivar]*ivar
        Ndic[refv]=parse_expr(str(sc))
        return Ndic 
    def reversetuple(self,tup1):
        if len(tup1)==2 :
            tup2=(tup1[1],tup1[0])
        else :
            print ("Check this case")
        return tup2

    def symplify_results(self,vec) :
#        w = symbols('w:{0}_:{1}'.format(self.nrows, self.ncols), cls=Dummy)
#        self.write_to_exfile("temp_w.py",w)

        global wlist
        Lvec=vec.tolist()
        varble=[]
        maxvar=0
        for ieq in range(len(Lvec)):
            maxvar=max(maxvar,len(Lvec[ieq][0].free_symbols))
            temv=list(Lvec[ieq][0].free_symbols)
            varble=list(set().union(varble,temv))
        if len(varble) == 0 :
            return Matrix(vec)
        Svarble=[]
        for ivar in range(len(varble)) :
            Svarble.append(Symbol(str(varble[ivar])))
        Jac=vec.jacobian(varble)
#        npJac=np.array(Jac).astype(np.float64)
#        nullspace=Jac.nullspace()
        Tnullspace=(Jac.T).nullspace()

        return Tnullspace

    def used_sym_eqR(self) :
        return [self.res_S_l_str,self.res_S_l_Const,self.res_S_Xw_str,self.res_S_Xw_Const] 

    def sym_relations(self) :
        return [self.res_S_l,self.res_S_Xw] 


    def symplify_depen(self,vec,usesolvesub=False) :
#        w = symbols('w:{0}_:{1}'.format(self.nrows, self.ncols), cls=Dummy)
#        self.write_to_exfile("temp_w.py",w)
        if False :
            Nvec=[]
            tNvec=cse(vec.expand())
            for iv in range(len(tNvec[0])):
                rtupl = self.reversetuple(tNvec[0][iv])
                Jac0=list(rtupl[0]).jacobian(rtupl[0].free_symbols)
                for ive in range(len(vec)):
                    scrt = str(tNvec[1][0][ive])
                    Nvec.append(parse_expr(scrt).subs(rtupl))

        global wlist
        Lvec=vec.tolist()
        varble=[]
        maxvar=0
        for ieq in range(len(Lvec)):
            maxvar=max(maxvar,len(Lvec[ieq][0].free_symbols))
            temv=list(Lvec[ieq][0].free_symbols)
            varble=list(set().union(varble,temv))
        if len(varble) == 0 :
            return Matrix(vec)
        Svarble=[]
        for ivar in range(len(varble)) :
            Svarble.append(Symbol(str(varble[ivar])))
        Jac=vec.jacobian(varble)
 #       nullspace=Jac.nullspace()
        npJac=np.array(Jac).astype(np.float64)
#        orthnorm = orth(npJac)

        orthnorm2 = orth(npJac.T)
        if type(orthnorm2[0]) is np.ndarray :
            n2deg = len(orthnorm2)
            ndeg = len(orthnorm2[0])
        else :
            ndeg = 1
            n2deg = len(orthnorm2)            
        x = symbols('x:{0}'.format(ndeg), cls=Dummy)
        self.write_to_exfile("temp_x.py",x,wl="xlist")
        if not usesolvesub :
            subexpress={}
            for ideg in range(ndeg):
                dic1={}
                for iv in  range(n2deg):
                    dic1[Svarble[iv]]=orthnorm2[iv][ideg]
                subexpress=self.subvareigenv(xlist[ideg],dic1,subexpress)

        Nvec=[]
        for ive in vec:
            Nvec.append(ive)
        sc=[]
        for ideg in range(ndeg):
            sc.append(0.)
        for ideg in range(ndeg):
            for iv in  range(n2deg):            
                if abs(orthnorm2[iv][ideg])> 1e-14 :
                    sc[ideg]=sc[ideg]+ orthnorm2[iv][ideg]*Svarble[iv]
        if usesolvesub :
            pstrsc =[]
            for ideg in range(ndeg):                
                strsc = "("+ str(sc[ideg])+"-"+str(xlist[ideg])+")"
                pstrsc.append(parse_expr(strsc))
            resol=solve(pstrsc)

            sub=[]
            for ikey in resol.keys():
                sub.append((ikey,resol[ikey]))
            for ive in range(len(Nvec)):
                Nvec[ive]=parse_expr(str(Nvec[ive])).subs(sub)
        else :
            for ideg in range(ndeg):
                sub=[]
                strsc = "("+ str(sc[ideg])+"-"+str(xlist[ideg])+")"
                pstrsc = parse_expr(strsc)
                resol=subexpress
                for ikey in resol.keys():
                    sub.append((ikey,resol[ikey]))
                for ive in range(len(Nvec)):
                    Nvec[ive]=parse_expr(str(Nvec[ive])).subs(sub)
        if ndeg >0 :
            tvarble=xlist
            for ive in range(len(Nvec)):
                temw=list(Nvec[ive].free_symbols)
                tvarble=list(set().union(tvarble,temw))
            global wlist
            wlist = tvarble
        return Matrix(Nvec)

    ''' setup the set of non linear equations nased on Xw,lw , self.L'''
    def set_non_lin_eqSelfCons(self,Xw,lw):        
        self.used_non_eq=[]
        self.min_used_non_eq=[]        
        min_ttf=[]
        ttf=[]
# estimation of phenotypes based on genetic equlibrium
        Xw_est=self.est_geno_from_alle(len(Xw.tolist()),lw,True)
        ph_eq=self.AnoSim*Xw_est
        strtext=""
        strtextEQ=""
        used_gen=[]
        for ix in range(len(Xw.tolist())):
            strloc= self.gen_to_string(self.genotypes[ix])
            findused=False
            for iss in range(len(self.symetry)) :
                if findused :
                    break
                s0=self.symetry[iss][0]
                s1=self.symetry[iss][1]
                if s0 == ix or s1 == ix :
                    for iused in used_gen :
                        if iused == s1 or iused == s0:
                            findused=True
                            break 
            if findused :
                continue
#                pass
            refalel=0
            strtext=""
            strtextEQ=""
            for iph_ in range(len(ph_eq)) :
                if self.AnoSim.tolist()[iph_][ix] > 0 : 
                    strtext="("+str(ph_eq[iph_])+")*("+str(Xw[ix])+")-(("+str(self.BnoSim.tolist()[iph_][0])+")*("+str(Xw_est[ix])+"))"
                    strtextEQ = "Pheno_Eq*"+self.aplot_genot_eqD[strloc][2]+"-[Pheno*"+strloc+"]"
                    tempexpr=parse_expr(strtext)
                    if tempexpr.is_number :
                        self.used_non_eq.append(strtextEQ)            
                        ttf.append(tempexpr)
                    else :
                        used_gen.append(ix)
                        min_ttf.append(tempexpr)                        
                        self.min_used_non_eq.append(strtextEQ)
                        self.used_non_eq.append(strtextEQ)            
                        ttf.append(tempexpr)
        return min_ttf,ttf   


    def set_non_lin_eq_Matrix(self,Xw,lw):        
        non_eq_matrix=[]
        strtext=""
        for i_ in range(len(lw)) :
            if self.L.tolist()[i_][ix] > 0 :
                if self.L.tolist()[i_][ix] > 0.6 :
                    if strtext=="":
                        strtext=strtext+"("+str(lw[i_])+")*("+str(lw[i_])+")"
                    else :    
                        strtext=strtext+"*("+str(lw[i_])+")*("+str(lw[i_])+")"
                else :
                    if strtext=="":
                        strtext=strtext+"("+str(lw[i_])+")"
                    else :    
                        strtext=strtext+"*("+str(lw[i_])+")"
#            strtext=strtext+"+("+str(-Xw[ix])+")"
            non_eq_matrix.append(parse_expr(strtext))
        return non_eq_matrix   

    def addedDofeq(self) :
        return self.nfreedom

    def used_sym(self):
        return self.symetryeq

    def used_nlin_eq(self):
        return self.used_non_eq

    def gen_and_aplot(self,B,SelfCons=True,Numerical=False,onlysymb=False) :
#    ''' Analitical allele frequency estimation based on phynotipic frequency '''         
        self.B = B
#   Lamda ij_kl ... genotype        
        Xw,X=self.solve()

        w = symbols('w:{0}_:{1}'.format(self.nrows, self.ncols), cls=Dummy)
        self.write_to_exfile("temp_w.py",w)

        if self.sympifyexpr :
            Xw = self.symplify_depen(Xw)
            self.res_S_Xw=self.symplify_results(Xw)
            self.res_S_Xw_Const=[]
            self.res_S_Xw_str=[]
            if len(self.res_S_Xw) > 0 :
                if type(self.res_S_Xw[0])==Matrix :
                    for iv in range(len(self.res_S_Xw)) :
                        self.res_S_Xw_str.append("")
                        self.res_S_Xw_Const.append(Xw.T*self.res_S_Xw[iv])
                        for igen in range(len(Xw)) :
                            strloc= self.gen_to_string(self.genotypes[igen])
                            self.res_S_Xw_str[iv]+="+"+str(round(self.res_S_Xw[iv][igen],3))+"*"+strloc
                else :
                    self.res_S_Xw_str.append("")
                    self.res_S_Xw_Const.append(Xw.T*self.res_S_Xw)
                    for igen in range(len(Xw)) :
                        strloc= self.gen_to_string(self.genotypes[igen])
                        self.res_S_Xw_str[0]+="+"+str(round(self.res_S_Xw[igen],3))+"*"+strloc


        lw=self.L*Xw
        l=self.L*X
        self.res_S_l=self.symplify_results(lw)
        self.res_S_l_Const=[]
        self.res_S_l_str=[]
        if len(self.res_S_l) > 0 :
            if type(self.res_S_l[0])==Matrix :
                for iv in range(len(self.res_S_l)) :
                    self.res_S_l_str.append("")
                    self.res_S_l_Const.append(lw.T*self.res_S_l[iv])
                    for igen in range(len(lw)) :                
                        strloc= self.aplotypeList[igen]
                        self.res_S_l_str[iv]+="+"+str(round(self.res_S_l[iv][igen],3))+"*["+strloc+"]"
            else :
                self.res_S_l_str.append("")
                self.res_S_l_Const.append(lw.T*self.res_S_l)
                for igen in range(len(lw)) :                
                    strloc= self.aplotypeList[igen]
                    self.res_S_l_str[0]+="+"+str(round(self.res_S_l[igen],3))+"*["+strloc+"]"

#    lamda i.... (allele frequensy in the population)
        n=lw.rows
#        one=Matrix.ones(1,n)
#        w = symbols('w:{0}_:{1}'.format(self.nrows, self.ncols), cls=Dummy)
#        self.write_to_exfile("temp_w.py",w)
#   only for 1 locus
#        f=[]
#        for i_ in range(len(lw)) :
#            for j_ in range(len(lw)) :
#                f.append(parse_expr(str(lw[i_]*lw[j_]-Xw[i_*(n-1)+j_

        self.nfreedom=self.A.cols-self.A.rows
        if self.nfreedom <=0 or onlysymb: 
            self.used_non_eq=[]
            self.used_non_eq.append("None Used # freedom = 0 ")
# solution comes only from the liniear part 
            S_lamdas=[]
            S_lamdasAN=[]
            S_lamdasw=[]
            S_lamdas_gen=[]
            S_lamdasAN_gen=[]
            S_lamdasw_gen=[]

            if onlysymb or ((sum(Xw)-1.)**2<0.000000001 and self.allpos(Xw)) :             
                S_lamdas.append(self.L*Matrix(Xw))
                S_lamdasw.append(self.L*Matrix(Xw))            

                S_lamdas_gen.append(Xw)
                S_lamdasw_gen.append(Xw)            
            else :
                S_lamdasAN.append(self.L*Matrix(Xw))
                S_lamdasAN_gen.append(Xw)                
        else :    
# set up set of nonlinear equations
            if SelfCons :
                f,fall=self.set_non_lin_eqSelfCons(Xw,lw)
            else :             
                f,fall=self.set_non_lin_eq(Xw,lw)

    #        testeq=self.set_non_lin_eq_Matrix(Xw,lw)
    # do not tray to solve all
            fs=tuple(f)
            Sep_solutions=[]
    #        solutions=solve(fs)
    #        if len(solutions)==0 :
    #            for ifs in fs :
    #                Sep_solutions.append(solve(ifs))
    #            solutions=Sep_solutions[0]
    #        for i in range(len(fs)-1) :
    #       tray only one
            isol=-1
            for i in [0] :
#                Sep_solutions.append(solve(fs[i:i+self.nfreedom]))
                if Numerical :

                    maxvar=0
                    varble=[]
                    for ieq in range(len(fs)):
                        maxvar=max(maxvar,len(fs[ieq].free_symbols))
                        temv=list(fs[ieq].free_symbols)
                        varble=list(set().union(varble,temv))
                    inval=[]
                    for ival in range(maxvar):
                        inval.append(0.)
                    if len(varble)>len(fs) :
                       for ieq in range(len(varble)-len(fs)) :
                           f.append(parse_expr("0"))
                       fs=tuple(f) 
                    neq=min(max(self.nfreedom,len(varble)),len(fs))
#                    neq=self.nfreedom
                    isol = isol+1                       
                    Sep_solutions.append(nsolve(fs[i:i+neq],varble,inval))
                    eqsol=[]
                    for ifs in fs:
                        tempxw=parse_expr(str(ifs))
                        ivar_i=-1
                        for ivar in varble :
                            ivar_i = ivar_i +1
                            tempxw=tempxw.subs(ivar,Sep_solutions[isol][ivar_i]) 
                        eqsol.append(tempxw)
#
                    S_lamdas=[]
                    S_lamdasAN=[]
                    S_lamdasw=[]
                    S_lamdas_gen=[]
                    S_lamdasAN_gen=[]
                    S_lamdasw_gen=[]

                    lamda=[]
                    lamdaw=[]
                    Xw_num=[]                    
                    for iXw in Xw:
                        tempxw=parse_expr(str(iXw))
                        try :
                            ivar_i=-1
                            for ivar in varble :
                                ivar_i = ivar_i +1
                                tempxw=tempxw.subs(ivar,Sep_solutions[isol][ivar_i]) 
                            lamda.append(tempxw)
#                            lamda.append(parse_expr(str(iXw)).subs(sub))
                            lamdaw.append(parse_expr(str(iXw)))
                        except :
                            pass
                    if (sum(lamda)-1.)**2<0.000000001 and self.allpos(lamda) :             
                        S_lamdas.append(self.L*Matrix(lamda))
                        S_lamdasw.append(self.L*Matrix(lamdaw))            

                        S_lamdas_gen.append(lamda)
                        S_lamdasw_gen.append(lamdaw)            
                    else :
                        S_lamdasAN.append(self.L*Matrix(lamda))
                        S_lamdasAN_gen.append(lamda)                
                    self.solved=True
                    return S_lamdas,S_lamdasw,S_lamdasAN,S_lamdas_gen,S_lamdasw_gen,S_lamdasAN_gen

#                    Sep_solutions.append(nsolve(fs,varble,inval))
                else :
                    # maxvar=0
                    # varble=[]
                    # for ieq in range(len(fs)):
                    #     maxvar=max(maxvar,len(fs[ieq].free_symbols))
                    #     temv=list(fs[ieq].free_symbols)
                    #     varble=list(set().union(varble,temv))
                    # inval=[]
                    # for ival in range(maxvar):
                    #     inval.append(0.)
                    # if len(varble)>len(fs) :
                    #    for ieq in range(len(varble)-len(fs)) :
                    #        f.append(parse_expr("0"))
                    #    fs=tuple(f) 
                    # neq=min(max(self.nfreedom,len(varble)),len(fs))                       
#                    neq=self.nfreedom
                    if self.sympifyexpr :
                        Sep_solutions.append(solve(fs))
                    else :
                        Sep_solutions.append(solve(fs[i:i+self.nfreedom]))                    
 
            solutions=Sep_solutions[0]

            global wlist
            S_lamdas=[]
            S_lamdasAN=[]
            S_lamdasw=[]
            S_lamdas_gen=[]
            S_lamdasAN_gen=[]
            S_lamdasw_gen=[]

            for isol in solutions :

                lamda=[]
                lamdaw=[]
                sub=[]
                for ilist in wlist :
                    if type(isol) is dict : 
                        if ilist in isol.keys():
                            iso=isol[ilist]
                            sub.append((ilist,iso))
                    else :
                        if ilist == isol:
                            iso=isol
                            sub.append((ilist,iso))
                for iXw in Xw:
                        try :
                            lamda.append(parse_expr(str(iXw)).subs(sub))
                            lamdaw.append(parse_expr(str(iXw)))
                        except :
                            pass

                eqsol=[]
                for ifs in fs:
                    tempxw=parse_expr(str(ifs))
                    tempxw=tempxw.subs(sub) 
                    eqsol.append(tempxw)

                if (sum(lamda)-1.)**2<0.000000001 and self.allpos(lamda) :             
                    S_lamdas.append(self.L*Matrix(lamda))
                    S_lamdasw.append(self.L*Matrix(lamdaw))            

                    S_lamdas_gen.append(lamda)
                    S_lamdasw_gen.append(lamdaw)            
                else :
                    S_lamdasAN.append(self.L*Matrix(lamda))
                    S_lamdasAN_gen.append(lamda)                
        self.solved=True
        return S_lamdas,S_lamdasw,S_lamdasAN,S_lamdas_gen,S_lamdasw_gen,S_lamdasAN_gen

    def genotypepop(self,B,SelfCons=True) :
#    ''' Analitical genotype frequency estimation based on phynotipic frequency '''
        self.B = B
        Xw,X=self.solve()
        lw=self.L*Xw
        l=self.L*X
        n=lw.rows
        one=Matrix.ones(1,n)
        w = symbols('w:{0}_:{1}'.format(self.nrows, self.ncols), cls=Dummy)
        self.write_to_exfile("temp_w.py",w)
#        f=[]
#        for iS in self.Sim :
#            i_,j_=iS
#            f.append(parse_expr(str(Xw[i_]-Xw[j_])))
#        for i_ in range(len(lw)) :
#            for j_ in range(len(lw)) :
#                f.append(parse_expr(str(lw[i_]*lw[j_]-Xw[i_*(n-1)+j_])))  
        if SelfCons :
            f,fall=self.set_non_lin_eqSelfCons(Xw,lw)
        else :              
            f,fall=self.set_non_lin_eq(Xw,lw)                  
#        fs=str(lw[0]*lw[0]-Xw[0])
#        f=parse_expr(fs)
#        print f
        fs=tuple(f)
        Sep_solutions=[]
        solutions=solve(fs)
        self.nfreedom=self.A.cols-self.A.rows
        if len(solutions)==0 :
#            for ifs in fs :
#                Sep_solutions.append(solve(ifs))
#            solutions=Sep_solutions[0]
            for i in range(len(fs)-1) :
                Sep_solutions.append(solve(fs[i:i+self.nfreedom]))
            solutions=Sep_solutions[0]

        global wlist

        S_lamdas=[]
        S_lamdasAN=[]
        S_lamdasw=[]
        for isol in solutions :
            lamda=[]
            lamdaw=[]
            sub=[]
            for ilist in wlist :
                if ilist in isol.keys():
                    iso=isol[ilist]
                    sub.append((ilist,iso))
            for iXw in Xw:
                    try :
                        lamda.append(parse_expr(str(iXw)).subs(sub))
                        lamdaw.append(parse_expr(str(iXw)))
                    except :
                        pass
            if (sum(lamda)-1.)**2<0.000000001 and self.allpos(lamda) :             
                S_lamdas.append(lamda)
                S_lamdasw.append(lamdaw)            
            else :
                S_lamdasAN.append(lamda)
        self.solved=True                
        return S_lamdas,S_lamdasw,S_lamdasAN

    def setupgen(self,alleletypes) : 
        self.aplotypes=[]
        for iloc in  alleletypes.split(";") :
#        if len(iloc.split(",")!=2) :
#            print "len(iloc.split(",")!=2" 
#            exit(1)
            self.aplotypes.append(iloc.split(","))
        self.aplotypes_dic={}
        iapl_int=-1
        self.aplotypeList=[]
        for iloc in self.aplotypes :
            for iapl in iloc :
                iapl_int=iapl_int+1
                self.aplotypes_dic[iapl]=iapl_int
                self.aplotypeList.append(iapl)
        

#    aplotypes = [[["A"],["a"]],[["R"],["r"]]]
        nloci = len(self.aplotypes)
        l_ij_dic={}
        ij=-1
        self.l_ij_kl=[]
        for li in range(nloci):
            charij=[]
            i_aplot=len(self.aplotypes[li])
            for i in range(i_aplot) :
                for j in range(i_aplot) :
                    charij.append(str(self.aplotypes[li][i])+","+str(self.aplotypes[li][j]))
            self.l_ij_kl.append(charij)

    def set_gent(self):
        icount=0
        self.genotypes=self.l_ij_kl[0]
        for second_list in self.l_ij_kl[1:] :
            self.genotypes = [list((a,b)) for a in self.genotypes for b in second_list]
            icount=icount+1
            if icount > 1 :
                for ilist in range(len(self.genotypes)) :
                    temp=[]
                    for iel in self.genotypes[ilist][0] :
                        temp.append(iel)
                    self.genotypes[ilist] = temp+self.genotypes[ilist][1:]    




#    def rec_gent(self,loclevel=0,alellevel=0,Cret=""):
#        if loclevel == 0 :
#            Cret=''
#        Cretold=Cret
#        for iloci in range(loclevel,len(self.l_ij_kl)) :           
#            for ialel in range(len(self.l_ij_kl[iloci])) :
#                Cret=Cretold
#                if loclevel == 0 :
#                    Cret=''
#                    Cret+=self.l_ij_kl[iloci][ialel]
#                else :
#                    Cret+=";"+self.l_ij_kl[iloci][ialel]
#                if loclevel == len(self.l_ij_kl)-1 :
#                    self.genotypes.append(Cret)
#                    Cret=''
#                    Cret=Cretold
#                else :    
#                    Cret=self.rec_gent(iloci+1,ialel,Cret)
#            else :
#                Cret=Cretold
#        return Cret
    ''' Estimate Symetries Assumimg either Gen-eq or loci in diferent chromosome  '''
    def Gen_symetries_assum(self):
        self.symetry=[]
        self.symetryeq=[]
        for intxy in range(len(self.genotypes)-1) :
            for intkl in range(intxy+1,len(self.genotypes)) :
                xy = self.genotypes[intxy]
                kl =self.genotypes[intkl]
                strxy = self.gen_to_string(xy)
                ifm=[]
                if type(xy)==list:
                    for ixy in xy :
                        ifm.append(ixy.split(","))
                else :
                    ifm.append(xy.split(","))
                if xy == kl : 
                    continue 
                strkl = self.gen_to_string(kl)
                jfm=[]
                if type(kl)==list:
                    for jkl in kl :
                        jfm.append(jkl.split(","))
                else :
                    jfm.append(kl.split(","))

                test = True
                for itest in range(len(ifm)):
                    test=test and  ( (ifm[itest][0]==jfm[itest][1] and ifm[itest][1]==jfm[itest][0]) or (ifm[itest][0]==jfm[itest][0] and ifm[itest][1]==jfm[itest][1]))  
                if test:   # diploid
                    self.symetry.append([self.genotypes_dic[strxy],self.genotypes_dic[strkl]])
                    self.symetryeq.append(strxy+" = "+strkl)
        self.redsymetries()
        return self.symetry

    ''' find redundant symetries '''
    def redsymetries(self):
        self.use_symetry=[]
        for s in self.symetry :
            self.use_symetry.append(True)

        iss=-1
        for s in self.symetry :
            iss=iss+1
            if self.use_symetry[iss]:
                itt = -1
                for t in self.symetry :
                    itt=itt+1
                    iww=-1
                    for w in self.symetry :
                        iww=iww+1
                        if s!=t and s!=w and t!=w :
                                if self.use_symetry[itt] and self.use_symetry[iww]:                            
                                    if s[0]==t[0] and ((s[1]==w[0] and t[1] ==w[1]) or (s[1]==w[1] and t[1] ==w[0]) ) :
                                        self.use_symetry[iww]=False     
                                    if s[1]==t[0] and ((s[0]==w[0] and t[1] ==w[1]) or (s[0]==w[1] and t[1] ==w[0]) ) :
                                        self.use_symetry[iww]=False     
                                    if s[0]==t[1] and ((s[1]==w[0] and t[0] ==w[1]) or (s[1]==w[1] and t[0] ==w[0]) ) :
                                        self.use_symetry[iww]=False     
                                    if s[1]==t[1] and ((s[0]==w[0] and t[0] ==w[1]) or (s[0]==w[1] and t[0] ==w[0]) ) :
                                        self.use_symetry[iww]=False     



    ''' Estimate standar Symetries with out Assumimg anything'''
    def Gen_symetries(self):
        self.symetry=[]
        self.symetryeq=[]
        for intxy in range(len(self.genotypes)-1) :
            for intkl in range(intxy+1,len(self.genotypes)) :
                xy = self.genotypes[intxy]
                kl =self.genotypes[intkl]
                strxy = self.gen_to_string(xy)
                ifm=[]
                if type(xy)==list:
                    for ixy in xy :
                        ifm.append(ixy.split(","))
                else :
                    ifm.append(xy.split(","))
                if xy == kl : 
                    continue 
                strkl = self.gen_to_string(kl)
                jfm=[]
                if type(kl)==list:
                    for jkl in kl :
                        jfm.append(jkl.split(","))
                else :
                    jfm.append(kl.split(","))

                test = True
                for itest in range(len(ifm)):
                    test=test and ifm[itest][0]==jfm[itest][1] and ifm[itest][1]==jfm[itest][0]  
                if test:   # diploid
                    self.symetry.append([self.genotypes_dic[strxy],self.genotypes_dic[strkl]])
                    self.symetryeq.append(strxy+" = "+strkl)
        self.redsymetries()
        return self.symetry

    def updateA(self,A) :
        self.A=A

    def updateB(self,B) :
        self.B=B


    def remove_simetries(self):
        self.A=self.AnoSim
        self.B=self.BnoSim
        return self.A,self.B

    ''' not used since it  creates Rank-deficient matrice that is not suported by sympy '''    
    def addnormeq_gen(self):
        ncols=self.A.cols
        addrow=[]
        for ic in range(ncols):
            addrow.append(1)
        nrows=self.A.rows
        self.A=self.A.row_insert(nrows,Matrix([addrow]))
        self.B=self.B.row_insert(nrows,Matrix([1]))

    ''' can be  used only to check since are 0 by constractions, add normalised equations in the solution of nonlinear problem '''    
    def addnorm_eqs(self,Xw,lw):
        self.norm_eq=[]
        self.norm_eq_STR=[]

        strtext=""
        strtextEQ=""
        for ix in range(len(Xw.tolist())):
            strloc= self.gen_to_string(self.genotypes[ix])
            strtext=strtext+"+("+str(+Xw[ix])+")"
            strtextEQ=strtextEQ+"+["+strloc+"]"                            
        strtext=strtext+"+("+str(-1)+")"
        strtextEQ=strtextEQ+"-1]"                            
        self.norm_eq.append(parse_expr(strtext))
        self.norm_eq_STR.append(strtextEQ)            
        strtext=""
        strtextEQ=""
        for i_ in range(len(lw)) :
            strtext=strtext+"+("+str(lw[i_])+")"
            strtextEQ=strtextEQ+"("+self.aplotypeList[i_]+")"
        strtext=strtext+"-1"
        strtextEQ=strtextEQ+"-1"                            
        self.norm_eq.append(parse_expr(strtext))
        self.norm_eq_STR.append(strtextEQ)            



    def addsymetries(self):
        self.AnoSim=self.A
        self.BnoSim=self.B
        ncols=self.A.cols
        simin =-1
        for isim in self.symetry:
            simin = simin+1
            if self.use_symetry[simin] :
#            if True :
                addrow=[]
                for ic in range(ncols):
                    addrow.append(0)
                addrow[isim[0]]=1
                addrow[isim[1]]=-1
                nrows=self.A.rows
                self.A=self.A.row_insert(nrows,Matrix([addrow]))
                self.B=self.B.row_insert(nrows,Matrix([0]))
#        self.addnormeq_gen()
#         creates Rank-deficient matrice that is not suported by sympy
        return self.A,self.B

    def nosymetries(self):
        if self.AnoSim != None :
            self.A=self.AnoSim
            self.B=self.BnoSim
            return self.A,self.B
        else :
            return self.A,self.B

    
if __name__ == "__main__":
#    genotypes.append("A,B")
    alleletypes="A,a;B,b"
    loci = len((alleletypes).split(";"))
#   A from gemotypes to phynotypes 
#    A = Matrix([[1,1,1,0], [0, 0, 0, 1], [0, 1, -1, 0]])
    A = Matrix([[1,1,1,0], [0, 0, 0, 1]])
#    A = Matrix([[1.,1.,0], [0, 0, 1.]])
#   from genotypes -> allele frequances       
#    L=Matrix([[2./2., 1./2., 1./2., 0], [0, 1./2, 1./2., 2./2.]])
#    L=Matrix([[2./2., 1./2., 0], [0, 1./2, 2./2.]])
#   imposed simetries
#    Sim=[(2,3)]
#   expetimenal values of the phenotypes plus symetries constrains 
#    Phenotype = Matrix([0.789473684,0.210526316,0])
    Phenotype = Matrix([0.789473684,0.210526316])

#   initialize the solver     
  
    Popsolver=phen_to_gen(alleletypes)
    Popsolver.initialize(A,Phenotype)
    print (Popsolver.genotypesD())
    symetries =Popsolver.Gen_symetries()
    A,Phenotype=Popsolver.addsymetries()        
#   solve for the messured phenotypes    
    al_pop,al_popw,sol_un=Popsolver.genotypepop(Phenotype)

    print (al_pop)
    print (al_popw)
    print (sol_un)

    al_pop,al_popw,sol_un=Popsolver.apop(Phenotype)

    print (al_pop)
    print (al_popw)
    print (sol_un)

