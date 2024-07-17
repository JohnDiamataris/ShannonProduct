import random
import numpy as np
import scipy
import csv
from scipy.linalg import null_space

class Phenomenon:
    
    def __init__(self,GenType,PhenTypeL=None,Weights=None,SampleSize=10,PhenDict=None):
        self.GenType=GenType
        self.PhenList=self.MakePhenList()
        self.Weights=Weights 
        self.PhenDict=PhenDict
        self.PopSampleSize=SampleSize
        self.PhenTypelist=[]
        self.PopSampledict={}
        self.PopSampleLista=[]
    

    def MakePhenList(self):
        PhenLista=[]
        for genes in self.GenType:
            phenelement=''
            for gene in genes:
                q=gene.split(',')
                if gene.islower():
                    if (q[0] == q[1]):
                        phenelement= phenelement+gene[0]
                    else: phenelement = phenelement + q[0] +q[1]
                elif gene.isupper():
                    if (q[0] == q[1]):
                        phenelement =phenelement + gene[0]
                    else: phenelement = phenelement + q[0] +q[1]
                else: 
                    if q[0].isupper():
                        phenelement += gene[0]
                    else:
                        if q[1].isupper():
                            phenelement = phenelement + gene[0]
            
            if phenelement not in PhenLista and phenelement[::-1] not in PhenLista:
                PhenLista.append(phenelement)                
                
        return PhenLista
    
   
    def PopSampleListExact(self):
        size={}
        if not self.Weights:
            self.Weights={}
            for i in self.PhenList:
                self.Weights[i]=1/len(self.PhenList)
        for i,j in self.Weights.items():
            size[i]= round(j*self.PopSampleSize) 
        
        if sum(size.values()) < self.PopSampleSize:
            rnd=('',0)
            for l,m in self.Weights.items():
                if (m*self.PopSampleSize) % 1-0.5 <= 0:
                    if m*self.PopSampleSize % 1 > rnd[1] % 1 :
                        rnd = (l,m*self.PopSampleSize )
            size[rnd[0]]=round(rnd[1]) + 1

        elif sum(size.values()) > self.PopSampleSize:
            rnd=('',0.9999)
            for l,m in self.Weights.items():
                if (m*self.PopSampleSize) % 1-0.5 >= 0:
                    if m*self.PopSampleSize % 1 < rnd[1] % 1 :
                        rnd = (l,m*self.PopSampleSize )
            size[rnd[0]]=round(rnd[1]) - 1
    
        pe={}
        for i in self.PhenDict.keys():
            
            choicesize=size[i]
            pe[i]=np.random.choice(np.arange(len(self.PhenDict[i])),choicesize)    
            pe[i]=[self.PhenDict[i][x] for x in pe[i]]

        # with open ('PopSample.csv','w') as PsF:
        #     writer=csv.writer(PsF)
        #     writer.writerow(PopSampleList) 
        self.PopSampledict=pe
        return pe

    def MakePhenTypeDict(self):
        if not self.PhenList==None:
            PhenTypeDict={}
            for ph in self.PhenList:
                if not ph in PhenTypeDict.keys():

                    PhenTypeDict[ph]=[]
            for g in self.GenType:
                G2Ph=''
                for k in g:
                    singlek=k.split(',')
                    if k.isupper():
                        if singlek[0] == singlek[1]:
                            G2Ph += singlek[0]
                        else:
                            G2Ph += singlek[0] + singlek[1]
                    elif k.islower():
                        if singlek[0] == singlek[1]:
                            G2Ph += singlek[0]
                        else:
                            G2Ph += singlek[0] + singlek[1]
                    else:
                        for i in singlek:
                            if i.isupper():
                                G2Ph += i
                if  G2Ph in PhenTypeDict.keys():
                    PhenTypeDict[G2Ph].append(g)
                else:
                    if G2Ph[::-1] in PhenTypeDict.keys():
                        PhenTypeDict[G2Ph[::-1]].append(g)
        
        self.PhenDict=PhenTypeDict
        return PhenTypeDict

        

    def PhenType(self):
        SampleGenPhenList=[]
        Sampldict=self.PopSampleListExact()
        SamplList=[]
        for i in Sampldict:
            SamplList=SamplList+Sampldict[i]
        self.PopSampleLista=SamplList
        if self.PhenDict==None:
            self.PhenDict=self.MakePhenTypeDict()
        PhenTypeList=[]
        PhDictTest={}
        for gen in SamplList:
            for Phen in self.PhenDict.items():
                if gen in Phen[1]:
                    SampleGenPhenList.append([gen,Phen[0]])
                    PhenTypeList.append(Phen[0])
                    if Phen[0] in PhDictTest.keys():
                        PhDictTest[Phen[0]].append(gen)
                    else:
                        PhDictTest[Phen[0]]=[gen]
                self.PhenTypelist=PhenTypeList
        return PhenTypeList,PhDictTest,SampleGenPhenList

    def PhenProb(self,prob=False):
        if self.PhenDict==None:
            self.PhenDict=self.MakePhenTypeDict()
        PhenTypeL=self.PhenType()
        PhenProb={}    
        for i in self.PhenDict.keys():
            c=PhenTypeL[0].count(i)
            if prob:
                PhenProb[i]=c/float(len(PhenTypeL[0]))
            else:
                PhenProb[i]=c
        # with open ('PhenProb.csv','w') as PhF:
        #     writer=csv.writer(PhF)
        #     for key, value in PhenProb.items(): 
        #         writer.writerow([key, value]) 
             
        return PhenProb,PhenTypeL[1] 

if __name__ == "__main__":
    GenTypeA=[['A,A'],['A,a'],['a,A'],['a,a']]
    PhenDictyA={'A':[['A,A'],['A,a'],['a,A']],'a':[['a,a']]}
    #GenTypeAB=[['A,A','B,B'],['A,A','B,b'],['A,A','b,B'],['A,A','b,b'],['A,a','B,B'],['A,a','B,b'],['A,a','b,B'],['A,a','b,b'],['a,A','B,B'],['a,A','B,b'],['a,A','b,B'],['a,A','b,b'],['a,a','B,B'],['a,a','B,b'],['a,a','b,B'],['a,a','b,b']]#from Popsolver.genotypes
    #PhenDicty={'AB':[['A,A','B,B'],['A,A','B,b'],['A,A','b,B'],['A,a','B,B'],['A,a','B,b'],['A,a','b,B'],['a,A','B,B'],['a,A','B,b'],['a,A','b,B']],'Ab':[['A,A','b,b'],['A,a','b,b'],['a,A','b,b']],'aB':[['a,a','B,B'],['a,a','B,b'],['a,a','b,B']],'ab':[['a,a','b,b']]}
    #GenTBomb = [['A,A','H,H'],['A,A','H,h'],['A,A','h,H'],['A,A','h,h'],['A,i','H,H'],['A,i','H,h'],['A,i','h,H'],['A,i','h,h'],['B,B','H,H'],['B,B','H,h'],['B,B','h,H'],['B,B','h,h'],['B,i','H,H'],['B,i','H,h'],['B,i','h,H'],['B,i','h,h'],['A,B','H,H'],['A,B','H,h'],['A,B','h,H'],['A,B','h,h'],['i,i','H,H'],['i,i','H,h'],['i,i','h,H'],['i,i','h,h'],]
    #PhenDicty={'A':[['A,A','H,H'],['A,A','H,h'],['A,A','h,H'],['A,i','H,H'],['A,i','H,h'],['A,i','h,H']],'B':[['B,B','H,H'],['B,B','H,h'],['B,B','h,H'],['B,i','H,H'],['B,i','H,h'],['B,i','h,H']],'AB':[['A,B','H,H'],['A,B','H,h'],['A,B','h,H']],'i':[['i,i','H,H'],['i,i','H,h'],['i,i','h,H']],'h':[['i,i','h,h'],['A,B','h,h'],['B,i','h,h'],['B,B','h,h'],['A,A','h,h'],['A,i','h,h']]}
    #PhenTBomb =['AH','Ah','BH','Bh','ABH','ABh','iH','ih']
    #Weights=[0.057,0.00025,0.00025,0.0001,0.316,0.0015,0.0015,0.0002,0.008,0.00045,0.00045,0.0001,0.119,0.001,0.001,0.0002,0.044,0.0005,0.0005,0.0001,0.435,0.0025,0.0025,0.0079]
    #PhenTypeA=['A','a']
    #PhenTypeAB=['AB','Ab','aB','ab']
    a=Phenomenon(GenTypeA,PhenDict=PhenDictyA)
    b=a.PhenProb()
    while not all(i>0 for i in b[0].values()):
         b=a.PhenProb()
    with open ('TestPhDict.csv','w') as PhT:
            writer=csv.writer(PhT)
            for key, value in b[1].items(): 
                writer.writerow([key, value])
    with open ('PhenProb.csv','w') as PhF:
            writer=csv.writer(PhF)
            for key, value in b[0].items(): 
                writer.writerow([key, value])
pass   