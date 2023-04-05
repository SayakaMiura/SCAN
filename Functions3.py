import os
from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
import numpy as np
import glob

def root_tree(OriNwk,Root):
    Out=OriNwk[:-4]+'_rooted.nwk'
    trees = list(Phylo.parse(OriNwk, 'newick'))
    for tree in trees:
        tree = tree.root_with_outgroup({'name': Root})
    Phylo.write(trees, Out, "newick")
def SumMutTree(File):
    File=open(File,'r').readlines()[1:]
    Dec2Anc={}
    Dec2Mut={}
    for i in File:
        i=i.strip().split('\t')
        Dec2Anc[i[0]]=i[1]
        Dec2Mut[i[0]]=i[2].split(';')
    Edge2Type={}
    for Dec in Dec2Anc:
        Anc=Dec2Anc[Dec]
        DecM=Dec2Mut.get(Dec,'?')[0]
        if DecM=='' and Anc=='Root': 
           DecM='Root'
           print (Dec,Anc,DecM)
           	   
        if DecM=='?':
            print (Dec)
            open('A','r').readlines()
        if Anc=='Root': AncM='Root'			
        else: AncM=Dec2Mut.get(Anc,'N')[0]
        if AncM=='' and Dec2Anc[Anc]=='Root':
           AncM='Root'		
          		   
        if AncM!=DecM:Edge2Type[AncM+'->'+DecM]='Order'
        if len(Dec2Mut[Dec])>1:
           Clu=Dec2Mut[Dec][1:]
           for i in Clu:
               if i!=DecM:Edge2Type[DecM+'->'+i]='Clu'
    return Edge2Type 
def MutTree2dot(File):
   Out=File[:-4]+'.gv'

   File=open(File,'r').readlines()[1:]
   Nodes=''
   Paths=''
   for i in File:
       i=i.split('\t')
       
       if len(i)>5: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4]+'\\nCellC:'+i[5].strip()+'\"]\n'	
       else: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4].strip()+'\"]\n'	   
       Paths+=i[1]+'->'+i[0]+'\n'	   

   out='digraph D {\n'+Nodes+Paths+'}\n'
   GetOut(Out,out)	
def Clean(Target):
	filelist = glob.glob(Target)
	for f in filelist:
		os.remove(f)


def GetHead(Head):
    Head=Head.strip().split('\t')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col

def ReadFasSeq(Fas):
 ID2Seq={}
 Fas=open(Fas,'r').readlines()
 Seq='start'
 for i in Fas:
  
   if i[0]=='>':
       if Seq!='start': ID2Seq[ID]=Seq
       ID=i[1:].strip()
       Seq=''

   else: Seq+=i.strip()
 ID2Seq[ID]=Seq
 return ID2Seq

def ReadMegSeq(Meg): #input is mega alignment file. out is name2seq dictionary and mega head
  Meg=open(Meg,'r').readlines()
  Read='s'
  out2=''
  NameOrder=[]
  Name2Seq={}
  SeqNum=0
  for i in Meg:
    if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
      	
        Read='n'
        Name=i.strip()
        NameOrder.append(Name)
        Name2Seq[Name]=''
        SeqNum+=1
    elif Read=='n': Name2Seq[Name]+=i.strip()
    elif Read=='s': out2+=i
  return NameOrder, Name2Seq
def ReadGV(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	
   
    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph D {','')	
        if i.find('->')!=-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge=i0[0].strip()+'->'+i0[1].split('[')[0].strip()			 
             Edge2In[Edge]=i.strip()
             if i.find('];')==-1: Read='E'
        elif Read=='E': 
             Edge2In[Edge]+=i.strip()
             if i.find('];')!=-1: Read='N'			 
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split('[')
           			 
             Node=i0[0].strip()
             MutIn=GetMutIn(i0[1])			 
     
             NodeMut[Node]=MutIn	
             Node2In[Node]=i
             if i.find('];')==-1: Read='NM'
        elif Read=='NM': 
             Node2In[Node]+=i
             NodeMut[Node]+=GetMutIn(i)		
         		 
             if i.find('];')!=-1: Read='N'				 

	 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In
def GetCla(Node,Dec2Anc):
   Cla=[]
   for D in Dec2Anc:
       In=[D]
       A=Dec2Anc[D]
       if A==Node: Find='y'
       else:	   
          Find='n'
          In.append(A)	
       	  
          while A in Dec2Anc:
                A=Dec2Anc[A]
                if A!=Node and Find=='n': In.append(A)
                elif A==Node: Find='y'
         			
       if Find=='y': Cla+=In	
   Cla+=[Node]		  
   return Cla 		  
def GetMutIn(i0):			 
             if i0.find('Mut:')!=-1:			 
                 MutLs=i0.split('Mut:')[-1].split('\"')[0]
                 if MutLs.find(',')!=-1: MutLs=MutLs.split(',')
                 else: MutLs=MutLs.split(';')
             else: MutLs=[]				 
             MutIn=[]
             for M in MutLs:
                 M=M.strip().split('\\n')			   
                 MutIn+=M	
             return MutIn				 
def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif

def UpMeg(Name2Seq,NameLs,Out):
    out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
    for Name in NameLs:
        if Name[0]!='#': Name='#'+Name
        if Name not in Name2Seq: Seq=Name2Seq[Name.replace('#','')]
        else: Seq=Name2Seq[Name]		
        out+=Name+'\n'+Seq+'\n'
    OutF=open(Out,'w')
    OutF.write(out)
    OutF.close()	
def prune_tree(OriNwk,ExtraLs):
   TreeLs=[open(OriNwk,'r').readlines()[0]]
   Cou=1 
   for Tst in TreeLs:
    tree=Phylo.read(StringIO(Tst.strip()), "newick")   
  
    for tip in ExtraLs:	
         tip=tip.replace('#','')
         if Tst.find(','+tip+':')!=-1 or Tst.find('('+tip+':')!=-1:	
     		 
            tree.prune(tip)
    Phylo.write(tree, OriNwk[:-4]+'_prune.nwk','newick')
    Cou+=1	
def BooSum(OriTreFile,OutFname,Rpath,cwd,Pattern):
 
    cwd=cwd.replace('\\','/')
    Rout=''
    Rout+='library(ape)\n'
    Rout+='true_tree <- ape::read.tree(\"'+OriTreFile+'\")\n'
    Rout+='lf <- list.files(\"'+cwd+'/\", pattern = \"'+Pattern+'\", full.names = TRUE)\n'
    Rout+='x<- ape::rmtree(length(lf), true_tree$Nnode)  ## length(lf) = num of replicate trees\n'
    Rout+='for(j in 1:length(lf)){\n'
    Rout+='  x[j] <- list(ape::read.tree(lf[j]))\n'
    Rout+='}\n'
    Rout+='b <- phangorn::plotBS(true_tree, x, p =10,  \"phylogram\")\n'
    Rout+='true_boot_support <- as.numeric(b$node.label)\n'
    Rout+='ape::write.tree(b, file = \"'+OutFname+'\")\n'
    GetOut('RunR.r',Rout)
    os.system(Rpath+' RunR.r')
  	
def SumClo(Cell2Seq,OutFasName,CloID,Cut):
    Seq2CellLs={}
    for Cell in Cell2Seq:
        Seq=Cell2Seq[Cell]
        Seq2CellLs[Seq]=Seq2CellLs.get(Seq,[])+[Cell]		
    Clo2Seq={}
    Clo2CellLs={}  
    ID=1
    out=''	
    outIn='Clone\tCellLs\n'	
    for Seq in Seq2CellLs:
        Clo=CloID+str(ID)

        if len(Seq2CellLs[Seq])>=Cut:		
            out+='>'+Clo+'\n'+Seq+'\n'
            Clo2Seq[Clo]=Seq
            Clo2CellLs[Clo]=Seq2CellLs[Seq]			
            outIn+=Clo+'\t'+','.join(Seq2CellLs[Seq])+'\n'		
        ID+=1		
    GetOut(OutFasName,out)
    GetOut(OutFasName[:-6]+'_CloneAnno.txt',outIn)	
    return Clo2Seq,Clo2CellLs	
def PairClone_by_cellLs(Ori2CellLs,Inf2CellLs,IClo2Seq,OutFasName):
    
    Miss=[]
    Extra=[]	
    Old2New={}	
    outInf='ReNameID\tOldID\n'	
    for Ori in Ori2CellLs:
       OCellLs=Ori2CellLs[Ori]
       Best=0
       BestInf=''
       for Inf in Inf2CellLs:
           ICellLs=Inf2CellLs[Inf]
           Count=0
           for I in ICellLs:
               if OCellLs.count(I)!=0: Count+=1
           if Count>Best:
               Best=Count
               BestInf=Inf	
       	   
       if Best==0: 
           Miss.append(Ori)
           outInf+=Ori+'\tMiss\n'		   
       else: 
           Old2New[BestInf]=Old2New.get(BestInf,[])+[Ori]	
           outInf+=Ori+'\t'+BestInf+'\t'+str(Best)+'\t'+','.join(Inf2CellLs[BestInf])+'\n'
       
    out=''
	
    Clo2Seq={}	
    for Inf in Inf2CellLs:
        if Inf not in Old2New:
            Extra.append(Inf)
            out+='>'+Inf+'\n'+IClo2Seq[Inf]+'\n'
            Clo2Seq[Inf]=IClo2Seq[Inf]	
            outInf+='Extra\t'+Inf+'\n'			
        else:
          NewLs=Old2New[Inf]
          for New in NewLs:		  
            out+='>'+New+'\n'+IClo2Seq[Inf]+'\n'	
            Clo2Seq[New]=IClo2Seq[Inf]			
    GetOut(OutFasName,out)
    GetOut(OutFasName[:-4]+'_inf.txt',outInf)	
    return Extra,Miss,Clo2Seq	
def GetOut(File,In):
    OutF=open(File,'w')
    OutF.write(In)
    OutF.close()	
