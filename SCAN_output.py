import sys
import os
import Functions3

Dir=sys.argv[1]
ID=sys.argv[2]
SCITEoutID=sys.argv[3] #G7snv_100_2_BEAMin_SCITE_ml0
COIcut=float(sys.argv[4])
Dot=Dir+os.sep+'SCAN_output.dot'
CellAnno=Dir+os.sep+SCITEoutID+'_CloneAnno.txt'
CellAnnoOut=CellAnno[:-4]+'1.txt'
Meg=Dir+os.sep+SCITEoutID+'_Clone.meg'
MegOut=Meg[:-4]+'1.meg'
Dou=Dir+os.sep+'Doublet_Detection_Result.txt'

Dec2Anc,NodeMut,Node2In,Edge2In	=Functions3.ReadGV(Dot)

PruneNode=[]
for E in Edge2In:
    In=Edge2In[E]
    COI=In.split('COI:')[1].split('\\n')[0]
    COI=float(COI)
    if COI<COIcut: 
        AncDec=E.split('->')
        Node=AncDec[1].strip()
        Cla=Functions3.GetCla(Node,Dec2Anc)		
        PruneNode+=Cla
        PruneNode+=Cla
PruneNode=list(set(PruneNode))

PruneMutCoi=[]
for N in PruneNode:
   	PruneMutCoi+=NodeMut[N]
	
PruneMutCon=[]
for N in NodeMut:
    Mut=NodeMut[N]
    if len(Mut)>1:
       PruneMutCon+=Mut[:-1]
	
Dou=open(Dou,'r').readlines()[1:]
DouCell=[]
for i in Dou:
    i=i.split('\t')
    if i[1]=='Doublet': DouCell.append(i[0])
print (len(DouCell))
CellAnno=open(CellAnno,'r').readlines()
out=CellAnno[0].strip()+'\tCon\tChild of <'+str(COIcut)+'\tCellType\n'
CellAnno=CellAnno[1:]
NonDoubClo=[]
CloLs=[]
for i in CellAnno:
    i=i.strip()
    if i !='': 
        i1=i.split('\t')
        Muts=i1[0].replace('C','').split(';')		
        CloLs+=i1[0].split(';')

        Conin=[]
        Plain=[]		
        for M in Muts: 
           if PruneMutCoi.count(M)==0: Conin.append('N')
           else: Conin.append('Y')
           if PruneMutCon.count(M)==0: Plain.append('N')
           else: Plain.append('Y')	
        Cell=i1[1]
        if DouCell.count(Cell)==0: NonDoubClo.append(i1[0])		
        if DouCell.count(Cell)==0: D='Singlet'
        else: D='Doublet'		
        out+=i+'\t'+';'.join(Plain)+'\t'+';'.join(Conin)+'\t'+D+'\n'
Functions3.GetOut(CellAnnoOut,out)
CloLs,Clo2Seq=Functions3.ReadMegSeq(Meg)
out='#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'
NonDoubClo=list(set(NonDoubClo))
print (len(NonDoubClo))
for Clo in CloLs:
    NID= Clo.replace('#','').replace('C','')
    if PruneMutCoi.count(NID)==0 and PruneMutCon.count(NID)==0 and NonDoubClo.count(Clo.replace('#',''))!=0:
        out+='#C'+NID+'\n'+Clo2Seq[Clo]+'\n'
Functions3.GetOut(MegOut,out)		
		
	