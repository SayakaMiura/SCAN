
import Functions3
import sys
import os
import glob

GV=sys.argv[1]
CellOrderFile=sys.argv[2]
OutID=sys.argv[3]
if os.path.exists(GV)==True and os.path.exists(CellOrderFile)==True:
    print (GV)
   
    ID=GV[:-3]
    CellOrder=open(CellOrderFile,'r').readlines()
   

    CellNum=len(CellOrder)
    print (CellNum)
    ScID2Cell={}
    c=0
    while c<CellNum:
        ScID2Cell[str(c)]=CellOrder[c].replace('#','')
        c+=1	
    
    Clones=open(GV,'r').readlines()
    Cell2Clo={}
    Clo2Cell={}
    CloLs=[]
    Dec2Anc={}
    Anc2Dec={}
    out='Clone\tCell\n'
    outSeq='#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'
    outCloSeq='#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'
    MutPosi=[]
    AkkMutPos=[]	
    for i in Clones:
      if i.find('-> s')!=-1:
       i=i.strip().split('->')
       Clo='C'+i[0].strip()
       Cell=ScID2Cell[i[1].strip()[:-1][1:]]#CellOrder[int(i[1].strip()[:-1][1:])]
    
       if Clo not in Clo2Cell:Clo2Cell[Clo]=[]
       Clo2Cell[Clo]+=[Cell]   
       Code=Cell in Cell2Clo
       if Code!=True: Cell2Clo[Cell]=[]
       Cell2Clo[Cell].append(Clo)
       CloLs.append(Clo)
      elif i.find('-> ')!=-1:
       i=i.strip().split('->')
       Dec='C'+i[1].strip().replace(';','')
       Anc='C'+i[0].strip()
       AkkMutPos+=[int(i[1].replace(';','')),int(i[0])]	   
      
       Dec2Anc[Dec]=Anc
       if Anc not in Anc2Dec: Anc2Dec[Anc]=[]
       Anc2Dec[Anc].append(Dec)
    NodeLs=[]
    for Anc in Anc2Dec:
        if len(Anc2Dec[Anc])>1: NodeLs.append(Anc)
    	
    
    CloLs=list(set(CloLs))
    NodeLs+=CloLs
    AkkMutPos=list(set(AkkMutPos))
    AkkMutPos.sort()
    SNVnum=AkkMutPos[-2]
    print ('SNVnum',SNVnum)	
 
    for Clo in CloLs:
         MutLs=[Clo[1:]]
         Code=Clo in Dec2Anc
         Node1=Clo
         while Code==True: 
             Node1=Dec2Anc[Node1]	 
             MutLs.append(Node1[1:])	 
             Code=Node1 in Dec2Anc
      
         c=0
         
         while c<SNVnum:
        
             Code=str(c+1) in MutLs
             if Code==True: 
                
                 MutPosi.append(c)			 
       
             c+=1
    for Clo in CloLs:
         MutLs=[Clo[1:]]
         Code=Clo in Dec2Anc
         Node1=Clo
         while Code==True: 
             Node1=Dec2Anc[Node1]	 
             MutLs.append(Node1[1:])	 
             Code=Node1 in Dec2Anc
    
         c=0
         Seq=''	 
         while c<SNVnum:
      	 
           if MutPosi.count(c)==0: Seq+='A'
           else:	   
             Code=str(c+1) in MutLs
             if Code==True: 
                 Seq+='T'
            		 
             else: Seq+='A'
           c+=1		 
         outCloSeq+='#'+Clo+'\n'+Seq+'\n'		 
         CellLs=Clo2Cell[Clo]
    	 
         for Cell in CellLs:	 
        
            out+=';'.join(Cell2Clo[Cell])+'\t'+Cell+'\n'		
            if len(Cell2Clo[Cell])==1:		
                outSeq+='#'+Cell.strip()+'\n'+Seq+'\n'   
    outSeq+='#Normal\n'+('A'*SNVnum)+'\n'   
    Functions3.GetOut(OutID+'_Clone.meg',outCloSeq)
    Functions3.GetOut(OutID+'_Cell.meg',outSeq)	
    Functions3.GetOut(OutID+'_CloneAnno.txt',out)
else:
    print (GV,CellOrderFile,'do not exist')	
