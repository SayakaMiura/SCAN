import os
import sys
import Functions3
import glob
import shutil
from scipy import stats

SIgv=sys.argv[1]
SIcellOr=sys.argv[2]
OutDir=sys.argv[3]
InFas=sys.argv[4]
CoiCut='0.4'#'0.2'#'0.4'#sys.argv[5]
ID=SIgv.split(os.sep)[-1][:-3]
Python='python'	

def GetFNR(OriFas,PreFas):
    if OriFas[-4:]=='.meg':
        OCellLs,OCell2Seq=Functions3.ReadMegSeq(OriFas)
    else: OCell2Seq=Functions3.ReadFasSeq(OriFas)
    PCellLs,PCell2Seq=Functions3.ReadMegSeq(PreFas)	
    Tc=0	
    FNc=0	
    Ac=0
    FPc=0	
    Cell2FNFP={}
    out='Cell\tType\tFNR\tFPR\tp_FNR\tp_FPR\t'	
    outIn=''	
    AllFP=[]
    AllFN=[]	
    for Ocell in OCell2Seq:
	   
        if '#'+Ocell in PCell2Seq:
            Oseq=OCell2Seq[Ocell]
            Pseq=PCell2Seq['#'+Ocell]
            Len=len(Oseq)
            c=0
            CellTc=0	
            CellFNc=0	
            CellAc=0
            CellFPc=0			
            while c<Len:
               if Pseq[c]=='T': 
                    CellTc+=1
                    if Oseq[c]=='A': CellFNc+=1
               if Pseq[c]=='A': 
                    CellAc+=1
                    if Oseq[c]=='T': CellFPc+=1					
               c+=1	
            Tc+=CellTc	
            FNc+=CellFNc	
            Ac+=CellAc
            FPc+=CellFPc
                  
            if CellTc==0:
                print ('normal predicted cell',Ocell,CellTc,CellAc)     
                Cell2FNFP[Ocell]={'FN':0,'FP':1.0*CellFPc/CellAc}	                
            else: 
                Cell2FNFP[Ocell]={'FN':1.0*CellFNc/CellTc,'FP':1.0*CellFPc/CellAc}	
                AllFP.append(1.0*CellFPc/CellAc)
                AllFN.append(1.0*CellFNc/CellTc)			
       	else: outIn+=Ocell+'\tDoublet\tNA\tNA\tNA\tNA\n'		
    Rate=1.0*FNc/Tc	
    FPrate=1.0*FPc/Ac	
    
    out+='FNR:'+str(Rate)+' FPR:'+str(FPrate)+'\n'+outIn
    for Cell in Cell2FNFP:
        AnnoFN=Test(Cell2FNFP[Cell]['FN'],AllFN)
        AnnoFP=Test(Cell2FNFP[Cell]['FP'],AllFP) 
    
        if AnnoFN.split('\t')[0]=='Singleton' and AnnoFP.split('\t')[0]=='Singleton': out+=Cell+'\tSingleton\t'
        else: out+=Cell+'\tDoublet\t'
        out+=str(Cell2FNFP[Cell]['FN'])+'\t'+str(Cell2FNFP[Cell]['FP'])+'\t'+AnnoFN.split('\t')[1]+'\t'+AnnoFP.split('\t')[1]+'\n'
        	
	
    return Rate,FPrate,out,Cell2FNFP	
def Test(Rate,Ls):
    Ave=1.0*sum(Ls)/len(Ls)
    if Ave>=Rate: return ('Singleton\tNA')
    else: 
       Res=stats.ttest_1samp(Ls, popmean=Rate)
       
       if Res[1]<0.05: return ("Doublet\t"+str(Res[1]))	
       else: return ("Singleton\t"+str(Res[1]))	   

if os.path.exists(OutDir)!=True: os.mkdir(OutDir)
print (OutDir)
os.system(Python+' FormatSCITEoutput.py '+SIgv+' ' +SIcellOr+' '+OutDir+' '+Python)	
if InFas[-4:]=='.meg':
    CellLs,Cell2Seq=Functions3.ReadMegSeq(InFas)
    out=''
    for i in CellLs:
        out+='>'+i.replace('#','')+'\n'+Cell2Seq[i]+'\n'
    Functions3.GetOut(InFas[:-4]+'.fasta',out)
    InFas=InFas[:-4]+'.fasta'
FNR,FPR,DoubLetoutIn,Cell2FNFP=GetFNR(InFas,OutDir+os.sep+ID+'_Cell.meg')
Functions3.GetOut(OutDir+os.sep+'Doublet_Detection_Result.txt',DoubLetoutIn)
print ('FNR',FNR,'FPR',FPR,len(Cell2FNFP))
os.system(Python+' SCAN_Pipeline_FIX.py '+InFas+' '+OutDir+os.sep+ID+'_dot.txt '+OutDir+' false '+str(FNR))
shutil.copy('DD_FastMOA_Result'+os.sep+'co_counts_matrix.txt',OutDir+os.sep+'co_counts_matrix.txt')
shutil.copy('DD_FastMOA_Result'+os.sep+'COI_matrix.txt',OutDir+os.sep+'COI_matrix.txt')
MOAfi=glob.glob('DD_FastMOA_Result'+os.sep+'*')
for i in MOAfi: os.remove(i)
os.system(Python+' SCAN_output.py '+OutDir+' '+ID+' '+SIgv[:-3].split(os.sep)[-1]+' '+CoiCut)
os.remove(OutDir+os.sep+'Labels.txt')
os.remove(OutDir+os.sep+'mutationConcurrentSequential_Result.txt')
os.remove(OutDir+os.sep+'ExpectedSequences.fasta')
os.remove(OutDir+os.sep+'Doublet_Detection_Result.txt')
os.remove(OutDir+os.sep+'co_counts_matrix.txt')
os.remove(OutDir+os.sep+'COI_matrix.txt')
os.remove(OutDir+os.sep+'MutTree_Cell.meg')
os.remove(OutDir+os.sep+'MutTree_Clone.meg')
os.remove(OutDir+os.sep+'MutTree_CloneAnno.txt')
os.remove(OutDir+os.sep+'MutTree_dot.txt')
