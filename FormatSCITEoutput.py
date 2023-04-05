import os
import sys
import Functions3

SIgv=sys.argv[1]
SIcellOr=sys.argv[2]
OutDir=sys.argv[3]
Python=sys.argv[4]
ID=SIgv.split(os.sep)[-1][:-3]
GVin=open(SIgv,'r').readlines()
out=''
for i in GVin:
   if i[:4]=='node' and i.find('color=lightgrey')!=-1: break
   else: out+=i
out+='}\n'

Functions3.GetOut(OutDir+os.sep+ID+'_dot.txt',out)

os.system(Python+' SCITE_Sum.py '+SIgv+' '+SIcellOr+' '+OutDir+os.sep+ID)   