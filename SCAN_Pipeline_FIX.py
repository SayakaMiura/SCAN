import os
import sys
import time
from Bio import SeqIO
from statistics import mean
from matplotlib.pyplot import get
import numpy as np
import pydot
import graphviz
from scipy.stats import fisher_exact
from scipy.stats.distributions import chi2

path = ''#os.path.join(os.path.expanduser('~'), os.getcwd())+os.sep
Python='python'
inputFasta = str(sys.argv[1])
inputDot  = str(sys.argv[2])
outputDir = str(sys.argv[3])
Remove_Doublet = str(sys.argv[4])
beta=float(sys.argv[5]) #FNR


fastaFilename =  inputFasta
gvMutTreeFilename =  inputDot


if os.path.exists(outputDir) != True:
   
    os.mkdir(outputDir)


def get_Seq_ID(fasta_object):
    get_listSEQ = []
    for fasta in fasta_object:
        name, sequence = fasta.id, str(fasta.seq)
        get_listSEQ.append(str(name)+':'+str(sequence))
       
    return get_listSEQ #, get_listID

def get_best_node(ind, new_pos):
    bn = new_pos[ind]  ## 1. CHANGE variable name for if using a different position list than 'new_pos'
    return bn

def get_match_count(seq1, seq2):
    match = []
    extra = []
    for i in range(len(seq1[1])):
        if seq1[1][i] == 'T' and seq2[1][i] == 'T':
            match.append(i)
            
        if seq1[1][i] == 'T' and seq2[1][i] == 'A':# or seq1[1][i] == '?' and seq2[1][i] == 'T':
            extra.append(i)

    if len(match) != 0:
        return len(match), match[-1], extra, seq2[0]
    else:
        return 0, 0, 0, 0

def get_variant_count(list_of_sequences, condition): # condition referes to missing data counted or not counted
    count = []
    for i in range(len(list_of_sequences[1])):
        if condition == True:
            count_of_variants = [seq[i] for seq in list_of_sequences if seq[i] != 'A'] #change to seq[i] != 'A'
            count.append(len(count_of_variants))
        if condition == False:
            count_of_variants = [seq[i] for seq in list_of_sequences if seq[i] == 'T']
            count.append(len(count_of_variants))
    return count

def get_sorted_position(list_of_postitions, dictionary1):
    sort_variant = []
    for cnt in list_of_postitions:
        sort_v = [ i+1 for i in range(len(list_of_postitions)) if dictionary1[i+1] == cnt]
        sort_variant.extend(sort_v)
    return list(dict.fromkeys(sort_variant))

def get_sorted_fasta(seqlist_type, sorted_pos):
    new_list = []
    seqlist = [seqlist_type[i].split(':')[1] for i in range(0, len(seqlist_type))]
    for i in range(len(seqlist)):
        srt = [seqlist[i][pos-1] for pos in sorted_pos]
        new_srt = seqlist_type[i].replace(seqlist[i],''.join(srt) )
        new_list.append(new_srt)
    return new_list

def make_hybrid(exp_1, exp_2):
    hybrid =[]
    for i in range(len(exp_1)):
        if exp_1[i] == 'T' and exp_2[i] != 'T':
            hybrid.append(exp_1[i])
        else:
            hybrid.append(exp_2[i])
    return ''.join(hybrid)

def get_best_mtch(result_list):
    h_mtch_list = [result_list[i][0] for i in range(len(result_list))]
    h_mtch = 0
    for i in range(len(h_mtch_list)):
        if h_mtch_list[i] >= h_mtch:
            h_mtch = h_mtch_list[i]
            bn = result_list[i]
            best_seq_ind = i
    return bn,best_seq_ind

def get_ML(seq1, seq2):
    alpha = 0.01
    beta = 0.2
    total_prob = []
    for i in range(len(seq1)):
        if seq1[i] == 'A' and seq2[i] == 'A':
            prob = 1 - alpha
        if seq1[i] == 'T' and seq2[i] == 'A':
            prob = alpha
        if seq1[i] == 'T' and seq2[i] == 'T':
            prob = 1 - beta
        if seq1[i] == 'A' and seq2[i] == 'T':
            prob = beta
        if seq1[i] == '?' and seq2[i] == 'T':
            prob = 1 - beta
        if seq1[i] == '?' and seq2[i] == 'A':
            prob = 1 - alpha

        total_prob.append(prob)
    result = np.prod(total_prob)
    return result

def LR_pval(L1,L2):
    LR = 2 * abs((np.log(L1) - np.log(L2)))
    p = chi2.sf(LR, 1)
    return ' %.30f' % p


def get_ancestor(lis1, lis2, tp):
    for i in range(0, len(lis1)):
        if lis1[i] == tp:
            ancs = lis2[i]
    return ancs

def get_expHaplotype(lineage, expS):
  
    replacement = {}
    for i in lineage:
        replacement[i-1] = 'T'
    replacer = replacement.get
    return ''.join([replacer(n,expS[n]) for n in range(len(expS))])

def get_fasta_format(reference, seq):
    return (">" + str(reference) + "\n" + str(seq) + "\n")

def write_file(fastalist, fname):
    file = open(path+ outputDir+ os.sep+fname, 'w')
    for x in range(0, len(fastalist)):
        file.write(str(fastalist[x]))
    file.close()
    return

def get_pairwiseVariants (ex):
    pairwise =[]
    for i in range (len( ex)):
        if i == len(ex)-1:
            ex = ex[i].strip()
            pairwise.append(int(ex[:-1]))
        else:
            pairwise.append(int(ex[i]))
    return pairwise

def get_observedFrequency(seq_list, target_v1, target_v2):
    observedMM = []
    observedMW = []	
    for sequence in seq_list:
      	
        if sequence[target_v1-1] =='T' and sequence[target_v2-1] == 'A':
            observedMW.append(sequence)

        if sequence[target_v1-1] =='T' and sequence[target_v2-1] == 'T':
            observedMM.append(sequence)
    return  observedMM, observedMW

def get_expectedFrequency(observedMM, observedMW, beta):
    total = observedMM + observedMW
    expectedMW = total * beta
    expectedMM = total - expectedMW
    return expectedMM, expectedMW

def get_contingencyTable(oMM, oMW, exMM, exMW):
    table = np.array([[oMM, oMW],[exMM, exMW]])
    return table

def fix_edges(alledges):
    for i in range(len(alledges)-1):
        if alledges[i][1] != alledges[i+1][0]:
            alledges[i][1] = alledges[i+1][0]
        else:
            pass
    return alledges


def get_newEdges(ancs2decnList, lineagePath, rootindex, allparents, spltnode):
    nodelabel = str(ancs2decnList[rootindex][1])
    tedges = [[int(val[0]), int(val[1])] for val in ancs2decnList]
    nodeList = []
    N = 2
    for i in range(len(lineagePath)):
        if lineagePath[i] not in allparents:
            nodeList.append(nodelabel)
            break

        if [lineagePath[i], lineagePath[i+1]] in tedges:
            x = tedges.index([lineagePath[i], lineagePath[i+1]])

        if lineagePath[i] in spltnode:
            nodeList.append(nodelabel)
            nodelabel =  str(ancs2decnList[x][1])
            nodeList = nodeList + [ancs2decnList[x][1]]
            continue

        if ancs2decnList[x][3] == 'concurrent':
            nodelabel += ','+ str(ancs2decnList[x][1])

        if ancs2decnList[x][3] == 'sequential':
            nodeList.append(nodelabel)
            nodelabel =  str(ancs2decnList[x][1])
            nodeList = nodeList + [ancs2decnList[x][1]]
    edgeList = [nodeList[n:n+N] for n in range(0, len(nodeList), N)]

    fnlEdgeList = fix_edges(edgeList)

    return fnlEdgeList

def add_edges(tr, nd1, nd2, cs, pv, ind1, ind2):
    tr.add_node(pydot.Node(ind1, label = 'Mut: ' + str(nd1))) #+ '\n' + 'COI:' +str(cs)))
    tr.add_node(pydot.Node(ind2 , label = 'Mut: ' + str(nd2))) #+ '\n' + 'COI:' +str(cs)))
    tr.add_edge(pydot.Edge(ind1, ind2, label = 'COI:' +str(cs) + '\n'+ 'pval: '+ str(pv),))
    return tr

def get_COI(nd, CM):
    if len(nd) > 1:
        c2p = nd
        prnts = c2p[0].split(',')
        chld = c2p[1].split(',')

        if len(prnts) >= 1 and len(chld) == 1:
            coiAvg = mean([float(CM[int(chld[0])][int(x)]) for x in prnts])
            return '{:.3f}'.format(coiAvg)

        if len(prnts) >= 1 and len(chld) > 1:
            totalCOI = []
            for eachChild in chld:
                coiAvg = mean([float(CM[int(eachChild)][int(x)]) for x in prnts])
                totalCOI.append(coiAvg)
            return '{:.3f}'.format(mean(totalCOI))
    else:
        pass

def get_pval(nd, ancs2decnList):
    tedges = [[int(val[0]), int(val[1])] for val in ancs2decnList]
    if len(nd) > 1:
        c2p = nd
        prnts = c2p[0].split(',')
        chld = c2p[1].split(',')
        if [int(prnts[-1]) , int(chld[0])] in tedges:
            x = tedges.index([int(prnts[-1]) , int(chld[0])])
            return '{:.3f}'.format(float(ancs2decnList[x][2]))

##### Reads dot/.gv mutation tree to get tips, ancestors, descendants

starttime = time.time()

graph = pydot.graph_from_dot_file(inputDot)[0]
edges = [(x.get_source().replace("\"", ""), x.get_destination().replace("\"", "")) for x in graph.get_edges() if x.obj_dict['attributes'].get("constraint", "True") == "True"]

ancsVariants= [int(val[0]) for val in edges ]
descVariants= [int(val[1]) for val in edges]
treeTipVariants= list(set(descVariants) - set(ancsVariants))
treeRoot = list(set(ancsVariants) - set(descVariants))[0]

allLineages = []
for tip in treeTipVariants:
    path1 = [tip]
    while tip < len(ancsVariants) + 1:
        newtip = get_ancestor(descVariants, ancsVariants, tip)
        tip = newtip
        path1.append(tip)
        lineage = path1[::-1]
    allLineages.append(lineage[1:])


##### Read fasta file
fasta_obj = SeqIO.parse(open(path+ fastaFilename),'fasta')

fastaList  = get_Seq_ID(fasta_obj)

Original_allSequenceList = [fastaList[i].split(':')[1] for i in range(0, len(fastaList))]


###### Make expected haplotypes from the input mutation tree and make a fasta file with all the expected sequences

expectedSequenceDictionary = {}
for ind in range(len(allLineages)):
    outgrpSeq = ['A' for i in range(len(Original_allSequenceList[0]))]
    expSeq = get_expHaplotype(allLineages[ind], outgrpSeq)
    expectedSequenceDictionary['expSeq'+str(ind+1)] = str(expSeq)


all_fasta = []
for i in range(len(expectedSequenceDictionary)):
    key = 'expSeq'+str(i+1)
    fasta = get_fasta_format(key, expectedSequenceDictionary[key])
    all_fasta.append(fasta)

outputFilename1 ="ExpectedSequences.fasta"
create = write_file(all_fasta, outputFilename1)

fastaObj = SeqIO.parse(open(path+ outputDir +os.sep+outputFilename1),'fasta')
expFastaList  = get_Seq_ID(fastaObj)
expectedSequenceList  = [expFastaList[i].split(':')[1] for i in range(0, len(expFastaList))]


####### Pre-Processing data steps

missing_data = False
variant_count = get_variant_count(expectedSequenceList, missing_data)
count_pos_dict = {i+1 : variant_count[i] for i in range(len(variant_count))} # variant positions w/ mutation count stored as dictionary, dictionary1
variant_count.sort(reverse = True)
sorted_variant_position = get_sorted_position(variant_count, count_pos_dict)


sortedExpectedSequence = get_sorted_fasta(expFastaList, sorted_variant_position)
sortedFastaSequence = get_sorted_fasta(fastaList, sorted_variant_position)


## -------------------------DOUBLET DETECTION-------------------------------------

####### DD Main script beigns ########


CellID = []
CellType = []
Best_Node_list = []


exp_seq_list=[]
hybrid_list = []
exp_seq2_list=[]

doubletList = []
Go='n'
if Go!='n':
 for i in range(len(sortedFastaSequence)):
    test = i
    CellID.append(sortedFastaSequence[test].split(':')[0])



    match_info = [get_match_count(sortedFastaSequence[test].split(':'), sortedExpectedSequence[i].split(':')) for i in range(len(sortedExpectedSequence)) ]
    t,m  = get_best_mtch(match_info)
    print (match_info,t,m)	
    BestNode_expseq1 = sorted_variant_position[t[1]]
    expected_seq1 = sortedExpectedSequence[m].split(':')[1]
    exp_seq_list.append(expected_seq1)



    BestNode_expseq2 = BestNode_expseq1
    expected_seq2 = expected_seq1
    for ex_pos in t[2]:
        if len(t[2]) > 0:
            match_info2 = [ get_match_count(sortedFastaSequence[test].split(':'),sortedExpectedSequence[i].split(':')) for i in range(len(sortedExpectedSequence)) if sortedExpectedSequence[i].split(':')[1][ex_pos] == 'T' ]
            match_seq = [ sortedExpectedSequence[i].split(':') for i in range(len(sortedExpectedSequence)) if sortedExpectedSequence[i].split(':')[1][ex_pos] == 'T']
            if len(match_info2) > 0:
                t2,m2  = get_best_mtch(match_info2)
                BestNode_expseq2 = sorted_variant_position[t2[1]]
                expected_seq2 = match_seq[m2][1]
            else:
                pass
        else:
            pass


    if BestNode_expseq2 != BestNode_expseq1:
        hybrid = make_hybrid(expected_seq1, expected_seq2)
        hybrid_list.append(hybrid)


        l1 = get_ML(sortedFastaSequence[test].split(':')[1], expected_seq1)
        l2 = get_ML(sortedFastaSequence[test].split(':')[1], hybrid)
        p = LR_pval(l1,l2)

        if float(p) > 0.05:
            CellType.append("Singleton")
            Best_Node_list.append(BestNode_expseq1)
            

        if float(p) < 0.05:
            CellType.append("Doublet")
            Best_Node_list.append(str(BestNode_expseq1) + ',' +str(BestNode_expseq2))
            doubletList.append(sortedFastaSequence[test].split(':')[0])
            
    else:
        CellType.append("Singleton")
        Best_Node_list.append(BestNode_expseq1)


 DDResultFilename = 'Doublet_Detection_Result0.05.txt'
 file = open(path + outputDir+os.sep+ DDResultFilename, 'w')
 for i in range(len(CellType)):
    file.write(str(CellID[i]) + '\t' + str(CellType[i]) + '\t' + str(Best_Node_list[i]) + '\n') #+ '\t' + str(TopHap_ID_list[i])  + '\n')
 file.close()


##### remove doublet cells from downstream analysis
 doubletList=[]
 for cl in doubletList:
    for fs in fastaList:
        if fs.split(':')[0] == cl:
            fastaList.remove(fs)
          
 removedDoubletlist = []
 for i in range(len(fastaList)):
    fasta = get_fasta_format(fastaList[i].split(':')[0], fastaList[i].split(':')[1])
    removedDoubletlist.append(fasta)

 outputFilename2 ="postDoubletDetection.fasta"
 create = write_file(removedDoubletlist, outputFilename2)

 fasta_obj = SeqIO.parse(open(path+ outputDir+ os.sep+ outputFilename2),'fasta')
 fastaList  = get_Seq_ID(fasta_obj)
 PostDoublet_allSequenceList = [fastaList[i].split(':')[1] for i in range(0, len(fastaList))]


## -------------------------MUTATION ORDER-------------------------------------

######## Mutation Order Main script begins

ancsToDescVariants = edges
treeRoot = list(set(ancsVariants) - set(descVariants))[0]
trInd = ancsVariants.index(treeRoot)

significance_alpha = 0.01


final_result = []

for x in range(len(ancsToDescVariants)):
    
    target_variant1 = int(ancsToDescVariants[x][0])
    target_variant2 = int(ancsToDescVariants[x][1])
   	

    if x == trInd or target_variant1==treeRoot or target_variant2==treeRoot:
        mutationOrder = 'Root'
        final_result.append([target_variant1, target_variant2,'nan', mutationOrder])

    else:
        if Remove_Doublet == "true":
            result = get_observedFrequency(Original_allSequenceList, target_variant1, target_variant2)
        else:
           
            result = get_observedFrequency(Original_allSequenceList, target_variant1, target_variant2) # before Doublet
        

        observed_mut_mut = len(result[0])
        observed_mut_wild = len(result[1])

        result2 = get_expectedFrequency(observed_mut_mut,observed_mut_wild, beta )

        expected_mut_mut = result2[0]
        expected_mut_wild = result2[1]

        ct = get_contingencyTable(observed_mut_mut, observed_mut_wild, expected_mut_mut, expected_mut_wild)

        oddsr, p = fisher_exact(ct)

        if p > significance_alpha:
            mutationOrder = 'concurrent'
        if p < significance_alpha:
            mutationOrder = 'sequential'
        final_result.append([target_variant1, target_variant2,p, mutationOrder])

mutOrderResultFilename = 'mutationConcurrentSequential_Result.txt'
file = open(path + outputDir + os.sep +mutOrderResultFilename, 'w')
for i in range(len(final_result)):
    file.write(str(final_result[i][0]) + '\t' + str(final_result[i][1]) + '\t' + str(final_result[i][2]) + '\t' + str(final_result[i][3]) + '\n') #+ '\t' + str(TopHap_ID_list[i])  + '\n')
file.close()

print(sys.platform)
## -------------------------COI SCORING & RUN FastMOA-------------------------------------

if sys.platform == "win32":
    if Remove_Doublet == 'true':
        os.system(Python+' FastMOA'+os.sep+'Generate_Label.py ' + str(outputFilename2) + ' '+ str(outputDir) +' '+  str(Remove_Doublet))
        os.system(Python+' FastMOA'+os.sep+'FastMOA_Update.py ' + outputDir + os.sep +str(outputFilename2) + ' ' +outputDir +os.sep+'Labels.txt -o DD_FastMOA_Result ' + '--flip_pass_thresh 1.0 --disable_graph_flipping')
    if Remove_Doublet == 'false':
        os.system(Python+' FastMOA'+os.sep+'Generate_Label.py ' +str(inputFasta) + ' '+ str(outputDir) +' '+  str(Remove_Doublet))
        os.system(Python+' FastMOA'+os.sep+'FastMOA_Update.py ' + str(inputFasta) + ' ' +outputDir +os.sep+'Labels.txt -o DD_FastMOA_Result ' + '--flip_pass_thresh 1.0 --disable_graph_flipping')
  
else:
    if Remove_Doublet == 'true':
        os.system(Python+' FastMOA'+os.sep+'Generate_Label.py ' +str(outputFilename2) + ' '+ str(outputDir) + ' '+  str(Remove_Doublet))
        os.system(Python+' FastMOA'+os.sep+'FastMOA_Update.py ' + outputDir + os.sep +str(outputFilename2) + ' ' +outputDir +os.sep+'Labels.txt -o DD_FastMOA_Result ' + '--flip_pass_thresh 1.0 --disable_graph_flipping')
    if Remove_Doublet == 'false':
        os.system(Python+' FastMOA'+os.sep+'Generate_Label.py ' +str(inputFasta) + ' '+ str(outputDir) +' '+  str(Remove_Doublet))
        os.system(Python+' FastMOA'+os.sep+'FastMOA_Update.py ' + str(inputFasta) + ' ' +outputDir +os.sep+'Labels.txt -o DD_FastMOA_Result ' + '--flip_pass_thresh 1.0 --disable_graph_flipping')

coiMatrixFile = open(path + 'DD_FastMOA_Result'+ os.sep +'COI_matrix.txt', 'r')
coiMatrix = coiMatrixFile.readlines()

coiMatrix = [coi.split() for coi in coiMatrix] # coiMatrix[row][column], coiMatrix[child][parent]

muInputFile = open(path + outputDir +os.sep + mutOrderResultFilename, 'r')
muRes = [val.split() for val in muInputFile] # muRes[0] = ['78', '1', '0.8503402540432693', 'concurrent'], muRes[0][0] = '78'

tedges = [[int(val[0]), int(val[1])] for val in muRes]



d = pydot.Dot(graph_type = 'digraph')

trueNode_count = []
trueNode = []
for val in ancsVariants:
    count = [i for i in ancsVariants if i == val]
    if len(count) > 1:
        trueNode_count.append([val, len(count)])
        trueNode.append(val)


allEdgeList = []
for lin in allLineages:
    edgls  = get_newEdges(muRes, lin, trInd, ancsVariants, trueNode)
    allEdgeList = allEdgeList + edgls



outTree = [['Root', allEdgeList[0][0]]]

[outTree.append(x) for x in allEdgeList if x not in outTree]

flat_list = [item for sublist in outTree for item in sublist]
All_node_outTree = []
[All_node_outTree.append(val) for val in flat_list if val not in All_node_outTree]
node_dict = dict(zip(All_node_outTree, [x for x in range(len(All_node_outTree))]))

prnt_list = []
chld_list = []
for z in range(len(outTree)-1):
    coiscore = get_COI(outTree[1:][z],coiMatrix)
    pvalues = get_pval(outTree[1:][z], muRes)
    if len(outTree[1:][z]) == 2:
        prnt_lbl = outTree[1:][z][0]
        chld_lbl = outTree[1:][z][1]
        add_edges(d,prnt_lbl, chld_lbl, coiscore, pvalues, node_dict[prnt_lbl], node_dict[chld_lbl])
    else:
        pass

ouput_graphviz_dot = d.create_dot()
d.write_dot(path + outputDir+ os.sep +'SCAN_output.dot')
if Remove_Doublet == 'true':
    d.write_png(path + outputDir+ os.sep +"SCAN_output.png")
if Remove_Doublet == 'false':
    d.write_png(path + outputDir+ os.sep +"SCAN_output.png")

endtime = time.time()
print('Total time in seconds: ' + str(endtime-starttime))
