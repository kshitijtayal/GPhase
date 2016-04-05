import numpy as np
import math
import sys
import scipy.stats.distributions
from scipy.stats import mode
from numpy import copy
import pandas as pd
from multiprocessing import Process, Manager
import multiprocessing
from collections import defaultdict
import MinHeap


WEIGHTS = [0.85355339,0.1464466]
TIMES = [0.58578643,3.41421356]


def prob(numberIndividual):
    theta = 1/math.log(2*numberIndividual)                   #theta = 1/log(2n)
    lamda1 = (theta * TIMES[0])/numberIndividual
    lamda2 = (theta * TIMES[1])/numberIndividual
    F00L1 = F11L1 = probAlphaiBetai(lamda1,2)
    F01L1 = F10L1 = probAlphaiBetai(lamda1,1)
    F00L2 = F11L2 = probAlphaiBetai(lamda2,2)
    F01L2 = F10L2 = probAlphaiBetai(lamda2,1)
    return F00L1,F01L1, F00L2,F01L2



def probAlphaiBetai(lamda,oddEven):
    proMutSingPos = []
    if oddEven % 2 == 1:
            for m in range(1, 50, 2):                                    #run for odd value of m 
                proMutSingPos.append(poisson_probability(m,lamda))        #list appending poission probabilities 
    else:   
            for m in range(0, 50, 2):                                    #run for even value of m 
                proMutSingPos.append(poisson_probability(m,lamda))        #list appending poission probabilities
    return sum(proMutSingPos)


def poisson_probability(actual, mean):
    # naive:   math.exp(-mean) * mean**actual / factorial(actual)

    # iterative, to keep the components from getting too large or small:
    p = math.exp(-mean)
    for i in xrange(actual):
        p *= mean
        p /= i+1
    return p



def probPair(beta,alpha,lamb):                # probability of beta given alpha
            
            if len(beta) != len(alpha):
                print 'Length not same'
            if lamb == 0:
                prob =  WEIGHTS[0]
                M00 = F00l1
                M01 = F01l1
            else: 
                prob = WEIGHTS[1]
                M00 = F00l2
                M01 = F01l2
            for x in range(len(beta)):
                if beta[x] == alpha[x]:
                    prob = prob * M00
                else :
                    prob = prob * M01
            return prob

def objective(m,n):
    return m*n
    
def skewness(Array):
    mydict = {}
    total = Array.shape[0]
    index = 0
    for column in range(Array.shape[1]):
        mydict[index] = abs(2 * (Array[:,column]==0).astype(int).sum() - total)
        index = index + 1
    #print mydict
    sortedcolumn = []
    for key, value in sorted(mydict.iteritems(), key=lambda (k,v): (v,k)):
        sortedcolumn.append(key)
    #print sortedcolumn
    return Array[:,sortedcolumn] , sortedcolumn



          
def SwitchAccuracy(dictAccuracy): 
              switchacc = []
              for key in dictAccuracy.keys():
                  total = 0
                  accurate = 0
                  total = dictAccuracy[key][0] - 1 
                  accurate = dictAccuracy[key][0] - 1 -  dictAccuracy[key][2]
                  if total == 0:
                        switchacc.append(1)
                        continue
                  switchacc.append(accurate/float(total))
              
              return sum(switchacc)/float(len(switchacc))	
              

              
def NewAccuracy(dictAccuracy): 
      
              acc = []
              for key in dictAccuracy.keys():
                  total = 0
                  accurate = 0
                  total = dictAccuracy[key][0] - 1 
                  accurate = dictAccuracy[key][0] - 1 -  dictAccuracy[key][3]
                  if total == 0:
                        acc.append(1)
                        continue
                  acc.append(accurate/float(total))
              return sum(acc)/float(len(acc))
              
def list_split(listsize, splitsize):
          
          if listsize%2 != 0:
               print 'Error... No. of haplotypes should be even'
               return 0
           
          splitrange = []
          
          s = listsize / splitsize
          
          if s%2 !=0:
              s = s + 1
          temp = listsize - (listsize%s)
          
          for j in range(splitsize):
              if j == 0 :
                  splitrange.append([j,j + s])
                  p = j + s
              elif j == splitsize-1:
                  splitrange.append([p,listsize])
                  break
              else: 
                  splitrange.append([p,p + s])
                  p = p +s
                   
          return splitrange

def matchscore(first_word, second_word):
    score = 0
    if len(first_word)!=len(second_word):
        print 'error'
    for x in range(len(first_word)):
        if first_word[x] == second_word[x]:
            score = score + 1
    return score


def mismatchscore(first_word, second_word):
    score = 0
    if len(first_word)!=len(second_word):
        print 'error'
    for x in range(len(first_word)):
        if first_word[x] != second_word[x]:
            score = score + 1
    return score

def NewSwitcherror(inferHaplotype,trueHaplotype): 
          
          switchrate = 0
          p = 0
          Compleinferhaplotype = [0 if x == 1 else 1 for x in inferHaplotype]
          for x in range(len(inferHaplotype)):
              if p == 0: 
                  if trueHaplotype[x] == inferHaplotype[x] :
                      continue
                  else:
                      p = 1
                      switchrate = switchrate + 1
              elif p ==1:
                  if trueHaplotype[x] != inferHaplotype[x]:
                      continue
                  else: 
                      p = 0
          return switchrate


def probPartialHaplotype(partialSolution,sortedArraybySkewness,dictColumnPosition):

    dictpreComputation = defaultdict(list)
    my_dicts1 = dict()
    objtoPair = {}                      # store pair to its objective function

    finalSolution = []
    heap = MinHeap.MinHeap()
    initialise = 0
    
    index = 0      
    for row in sortedArraybySkewness[:,-len(partialSolution):]:                # Precomputation of probability given partial Solution
            dictpreComputation[index] = [probPair(partialSolution, row, 0),probPair(partialSolution, row, 1)]
            index = index + 1

    optimalSolution = [[partialSolution,partialSolution]]
    my_dicts1["dict" + str(partialSolution)] = dictpreComputation

    for coloumnIndex in range(-(1+len(partialSolution)),-sortedArraybySkewness.shape[1]-1,-1): #ignore last coloumn
        
        
        tempArray =  sortedArraybySkewness[:,coloumnIndex]
        
        my_dicts = my_dicts1
        my_dicts1 = {}
        
        
        if initialise == 1 :
            optimalSolution = []
            for key,pair in objtoPair.items():
                optimalSolution.append(pair)
            heap = MinHeap.MinHeap()
            objtoPair = {}
            
        for element1,element2 in optimalSolution:
            
            dictpreComputationelm1 = my_dicts["dict" + str(element1)]
            dictpreComputationelm2 = my_dicts["dict" + str(element2)]
            dictpreComputation10 = defaultdict(list)
            dictpreComputation11 = defaultdict(list)
            dictpreComputation20 = defaultdict(list)
            dictpreComputation21 = defaultdict(list)
            
            for rownumber in range(tempArray.shape[0]):
                if tempArray[rownumber] == 0:
                    dictpreComputation10[rownumber] = [F00l1 * dictpreComputationelm1[rownumber][0] , F00l2 * dictpreComputationelm1[rownumber][1]]
                    dictpreComputation11[rownumber] = [F01l1 * dictpreComputationelm1[rownumber][0] , F01l2 * dictpreComputationelm1[rownumber][1]]
                    dictpreComputation20[rownumber] = [F00l1 * dictpreComputationelm2[rownumber][0] , F00l2 * dictpreComputationelm2[rownumber][1]]
                    dictpreComputation21[rownumber] = [F01l1 * dictpreComputationelm2[rownumber][0] , F01l2 * dictpreComputationelm2[rownumber][1]]

                else:
                    dictpreComputation10[rownumber] = [F01l1 * dictpreComputationelm1[rownumber][0] , F01l2 * dictpreComputationelm1[rownumber][1]]
                    dictpreComputation11[rownumber] = [F00l1 * dictpreComputationelm1[rownumber][0] , F00l2 * dictpreComputationelm1[rownumber][1]]
                    dictpreComputation20[rownumber] = [F01l1 * dictpreComputationelm2[rownumber][0] , F01l2 * dictpreComputationelm2[rownumber][1]]
                    dictpreComputation21[rownumber] = [F00l1 * dictpreComputationelm2[rownumber][0] , F00l2 * dictpreComputationelm2[rownumber][1]]

            tempProb10 = sum(sum(summ) for summ in dictpreComputation10.values())
            tempProb11 = sum(sum(summ) for summ in dictpreComputation11.values())
            tempProb20 = sum(sum(summ) for summ in dictpreComputation20.values())
            tempProb21 = sum(sum(summ) for summ in dictpreComputation21.values())
    
            solutionspace = []
            if tempProb10 > tempProb11 :
                solutionspace.append([[0] + element1,[1] + element2])
             
            else : solutionspace.append([[1] + element1,[0] + element2])
              
                
            if tempProb20 > tempProb21 and tempProb10 > tempProb11 :
                solutionspace.append([[1] + element1,[0] + element2])
               
            elif tempProb20 < tempProb21 and tempProb10 < tempProb11 :
                solutionspace.append([[0] + element1,[1] + element2])
               
            
            for inferhap11, inferhap22 in solutionspace: 
                
                if      inferhap11 == [0] + element1 and inferhap22 == [1] + element2 :
                        dictpreComputation0 =  dictpreComputation10
                        dictpreComputation1 =  dictpreComputation21
                        objectValue = objective(tempProb10,tempProb21)
                
                elif    inferhap11 == [1] + element1 and inferhap22 == [0] + element2 :
                        dictpreComputation0 =  dictpreComputation11
                        dictpreComputation1 =  dictpreComputation20
                        objectValue = objective(tempProb11,tempProb20)
                
                if heap.currentSize < 20:
                
                    heap.insert(objectValue)
                    objtoPair[objectValue] = [inferhap11,inferhap22]
                    my_dicts1["dict" + str(inferhap11)] = dictpreComputation0
                    my_dicts1["dict" + str(inferhap22)] = dictpreComputation1
                   
            
                else :
                    if  objectValue <  heap.heapList[1]:   #mininmum value check condition
                            continue
                    else:
                        minpop = heap.delMin()             #Remove mininmum from heap
                        objtoPair.pop(minpop, None)        #Remove minimum from dictionary
                    
                        heap.insert(objectValue)
                        objtoPair[objectValue] = [inferhap11,inferhap22]
                        my_dicts1["dict" + str(inferhap11)] = dictpreComputation0
                        my_dicts1["dict" + str(inferhap22)] = dictpreComputation1
            
            if initialise == 0:
                initialise = 1
                
    try:
            pairSolution = objtoPair[max(heap.heapList)]
    
    except  KeyError:
            print 'Key Error'
            return sortedArraybySkewness[1,]
    
    
    try:
            partialSolution = pairSolution[1]
    
    except  UnboundLocalError:
            print 'Unbound Local Error'
            return sortedArraybySkewness[1,]
    

            

    for key, value in sorted(dictColumnPosition.iteritems(), key=lambda (k,v): (v,k)):
        finalSolution.append(int(partialSolution[key]))
    return finalSolution


def Computation(hap_range,haplotypeMatrix,haplotype_TestMatrix,dictAccuracy) :
   

    for value in range(hap_range[0],hap_range[1] ,2):
            #print 'Haplotype number:' ,value
            #trainHaplotype = np.delete(haplotypeMatrix,[value,value+1],axis = 0)     #training haplotype , removed test(inferenced) haplotype\n",
            
	    hapwoutDot1 = haplotype_TestMatrix[value]
            hapwoutDot2 = haplotype_TestMatrix[value + 1]
	    #print hapwoutDot1
	    #print hapwoutDot2		
            #trainHaplotype = np.delete(haplotypeMatrix,dotposition,axis = 1)         #training haplotype , removed coloumn with dot in test dataset
            trainHaplotype = haplotypeMatrix
	    partialSolution = []
            for f, b in zip(hapwoutDot1,  hapwoutDot2):
                if str(f)+'|'+str(b) == '0|0':
                    partialSolution.append(0)
                elif str(f)+'|'+str(b) == '1|1':
                    partialSolution.append(1)
                else:
                    partialSolution.append('')
            #print 'Partial Solution is' +  str(partialSolution )
            dictColumnPosition = {}
            ambiguousPosition = [x for x in range(len(partialSolution))if partialSolution[x] == '']  # ambiguous position in array,
	    if len(ambiguousPosition) <=1:
		continue
            knownPosition = [x for x in range(len(partialSolution))if partialSolution[x] != '']  # ambiguous position in array,
            partialAmbiguousArrayNoS = trainHaplotype[:,ambiguousPosition]
            partialKnownArray = trainHaplotype[:,knownPosition]
	    posit = []
            partialAmbiguousArray , posit = skewness(partialAmbiguousArrayNoS)
            reaarangedArray = np.concatenate((partialAmbiguousArray,partialKnownArray),axis=1) #concatenate partialAmbiguousArray,partialKnownArray,
            #print reaarangedArray
            partialSolutionwithoutSpace = [partialSolution[x] for x in range(len(partialSolution))if partialSolution[x] != '']  #partial solution,
            #print partialSolutionwithoutSpace
            pos = 0
            for element in posit:            #0th coloumn of skewed coloum contain 1st coloumn of original array\n",
                dictColumnPosition[pos] = ambiguousPosition[element]
                pos += 1
            for element in knownPosition:
                dictColumnPosition[pos] = element
                pos += 1
            
            #print dictColumnPosition
            finalSolution = probPartialHaplotype(partialSolutionwithoutSpace,reaarangedArray,dictColumnPosition)
       
            
            inferSolAmbigious = [ finalSolution[i] for i in ambiguousPosition]
            
	    hap1sol = [ hapwoutDot1[i] for i in ambiguousPosition]   #Real ambiguous position sol1 for switch-error
            hap2sol = [ hapwoutDot2[i] for i in ambiguousPosition]   #Real ambiguous position sol2 for switch-error
            
            
            no_of_switch = min(NewSwitcherror(inferSolAmbigious,hap1sol),NewSwitcherror(inferSolAmbigious,hap2sol))
            
            no_of_match = max(matchscore(inferSolAmbigious,hap1sol),matchscore(inferSolAmbigious,hap2sol))
            
            no_of_mismatch = min(mismatchscore(inferSolAmbigious,hap1sol),mismatchscore(inferSolAmbigious,hap2sol))
            
            dictAccuracy[value] = [len(ambiguousPosition),no_of_match,no_of_switch,no_of_mismatch]
    


def main(haplotypeMatrix,haplotype_TestMatrix):
        
        #dictAccuracy = {}   # {0 : [no. of snp , correctly identified]}\n",
        manager = Manager()
        dictAccuracy = manager.dict()   # synchronize dictionary for multiprocessing	

        nprocs = []   # saves the process
        
        #list_split1 = [[0,100],[100,202]]
        
        #for item in list_split(haplotype_TestMatrix.shape[0],4): 
        #	print 'range of haplotype given to each thread: ' ,item
        
        for item in list_split(haplotype_TestMatrix.shape[0],arg3):                      #Specify number of thread
            n = multiprocessing.Process(target=Computation, args=(item,haplotypeMatrix,haplotype_TestMatrix, dictAccuracy))   # multiprocessing 
            nprocs.append(n)
            n.start()
        
        for i in nprocs:   
             i.join()        # waiting for all the process to finish
        
        if     bool(dictAccuracy) == False:
                 print        
        else : 
		#print dictAccuracy	
		print 'Switch Accuracy is : ', SwitchAccuracy(dictAccuracy)
        	print 'Accuracy is : ', NewAccuracy(dictAccuracy)
        
    
if __name__ == '__main__':
    
    try:
        arg1 = sys.argv[1]
	arg2 = sys.argv[2]
	arg3 = int(sys.argv[3])
    except IndexError:
        print "Usage: GPhase.py <arg1> <arg2> <arg3>"
        sys.exit(1)
  
    df=pd.read_csv(arg1, sep=',',header=None)    #get training haplotype file
    haplotypeMatrix = df.as_matrix()
    csv = np.genfromtxt (arg2, delimiter=",",dtype=int) #test
    haplotype_TestMatrix = copy(csv)
    haplotype_TestMatrix[csv==-1] = 0  # insert 0 for missing values
    global F00l1,F01l1,F00l2,F01l2 
    print haplotypeMatrix.shape[0]
    F00l1,F01l1,F00l2,F01l2 = prob(haplotypeMatrix.shape[0])
    if haplotypeMatrix.shape[1] == haplotype_TestMatrix.shape[1] and haplotype_TestMatrix.shape[0] > 1: # Number of Markers should be same for test and train data
    	main(haplotypeMatrix,haplotype_TestMatrix)


        
