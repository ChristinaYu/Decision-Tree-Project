import csv
import math
import random
import sys

## output file, with tree node index and attribute 
outcsv = open('./id3.csv','wb')
output = csv.writer(outcsv)

def Main():
  num_seq, gene_seq, class_seq = [],[], []

  # extract data from training file
  csvfile = open('./training.csv', 'rb')
  csvreader = csv.reader(csvfile)
  for row in csvreader:
    num_seq.append(row[0])
    gene_seq.append(row[1])
    class_seq.append(row[2])
  num_seq = num_seq[0:-1]
  gene_seq= gene_seq[0:-1]   
  class_seq = class_seq[0:-1]   
  
  # start finding proper attributes
  Find_Attr(gene_seq, class_seq,0,'NULL')
  
  
## set a certain position of a gene sequence as an attribute
## split the genes into to 8 groups based on the character on the position
def Attr_Position(position, gene_seq, class_seq):
  gene_dic = {'A':[], 'T':[], 'C':[],'G':[], 'D':[], 'N':[], 'S':[], 'R':[]}
  class_dic = {'A':[], 'T':[], 'C':[],'G':[], 'D':[], 'N':[], 'S':[], 'R':[]}
  for i in range(len(gene_seq)):
      gene_dic[gene_seq[i][position]].append(gene_seq[i])
      class_dic[gene_seq[i][position]].append(class_seq[i])
  return gene_dic, class_dic

            
## go through all attributes to find attribute with maximum gain recursively
## if one branch does not pass chisquare test, stop it
## export the result to the output csv file
def Find_Attr(gene_seq, class_seq, index, branch_char):
  char_set = ['A', 'T', 'C', 'G', 'D', 'N', 'S', 'R']
  gain = []
  method = 0

  if(sys.argv[1] == 'G'): method = 1
  elif(sys.argv[1] == 'E'): method = 2
  else: 
    print("Error! Input 'arg1=E' as Entropy, 'arg1=G' as Gini index. Please retry!")
    exit()
  for pos in range(60):
    gene_dic, class_dic = Attr_Position(pos,gene_seq,class_seq)
    if(method == 1):
      gain.append(GiniGain(class_seq,class_dic))
    elif(method == 2):
      gain.append(Gain(class_seq,class_dic))

  max_index = gain.index(max(gain))
  gene_dic, class_dic = Attr_Position(max_index,gene_seq,class_seq)

  ## chi-square test
  gene_seq_chi = []
  for gene in gene_seq:
    gene_seq_chi.append(gene[max_index])
  pass_chisquare, critical_value = Chisquare(gene_seq_chi,class_seq)
  if(sys.argv[2] != '0'):
    output.writerow([index,max_index,pass_chisquare, Majority(class_seq)])
  else:
    output.writerow([index,max_index,index<8**60, Majority(class_seq)])

  #output.writerow([index,max_index,max(gain),len(gene_seq),pass_chisquare, critical_value,branch_char,Majority(class_seq)])
  if (not pass_chisquare):
    return 0
  if (sys.argv[2] == '0'):
    if (index > 8**60): ## tree depth = 60
      return 0
  for k in range(8):
    if (sys.argv[2] == '0'):
      if (Entropy(class_dic[char_set[k]])!=0):
      #if (len(class_dic[char_set[k]])!=0):
        Find_Attr(gene_dic[char_set[k]],class_dic[char_set[k]], index*8+k+1,char_set[k])
      else:
        output.writerow([index*8+k+1,max_index,False, Majority(class_seq)])
    else:
      Find_Attr(gene_dic[char_set[k]],class_dic[char_set[k]], index*8+k+1,char_set[k])
  return 0
     
## once a branch is stopped, find the majority of the classes of the 
## gene sequences in the branch, set it as the final answer of the branch
def Majority(seq):
  count = [seq.count("N"), seq.count("IE"), seq.count("EI")]
  class_name = ['N','IE','EI']
  ## handle situation of a tie
  if(count[0] == count[1] and count[0]==count[2]):
    return class_name[random.randint(0,2)]
  elif(count[0] == count[1] and count[0]>count[2]):
    return class_name[random.randint(0,1)]
  elif(count[1] == count[2] and count[1]>count[0]):
    return class_name[random.randint(1,2)]
  elif(count[0] == count[2] and count[0]>count[1]):
    return class_name[random.randint(0,1)*2]
  else:  
    max_count = count.index(max(count))
    return class_name[max_count]

# calculate gini index gain
def GiniGain(dataset, class_dic):
  total = len(dataset)
  giniGain = Gini(dataset)
  for k in class_dic:
    if (Entropy(class_dic[k])!=0):
      giniGain = giniGain - float(len(class_dic[k]))/total*Gini(class_dic[k])
  return giniGain

# compute the gini value for attribute
def Gini(seq):
  if (len(seq)==0):
    return 0      
  n = float(seq.count("N"))
  ie = float(seq.count("IE"))
  ei = float(seq.count("EI"))
  total = n+ei+ie
  n = n / total
  ie = ie / total
  ei = ei / total
  gini = 1 - (n*n + ie*ie + ei*ei)
  return gini
  
## with S split into subsets and given entropy of S, calculate gain
def Gain(dataset, class_dic):
  total = len(dataset)
  gain = Entropy(dataset)
  for k in class_dic:
    if (Entropy(class_dic[k])!=0):
      gain = gain - float(len(class_dic[k]))/total*Entropy(class_dic[k])
  return gain 
  
## given a sequence of n/ie/ei, calculate entropy
def Entropy(seq):   
  if (len(seq)==0):
    return 0      
  n = float(seq.count("N"))
  ie = float(seq.count("IE"))
  ei = float(seq.count("EI"))
  total = n+ei+ie
  n = n / total
  ie = ie / total
  ei = ei / total
  entropy = - entro_log(n) - entro_log(ie) - entro_log(ei)
  return entropy

## appendix to Entropy function
def entro_log(x):
  if (x==0):
    return 0
  else:
    return x * math.log(x,2)

## chi-square test
def Chisquare(char_seq, class_seq):#real count, expect
  ei_dic = {'A':[0,0], 'T':[0,0], 'C':[0,0],'G':[0,0], 'D':[0,0], 'N':[0,0], 'S':[0,0], 'R':[0,0]} 
  ie_dic = {'A':[0,0], 'T':[0,0], 'C':[0,0],'G':[0,0], 'D':[0,0], 'N':[0,0], 'S':[0,0], 'R':[0,0]}
  n_dic = {'A':[0,0], 'T':[0,0], 'C':[0,0],'G':[0,0], 'D':[0,0], 'N':[0,0], 'S':[0,0], 'R':[0,0]}
  
  char_set = ['A', 'T', 'C', 'G', 'D', 'N', 'S', 'R']
  for x in range(len(class_seq)):
    if (class_seq[x]=='IE'):
      ie_dic[char_seq[x]][0]+=1
    elif (class_seq[x]=='EI'):
      ei_dic[char_seq[x]][0]+=1
    elif (class_seq[x]=='N'):
      n_dic[char_seq[x]][0]+=1
  total = len(char_seq)
  if (total == 0):
    return False,0; 
  class_total = {'EI':0,'IE':0,'N':0}
  char_total = {'A':0, 'T':0, 'C':0, 'G':0, 'D':0, 'N':0, 'S':0, 'R':0}
  critical_value = 0
  for c in char_set:
    class_total['EI'] += ei_dic[c][0] 
    class_total['IE'] += ie_dic[c][0]
    class_total['N'] += n_dic[c][0]
    char_total[c] = ei_dic[c][0] + ie_dic[c][0] + n_dic[c][0]  
  for c in char_set:
    ei_dic[c][1] = float(char_total[c] * class_total['EI']) / total
    ie_dic[c][1] = float(char_total[c] * class_total['IE']) / total
    n_dic[c][1] = float(char_total[c] * class_total['N']) / total
  for c in char_set:
    critical_value += cv(ei_dic[c]) + cv(ie_dic[c]) + cv(n_dic[c])

  if(sys.argv[2] == '95'): x = 23.685
  elif(sys.argv[2] == '99'): x = 29.141
  # when x = 0, chi-square test is not taken into effect
  # elif(sys.argv[2] == '0'): x = 0 
  else: 
    print("Error! Input arg2=95 or arg2=99. Please retry!")
    exit()
    
  if (critical_value > x):
    return True, critical_value
  else:
    return False, critical_value  

## appendix to critical value calculation
def cv(values):
  real = values[0]
  expect = values[1]
  if (expect == 0):
    return 0
  return float((real-expect)**2) / expect

Main()