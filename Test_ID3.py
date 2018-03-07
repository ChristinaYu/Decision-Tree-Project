import csv

## csv reader for tree nodes
treecsv = open('./id3.csv', 'rb')
treereader = csv.reader(treecsv)

## csv writer for final answer
outcsv = open('./answer.csv','wb')
output = csv.writer(outcsv)
output.writerow(["id", "class"])

## read tree nodes from file
tree = {}
for row in treereader:
  tree[int(row[0])] = [row[1], row[2], row[3]]
  
## read testing dataset and do classification 
## according to the tree nodes
def Main():
  csvfile = open('./testing.csv', 'rb')
  csvreader = csv.reader(csvfile)
  for row in csvreader:
    answer = Classification(row[1])
    output.writerow([row[0], answer])

## if a tree node shows "false" means this branch
## is the end of the classification, then return the 
## class of this brach as the final answer
def Classification(gene_seq):
  node = 0
  char_set = ['A', 'T', 'C', 'G', 'D', 'N', 'S', 'R']
  while(1):
    if (tree[node][1] == 'False'):
      return tree[node][2]
      break
    else: 
      node = node*8 + char_set.index(gene_seq[int(tree[node][0])]) + 1
      
Main()