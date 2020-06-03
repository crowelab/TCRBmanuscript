import gzip 
import os
import errno
import sys
import csv
import time
import argparse
import numpy as np
import random 
from itertools import combinations
from scipy.stats import rankdata

def read_v3j_clonotypes(clonotypes,clonotype_file,lowercdr3length,uppercdr3length,blacklist,vdjcount=False):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''

   if  os.path.splitext(clonotype_file)[-1] in [ '.gz', '.gzip']:
        with gzip.open(clonotype_file,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                       triple='%s %s %s'%(row[0],row[1],row[2])
                       #if row[0] not in BADGERMLINES:
                       if triple in clonotypes:
                           if vdjcount:
                              clonotypes[triple]=clonotypes[triple]+float(row[4])
                           else:
                              clonotypes[triple]=clonotypes[triple]+float(row[3])
                       else:
                           if float(row[3])>0: 
                              clonotypes[triple]={}
                              if vdjcount:
                                 clonotypes[triple]=float(row[4])
                              else:
                                 clonotypes[triple]=float(row[3])
            csv_file.close()
   else:
        with open(clonotype_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                    triple='%s %s %s'%(row[0],row[1],row[2])
                    if tripe in clonotypes:
                        if vdjcount:
                           clonotypes[triple]=clonotypes[triple]+float(row[4])
                        else:
                           clonotypes[triple]=clonotypes[triple]+float(row[3])
                    else:
                        if float(row[3])>0:
                           clonotypes[triple]={}
                           if vdjcount:
                              clonotypes[triple]=float(row[4])
                           else:
                              clonotypes[triple]=float(row[3])
            csv_file.close()
   return clonotypes

start = time.time()
CLI=argparse.ArgumentParser()
CLI.add_argument(
     '--clonotype-files',
     nargs="*",
     type=str,
     default=[],
)
CLI.add_argument(
     '--lower-bound',
     type=int,
     default=2,
)
CLI.add_argument(
     '--upper-bound',
     type=int,
     default=50,
)
CLI.add_argument(
     '--outfile',
     type=str
)
CLI.add_argument(
     '--group-size',
     type=int,
     default=2
)
CLI.add_argument("--dump-clonotypes",
     action='store_true',
     help='consider abundance calculation',
     default=False)
CLI.add_argument("--common-clonotypes-file",
                  type=str,
                  help='consider abundance calculation',
                  default='common-clonotypes.dat.gz')


args=CLI.parse_args()
sizes=[]
table=[]
vj_clonotypes={}
rank_table=[]
if args.group_size > len(args.clonotype_files):
   sys.exit(1)

'''------------------------------'''
''' Read in the number of papers '''
'''------------------------------'''
for clonotype_file in args.clonotype_files:
    vj_clonotypes=read_v3j_clonotypes(vj_clonotypes,clonotype_file,args.lower_bound,args.upper_bound,[],False)
    table.append(vj_clonotypes.keys())
    rank_table.append(rankdata(vj_clonotypes.values()))
    sizes.append(len(vj_clonotypes.keys()))
    '''----------------------------------------'''
    ''' likely you should use a simple counter '''
    '''----------------------------------------'''
    vj_clonotypes.clear()

outfile=gzip.open(args.outfile,'wb')
indeces=range(0,len(table))

''' --------------------------- '''
'''  Check overlaps for n-wise  '''
''' --------------------------- '''
group_size=2
common={}
while group_size <= len(args.clonotype_files):
    for index in combinations(indeces, group_size):
        list_of_sets=[]
        clonotype_counts=[]
        file_names=[]
        for l in index:
            list_of_sets.append(set(table[l]))
            clonotype_counts.append(sizes[l])
            file_names.append(args.clonotype_files[l])
        outfile.write('%s %s %d\n'%(' '.join(file_names),' '.join(map(str,clonotype_counts)),len(set.intersection(*list_of_sets))))
        #print '%s %s %d\n'%(' '.join(file_names),' '.join(map(str,clonotype_counts)),len(set.intersection(*list_of_sets)))
        if group_size == len(args.clonotype_files) and args.dump_clonotypes:
           common_clonotypes=gzip.open(args.common_clonotypes_file,'wb')
           for clonotype in list(set.intersection(*list_of_sets)):
               common_clonotypes.write('%s\n'%clonotype)
               common[clonotype]=-1
           common_clonotypes.close()
        del list_of_sets
        del clonotype_counts
        del file_names
    group_size=group_size+1
outfile.close()
final_table={}
for group_listing,rank_listing in zip(table,rank_table):
    for i in range(0,len(group_listing)):
        if group_listing[i] in common:
            if group_listing[i] in final_table:
               final_table[group_listing[i]].append(rank_listing[i])
            else:
               final_table[group_listing[i]]=[]
               final_table[group_listing[i]].append(rank_listing[i])

#for key in final_table.keys():
#    print key,final_table[key][0],final_table[key][1]
