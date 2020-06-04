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

def read_vj_clonotypes(clonotypes,clonotype_file,lowercdr3length,uppercdr3length,blacklist=False):
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
                       if triple in clonotypes:
                           clonotypes[triple]=clonotypes[triple]+float(row[3])
                       else:
                           clonotypes[triple]=float(row[3])
            csv_file.close()
   else:
        with open(clonotype_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ') 
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                    triple='%s %s %s'%(row[0],row[1],row[2])
                    if triple in clonotypes: 
                        clonotypes[triple]=clonotypes[triple]+float(row[3])
                    else:
                        clonotypes[triple]=float(row[3])
            csv_file.close()
   return clonotypes

def create_subsample_list(clonotypes,abundance=False):
   '''
    This function will create a subsampled list of entries. It assumes 
    the input dictionary has only one value.
   '''
   subsample=[]
   if abundance==True:
      for k in clonotypes.keys():
          for i in range(0,int(clonotypes[k])):
              subsample.append(k)
      return subsample
   else:
      return (clonotypes.keys())

def common_clonotypes(clonotypesA,clonotypesB):
   return list(set(clonotypesA) & set(clonotypesB))

def species_count(clonotypes):
    return len(clonotypes)

'''------------------------------------------------------'''
''' You must look things up in your reference dictionary '''
'''------------------------------------------------------'''
def individual_count(dictionary_clonotypes,clonotypes):
    count=float(0)
    for key in clonotypes:
        count = dictionary_clonotypes[key]+count
    return count

def simpson(sample):
    v3j_clonotypes={}
    '''--------------------------------'''
    '''Create a dictionar for sampling '''
    '''--------------------------------'''
    for element in sample:
        if element in v3j_clonotypes:
            v3j_clonotypes[element]=+v3j_clonotypes[element]+1.0
        else:
            v3j_clonotypes[element]=float(1.0)

    '''------------------------------'''
    ''' Normalize items in the array '''
    '''------------------------------'''
    v3j_clonotypes_sorted=sorted(v3j_clonotypes.items(), key=lambda item: item[1],reverse=True)
    sum_total_individuals=0
    for entry in v3j_clonotypes_sorted:
       sum_total_individuals=sum_total_individuals+entry[1]

    '''--------------------------------'''
    ''' Compute the simpson index here '''
    '''--------------------------------'''
    simpson=0.0000
    for entry in v3j_clonotypes_sorted:
        simpson= simpson + float(entry[1]*(entry[1]-1.0)) / float(( sum_total_individuals*(sum_total_individuals-1)))
    return simpson

def berger_parker(sample):
    v3j_clonotypes={}
    '''--------------------------------'''
    '''Create a dictionar for sampling '''
    '''--------------------------------'''
    for element in sample:
        if element in v3j_clonotypes:
            v3j_clonotypes[element]=+v3j_clonotypes[element]+1.0
        else:
            v3j_clonotypes[element]=float(1.0)

    '''------------------------------'''
    ''' Normalize items in the array '''
    '''------------------------------'''
    v3j_clonotypes_sorted=sorted(v3j_clonotypes.items(), key=lambda item: item[1],reverse=True)
    sum_total_individuals=float(0)
    max_individual=float(0)
    for entry in v3j_clonotypes_sorted:
       sum_total_individuals=sum_total_individuals+entry[1]
       if entry[1] > max_individual:
           max_individual=entry[1]
    return max_individual/sum_total_individuals

def species_richnesss(sample):
    v3j_clonotypes={}
    '''--------------------------------'''
    '''Create a dictionar for sampling '''
    '''--------------------------------'''
    for element in sample:
        if element in v3j_clonotypes:
            v3j_clonotypes[element]=+v3j_clonotypes[element]+1.0
        else:
            v3j_clonotypes[element]=float(1.0)
    return len(v3j_clonotypes)

def sorensen(NA,NB):
    if NA>NB:
        return (( 2.0*NB)/(NA+NB) ) 
    elif NB>NA:
       return (( 2.0*NA)/(NA+NB) )
    else:
       return (( 2.0*NB)/(NA+NB) )

start = time.time()
CLI=argparse.ArgumentParser()
CLI.add_argument(
     '--clonotype-file1',
     type=str
)
CLI.add_argument(
     '--clonotype-file2',
     type=str
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
     '--number-sweeps',
     type=int,
     default=1000,
)
CLI.add_argument(
     '--initialize',
     type=int,
     default=100,
)
CLI.add_argument(
     '--outfile',
     type=str
)

CLI.add_argument(
     '--outfile-summary-statistics',
     type=str,
     default='statistics.txt'
)

CLI.add_argument("--abundance",
                  action='store_true',
                  help='consider abundance calculation',
                  default=False)

CLI.add_argument("--alpha",
                  action='store_true',
                  help='consider abundance calculation',
                  default=False)

CLI.add_argument("--beta",
                  action='store_true',
                  help='consider abundance calculation',
                  default=False)

CLI.add_argument("--alpha-metrics",
                  choices=('simpson', 'berger-parker', 'species-richness'),
                  help='consider abundance calculation')

CLI.add_argument("--beta-metrics",
                  choices=('overlap','percentage-overlap','morisita-horn'),
                  help='consider abundance calculation')

CLI.add_argument("--standardize-sample-size",
                  action='store_true',
                  help='use the smaller of the two population sizes',
                  default=False)

CLI.add_argument("--set-sample-size",
                  action='store_true',
                  help='manually set the sample size',
                  default=False)

CLI.add_argument("--sample-size",
                  type=int,
                  help='use the smaller of the two population sizes',
                  default=1000000)



args=CLI.parse_args()

#if args.subgroups > len(args.clonotype_files):
#   print 'Cannot have more subgroups than number of clonotype_files'
#   sys.exit(1)

sizes=[]
table=[]
vj1_clonotypes={}
vj2_clonotypes={}
'''------------------------------'''
''' Read in the number of papers '''
'''------------------------------'''
vj1_clonotypes=read_vj_clonotypes(vj1_clonotypes,args.clonotype_file1,args.lower_bound,args.upper_bound)
vj2_clonotypes=read_vj_clonotypes(vj2_clonotypes,args.clonotype_file2,args.lower_bound,args.upper_bound)

'''-----------------------------------------------------------------'''
''' Determine what individuals you would like to read into the file '''
'''-----------------------------------------------------------------'''
minimum_size=0
if args.abundance:
   individuals_a=int(individual_count(vj1_clonotypes,vj1_clonotypes.keys()))
   individuals_b=int(individual_count(vj2_clonotypes,vj2_clonotypes.keys()))
   minimum_size=min(individuals_a,individuals_b)
   if args.set_sample_size==True:
       minimum_size=args.sample_size
else:
   individuals_a=int(species_count(vj1_clonotypes.keys()))
   individuals_b=int(species_count(vj2_clonotypes.keys()))
   minimum_size=min(individuals_a,individuals_b)
   if args.set_sample_size==True:
       minimum_size=args.sample_size

'''-----------------------------------------'''
''' Create one subsampled list to work with '''
'''-----------------------------------------'''
subsample=create_subsample_list(vj1_clonotypes,abundance=args.abundance)+create_subsample_list(vj2_clonotypes,abundance=args.abundance)

'''--------------------------------------------------------------------------'''
''' Initialize large array for subsampling by randomly shuffling a few times '''
'''--------------------------------------------------------------------------'''
values=np.zeros(args.number_sweeps)
for i in range(1,args.initialize):
        random.shuffle(subsample)

outfile=gzip.open(args.outfile,'wb')
for i in range(1,args.number_sweeps):
       random.shuffle(subsample)
       '''------------------------------------'''
       ''' Compute the beta diversity measure '''
       '''------------------------------------'''
       if args.beta==True:
           if args.beta_metrics=='overlap':
               if args.standardize_sample_size:
                  values[i]=float(len(common_clonotypes(subsample[0:minimum_size],subsample[-minimum_size:]))) 
                  outfile.write('%d %d %f\n'%(minimum_size,minimum_size, values[i]))
               else:
                  values[i]=float(len(common_clonotypes(subsample[0:individuals_a],subsample[individuals_a:])))
                  outfile.write('%d %d %f\n'%(individuals_a,individuals_b, values[i]))
           elif args.beta_metrics=='percentage-overlap':
               if args.standardize_sample_size:
                  shared=common_clonotypes(subsample[0:minimum_size],subsample[-minimum_size:])
                  values[i]=100.0*float(len(shared))/float(minimum_size)
                  outfile.write('%d %d %d %f\n'%(minimum_size,minimum_size, len(shared), values[i]))
               else: 
                  shared=common_clonotypes(subsample[0:individuals_a],subsample[individuals_a:]) 
                  values[i]=100.0*float(len(shared))/float(min(individuals_a,individuals_b))
                  outfile.write('%d %d %f\n'%(individuals_a,individuals_b, len(shared), values[i]))
outfile.close()

'''---------------------'''
''' Dump out statistics '''
'''---------------------'''
outfile=gzip.open(args.outfile_summary_statistics,'wb')
outfile.write('%s %s %f %f %f\n'%(args.clonotype_file1,args.clonotype_file2,np.median(values),np.mean(values),np.std(values)))
outfile.close()
