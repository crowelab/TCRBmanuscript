#!/usr/bin/env python
import math
import gzip 
import os
import errno
import sys
import csv
import time
import operator
import itertools
import argparse
import numpy as np
from scipy.stats import t
from random import shuffle
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
font = {'family' : 'normal',
#        'weight' : 'bold',
        'size'   : 24}

matplotlib.rc('font', **font)


def read_v3j_clonotypes(clonotypes,clonotype_file,lowercdr3length,uppercdr3length,blacklist=False,vdjcount=False):
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
                           clonotypes[triple]['clonotype-count']=clonotypes[triple]['clonotype-count']+float(row[3])
                           if vdjcount:
                               clonotypes[triple]['vdj-count']=clonotypes[triple]['vdj-count']+float(row[4])
                       else:
                           if float(row[3])>0: 
                              clonotypes[triple]={}
                              clonotypes[triple]['clonotype-count']=float(row[3])
                              if vdjcount:
                                 clonotypes[triple]['vdj-count']=float(row[4])
            csv_file.close()
   else:
        with open(clonotype_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                    triple='%s %s %s'%(row[0],row[1],row[2])
                    if tripe in clonotypes:
                        clonotypes[triple]['clonotype-count']=clonotypes[triple]['clonotype-count']+float(row[3])
                        if vdjcount:
                           clonotypes[triple]['vdj-count']=clonotypes[triple]['vdj-count']+float(row[4])
                    else:
                        if float(row[3])>0:
                           clonotypes[triple]={}
                           clonotypes[triple]['clonotype-count']=float(row[3])
                           if vdjcount:
                              clonotypes[triple]['vdj-count']=float(row[4])
            csv_file.close()
   return clonotypes

'''----------------'''
'''read in vj only '''
'''----------------'''

def read_vj(clonotypes,clonotype_file,lowercdr3length,uppercdr3length,blacklist=False,vdjcount=False):
   '''
    This function will record the VJ frequency count.
   '''
   if  os.path.splitext(clonotype_file)[-1] in [ '.gz', '.gzip']:
        with gzip.open(clonotype_file,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                       triple='%s %s'%(row[0],row[1])
                       #if row[0] not in BADGERMLINES:
                       if triple in clonotypes:
                           clonotypes[triple]['clonotype-count']=clonotypes[triple]['clonotype-count']+float(row[3])
                           if vdjcount:
                               clonotypes[triple]['vdj-count']=clonotypes[triple]['vdj-count']+float(row[4])
                       else:
                           if float(row[3])>0: 
                              clonotypes[triple]={}
                              clonotypes[triple]['clonotype-count']=float(row[3])
                              if vdjcount:
                                 clonotypes[triple]['vdj-count']=float(row[4])
            csv_file.close()
   else:
        with open(clonotype_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                    triple='%s %s'%(row[0],row[1])
                    if tripe in clonotypes:
                        clonotypes[triple]['clonotype-count']=clonotypes[triple]['clonotype-count']+float(row[3])
                        if vdjcount:
                           clonotypes[triple]['vdj-count']=clonotypes[triple]['vdj-count']+float(row[4])
                    else:
                        if float(row[3])>0:
                           clonotypes[triple]={}
                           clonotypes[triple]['clonotype-count']=float(row[3])
                           if vdjcount:
                              clonotypes[triple]['vdj-count']=float(row[4])
            csv_file.close()
   return clonotypes

def pair_clonotype_entries(clonotypesA,clonotypesB):
    table=[]
    '''--------------------'''
    '''Create union of keys'''
    '''--------------------'''
    for element in list (set(clonotypesA.keys()) | set(clonotypesB.keys())):
        if element in clonotypesA and element in clonotypesB:
            table.append( (element,clonotypesA[element]['clonotype-count'],clonotypesB[element]['clonotype-count']) )
        elif element in clonotypesA:
            table.append( (element,clonotypesA[element]['clonotype-count'],float(0)) )
        elif element in clonotypesB:
            table.append( (element,float(0), clonotypesB[element]['clonotype-count']) )
        else:
            print 'Error with clonotype labels'
            sys.exit(1)
    return table

def species_count(clonotypes):
    return len(clonotypes)

def individual_count(clonotypes):
    count=float(0)
    for key in clonotypes.keys():
        count = clonotypes[key]['clonotype-count']+count
    return count


def morisitaHorn(data1, data2):

    if len(data1) != len(data2):
        raise ValueError("Error while calculating MorisitaHorn similarity index. The two input data lists must have the same length!\n")

    try: 
       x_norm=np.array(data1)/sum(np.array(data1))
       y_norm=np.array(data2)/sum(np.array(data2))
       sumXiYi = 0
       sumXiSqr = 0
       sumYiSqr = 0

       for x,y in zip(x_norm,y_norm): 
          sumXiYi += x*y
          sumXiSqr += x*x
          sumYiSqr += y*y

       return 2*sumXiYi / (sumXiSqr+sumYiSqr)
    except ZeroDivisionError:
       print("Cannot divide by zero")

def permutationTest(data1, data2,cutoff=0.90,count=1000):
   overlaps=[]
   for i in range(count):
      temp=[]
      x=[i for i in range(0,len(data1))]
      shuffle(x)
      for index in x:
          temp.append(data1[index])
      overlaps.append(morisitaHorn(temp,data2))
      del temp
   return float(len( sorted(i for i in overlaps if i >=cutoff))+1)/float(len(data1))

def scatter_colormap(masterlistA, masterlistB, colors):
    font = {'family': 'serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 16,
    }
    for i,j,color in zip(masterlistA,masterlistB,colors):
        plt.scatter(np.array(i)/sum(np.array(i)), np.array(j)/sum(np.array(j)),s=80,facecolors='none',edgecolors=color,marker='o')
    plt.plot([-1,0,1],[-1,0,1],'--',color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(10e-9, 1)
    plt.ylim(10e-9, 1)
    axes = plt.gca()
    #axes.set_yticks([10E-8,10E-7,10E-6, 10E-5, 10E-4,10E-3,10E-2,10E-1])
    #axes.set_xticks([10E-8,10E-7,10E-6, 10E-5, 10E-4,10E-3,10E-2,10E-1])
    axes.set_yticks([10E-8,10E-6,10E-4,10E-2])
    axes.set_xticks([10E-8,10E-6,10E-4,10E-2])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("MORISTA.pdf", bbox_inches='tight')
    plt.savefig("MORISTA.png", bbox_inches='tight')
    plt.show()

def main():

   CLI=argparse.ArgumentParser()
   CLI.add_argument(
     '--clonotype-files1',
     nargs="*",
     type=str,
     default=[],
   )
   CLI.add_argument(
     '--clonotype-files2',
     nargs="*",
     type=str,
     default=[],
   )
   CLI.add_argument(
     '--lower-bound-cdr3-length',
     type=int,
     default=2,
   )
   CLI.add_argument(
     '--upper-bound-cdr3-length',
     type=int,
     default=50,
   )
   CLI.add_argument("--vj-frequency-only",
                 action='store_true',
                 help='this will record the vj frequency only',
                 default=False)

   CLI.add_argument("--dump-ordered-table",
                 action='store_true',
                 help='dump ordered table showing corresponding clonotypes',
                 default=False)
   CLI.add_argument("--permutation-test",
                 action='store_true',
                 help='Perform a permutation test',
                 default=False)
   CLI.add_argument("--randomization-count",
                 type=int,
                 help='how many randomizations to do',
                 default=100)

   CLI.add_argument("--create-plot",
                 action='store_true',
                 help='generate a plot',
                 default=False)
   CLI.add_argument(
     '--hex-colors',
     nargs="*",
     type=str,
     default=[],
   )

   args=CLI.parse_args()

   masterlistA=[]
   masterlistB=[]
   clonotypesA={}
   clonotypesB={}
   if args.vj_frequency_only==True:
      for files1,files2 in zip(args.clonotype_files1,args.clonotype_files2):
         clonotypesA=read_vj(clonotypesA,files1,args.lower_bound_cdr3_length,args.upper_bound_cdr3_length,blacklist=False,vdjcount=False)
         clonotypesB=read_vj(clonotypesB,files2,args.lower_bound_cdr3_length,args.upper_bound_cdr3_length,blacklist=False,vdjcount=False)

         '''------------------------------'''
         ''' Store the intermediate table '''
         '''------------------------------'''
         if args.dump_ordered_table:
            sorted_table=[]
            sorted_table=sorted(pair_clonotype_entries(clonotypesA,clonotypesB), key=lambda item: item[1],reverse=True)
            for element in sorted_table:
               print '%s %f %f'%(element[0],element[1],element[2])
            del sorted_table

         '''------------------------------'''
         ''' Store the intermediate table '''
         '''------------------------------'''
         data1=[]
         data2=[]
         for element in pair_clonotype_entries(clonotypesA,clonotypesB):
            data1.append(element[1])
            data2.append(element[2])

         masterlistA.append(data1)
         masterlistB.append(data2)

         '''---------------------------'''
         ''' Compute the Morista-Index '''
         '''---------------------------'''
         m=morisitaHorn(data1, data2)
         if args.permutation_test:
               print '%s %s %d %d %d %d %f %f(%d)' % (files1,files2,int(species_count(clonotypesA)),\
                                                      int(individual_count(clonotypesA)),\
                                                      int(species_count(clonotypesB)),\
                                                      int(individual_count(clonotypesB)),m,
                                                      permutationTest(data1, data2,cutoff=m,count=args.randomization_count),args.randomization_count)
         else:
               print '%s %s %d %d %d %d %f' % (files1,files2,int(species_count(clonotypesA)),\
                                                     int(individual_count(clonotypesA)),\
                                                     int(species_count(clonotypesB)),\
                                                     int(individual_count(clonotypesB)),m)
         '''---------------'''
         ''' Clean up here '''
         '''---------------'''
         del data1
         del data2
         clonotypesA.clear()
         clonotypesB.clear()

   else:
      for files1,files2 in zip(args.clonotype_files1,args.clonotype_files2):
         clonotypesA=read_v3j_clonotypes(clonotypesA,files1,args.lower_bound_cdr3_length,args.upper_bound_cdr3_length,blacklist=False,vdjcount=False)
         clonotypesB=read_v3j_clonotypes(clonotypesB,files2,args.lower_bound_cdr3_length,args.upper_bound_cdr3_length,blacklist=False,vdjcount=False)
          
         '''------------------------------'''
         ''' Store the intermediate table '''
         '''------------------------------'''
         if args.dump_ordered_table:
            sorted_table=[]
            sorted_table=sorted(pair_clonotype_entries(clonotypesA,clonotypesB), key=lambda item: item[1],reverse=True)
            for element in sorted_table:
               print '%s %f %f'%(element[0],element[1],element[2])
            del sorted_table

         '''---------------------------------------'''
         ''' Store arrays to comptue Morisita-Horn '''
         '''---------------------------------------'''
         data1=[]
         data2=[]
         for element in pair_clonotype_entries(clonotypesA,clonotypesB):
            data1.append(element[1])
            data2.append(element[2])

         '''---------------------------'''
         ''' Compute the Morista-Index '''
         '''---------------------------'''
         m=morisitaHorn(data1, data2)

         if args.permutation_test:
               print '%s %s %d %d %d %d %f %f(%d)' % (files1,files2,int(species_count(clonotypesA)),\
                                                      int(individual_count(clonotypesA)),\
                                                      int(species_count(clonotypesB)),\
                                                      int(individual_count(clonotypesB)),m,
                                                      permutationTest(data1, data2,cutoff=m,count=args.randomization_count),args.randomization_count)
         else:
               print '%s %s %d %d %d %d %f' % (files1,files2,int(species_count(clonotypesA)),\
                                                     int(individual_count(clonotypesA)),\
                                                     int(species_count(clonotypesB)),\
                                                     int(individual_count(clonotypesB)),m)
         '''---------------'''
         ''' Clean up here '''
         '''---------------'''
         del data1
         del data2
         clonotypesA.clear()
         clonotypesB.clear()

   if args.create_plot==True:
       #print np.asarray(masterlistA)
       #print np.asarray(masterlistB)
       scatter_colormap(masterlistA, masterlistB, args.hex_colors)
if __name__ == '__main__':
    main()
