
'''-------------------------------'''
''' Script for creating a heatmap '''
'''-------------------------------'''
'''C.SOTO'''

import os
import sys
import re
import gzip
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

'''-----------------------------------'''
''' Use Seaborn to create the heatmap '''
'''-----------------------------------'''
import seaborn as sns



'''---------------'''
''' Constant data '''
'''---------------'''
JCHAINS={ 'IGHJ' :{'IGHJ1':0,'IGHJ2':1,'IGHJ3':2,'IGHJ4':3,'IGHJ5':4,'IGHJ6':5},
          'TRBJ': {'TRBJ1-1':0, 'TRBJ1-2':1,'TRBJ1-3':2,'TRBJ1-4':3,'TRBJ1-5':4,'TRBJ1-6':5,'TRBJ2-1':6, 'TRBJ2-2':7,'TRBJ2-3':8,'TRBJ2-4':9,'TRBJ2-5':10,'TRBJ2-6':11,'TRBJ2-7':12 }
          }


def read_vj(clonotype_file,lowercdr3length,uppercdr3length,chaintype='IGHV',jlist={'IGHJ1':0, 'IGHJ2':1, 'IGHJ3':2,'IGHJ4':3,'IGHJ5':4,'IGHJ6':5},blacklist=['IGHV4/OR15-8']):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''
   clonotypes={}

   if  os.path.splitext(clonotype_file)[-1] in [ '.gz', '.gzip']:
        with gzip.open(clonotype_file,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
               if row[0] not in blacklist and chaintype in row[0]:
                   if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                       vgerm='%s'%(row[0])
                       if vgerm in clonotypes:
                           clonotypes[vgerm][jlist[row[1]]]=clonotypes[vgerm][jlist[row[1]]]+1
                       else:
                           clonotypes[vgerm]=[0]*len(jlist)
                           clonotypes[vgerm][jlist[row[1]]]=1
            csv_file.close()
   else:
        with open(clonotype_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ') 
            for row in csv_reader:
                if row[0] not in blacklist and chaintype in row[0]:
                   if len(row[3]) > lowercdr3length and len(row[3]) < uppercdr3length:
                       vgerm='%s'%(row[0])
                       if vgerm in clonotypes:
                           clonotypes[vgerm][jlist[row[1]]]=clonotypes[vgerm][jlist[row[1]]]+1
                       else:
                           clonotypes[vgerm]=[0]*len(jlist)
                           clonotypes[vgerm][jlist[row[1]]]=1
            csv_file.close()
   return clonotypes

def sum_frequency_heatmap(vj_table):
    '''
    Normalize the heatmap using the sum of the frequencies
    '''
    total_sum=0
    for k in vj_table.keys():
       total_sum=total_sum+sum(vj_table[k])
    return total_sum


def sum_frequency_heatmaps(vj_tables):
    '''
    Normalize the heatmap using the sum of the frequencies
    '''
    total_sum=0
    for table in vj_tables:
       for k in table.keys():
          total_sum=total_sum+sum(table[k])
    return total_sum


def normalize_heatmap(vj_table):
    '''
    Normalize the heatmap using the sum of the frequencies
    '''
    vj_usage={}
    total_sum=0
    for k in vj_table.keys():
       total_sum=total_sum+sum(vj_table[k])
    return total_sum


def remove_sparse_rows_heatmap(vj_table,cut=10):
    '''
    Normalize the heatmap using the sum of the frequencies
    '''
    vj_usage={}
    for k in vj_table.keys():
       if sum(vj_table[k]) >= cut:
            vj_usage[k]=[]
            vj_usage[k]=vj_table[k]
    return vj_usage


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

'''-----------------------'''
''' Commandline arguments '''
'''-----------------------'''

CLI=argparse.ArgumentParser()
CLI.add_argument(
     '--clonotype-files',
     nargs="*",
     type=str,
     help='clonotype files containing Vgermline Jgermline CDR3aa SomaticCount',
     default=[],
)
CLI.add_argument(
     '--output-csv-files',
     nargs="*",
     help='CSV file for VJ frequency counts',
     type=str,
     default=[],
)
CLI.add_argument("--dump-csv",
                  action='store_true',
                  help='will read in data files',
                  default=False)
CLI.add_argument(
     '--output-aggregate-csv-file',
      help='aggregate csv file containing normalized data',
     type=str,
     default='aggregate-heatmap.csv',
)
CLI.add_argument(
     '--v-chain-type',
     help='chain types VH or VB',
     choices=['IGHV', 'TRBV'],
     default='TRBV'
)
CLI.add_argument(
     '--j-chain-type',
     help='chain types VH or VB',
     choices=['IGHJ', 'TRBJ'],
     default='TRBJ'
)

'''--------------------------'''
''' Heatmap specific options '''
'''--------------------------'''
CLI.add_argument(
     '--heatmap-vmin-vmax',
     nargs="*",
     type=float,
     help='min and max values for colorbar',
     default=[-1.0,1.0],
)
CLI.add_argument(
     '--heatmap-figsize',
     nargs="*",
     type=float,
     help='figure size in inches (width, height)',
     default=[10.5,10.5],
)

CLI.add_argument(
     '--heatmap-colobar-ticks',
     nargs="*",
     type=float,
     help='tick locations on colorbar',
     default=[-1.0, -0.5, 0, 0.5,1.0],
)
CLI.add_argument(
     '--heatmap-output-pdf',
     type=str,
      help='output filename for PDF image',
     default='Plot.pdf',
)
CLI.add_argument(
     '--heatmap-output-png',
     type=str,
     help='output filename for PNG image',
     default='Plot.png',
)
CLI.add_argument(
     '--heatmap-diverging-cmap',
     choices=['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'],
     help='possible diverging heatmaps',
     default='coolwarm',
)
CLI.add_argument("--heatmap-display",
                  action='store_true',
                  help='display plot on screen',
                  default=False)


'''--------------------------'''
''' Parse the arguments here '''
'''--------------------------'''
args=CLI.parse_args()

'''--------------------'''
''' Initialize storage '''
'''--------------------'''
datasets=[]
orders=[]
counter=0
for i in range(0,len(args.clonotype_files)):
    datasets.append({})

'''--------------------------------------------------------'''
''' Read in files and generate statisitcs on all data sets '''
'''--------------------------------------------------------'''
counter=0
for clonotypefile in args.clonotype_files:
    '''--------------------------------------------'''
    ''' Read in a list of VJ heatmaps for analysis '''
    '''--------------------------------------------'''
    vj_usage=read_vj(clonotypefile,2,40,chaintype=args.v_chain_type,jlist=JCHAINS[args.j_chain_type], blacklist=['TRBV29/OR9-2','TRBV20/OR9-2'])
    ordered=vj_usage.keys()
    ordered.sort(key=natural_keys)

    '''----------------------------------'''
    ''' Remove sparse rows from heatmap  '''
    '''----------------------------------'''
    clean=remove_sparse_rows_heatmap(vj_usage,cut=10)

    '''--------------------------------'''
    '''Dump out the normalized heatmap '''
    '''--------------------------------'''
    if args.dump_csv and len(args.output_csv_files)==len(args.clonotype_files):
       outfile=gzip.open(args.output_csv_files[counter],'wb')
       for k in ordered:
          if k in clean:
              outfile.write('%s,%s\n'%(k,','.join(map(str,clean[k]))))
       outfile.close()

    '''------------------------------------------------------'''
    ''' Cycle through the clean list and generate data array '''
    '''------------------------------------------------------'''
    non_flattened_list=[]
    for k in ordered:
        if k in clean:
           non_flattened_list.append(vj_usage[k])
    counter=counter+1
    clean.clear()
    del ordered
    vj_usage.clear()

'''------------------------'''
''' Compute stats for this '''
'''------------------------'''
counts=np.array([item for sublist in non_flattened_list for item in sublist])
mean  = np.mean(counts,axis=0)
total = np.sum(counts,axis=0)
std   = np.std(counts,axis=0)

'''---------------'''
'''  Clean up     '''
'''---------------'''
del non_flattened_list

counter=0
for clonotypefile in args.clonotype_files:
    '''--------------------------------------------'''
    ''' Read in a list of VJ heatmaps for analysis '''
    '''--------------------------------------------'''
    vj_usage=read_vj(clonotypefile,2,40,chaintype=args.v_chain_type,jlist=JCHAINS[args.j_chain_type], blacklist=['IGHV4/OR15-8'])
    ordered=vj_usage.keys()
    ordered.sort(key=natural_keys)

    '''----------------------------------'''
    ''' Remove sparse rows from heatmap  '''
    '''----------------------------------'''
    clean=remove_sparse_rows_heatmap(vj_usage,cut=10)

    '''--------------------------------------------'''
    ''' Create new heatmap with normalized values  '''
    '''--------------------------------------------'''
    for k in ordered:
        if k in clean:
            if k in datasets[counter]:
               pass
            else:
               datasets[counter][k]=[]
               datasets[counter][k] = [(x -mean)/std for x in map(float,vj_usage[k])]
    '''-----------------------------------------'''
    ''' Store the data sets for later analysis  '''
    '''-----------------------------------------'''
    orders.append(set(clean.keys()))
    
    '''---------------'''
    '''  Clean up     '''
    '''---------------'''
    vj_usage.clear()
    clean.clear()
    
    '''----------------------------'''
    ''' Increment counter variable '''
    '''----------------------------'''
    counter=counter+1


'''----------------------------------------'''
''' Generate a complete set of unique keys '''
'''-----------------------------------------'''
union=list(set.union(*orders))
union.sort(key=natural_keys) 

'''----------------------------'''
''' Append all tables together '''
'''----------------------------'''
if args.v_chain_type=='IGHV':
   size=6
elif args.v_chain_type=='TRBV':
   size=13
else:
   print 'Lights have not been worked out'
   sys.exit(1)

counter=0
temp={}
for counter in range(0,len(args.clonotype_files)):
    for k in union:
        if k in datasets[counter]:
            if k in temp:
               temp[k].extend(datasets[counter][k])
               temp[k].append(-99999.99)
               temp[k].append(-99999.99)
            else:
               temp[k]=[]
               temp[k].extend(datasets[counter][k])
               temp[k].append(-99999.99)
               temp[k].append(-99999.99)
        else: 
           if k in temp: 
               temp[k].extend([0]*size)
               temp[k].append(-99999.99)
               temp[k].append(-99999.99)
           else:
               temp[k]=[]
               temp[k]=[0]*size
               temp[k].append(-99999.99)
               temp[k].append(-99999.99)

    counter=counter+1  

'''---------------------------'''
''' Dump out combined heatmaps'''
'''---------------------------'''
outfile=gzip.open(args.output_aggregate_csv_file,'wb')
for k in union:
    outfile.write('%s,%s\n'%(k,','.join(map(str,temp[k]))))
outfile.close()

fig, ax = plt.subplots()
fig.set_size_inches(args.heatmap_figsize[0], args.heatmap_figsize[1])
df = pd.read_csv(args.output_aggregate_csv_file,header=None)
sns.set()
sns.set_style('darkgrid', {'axes.linewidth': 2, 'axes.edgecolor':'black'})
sns.heatmap(df.as_matrix(columns=np.arange(1,len(df.columns))), 
            ax=ax,linewidths=1.5,
            yticklabels=[item for sublist in df.as_matrix(columns=[0]) for item in sublist],
            xticklabels=[],
            square=True,
            cmap=args.heatmap_diverging_cmap,
            vmin=args.heatmap_vmin_vmax[0],vmax=args.heatmap_vmin_vmax[1],
            mask=df.as_matrix(columns=np.arange(1,len(df.columns))) <=-9999.0)

'''------------------------'''
''' Set heatmap attributes '''
'''------------------------'''
ax.set_yticklabels(ax.get_yticklabels(),rotation=0,fontsize=10)

'''-------------------------'''
''' Set colorbar atrributes '''
'''-------------------------'''
cbar = ax.collections[0].colorbar
cbar.set_ticks(args.heatmap_colobar_ticks)
#cbar.set_ticklabels(['0.0','1.0', '2.0', '3.0', '4.0','5.0','6.0'])
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_aspect(75)
cbar.ax.set_position((0.8, 0.175, 0.175, 0.64))

'''------------'''
''' Save files '''
'''------------'''
plt.savefig(args.heatmap_output_pdf)
plt.savefig(args.heatmap_output_png)

if args.heatmap_display:
   plt.show()
