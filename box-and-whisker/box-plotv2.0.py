import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import csv
import gzip

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
                       #if row[0] not in BADGERMLINES:
                       if clonotypes.has_key(triple):
                           clonotypes[triple]['clonotype-count']=clonotypes[triple]['clonotype-count']+int(row[3])
                       else:
                           clonotypes[triple]={}
                           clonotypes[triple]['clonotype-count']=int(row[3])
            csv_file.close()
   else:
        with open(clonotype_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if len(row[2]) > lowercdr3length and len(row[2]) < uppercdr3length:
                    triple='%s %s %s'%(row[0],row[1],row[2])
                    if clonotypes.has_key(triple):
                        clonotypes[triple]['clonotype-count']=clonotypes[triple]['clonotype-count']+int(row[3])
                    else:
                        clonotypes[triple]={}
                        clonotypes[triple]['clonotype-count']=int(row[3])
            csv_file.close()
   return clonotypes

def get_unique_cdr3s(vj_clonotypes):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''
   unique_cdr3={}
   for key in vj_clonotypes.keys():
       cdr3=key.split(' ')[-1]
       if cdr3 in unique_cdr3:
           unique_cdr3[cdr3]=unique_cdr3[cdr3]+1
       else:
           unique_cdr3[cdr3]=0
           unique_cdr3[cdr3]=unique_cdr3[cdr3]+1
   return unique_cdr3


CLI=argparse.ArgumentParser()
CLI.add_argument("--clonotype-file",
                  action='store_true',
                  help='will read in clonotype files',
                  default=False)

CLI.add_argument("--data-file",
                  action='store_true',
                  help='will read in data files',
                  default=False)
CLI.add_argument(
     '--files',
     nargs="*",
     type=str,
     default=[],
)
CLI.add_argument(
     '--hex-colors',
     nargs="*",
     type=str,
     default=[],
)
CLI.add_argument(
     '--y-tick-markers',
     nargs="*",
     type=int,
     default=[],
)
CLI.add_argument(
     '--y-limits',
     nargs="*",
     type=int,
     default=[],
)
CLI.add_argument(
     '--font',
     type=str,
     default='arial',
)
CLI.add_argument(
     '--font-size',
     type=int,
     default=28,
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
     '--summary-statistics-outfile',
     type=str,
     default='summary-statistics.dat.gz',
)
CLI.add_argument(
     '--png-file',
     type=str,
     default='BoxPlot.png',
)
CLI.add_argument(
     '--pdf-file',
     type=str,
     default='BoxPlot.pdf',
)

CLI.add_argument("--scatter",
                  action='store_true',
                  help='include scattered points too',
                  default=False)

CLI.add_argument("--scatter-point-size",
                  type=int,
                  default=5)

CLI.add_argument("--scatter-point-alpha",
                  type=float,
                  default=0.3)
CLI.add_argument("--scatter-point-marker",
                  type=str,
                  default='.')


args=CLI.parse_args()
if args.clonotype_file==True and args.data_file==True:
    print('You can only have one file type')
    sys.exit(1)

data=[]
colors=[]
outfile=gzip.open(args.summary_statistics_outfile,'wb')

for inputfile,color in zip(args.files, args.hex_colors):
   '''-----------------'''
   ''' Store the color '''
   '''-----------------'''
   colors.append(color)
   '''------------------------'''
   ''' Read in the first file '''
   '''------------------------'''
   if args.clonotype_file:
      vj_clonotypes={}
      cdr3_aa={}
      vj_clonotypes=read_vj_clonotypes(vj_clonotypes,inputfile,args.lower_bound,args.upper_bound)

      '''----------------------------'''
      ''' Grab only the CDR3 lengths '''
      '''----------------------------'''
      cdr3_aa=get_unique_cdr3s(vj_clonotypes)

      '''---------------------------------'''
      ''' Create an array of CDR3 lengths '''
      '''---------------------------------'''
      cdr3=np.array(map(len, cdr3_aa.keys()))
      vj_clonotypes.clear()
      cdr3_aa.clear()

      #print np.ptp(cdr3,axis=0)
      #print np.percentile(cdr3,25,axis=0)
      #print np.mean(cdr3,axis=0)
      #print np.median(cdr3,axis=0)
      #print np.std(cdr3,axis=0)
   elif args.data_file:
      cdr3=np.loadtxt(inputfile,)

   data.append(cdr3)
   outfile.write('%s %d %f %f %f\n'%(inputfile,np.size(cdr3),np.mean(cdr3,axis=0),np.median(cdr3,axis=0),np.std(cdr3,axis=0)))
   del cdr3
outfile.close()
'''----------------------------------------------------------------------'''
''' I have not generatlized this script. At the moment it requires five  '''
''' data sets. Each data set should be a file with a single column       '''
''' that contains the CDR3 lengths.                                      '''
'''----------------------------------------------------------------------'''


# basic plot
axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
'''---------------------------------------'''
''' Set the proper spacing for the xticks '''
'''---------------------------------------'''
axes.set_xticks(axes.get_xticks()[::5])

'''-----------------------------'''
''' Set the min-max CDR3 length '''
'''-----------------------------'''
axes.set_ylim(args.y_limits[0],args.y_limits[1])

'''------------------------------------'''
''' Set the colorbars for the BOX plot '''
'''------------------------------------'''
bplot1=plt.boxplot(data,0, '',patch_artist=True)
#colors = ['#009999', '#ff9933', '#ff44ab', '#ad9de7', '#abd9ac']
i=0
for box in bplot1['boxes']:
   box.set(facecolor = colors[i])
   box.set(linewidth=2)
   box.set(color=colors[i])
   i=i+1

'''-------------------------------------'''
''' Set the median bar for the BOX plot '''
'''-------------------------------------'''
i=0
for median in bplot1['medians']:
   median.set(color='black')
   median.set(linewidth=2)
   i=i+1


'''---------------------------------------------'''
''' Color the whiskers and caps in the BOX plot '''
'''---------------------------------------------'''
i=0
#colors = ['#009999','#009999','#ff9933','#ff9933','#ff44ab','#ff44ab','#ad9de7','#ad9de7','#abd9ac','#abd9ac']
#colors = [val for val in colors for _ in (0, 1)]
for whisker,cap in zip(bplot1['whiskers'],bplot1['caps']):
   whisker.set(linewidth=2)
   #whisker.set(color=colors[i])
   whisker.set(color=[val for val in colors for _ in (0, 1)][i])
   cap.set(linewidth=2)
   #cap.set(color=colors[i])
   cap.set(color=[val for val in colors for _ in (0, 1)][i])
   i=i+1

'''------------------------------------------------'''
''' This will set the font size in the tick labels '''
'''------------------------------------------------'''
for label in (axes.get_xticklabels() + axes.get_yticklabels()):
    label.set_fontname(args.font)
    label.set_fontsize(args.font_size)

'''-------------------'''
''' What does this do '''
'''-------------------'''
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

#axes.set_yticks([0, 5, 10, 15, 20])
axes.set_yticks(args.y_tick_markers)

if args.scatter:
   indeces = np.arange(1,len(data)+1) 
   #marker_style = dict(color='black',  marker='o')
   for i,d1 in enumerate(data):
       marker_style = dict(color=colors[i],  marker=args.scatter_point_marker)
       plt.scatter(np.random.normal(indeces[i], 0.04, size=len(data[i])),data[i],s=args.scatter_point_size,alpha=args.scatter_point_alpha,zorder=10, facecolors='none',linewidth=1.5,**marker_style)

#if args.display:
plt.savefig(args.png_file, bbox_inches='tight')
plt.savefig(args.pdf_file, bbox_inches='tight')
#plt.show()
