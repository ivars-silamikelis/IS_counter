import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import Levenshtein as lev
import sys

prefix=sys.argv[1]
for end in ["starts","ends"]:
  cont=[]
  print "processing "+prefix+"."+end+".fa"
  for rec in SeqIO.parse(prefix+"."+end+".fa","fasta"):
      cont.append(rec)
  distances=[]
  cont_string=[]
  if end == "ends":
    cont_string=[x.seq.tostring() for x in cont]
  elif end == "starts":		#If reads belong to start, the these reads are reversed to do pairwise comparison correctly
    cont_string=[x.seq.tostring()[::-1] for x in cont]
  cont_sorted=sorted(cont_string)
  #Pairwise calculation of levenshtein distances 
  #Note that comparison is done on equal length sequences 
  for rec_q in cont_sorted:
    small_dist=[]
    for rec_s in cont_sorted:
      if (len(rec_s)>=len(rec_q)):
        small_dist.append(lev.distance(rec_q,rec_s[0:len(rec_q)]))
      else:
        small_dist.append(lev.distance(rec_s,rec_q[0:len(rec_s)]))

    distances.append(small_dist)
  mat=pd.DataFrame(np.array(distances))
  plt.matshow(mat)
  plt.title("Pairwise Levenshtein distances for "+prefix+" "+end) 
  plt.savefig(prefix+"."+end+".png")
