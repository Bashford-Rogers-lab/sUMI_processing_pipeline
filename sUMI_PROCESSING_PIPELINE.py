# sUMI_PROCESSING_PIPELINE developed by Rachael Bashford-Rogers (2020)
# at the University of Oxford and Univeristy of Cambridge
# E-mail: rbr1@well.ox.ac.uk

# If you find the methods in sUMI_PROCESSING_PIPELINE, please cite the following reference:
# Ahmed, R. et al. 2020 ()

# Copyright (C) 2020  Dr Rachael Bashford-Rogers

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/python
import math
import sys
from collections import defaultdict
import os
import commands
from operator import itemgetter
import re

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break	
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def fuzzy_substring(needle, haystack):
    """Calculates the fuzzy match of needle in haystack,
    using a modified version of the Levenshtein distance
    algorithm.
    The function is modified from the levenshtein function
    in the bktree module by Adam Hupp"""
    m, n = len(needle), len(haystack)
    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m
    row1 = [0] * (n+1)
    for i in range(0,m):
        row2 = [i+1]
        for j in range(0,n):
            cost = ( needle[i] != haystack[j] )

            row2.append( min(row1[j+1]+1, # deletion
                               row2[j]+1, #insertion
                               row1[j]+cost) #substitution
                           )
        row1 = row2
    return min(row1)

def Get_barcodes(barcode, barcode_file,rc):
  barcode = barcode.split(",")
  fh=open(barcode_file,"r")
  barcodes,allbarcodes,barcode_stem,barcode_only,barcode_only_rc = [],[],[],[],[]
  rc_barcodes = []
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      allbarcodes.append(l[1])
      if(l[0] in barcode):
        barcodes.append(l[1])
        rc_barcodes.append(Reverse_comp(l[1], rc))
        barcode_only.append(l[2])
        barcode_only.append(Reverse_comp(l[2], rc))
      elif(l[0]=="STEM"):barcode_stem.append([l[1],l[1][0:len(l[1])/2], l[1][len(l[1])/2:len(l[1])]])
  fh.close()
  return(barcodes,allbarcodes,barcode_stem, rc_barcodes,barcode_only,barcode_only_rc)

def Init_rc():
  b1 = ["A","T","G","C","N","."]
  b2 = ["T","A","C","G","N","."]
  rc = {}
  for i in range(0,6):
    rc[b1[i]]=b2[i]
  return(rc)

def Reverse_comp(seq, rc):
  s = ''
  l=len(seq)
  for i in range(0,l):
    j=l-i-1
    if(seq[j] not in rc):
      print seq
    else:
      s=s+rc[seq[j]]
  return(s)

def Get_sequences(file):
  command = "gunzip "+file
  commands.getoutput(command)
  fh = open (file, "r")
  seqs = {}
  for header,sequence in fasta_iterator(fh):
    header = header.replace(":","").split()[0].replace("-","")
    header = header.split("#")[0].split("/")[0]
    seqs[header]=sequence
  fh.close()
  return(seqs)

def Split_reads_and_join(Seq_file1, Seq_file2, Tmp_file,id,gene_specific_primer,dir,sample,sUMI_sample_barcode_id, sUMI_sample_barcode_file):
  Split_reads_and_join_pre(Seq_file1, Seq_file2, Tmp_file,id,gene_specific_primer,dir,sample,sUMI_sample_barcode_id, sUMI_sample_barcode_file)
  Check_split(Tmp_file, sUMI_sample_barcode_id,sUMI_sample_barcode_file)
  return()

def Check_split(outfile, sUMI_sample_barcode_id,UMI_sample_barcode_file):
  rc = Init_rc()
  barcodes,allbarcodes,barcode_stem,rc_barcodes,barcode_only,barcode_only_rc = Get_barcodes(sUMI_sample_barcode_id,UMI_sample_barcode_file,rc)
  allbarcodes_rc = []
  for b in allbarcodes:
    allbarcodes_rc.append(Reverse_comp(b,rc))
  out,ind = '',0
  outfile1 = outfile.replace(".fasta","_CHECKED.fasta")
  fh=open(outfile1,"w")
  fh.close()
  fh = open(outfile,"r")
  seqs = {}
  for header,sequence in fasta_iterator(fh):
    seqs[header] = sequence
  fh.close()
  for header in seqs:
    sequence = seqs[header]
    score = []
    for i in range(0,len(allbarcodes)):
      p1 = len(Get_match(allbarcodes[i], sequence))
      p2 = len(Get_match(allbarcodes_rc[i], sequence))
      score.append(p1+p2)
    if(sum(score)<=1):
      out=out+">"+header+"\n"+sequence+"\n"
      ind = ind+1
      if(ind>500):
        Write_out(out, outfile1)
        out, ind = '',0
  fh.close()
  Write_out(out, outfile1)
  return()

def Join_reads(s1, s2,rc,length):
  seq=''
  (s2)=Reverse_comp(s2,rc)
  (l1,l2)=(len(s1),len(s2))
  failed = 1
  for i in range(0,100):
    ind=(i*5)+5
    if(ind<(l1-length)):
      (s11,s22,p,overlap)=Trim(s1,s2,l1,l2,ind,length)
      if(p==1):
        seq = s1[0:s1.index(overlap)]+s2[(s2.index(overlap)):l2].lower()
        if(len(seq)>120):
          failed = 0
          break
  return(seq,failed)

def Trim(s1,s2,l1,l2,indent,length):
  (i,p)=(indent,0)
  sample1=s1[i:i+length]
  index=s2.find(sample1)
  s1a, s2a = s1, s2
  if (index!=-1):
    if(index>i):
      s2a=s2a[index-i:l2]
    else:
      s1a=s1a[i-index:l1]
    min_len=min([len(s1),len(s2)])
    s1a=s1a[0:min_len]
    s2a=s2a[0:min_len]
    p=1
  return (s1a, s2a, p,sample1)

def Get_match(primer, seq):
  loc = []
  if (seq.count(primer)!=0):
    for m in re.finditer(primer, seq):
      loc = loc+[m.start()]
  return(loc)

def Split_reads_and_join_pre(file1, file2, outfile,id,primer_file,dir,sample,sUMI_sample_barcode_id, sUMI_sample_barcode_file):
  fh=open(outfile,"w")
  fh.close()
  rc = Init_rc()
  file1 = dir+"FASTQ_FILES/Sequences_"+sample+"_1.fasta"
  file2 = dir+"FASTQ_FILES/Sequences_"+sample+"_2.fasta"
  barcodes,allbarcodes,barcode_stem,rc_barcodes,barcode_only,barcode_only_rc = Get_barcodes(sUMI_sample_barcode_id, sUMI_sample_barcode_file,rc)
  (rc)= Init_rc()
  (seqs1)=Get_sequences(file1)
  (seqs2)=Get_sequences(file2)
  print "number of raw sequences:\tRead1:",len(seqs1),"\tRead2:", len(seqs2)
  out,ind = '',0
  fh = open (outfile, "w")
  fh.close()
  tot,joining,total,length,primer_match=0,0,0,20,0
  count_joined, count_barcode, count_barcode_saved = 0,0,0
  for id in seqs1:
    if (id in seqs2):
      total=total+1
      s1,s2 = seqs1[id], seqs2[id]
      seq,failed =Join_reads(seqs1[id], seqs2[id],rc,length)
      if(failed==0 and len(seq)>150):
        seq,joining = seq.upper(),joining+1
        count_joined = count_joined+1
        for primer in barcodes:
          p = Get_match(primer, seq)
          if(len(p)==0):
            seq = Reverse_comp(seq,rc)
            p = Get_match(primer, seq)
            break
          else:break
        if(len(p)!=0):
          if(p[0]<len(seq)/2):
            count_barcode = count_barcode+1
            out=out+">"+id+"\n"+seq+"\n"
            ind,primer_match = ind+1,primer_match+1
            if(ind>200):
              Write_out(out, outfile)
              out,ind = '',0
        else:
          for prim in barcode_stem:
            p = Get_match(prim[0], seq)
            if(len(p)==0):
              seq = Reverse_comp(seq,rc)
              p = Get_match(prim[0], seq)
            if(len(p)!=0):
              start = p[0]+len(prim[0])-3
              start = p[0]+len(prim[0])-6
              for i in range(0,len(barcodes)):
                bc_halves = [barcodes[i][2:len(barcodes[i])/2], barcodes[i][len(barcodes[i])/2:len(barcodes[i])-2]]
                bc_shift = [2,len(barcodes[i])/2]
                test = seq[start: start+len(barcodes[i])+10]
                for j in range(0,len(bc_halves)):
                  if(test.count(bc_halves[j])!=0):
                    start = barcodes[i].index(barcode_only[i])
                    p=test.index(bc_halves[j])-bc_shift[j]
                    sub_test = test[p+start:p+len(barcode_only[i])+start]
                    a = fuzzy_substring(sub_test,barcode_only[i])
                    if(a<=1):
                      count_barcode, count_barcode_saved = count_barcode+1, count_barcode_saved+1
                      out=out+">"+id+"\n"+seq+"\n"
                      ind,primer_match = ind+1,primer_match+1
                      if(ind>200):
                        Write_out(out, outfile)
                        out,ind = '',0
  Write_out(out, outfile)
  del out,ind
  print id+"\tTotal sequences:",total,"\tNumber of joined sequences:",joining, "Number of barcodes matched:",primer_match,"% retained:",primer_match*100.0/total, "IDs saved",count_barcode_saved
  return()

def Check_barcodes_MALBAC(primer_tag_file_count,primer_tag_file,Fail_file,Output_trim,threshold,cd_hit_directory,tmp_file):
  outs = ['','#ID\tnumber_of_reads\ttotal_reads_with_BC\tJ_barcode\tV_barcode\tbp_mismatches_to_consensus\tBC_accepted\tconsensus\tsequence\n']
  files = [Output_trim,primer_tag_file]
  for i in range(0,len(files)):
    fh=open(files[i],"w")
    fh.write(outs[i])
    fh.close()
  fh=open(primer_tag_file_count,"r")
  seqs = Tree()
  for l in fh:
    if (l[0]!="#"):
      l=l.strip().split()
      v_tag, j_tag, sequence, header = l[2],l[1],l[3],l[0]
      seqs[j_tag+"\t"+v_tag][sequence][header].value = 1
  fh.close()
  print len(seqs), "Unique tags"
  min_depth = (1.0/(1-threshold))
  total_tags,ind,passed_seqs_total = 0,0,0
  fail_less_than_threshold = 0
  outp=''
  for t1 in seqs:
    t,u_seq, u_freq,u_header = 0,[],[],[]
    for s in seqs[t1]:
      f = 0 
      for h in seqs[t1][s]:
        f = f+int(h.split(":")[1])
      total_tags,u_seq, u_freq, u_header, ind =total_tags+f, u_seq+[s], u_freq+[f], u_header+[h], ind+1
    f = sum(u_freq)
    h = h.split(":")[0].split("__")[0]+"__"+str(sum(u_freq))
    if(sum(u_freq)>20):print "BC group size:\t",len(u_freq),t1.split() 
    if(len(u_freq)==1): 
      passed_seqs_total = passed_seqs_total+f
      outp=outp+h+"\t"+str(f)+"\t"+str(f)+"\t"+t1+"\t0\tYES\t"+s+"\t"+s+"\n"
    elif(len(u_freq)<500):
      #print "BC group size:\t",len(u_freq),t1.split(), u_freq#, u_seq
      if(max(u_freq)>sum(u_freq)*threshold):
        nz = [i for i in range(len(u_freq)) if u_freq[i]!=max(u_freq)]
        passed_seqs_total = passed_seqs_total+f
        consensus = u_seq[nz[0]]
        outp = outp+"\t".join(map(str, [h, len(u_freq), sum(u_freq),t1,0,"YES",consensus, consensus ]))+"\n"
      elif(len(u_freq)>15): ## clustering first then alignment
        out_cluster,ids = '',{}
        for i in range(0,len(u_seq)):
          out_cluster = out_cluster+">"+u_header[i]+"\n"+u_seq[i]+"\n"
          ids[u_header[i]] = [u_seq[i], u_freq[i]]
        consensus,pass_consensus = Get_consensus_sequence_large(out_cluster, Fail_file,len(u_seq[i]), ids,sum(u_freq),threshold,tmp_file,cd_hit_directory)
        if(consensus!='' and pass_consensus==1):
          outp = outp+"\t".join(map(str, [h, len(u_freq), sum(u_freq),t1,0,"YES",consensus, consensus ]))+"\n"
          passed_seqs_total = passed_seqs_total+f
        else:fail_less_than_threshold = fail_less_than_threshold+1
      else: 
        consensus, pass_consensus=Get_consensus_sequence_cluster(u_seq, u_freq,tmp_file,threshold)
        if(consensus.count("_")==0 and pass_consensus==1):
          outp = outp+"\t".join(map(str, [h, len(u_freq), sum(u_freq),t1,0,"YES",consensus, consensus ]))+"\n"
          passed_seqs_total = passed_seqs_total+f
        else:
          fail_less_than_threshold = fail_less_than_threshold+1
    if(ind>200):
       Write_out(outp, primer_tag_file)
       outp, ind = '',0
  Write_out(outp, primer_tag_file)
  outp, ind = '',0
  print total_tags, passed_seqs_total, fail_less_than_threshold
  return()

def Get_consensus_sequence_cluster(u_seq, u_freq,tmp_file,threshold):
  out=''
  for i in range(0,len(u_seq)):
    out=out+">"+str(i)+"\n"+u_seq[i]+"\n"
  fh=open(tmp_file+"txt", "w")
  fh.write(out)
  fh.close()
  insert = ''
  if( len(u_seq) > 2000):insert = '--parttree'
  mafft = "/apps/well/mafft/7.149/bin/mafft "
  command1 = mafft+" --retree 2 "+insert+" "+tmp_file+"txt > "+tmp_file+"aligned"
  commands.getoutput(command1)
  fh=open(tmp_file+"aligned","r")
  max_seqs={}
  for header,sequence in fasta_iterator(fh):
    max_seqs[sequence.upper()]=int(header)
  fh.close()
  bases=["A","T","G","C","-"]
  base_dict = {}
  for b in range(0,len(bases)):
    base_dict[bases[b]]=b
  consensus = ''
  start,pass_consensus = 0,1
  for i in range(0,len(sequence)):
    f = [0]*len(bases)
    for s in max_seqs:
      if(s[i] not in base_dict):
        print i, s[i], max_seqs
      f[base_dict[s[i]]]=f[base_dict[s[i]]]+u_freq[max_seqs[s]]
    if(f[4]==0):start==1
    if(max(f)*1.0/sum(f) >= threshold):
      if(bases[f.index(max(f))]!="-"):
        consensus=consensus+bases[f.index(max(f))]
    else:
      pass_consensus = 0
      if(start==1):
        for j in range(0,5):
          if(f[j]!=0):
            consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
        consensus=consensus+"_"
      else:
        f = f[0:4]
        if(bases[f.index(max(f))]!="-"):
          consensus=consensus+bases[f.index(max(f))]
        else:
          for j in range(0,4):
            if(f[j]!=0):
              consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
            consensus=consensus+"_"
  return(consensus,pass_consensus)

def Get_consensus_sequence_large(out_cluster, Fail_file,l1,ids,sum_u_freq,threshold,tmp_file,cd_hit_directory):
  fh=open(Fail_file,"w")
  fh.write(out_cluster)
  fh.close()
  Cluster_i(Fail_file,Fail_file+"cls",(l1-3.0)/l1,cd_hit_directory)
  cluster = Tree()
  fh=open(Fail_file+"cls.bak.clstr","r")
  for l in fh:
    l=l.strip().split()
    clust, id = l[0],l[2].replace("...","").replace(">","")
    cluster[clust][id].value = 1
  fh.close()
  max_s, max_clust = 0,''
  for c in cluster:
    f = 0
    for id in cluster[c]:
      f = f+ids[id][1]
    if(max_s<f):max_clust,max_s = c,f
  consensus,pass_consensus = '',0
  if(sum_u_freq*threshold<max_s):
    seq_align,s_freq = [],[]
    for id in cluster[max_clust]:
      seq_align.append(ids[id][0])
      s_freq.append(ids[id][1])
    consensus = Get_consensus_sequence(seq_align,s_freq,tmp_file,threshold)
    if(consensus.count("_")==0 and len(consensus)>3):pass_consensus=1
    else:pass_consensus = 0
  return(consensus, pass_consensus)

def Get_consensus_sequence(u_seq, u_freq,tmp_file,threshold):
  out=''
  for i in range(0,len(u_seq)):
    out=out+">"+str(i)+"\n"+u_seq[i]+"\n"
  fh=open(tmp_file+"txt", "w")
  fh.write(out)
  fh.close()
  insert = ''
  if( len(u_seq) > 2000):insert = '--parttree'
  command1 = "/software/pubseq/bin/mafft-6.857/bin/mafft --retree 2 "+insert+" "+tmp_file+"txt > "+tmp_file+"aligned"
  commands.getoutput(command1)
  fh=open(tmp_file+"aligned","r")
  max_seqs={}
  pfh= Check_fasta_not_empty(fh)
  if(pfh==1):
    for header,sequence in fasta_iterator(fh):
      max_seqs[sequence.upper()]=int(header)
    fh.close()
    bases=["A","T","G","C","-"]
    base_dict = {}
    for b in range(0,len(bases)):
      base_dict[bases[b]]=b
    consensus = ''
    start = 0
    for i in range(0,len(sequence)):
      f = [0]*len(bases)
      for s in max_seqs:
        f[base_dict[s[i]]]=f[base_dict[s[i]]]+u_freq[max_seqs[s]]
      if(f[4]==0):start==1
      if(max(f)*1.0/sum(f) >= threshold):
        if(bases[f.index(max(f))]!="-"):
          consensus=consensus+bases[f.index(max(f))]
      else:
        if(start==1):
          for j in range(0,5):
            if(f[j]!=0):
              consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
          consensus=consensus+"_"
        else:
          f = f[0:4]
          if(bases[f.index(max(f))]!="-"):
            consensus=consensus+bases[f.index(max(f))]
          else:
            for j in range(0,4):
              if(f[j]!=0):
                consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
              consensus=consensus+"_"
    else:
      print "ERROR", u_freq, u_seq
  else:consensus =""
  return(consensus)

def Trim_sequences(Tmp_file,Fail_file,Output_trim,gene_specific_primer,dir,sample,sUMI_sample_barcode_id, sUMI_sample_barcode_file, cd_hit_directory, id, primer_tag_file, primer_tag_file_count,Tmp_file2):
  (rc)= Init_rc()
  forward, reverse,  barcoded_j, barcoded_v,v_ref = Get_primers_split(rc, gene_specific_primer)
  fh_out = open(Output_trim,"w")
  fh_out.close()
  fh_out = open(Fail_file,"w")
  fh_out.close()
  threshold_BARCODE = 0.80 ### Certainty for accepting a barcode
  if(barcoded_j==1 and barcoded_v==1):
    Double_barcoded_trimming(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,gene_specific_primer,Tmp_file2,primer_tag_file_count,threshold_BARCODE,Tmp_file2)
  return()

def Double_barcoded_trimming(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_file,tmp_file,primer_tag_file_count,threshold_BARCODE,Tmp_file2):
  Read_untrimmed_file_double(rc,Tmp_file,Fail_file,Output_trim,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,forward, reverse)
  Check_barcodes_MALBAC(primer_tag_file_count,primer_tag_file,Tmp_file2,Output_trim,threshold_BARCODE,cd_hit_directory,Tmp_file)
  Print_trimmed_sequences(Output_trim, primer_tag_file, primer_tag_file_count)
  return()

def Print_trimmed_sequences(Output_trim, primer_tag_file, primer_tag_file_count):
  fh=open(primer_tag_file,"r")
  unique_seqs,seqs = Tree(),{}
  passed,failed,uniq_failed = 0,0,{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[6]=="YES"):
        passed = passed+1
        consensus = l[7]
        seqs[l[3]+":"+l[4]] = [l[0],consensus]
      else:
        failed = failed+1
        if(l[3] in uniq_failed):uniq_failed[l[3]]= uniq_failed[l[3]]+1
        else:uniq_failed[l[3]] = 1
  fh.close()
  for bc in seqs:
    id = seqs[bc][0].split(":")[0].split("__")[0]
    unique_seqs[seqs[bc][1]][id].value = 1
  del seqs
  out,ind = '',0
  for seq in unique_seqs:
    f = len(unique_seqs[seq])
    for id in unique_seqs[seq]:
      break
    out=out+">"+id.split("__")[0]+"__"+str(f)+"\n"+seq+"\n"
    ind = ind+1
    if(ind>100):
      Write_out(out, Output_trim)
      out, ind = '',0
  Write_out(out, Output_trim)
  return()

def Read_untrimmed_file_double(rc,Tmp_file,Fail_file,Output_trim,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,forward, reverse):
  for f in [primer_tag_file_count]:
    fh=open(f, "w")
    fh.close()
  fh = open (Tmp_file, "r")
  seqs = Tree() 
  maxl=1000000
  minl = 110
  J_found,v_found,tot,indexing=0,0,0,0
  seqs1,t=Tree(),0
  for header,seq in fasta_iterator(fh):
    seq=seq.upper()
    seqs1[seq][header].value=1
    t=t+1
  out,ind = "#ID\tJ_tag\tV_tag\tSequence\ttype_rev\ttype_for\tprimer_rev\tprimer_for\n", 0
  total,pass_r, pass_f, pass_all = 0,0,0,0
  for seq in seqs1:
    for header in seqs1[seq]:
      break
    header = header+":"+str(len(seqs1[seq]))
    number = len(seqs1[seq])
    passes,j_tag,v_tag=0,'',''
    type_rev, type_for = '',''
    primer_rev,primer_for = '',''
    total = total+1
    for i in range(0,len(reverse)):
      pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])==-1 or len(pj)==0):
        seq = Reverse_comp(seq, rc) 
        pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])!=-1):
        bc_len = len(reverse[i][3])+3#len(reverse[i][2])/5
        if(pj[0]>bc_len):
          j_tag = seq[pj[0]+len(reverse[i][2]):pj[0]+len(reverse[i][2])+bc_len]
          seq = seq[0:pj[0]+bc_len]
          type_rev,primer_rev = reverse[i][1], reverse[i][5]
          passes = 1
          break
    #if(passes!=1):
    if(passes ==1000):
      for i in range(0,len(reverse)):
        words =reverse[i][7]
        p=[]
        for w in words:
          pj=Get_match(w[0], seq)
          if(max(pj+[-1])!=-1 and len(pj)>0):
            if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
            else:pj[0]=-1
        if(len(p)>1):
          bc_len = len(reverse[i][3])+3
          pj = max(set(p), key=p.count)
          start = pj
          j_tag = seq[pj-bc_len:pj+1]
          #start = pj-bc_len
          primer_match_seq = seq[start:start+len(reverse[i][6])]
          seq = seq[pj+bc_len:len(seq)]
          seq = Reverse_comp(seq, rc)
          type_rev,primer_rev = reverse[i][1], reverse[i][5]
          passes = 1
          break
      if(passes!=1):
        seq = Reverse_comp(seq, rc) 
        for i in range(0,len(reverse)):
          words =reverse[i][6]
          p=[]
          for w in words:
            pj=Get_match(w[0], seq)
            if(max(pj+[-1])!=-1 and len(pj)>0):
              if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
              else:pj[0]=-1
          if(len(p)>1):
            bc_len = len(reverse[i][3])+3
            pj = max(set(p), key=p.count)
            j_tag = seq[pj-bc_len:pj+1]
            start = pj-bc_len 
            primer_match_seq = seq[start:start+len(reverse[i][6])]
            seq = seq[pj+bc_len:len(seq)]
            seq = Reverse_comp(seq, rc)
            type_rev,primer_rev = reverse[i][1], reverse[i][5]
            passes = 1
            break
    if(passes==1):
      pass_r = pass_r+1
      for i in range(0,len(forward)):
        pv=Get_match(forward[i][0], seq)
        if(max(pv+[-1])!=-1):
          v_tag='-'
          seq = seq[pv[0]+len(forward[i][0])-1 :len(seq)]
          type_for, primer_for = forward[i][1],forward[i][2]
          passes = 2
          break
      if(passes!=2):
        for i in range(0,len(forward)):
          words =forward[i][6]
          p=[] 
          for w in words:
            pv=Get_match(w[0], seq) 
            if(max(pv+[-1])!=-1 and len(pv)>0):
              if(pv[0]<len(seq)/2):p=p+[max([0,pv[0]-w[1]])]
              else:pv[0]=-1
          if(len(p)>3):
            pv = max(set(p), key=p.count)
            bc_len = len(forward[i][3])+1
            v_tag = seq[pv-bc_len:pv]
            type_for, primer_for = forward[i][1],forward[i][2]
            seq = seq[pv+len(forward[i][2])-1 :len(seq)]
            type_for, primer_for = forward[i][1],forward[i][2]
            passes = 2
            break
    if (passes ==2):
      pass_f = pass_f+1
      if(len(seq)>minl and len(seq)<maxl):
        if(seq.count("N")==0 and len(j_tag)>2 and len(v_tag)>2):
          out=out+"\t".join(map(str,[header, j_tag, v_tag, seq, type_rev,type_for,primer_rev,primer_for]))+"\n"
          pass_all = pass_all+1
          ind = ind+1
          if(ind>200):
            Write_out(out, primer_tag_file_count)
            out, ind = '',0
  fh.close()
  Write_out(out, primer_tag_file_count)
  out, ind = '',0
  print "\tTotal sequences:\t"+str(total)+"\tNumber of REV primers found:\t"+str(pass_r)+"\tNumber of FOR primers found:\t"+str(pass_f)+"\tPassed sequences:\t"+str(pass_all)
  return()

def Get_primers_split(rc, primer_file):
  fh=open(primer_file, "r")
  forward, reverse,v_ref = [],[],[]
  barcoded_j, barcoded_v = 0,0
  word_size = 8
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>2):
        header, sequence, uni, BC,gene_specific =l[0],l[1],l[2],l[3],l[4]
        gene_specific=gene_specific.upper()
        words,words1 = [],[]
        for i in range(word_size,len(gene_specific)):
          words.append([gene_specific[i-word_size:i],i-word_size])
        for i in range(word_size,len(uni)):
          words1.append([uni[i-word_size:i],i-word_size])
        if(header.count("J")!=0 or header.count("REV_CONST")!=0 or header.count("REVERSE")!=0):
          if(header.count("REV_CONST")!=0):inside=1
          if(sequence.count("N")!=0):barcoded_j=1
          sequence=sequence.upper()
          clas = '-'
          l=len(gene_specific)
          reverse =reverse+[[sequence, clas, uni, BC,gene_specific,header,words,words1]]
        else:
          if(sequence.count("N")!=0):barcoded_v = 1
          sequence=sequence.upper()
          clas = '-'
          l=len(gene_specific)
          forward = forward+[[sequence, clas, uni, BC,gene_specific,header,words,words1]]
      else:
        header, sequence=l[0],l[1]
        sequence=sequence.upper()
        words = []
        for i in range(word_size,len(sequence)):
          words.append([sequence[i-word_size:i],i-word_size])
        if(header.count("J")!=0 or header.count("REV_CONST")!=0 or header.count("REVERSE")!=0):
          if(header.count("CONST")!=0):inside=1
          if(sequence.count("N")!=0):barcoded_j=1
          sequence=sequence.upper()
          clas = "-"
          reverse =reverse+[[sequence,clas, sequence,header, words]]
        elif(header.count("REF")==0):
          if(sequence.count("N")!=0):barcoded_v = 1
          sequence=sequence.upper()
          clas = "-"
          forward = forward+[[sequence, clas,header,words]]
        else:
          sequence=sequence.upper()
          clas = "-"
          v_ref = v_ref+[[sequence, clas,header,words]]
  fh.close()
  return(forward, reverse,  barcoded_j, barcoded_v,v_ref)

def Check_fasta_not_empty(fh):
  pfh = 0
  for l in fh:
    if(len(l)!=0):pfh =1
    break
  return(pfh)

def Get_freq (id):
  id = id.split("__")
  return(int(id[1]))

def Cluster_i (Reduced_file,tmp_file,diff,cd_hit_directory):
  command= cd_hit_directory+"cd-hit -i "+Reduced_file+" -o "+tmp_file+" -c "+str(diff)+" -d 180 -T 10  -M 0 -AL 40 "
  os.system(command)
  return()

def Check_same(ids, seqs_sub, same,inv):
  l=len(ids)
  ind = 0
  for i in range(0,l):
    for j in range(i,l):
      if(seqs_sub[i]==seqs_sub[j]):
        same [ids[i]][ids[j]].value = 1
        inv[ids[j]]=ids[i]
        ind = ind+1
      else:
        s1 = seqs_sub[i]
        s2 = seqs_sub[j]
        (sa,sb,p)=Trim_sequences1(s1,s2,len(s1),len(s2))
        if(sa==sb):
          same [ids[i]][ids[j]].value = 1
          inv[ids[j]]=ids[i]
          ind = ind+1
    if(ind>l/2):break
  return(same,inv)

def Trim_sequences1(s1,s2,l1,l2):
  p=0
  i=15
  sample1=s1[i:i+25]
  index=s2.find(sample1)
  if (index!=-1):
    if(index>i):
      if(index-i <=20):s2=s2[index-i:l2]
    else:
      if(i-index <=20):s1=s1[i-index:l1]
    min_len=min([len(s1),len(s2)])
    if ((max([len(s1),len(s2)]) - min_len) <25):
      s1=s1[0:min_len]
      s2=s2[0:min_len]
      p=1
  else:
    i=l1-50
    sample1=s1[i:i+25]
    index=s2.find(sample1)
    if (index!=-1):
      if(index>i):
        if(index-i <=20):s2=s2[index-i:l2]
      else:
        if(i-index <=20):s1=s1[i-index:l1]
      min_len=min([len(s1),len(s2)])
      if ((max([len(s1),len(s2)]) - min_len) <25):
        s1=s1[0:min_len]
        s2=s2[0:min_len]
        p=1
      else:
        p=0
    else:
      p=0
  return (s1, s2, p)

def Get_clusters (file):
  fh=open(file,"r")
  cluster=Tree()
  for l in fh:
    l=l.strip()
    l=l.split()
    id = l[2].replace("...","")
    id = id.replace(">","")
    cluster[l[0]][id].value=1
  fh.close()
  return(cluster)

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def Write_out(out, file):
  fh=open(file, "a")
  fh.write(out)
  fh.close()
  return ()

def Get_locations(ref_locations):
  fh=open(ref_locations,"r")
  locations = {}
  for l in fh:
    l=l.strip().split()
    locations[l[0]] = l[1]
  fh.close()
  return(locations)

###########################
if(len(sys.argv)<8):
  print "\nUsage:\npython Read_processing_and_quality_PUBLIC.py <input directory> <sample ID> <sample ID of pooled library> <sUMI primer file> <sUMI_sample barcode id> <sample barcode file> <cd-hit location>"
  quit()

dir = sys.argv[1]
id = sys.argv[2]
sample = sys.argv[3]
gene_specific_primer = sys.argv[4]
sUMI_sample_barcode_id = sys.argv[5]
sUMI_sample_barcode_file =sys.argv[6]
cd_hit_directory = sys.argv[7]
########################### Files for QC and filtering
Seq_file1 = dir+"/Sequences_"+id+"_1.fasta"
Seq_file2 = dir+"/Sequences_"+id+"_2.fasta"
Tmp_file = dir+"TMP/Untrimmed_"+id+".fasta"
Tmp_file2 = dir+"TMP/Tmp_"+id+".fasta"
Trim1=dir+"Filtered_sequences_"+id+".fasta"
Fail_file = dir+"TMP/Fail_filtered_"+id+".fasta"
primer_tag_file = dir+"TMP/Barcode_filtering_information_"+id+".txt"
primer_tag_file_count = dir+"TMP/All_barcodes_"+id+".txt"
######################### Commands
commands.getoutput("mkdir "+ dir+"/TMP/")
Split_reads_and_join(Seq_file1, Seq_file2, Tmp_file,id,gene_specific_primer,dir,sample,sUMI_sample_barcode_id, sUMI_sample_barcode_file) 
Trim_sequences(Tmp_file,Fail_file,Trim1,gene_specific_primer,dir,sample,sUMI_sample_barcode_id, sUMI_sample_barcode_file, cd_hit_directory, id,primer_tag_file, primer_tag_file_count,Tmp_file2)

