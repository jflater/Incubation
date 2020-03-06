#!/usr/bin/python
"""
This script find representative sequence from list file
usage: python get_seq.py fastafile listfile
python get_seq.py stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list 
"""

import sys
import screed
import datetime
import multiprocessing
from multiprocessing import Pool

def get_diagonal(diagonal_dist, seq_name, seq,i):
    diagonal_dist[i] = []
    for j in range(i+1,len(seq_name)):
            dist = get_distance(seq[i], seq[j])
            diagonal_dist[i].append(dist)

def get_distance(seq1,seq2):
    seq1.upper()
    seq2.upper()
    dist = 0
    if not len(seq1) == len(seq2):
	print "sequences are not aligned"
        sys.exit()
    match = 0
    tot = 0
    nu = ["A","G","C","T"]
    for i in range(0,len(seq1)):
        if seq1[i] in nu or seq2[i] in nu:
            tot += 1
            if seq1[i] == seq2[i]:
                match += 1
    proportion =  float(match) / float(tot)
    dist = 1 - proportion
    return dist

def get_center(otu_name, li, fasta):
    sub_fasta = {}
    seq_name = []
    seq = []
    name_num = {}
    for n, one_list in enumerate(li):
        seq_name.append(one_list)
        seq.append(fasta[one_list])
        name_num[one_list] = n
    diagonal_dist = {}

    for i in range(0,len(seq_name)):
        diagonal_dist[i] = []
        for j in range(i+1,len(seq_name)):
            dist = get_distance(seq[i], seq[j])
            diagonal_dist[i].append(dist)

    #make full table
    full_dist ={}
    for n in range(0,len(li)):
        full_dist[n] = diagonal_dist[n]
        for i in range(n,-1,-1):
            if n == i:
                full_dist[n].insert(0,0)
                continue
            full_dist[n].insert(0,diagonal_dist[i][n])

    #calculate smallest
    smallest = 1000000000
    smallest_id = 0
    for item in full_dist.items():
        if sum(item[1]) < smallest:
            smallest = sum(item[1])
            smallest_id = item[0]
    print ">"+otu_name+" "+li[smallest_id]+'\n'+fasta[li[smallest_id]]

def main():
    #read fasta
    #print datetime.datetime.now()
    fasta = {}
    for record in screed.open(sys.argv[1]):
        fasta[record.name] = record.sequence

    #read list
    otus = []
    li = []
    for n,line in enumerate(open(sys.argv[2],'r')):
          spl = line.strip().split('\t')
          if n == 0:
              otus = spl[2:len(spl)]
          elif n == 1:
              for i in range(2,len(spl)):
                  li = spl[i].split(',')
                  get_center(otus[i-2], li, fasta)
    #print datetime.datetime.now()

if __name__ == '__main__':
    main()
