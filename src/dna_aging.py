# -*- coding: utf-8 -*-

from __future__ import division
from Bio import SeqIO
from Bio import Seq


import random
import math
import pysam
import sys

# for mutating the character on the given position
def change_char(seq, position, change): 
    return seq[:position]+change+seq[position+1:]

#clipping the sequence according to a given length distribution

def clipper(number_to_clip, sequence, random_sampled_binary, clip):

    # clipping the sequence upto the sampled length     
    # 0 and 1 plays a crucial role here because if the length of the distribution is 100 or 99, 
    # it clips everything, we don't want that
    if number_to_clip == 0:
        clip.append(0)
        clip.append(0)
        return sequence

    elif number_to_clip == 1:
        if random_sampled_binary[0] is 1:
            clip.append(0)
            clip.append(1)
            return sequence[:-1]
        else:
            clip.append(1)
            clip.append(0)
            return sequence[1:]
    else:
        #if the number can be divided by 0
        if(number_to_clip % 2 is 0):
            toClip = int(number_to_clip/2)
            clip.append(toClip)
            clip.append(toClip)
            return sequence[toClip:-toClip]
        #if the number cannot be divided by 0 
        else:
        
            if random_sampled_binary[0] is 1:
                toClip_beginning = int(math.ceil(number_to_clip/2))
                toClip_end = int(math.floor(number_to_clip/2))
                clip.append(toClip_beginning)
                clip.append(toClip_end)
                return sequence[toClip_beginning:-toClip_end]
            else:
                toClip_beginning = int(math.floor(number_to_clip/2))
                toClip_end = int(math.ceil(number_to_clip/2))
                clip.append(toClip_beginning)
                clip.append(toClip_end)
                return sequence[toClip_beginning:-toClip_end]
         #the list 'clip' is to define the clipped lengths. It is used later for phred quality changes. 

# applying the pmd changes at both ends of the sequence

def pmd_changes(temp_sequence, CtoT, GtoA, changes):

    # this is for C to T conversion 
    temp=temp_sequence
    #its length is 15. Thus, it does not affect the time complexity significantly(negligible)
    for number,item in enumerate(CtoT):  
        percentage = []
        # if the pmd data is reasonable, it will be at most ~35-40
        for i in range(0, item):
            percentage.append(1)
        for j in range(item, 100): # it fills the rest with 0s upto length 100
            percentage.append(0)
        # changing the character done here
        if (random.sample(percentage, 1)[0] is 1) and (temp[number] is 'C'):
            temp = change_char(temp,number,'T')
            changes.append(number)

    # the same rules above applies here, but applied to the changed sequences.
    # this is for G to A conversion

    for number,item in enumerate(GtoA):
        percentage = []
        for i in range(0,item):
            percentage.append(1)
        for j in range(item, 100):
            percentage.append(0)
        if (random.sample(percentage, 1)[0] is 1) and (temp[len(temp)-1-number] is 'G'):
            temp = change_char(temp,len(temp)-1-number,'A')
            changes.append(len(temp)-1-number)
    return temp


                  
if __name__ == '__main__':


    distribution_array=[]
    CtoT=[]
    GtoA=[]
    

    # if bam file should be read, this can be used
    '''bam_file = pysam.AlignmentFile('./ex.bam', 'rb')
    for read in bam_file:
        if (int(read.query_length) <=100):
            distribution_array.append(read.query_length)
    '''
    #'/mnt/NEOGENE1/projects/analysisProb_2019/wgcapture/kilinc16/bam_reads/Bon002_03051102_15.merged.hs37d5.fa.cons.90perc.bed'
    with open(sys.argv[1], 'r') as distribution:
        for line in distribution:
            chromosome, start, end, length, sequence  = line.split('\t')
            if (int(length) <=100):
                distribution_array.append(int(length))
    #'/mnt/NAS/projects/2019_epipaleomix/dna_aging/Alh139.pmd.txt'
    with open(sys.argv[2], 'r') as pmd:
        next(pmd)
        for line in pmd:
            position, ct, ca, cg, cc, ga, gt,gC, gg  = line.split(' ')
            CtoT.append(int(round(100*float(ct))))
            GtoA.append(int(round(100*float(ga))))
    
    with open('output.fq', 'w') as output:
        #"/mnt/NAS/projects/2019_arda/simu30x_lane0.read1.fq"
        for index,record in enumerate(SeqIO.parse(sys.argv[3], "fastq")):
            clip = []
            changes = []
            new_letter_annotations_seq=""

            Sequence = str(record.seq)
            letter_annotations = record.letter_annotations
            record.letter_annotations = {}

            number_to_clip = 100 - random.sample(distribution_array, 1)[0]
            
            # binary sampling done to decide from which end to clip the sequence and clipping
            temp_sequence = clipper(number_to_clip, Sequence, random.sample([0,1], 1), clip)
            record.seq = Seq.Seq(pmd_changes(temp_sequence, CtoT, GtoA, changes))
            
            # change phred quality score to 15(estimated for aDNA) according to changes done for pmd(if else for clipping the annotations)
            if clip[1] is not 0:
                new_letter_annotations_seq = letter_annotations['phred_quality'][clip[0]:-clip[1]]
            else:
                new_letter_annotations_seq = letter_annotations['phred_quality'][clip[0]:]

            for number in changes:
                new_letter_annotations_seq[number] = 15

            record.letter_annotations = {'phred_quality': new_letter_annotations_seq}
            SeqIO.write(record, output, "fastq")
