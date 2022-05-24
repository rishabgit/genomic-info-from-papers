#!/usr/bin/env python

# from __future__ import division
# from __future__ import print_function

import os.path
import argparse
import re

import gffutils
from pyfaidx import Fasta

from wbpreader.translation import translate_dict
from wbpreader.translation import vars_from_GFF
from wbpreader.translation import sub_is
from wbpreader.translation import var_abspos
from wbpreader.translation import calculate_flanks


"""

Script for generating variants Ace from mutation

"""


# TODO
# Make it work for variants on reverese strand too - DONE
# Implement sub-routine clone_from_gene in translation.py to get data for the field Mapping_target
# Mapping_target	 "Y48B6A" <-clone name
# Mapping_target	 "F13G3" <-clone name
# Calculate Source_location
# Source_location	 200 "CHROMOSOME_I" 4393280 4393280
# Source_location	 200 "CHROMOSOME_I" 6461401 6461401
# Compare Source_location to current ones - if there is overlap they may need to merge
# 


epi = ('\
    \n\
	Takes input variant and transcripts, and generates ace-format data\n\
    The inupt file should contain 3 columns: \n\
    1. Gene name eg WBGene00007063\n\
    2. Transcript name eg 2L52.1b\n\
    3. Variant code eg V600E K342* or W256amber\n\
    Example: \n\
    WBGene00007065	3R5.1b	    Q214L \n\
    WBGene00007064	2RSSE.1a	L10K \n\
    WBGene00007064	2RSSE.1b	A4ochre \n\
    WBGene00007064	2RSSE.1b	S7amber \n\
    \n\
    Note: The genome fasta, protein fasta and GFF needs to be from the same source, and pertain  \n\
    to the genome and annotation the variant originates from \n\
    Note: If the GFF contains variants with VEP results, a smaller ACE file will be created for  \n\
    those new ones overlapping with existing ones\n\
    Note: Wormbase GFFs have to be pre-processed using the script \n\
    \n\
    \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script generates ace files from mutations', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='inp', action='store', required=True, help="Variant eg K342* or W256amber")
parser.add_argument('-p', '--protein', default=None, dest='pro', action='store', required=True, help="Protein fasta file")
parser.add_argument('-gff', '--gff', default=None, dest='gff', action='store', required=True, help="GFF file with full path")
parser.add_argument('-g', '--genome', default=None, dest='genome', action='store', required=True, help="Genome fasta file")

args = parser.parse_args()


def is_valid_file(arg):
    if not os.path.isfile(arg):
        print("The file does not exist!", arg)
        exit(1)
    else:
        pass
        #print ("File exists", arg)
        #return open(arg, 'r')  # return an open file handle


# print ("Do", args, args.inp)

is_valid_file( args.inp)
is_valid_file( args.pro)
# is_valid_file( args.gff)
# is_valid_file( args.genome)


# read the input file
inp = open(args.inp, 'r')

# Get a translation table
prots=translate_dict('all')
#print (prots)

# Get all known variants from GFF
known_vars = vars_from_GFF(args.gff)
#print (known_vars)

#exit(0)


# create an object of new bed file and open in to write data.
output=args.inp+".ace"
out = open(output, 'w')

vars={}

# Save the variants and transcripts
for r in inp:
    a = r.split()
    #print (a[0],a[1],sep='\t')
    v=a[0]+':'+a[1]
    if v in vars:
        vars[v][a[2]]={}
        vars[v][a[2]]['PUB']=a[3]
    else:
        vars[v]={}
        vars[v][a[2]]={}
        vars[v][a[2]]['PUB']=a[3]


# Read in the gff coords + chromosome
fn = gffutils.example_filename(args.gff)
#print(open(fn).read()) 
db = gffutils.create_db(fn, ':memory:')



# Make a hash containing the clones all transcripts belong to
clone={}

for t in vars:
    a=t.split(':')
    a[0]='Gene:'+a[0]
    clone[a[1]]='n'

    clone=clone_from_ts(gff,db,clone)


    #for v in vars[t]:
    #    clone=clone_from_ts
    #    vars[t][v]['CLONE']=clone



quit()


# Calculate the absolute genome position and chr

for t in vars:
    a=t.split(':')
    a[0]='Gene:'+a[0]
    #print (a[0],a[1],sep='\t',end='\t')


    for v in vars[t]:

        vari=re.split('(\d+)',v)
        varelstart=int(vari[1])*int(3)
        varelend=int(varelstart)+int(3)
        #print (vari,varelstart,varelend)
        try:
            gene = db[a[0]]
        except:
            print("WARNING: Gene", a[0] , "not found in GFF file, skipped")
            break

        #print (v, gene['Name'],gene.start, gene.end, varelstart, sep='\t')

        # ABSPOS calculation for genes on the positive strand
        (chro,abspos,geneori)=var_abspos(gene,a[1],varelstart,db)
        vars[t][v]['ABSPOS']=abspos
        vars[t][v]['CHR']=chro
        vars[t][v]['ORI']=geneori


# Read in the genomic sequences
fas = Fasta(args.genome)



# Check if the extracted regions match and look up the flanking sequences
# Extrapolate the alt

for t in vars:
    a=t.split(':')
    #a[0]='Gene:'+a[0]
    #print (a[0],a[1],sep='\t',end='\t')

    for v in vars[t]:
        vari=re.split('(\d+)',v)
        (varseq,left,right)=calculate_flanks(vars[t][v],fas,prots,vari[0])
        (ref,alt)=sub_is(varseq, vari[2])
        vars[t][v]['REF']=varseq
        vars[t][v]['LEFT']=left
        vars[t][v]['RIGHT']=right
        vars[t][v]['ALT']=alt



# Produce the ace file
print ("produce_ace",'\n')

i=99000000

for t in vars:
    i=i+1
    a=t.split(':')

    for v in vars[t]:
        vari=re.split('(\d+)',v)


        # See if the variant has been found before
        varid="WBVar" + str(i)
        
        ts=a[1]+'.1'
        ts2=a[1]+'.2'

        # Make stop codons general
        vmod=v.replace('amber','*')
        vmod=vmod.replace('opal','*')
        vmod=vmod.replace('ochre','*')



        # Variant is known since before, so we only add minimal data
        if ts in known_vars:
            if v in known_vars[ts]:
                #print("Known vars", a[0], a[1],v , known_vars[a[1]][v])
                (varid,chr,start,end)=known_vars[ts][v].split(':')
                print ("KNOWN", varid,chr,start,end)
            elif vmod in known_vars[ts]:
                #print("Known vars", a[0], a[1],v , known_vars[a[1]][v])
                (varid,chr,start,end)=known_vars[ts][vmod].split(':')
                print ("KNOWN", varid,chr,start,end)
            else:
                print("KNOWN TS", a[0], ts,v, known_vars[ts])
                pass
        elif ts2 in known_vars:
            if v in known_vars[ts2]:
                #print("Known vars", a[0], a[1],v , known_vars[a[1]][v])
                (varid,chr,start,end)=known_vars[ts2][v].split(':')
                print ("KNOWN", varid,chr,start,end)
            elif vmod in known_vars[ts2]:
                #print("Known vars", a[0], a[1],v , known_vars[a[1]][v])
                (varid,chr,start,end)=known_vars[ts2][vmod].split(':')
                print ("KNOWN", varid,chr,start,end)
            else:
                print("KNOWN TS", a[0], ts2,v, known_vars[ts2])
                pass

            
        # This is a new variant, so we add as much information as we can
        else: 
            print("NOT KNOWN vars", ts,v)
            pass


        print ("Variation : \"", varid , "\"", sep='')
        print ("Evidence\t ", "Paper_evidence \"", vars[t][v]['PUB'], "\"", sep='')
        print("Transcript\t \"",a[1],".1\" Missense ", vari[1], " \"", vari[0] , " to ", vari[2],"\"",sep='')
        #Transcript       "Y51H7C.11.1" Missense 151 "G to E"
        print ("Substitution ", "\"",vars[t][v]['REF'] ,  "\" \"", vars[t][v]['ALT'] ,  "\" Paper_evidence \"", vars[t][v]['PUB']  ,"\"",sep='') # Sequence_details Type_of_mutation

        if vars[t][v]['ORI']=='+': 
            try:
                print("Source_location 200 \"CHROMOSOME_", vars[t][v]['CHR'] , "\"  " , vars[t][v]['ABSPOS'] , " ", vars[t][v]['ABSPOS']+2 , sep='')
            except:
                pass
        elif vars[t][v]['ORI']=='-': 
            try:
                print("Source_location 200 \"CHROMOSOME_", vars[t][v]['CHR'] , "\"  " , vars[t][v]['ABSPOS'] , " ", vars[t][v]['ABSPOS']+2 , sep='')
            except:
                pass
        else:
            print("WARNING: this should not happen")

        try:
            print("HGVSg\t \"CHROMOSOME_",vars[t][v]['CHR'] ,':g.',vars[t][v]['ABSPOS'],ref,'>',alt, "\"", sep='')
            #HGVSg	 "CHROMOSOME_V:g.20824795C>T"
            #HGVSg	 "CHROMOSOME_V:g.20821901C>T"
        except:
            pass

        try:
            print ("Sequence_details Flanking_sequences \"", vars[t][v]['LEFT'] , "\" \"" , vars[t][v]['RIGHT'] ,"\"",sep='')
        except:
            pass
        print("Sequence_details SeqStatus Sequenced")
        print("Variation_type Allele")
        print("Variation_type Predicted_SNP")
        print("Sequence_details  Mapping_target \"###\"")
        print("Species \"Caenorhabditis elegans\"")
        print("Origin\t Status Live")
        print("Method\t \"Substitution_allele\"")
        print("Reference\t \"", vars[t][v]['PUB'], "\"" ,sep='')
        print("Affects Gene \"", a[0], "\"", sep='')
        print("Delete Not_sequenced or Pending_curation and replace with Sequenced")


       
        print()

'''
Variation : "WBVar00000002" -O "2013-08-13_09:50:15_pad"
Evidence	 -O "2010-02-05_10:07:16_mt3" Paper_evidence -O "2010-02-05_10:07:16_mt3" "WBPaper00035215" -O "2010-02-05_10:07:16_mt3"
Name	 -O "2010-02-05_10:07:16_mt3" Public_name -O "2010-02-05_10:07:16_mt3" "ac1" -O "2010-02-05_10:07:16_mt3"
Sequence_details	 -O "2010-02-05_10:07:16_mt3" Flanking_sequences -O "2010-02-05_10:07:16_mt3" "aattcgaaaaagtcgagttcacagccggtg" -O "2010-02-05_10:07:16_mt3" "agcccacagagatgacccaattttcgcgg
a" -O "2010-02-05_10:07:16_mt3"
Sequence_details	 -O "2010-02-05_10:07:16_mt3" Mapping_target -O "2012-12-14_16:52:32_mt3" "Y51H7C" -O "2012-12-14_16:52:32_mt3"
Sequence_details	 -O "2010-02-05_10:07:16_mt3" Type_of_mutation -O "2010-02-05_10:07:16_mt3" Substitution -O "2010-02-05_10:07:16_mt3" "g" -O "2010-02-05_10:07:16_mt3" "a" -O "2010-02-05_10:07:16_mt3" Paper_evidence -O "2010-02-05_10:07:16_mt3" "WBPaper00035215" -O "2010-02-05_10:07:16_mt3"
Sequence_details	 -O "2010-02-05_10:07:16_mt3" SeqStatus -O "2010-02-05_10:07:16_mt3" Sequenced -O "2010-02-05_10:07:16_mt3"
Variation_type	 -O "2010-02-05_10:07:16_mt3" Allele -O "2010-02-05_10:07:16_mt3"
Origin	 -O "2009-11-19_10:37:46_CGC_strain_update" Species -O "2010-02-05_10:07:16_mt3" "Caenorhabditis elegans" -O "2010-02-05_10:07:16_mt3"
Origin	 -O "2009-11-19_10:37:46_CGC_strain_update" Strain -O "2018-03-09_15:35:45_CGC_strain_update" "WBStrain00000321" -O "2018-03-09_15:35:45_CGC_strain_update"
Origin	 -O "2009-11-19_10:37:46_CGC_strain_update" Laboratory -O "2010-02-05_10:07:16_mt3" "AY" -O "2010-02-05_10:07:16_mt3"
Origin	 -O "2009-11-19_10:37:46_CGC_strain_update" Status -O "2010-09-14_09:18:54_mt3" Live -O "2010-09-14_09:18:54_mt3"
Affects	 -O "2009-11-19_10:45:16_CGC_strain_update" Gene -O "2009-11-19_10:45:16_CGC_strain_update" "WBGene00021789" -O "2009-11-19_10:45:16_CGC_strain_update"
Reference	 -O "2010-05-04_10:16:36_mt3" "WBPaper00035215" -O "2010-05-04_10:16:36_mt3"
Reference	 -O "2010-05-04_10:16:36_mt3" "WBPaper00047105" -O "2015-10-21_16:22:47_caltech_Paper"
Method	 -O "2010-02-05_10:07:16_mt3" "Substitution_allele" -O "2010-02-05_10:07:16_mt3"


Variation : "WBVar00000007" -O "2013-08-13_09:45:32_pad"
Evidence	 -O "2004-03-19_17:26:22_ck1" Author_evidence -O "2004-03-19_17:26:22_ck1" "Curtis L" -O "2004-03-19_17:26:22_ck1"
Name	 -O "2005-02-11_17:13:21_mt3" Public_name -O "2006-08-08_10:18:36_Public_name_patch" "ad446" -O "2006-08-08_10:18:36_Public_name_patch"
Sequence_details	 -O "2004-03-19_17:26:22_ck1" Flanking_sequences -O "2004-03-19_17:26:22_ck1" "agactgggcaaaaatctttgacgatctcgaaaatgtggtggtgaa" -O "2004-03-24_12:40:04_ck1" "aaaatacaaaataa
tgaagtaggcacgtgtatgtaggcag" -O "2004-03-24_12:40:04_ck1"
Sequence_details	 -O "2004-03-19_17:26:22_ck1" Mapping_target -O "2012-12-14_16:52:32_mt3" "C05D2"C05D2" -O "2012-12-14_16:52:32_mt3"

Variation : "WBVar00000009" -O "2013-08-13_09:50:15_pad"
Evidence	 -O "2004-09-17_15:15:15_rem" Paper_evidence -O "2004-09-17_15:15:15_rem" "WBPaper00006439" -O "2004-09-17_15:15:15_rem"
Name	 -O "2005-02-11_17:13:21_mt3" Public_name -O "2006-08-08_10:18:36_Public_name_patch" "ad451" -O "2006-08-08_10:18:36_Public_name_patch"
Sequence_details	 -O "2004-09-17_15:15:15_rem" Flanking_sequences -O "2004-09-17_15:15:15_rem" "aatgactacaaactacgatggtctcctgag" -O "2004-09-17_15:15:15_rem" "agtacggtaatattactacgttgcaaatt
c" -O "2004-09-17_15:15:15_rem"
Sequence_details	 -O "2004-09-17_15:15:15_rem" Mapping_target -O "2012-12-14_16:52:32_mt3" "Y48B6A" -O "2012-12-14_16:52:32_mt3"
Sequence_details	 -O "2004-09-17_15:15:15_rem" Type_of_mutation -O "2004-09-17_15:15:15_rem" Substitution -O "2004-09-17_15:15:15_rem" "g" -O "2008-05-02_10:12:16_mt3" "a" -O "2004-09-17_
15:15:15_rem"
Sequence_details	 -O "2004-09-17_15:15:15_rem" SeqStatus -O "2006-02-23_16:08:40_mt3" Sequenced -O "2006-02-23_16:08:40_mt3"
Variation_type	 -O "2005-02-09_11:38:48_mt3" Allele -O "2005-02-09_11:38:48_mt3"
Origin	 -O "2005-01-31_16:47:30_ar2" Species -O "2003-07-21_12:08:26_ck1" "Caenorhabditis elegans" -O "2003-07-21_12:08:26_ck1"
Origin	 -O "2005-01-31_16:47:30_ar2" Laboratory -O "2004-05-17_13:53:37_ck1" "DA" -O "2004-05-17_13:53:37_ck1"
Origin	 -O "2005-01-31_16:47:30_ar2" Status -O "2005-05-25_17:33:42_mt3" Live -O "2005-05-25_17:33:42_mt3"
Affects	 -O "2004-04-07_11:23:34_wormpub" Gene -O "2004-04-07_11:23:34_wormpub" "WBGene00001133" -O "2004-04-07_11:23:34_wormpub"
Reference	 -O "2010-04-01_13:56:13_mt3" "WBPaper00017350" -O "2010-04-01_13:56:13_mt3"
Reference	 -O "2010-04-01_13:56:13_mt3" "WBPaper00015559" -O "2010-04-01_13:56:13_mt3"
Method	 -O "2003-01-06_10:01:54_ck1" "Substitution_allele" -O "2004-09-17_15:15:15_rem"

Variation : "WBVar02153522" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Evidence	 -O "2021-07-10_01:45:23_curation_Variation" Person_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPerson1928" -O "2021-07-10_01:45:23_curation_Variation"
Name	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Public_name -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "ix259" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" Flanking_sequences -O "2021-07-10_01:45:23_curation_Variation" "cattcttgttcggctttttcgagacgttct" -O "2021-07-10_01:45:23_curation_Variation" "tgaaaaataatttttttttggaaattttct" -O "2021-07-10_01:45:23_curation_Variation"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" Mapping_target -O "2021-07-10_01:45:23_curation_Variation" "F15D3" -O "2021-07-10_01:45:23_curation_Variation"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" Type_of_mutation -O "2021-07-10_01:45:23_curation_Variation" Substitution -O "2021-07-10_01:45:23_curation_Variation" "c" -O "2021-07-10_01:45:23_curation_Variation" "t" -O "2021-07-10_01:45:23_curation_Variation"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" SeqStatus -O "2021-07-10_01:45:23_curation_Variation" Sequenced -O "2021-07-10_01:45:23_curation_Variation"
Affects	 -O "2021-07-10_01:45:07_curation_Gene" Gene -O "2021-07-10_01:45:07_curation_Gene" "WBGene00001131" -O "2021-07-10_01:45:07_curation_Gene"
Possibly_affects	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBGene00001131" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Paper_evidence -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBPaper00060694" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Possibly_affects	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBGene00001131" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Remark -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "CGC_name dys-1" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Isolation	 -O "2021-07-10_01:45:23_curation_Variation" Mutagen -O "2021-07-10_01:45:23_curation_Variation" "EMS" -O "2021-07-10_01:45:23_curation_Variation"
Isolation	 -O "2021-07-10_01:45:23_curation_Variation" Forward_genetics -O "2021-07-10_01:45:23_curation_Variation" "standard phenotypic screen" -O "2021-07-10_01:45:23_curation_Variation"
Reference	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBPaper00060694" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "[2021-07-07T15:47:18.438Z WBPerson51134] New Variation: Reference WBPaper00060694\; Gene WBGene00001131 (dys-1)" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Curator_confirmed -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBPerson51134" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "alt_det = g to a mut_det = ag to aa" -O "2021-07-10_01:45:23_curation_Variation" Paper_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPaper00060694" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "alt_det = g to a mut_det = ag to aa" -O "2021-07-10_01:45:23_curation_Variation" Person_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPerson1928" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "alt_det = g to a mut_det = ag to aa" -O "2021-07-10_01:45:23_curation_Variation" Curator_confirmed -O "2021-07-10_01:45:23_curation_Variation" "WBPerson51134" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "Variation information submitted by WBPerson1928 on 2021-04-23_02:24:25 via the Allele submission form." -O "2021-07-10_01:45:23_curation_Variation" Curator_confirmed -O "2021-07-10_01:45:23_curation_Variation" "WBPerson51134" -O "2021-07-10_01:45:23_curation_Variation"
Method	 -O "2021-07-10_01:45:23_curation_Variation" "Substitution_allele" -O "2021-07-10_01:45:23_curation_Variation"


Variation : "WBVar02153521" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Evidence	 -O "2021-07-10_01:45:23_curation_Variation" Person_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPerson1928" -O "2021-07-10_01:45:23_curation_Variation"
Name	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Public_name -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "ix261" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" Flanking_sequences -O "2021-07-10_01:45:23_curation_Variation" "CACAAACCATAAGTAGTCAAATTGAATACT" -O "2021-07-10_01:45:23_curation_Variation" "CTAGCCGATCTCCTGCCAGTGAATCGGCAC" -O "2021-07-10_01:45:23_curation_Variation"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" Mapping_target -O "2021-07-10_01:45:23_curation_Variation" "R06C1" -O "2021-07-10_01:45:23_curation_Variation"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" Type_of_mutation -O "2021-07-10_01:45:23_curation_Variation" Substitution -O "2021-07-10_01:45:23_curation_Variation" "C" -O "2021-07-10_01:45:23_curation_Variation" "T" -O "2021-07-10_01:45:23_curation_Variation"
Sequence_details	 -O "2021-07-10_01:45:23_curation_Variation" SeqStatus -O "2021-07-10_01:45:23_curation_Variation" Sequenced -O "2021-07-10_01:45:23_curation_Variation"
Origin	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Species -O "2021-07-10_01:45:23_curation_Variation" "Caenorhabditis elegans" -O "2021-07-10_01:45:23_curation_Variation"
Origin	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Strain -O "2021-07-10_01:45:23_curation_Variation" "WBStrain00048692" -O "2021-07-10_01:45:23_curation_Variation"
Origin	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Production_method -O "2021-07-10_01:45:23_curation_Variation" CRISPR_Cas9 -O "2021-07-10_01:45:23_curation_Variation"
Origin	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Status -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Live -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Affects	 -O "2021-07-10_01:45:07_curation_Gene" Gene -O "2021-07-10_01:45:07_curation_Gene" "WBGene00001836" -O "2021-07-10_01:45:07_curation_Gene"
Possibly_affects	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBGene00001836" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Paper_evidence -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBPaper00060694" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Possibly_affects	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBGene00001836" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Remark -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "CGC_name hda-3" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Reference	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBPaper00060694" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "[2021-07-07T13:24:49.168Z WBPerson51134] New Variation: Reference WBPaper00060694\; allele of WBGene00001836 (hda-3)" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" Curator_confirmed -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "WBPerson51134" -O "2021-07-10_01:44:05_NewNS_data_ace_10072021"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "OF1355 strain also has a synonymous mutation (g to a) at R(269)." -O "2021-07-10_01:45:23_curation_Variation" Paper_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPaper00060694" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "OF1355 strain also has a synonymous mutation (g to a) at R(269)." -O "2021-07-10_01:45:23_curation_Variation" Person_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPerson1928" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "OF1355 strain also has a synonymous mutation (g to a) at R(269)." -O "2021-07-10_01:45:23_curation_Variation" Curator_confirmed -O "2021-07-10_01:45:23_curation_Variation" "WBPerson51134" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "alt_det = g to a mut_det = G(271)E" -O "2021-07-10_01:45:23_curation_Variation" Paper_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPaper00060694" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "alt_det = g to a mut_det = G(271)E" -O "2021-07-10_01:45:23_curation_Variation" Person_evidence -O "2021-07-10_01:45:23_curation_Variation" "WBPerson1928" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "alt_det = g to a mut_det = G(271)E" -O "2021-07-10_01:45:23_curation_Variation" Curator_confirmed -O "2021-07-10_01:45:23_curation_Variation" "WBPerson51134" -O "2021-07-10_01:45:23_curation_Variation"
Remark	 -O "2021-07-10_01:44:05_NewNS_data_ace_10072021" "Variation information submitted by WBPerson1928 on 2021-04-22_19:21:33 via the Allele submission form." -O "2021-07-10_01:45:23_curation_Variation" Curator_confirmed -O "2021-07-10_01:45:23_curation_Variation" "WBPerson51134" -O "2021-07-10_01:45:23_curation_Variation"
Method	 -O "2021-07-10_01:45:23_curation_Variation" "Engineered_allele" -O "2021-07-10_01:45:23_curation_Variation"

Variation : "WBVar02153513" -O "2021-06-23_17:11:25_pad"
Name	 -O "2021-07-12_00:35:14_Variation_public_name_ace" Public_name -O "2021-07-12_00:35:14_Variation_public_name_ace" "xf61" -O "2021-07-12_00:35:14_Variation_public_name_ace"
Sequence_details	 -O "2021-06-23_17:11:25_pad" SeqStatus -O "2021-06-23_17:11:25_pad" Pending_curation -O "2021-06-23_17:11:25_pad"
Variation_type	 -O "2021-06-23_17:11:25_pad" Engineered_allele -O "2021-06-23_17:11:25_pad"
Origin	 -O "2021-06-23_17:11:25_pad" Species -O "2021-06-23_17:11:25_pad" "Caenorhabditis elegans" -O "2021-06-23_17:11:25_pad"
Origin	 -O "2021-06-23_17:11:25_pad" Strain -O "2021-06-23_17:11:25_pad" "WBStrain00048727" -O "2021-06-23_17:11:25_pad"
Origin	 -O "2021-06-23_17:11:25_pad" Strain -O "2021-06-23_17:11:25_pad" "WBStrain00048721" -O "2021-06-23_17:11:25_pad"
Origin	 -O "2021-06-23_17:11:25_pad" Laboratory -O "2021-06-23_17:11:25_pad" "RFK" -O "2021-06-23_17:11:25_pad"
Origin	 -O "2021-06-23_17:11:25_pad" Production_method -O "2021-06-23_17:11:25_pad" CRISPR_Cas9 -O "2021-06-23_17:11:25_pad"
Origin	 -O "2021-06-23_17:11:25_pad" Status -O "2021-06-23_17:11:25_pad" Live -O "2021-06-23_17:11:25_pad"
Method	 -O "2021-06-23_17:11:25_pad" "Engineered_allele" -O "2021-06-23_17:11:25_pad"



'''


out.close()






quit()

'''

    #### FILTER OUT #####
    # Shared called total
    # Filter out sites which
    chr = r.chrom
    pos = r.pos
    id = str(r.id)
    varID=':'.join([id.split(":")[0],id.split(":")[1]])
    #altb = r.ref
    #altb = r.alts
    score = r.qual
    filter = r.filter
    info = r.info
    format = r.format
    samples = r.samples
    end = r.stop # r.info["END"]
    strand='.'
    svtype='NA'

    if 'SVTYPE' in r.info.keys():
        svtype = r.info.get('SVTYPE', "")

    #for key in r.info.keys():
    #    data = r.info.get(key, "")
    #    print (key,data)
    #    svtype=r.info['SVTYPE']

    # FORMAT
    #['PR', 'SR', 'RC', 'BC', 'CN', 'MCC']


    # Split out Manta calls

    if re.match(r'Manta', varID):
        pass
        #print("Manta",chr, pos, end, varID, score, strand, sep='\t')




    #print (list((r.header.filters)))
    #print(list((r.header.formats)))

    elif re.match(r'Canvas', varID):

        # Extract relevant information
        cn='NA'
        mcc='NA'
        filter='NA'
        rc='NA'

        if 'CN' in r.samples[0].keys():
            cn=r.samples[0]['CN']
        if 'MCC' in r.samples[0].keys():
            mcc=r.samples[0]['MCC']
        if 'RC' in r.samples[0].keys():
            rc=r.samples[0]['RC']
        for key in r.filter.keys():
            filter=key


        print (chr,pos,end,varID,score,strand,filter,svtype,rc,cn,mcc,sep='\t', file=out)

    else:

        print("Unknown",chr, pos, end, varID, score, strand, sep='\t')

out.close()
'''

"""

# open up file
with open(args.vcf, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter ='\t')
    for row in reader:
        tsvout =""
        # pass header
        if row[0].startswith("#"):
            pass
        # CANVAS
        # only count canvas if Gain or Loss
        elif "Canvas" in row[2] and "REF" not in row[2]:
            # add to tsv file
            varID= row[2].split(":")
            end = row[7].split('END=')[-1]
            # if fist line do not add line break
            if len(tsvout) < 1:
            # sample chr start end QUAL Filter tool change
                tsvout += sample + "\t"+ row[0] + "\t"+ row[1]+ "\t"+ end +"\t"+ row[5]+"\t"+ row[6] +\
                      "\t"+ varID[0]+ "\t"+ varID[1]
            else:
                tsvout += "\n"+ sample + "\t"+ row[0] + "\t"+ row[1]+ "\t"+ end +"\t"+ row[5]+"\t"+ row[6] + \
                          "\t" + varID[0]+ "\t"+ varID[1]
            # add to tallies
            if "PASS" in row:
                canvasT += 1
                if "LOSS" in row[2]:
                    canvasLoss += 1
                elif "GAIN" in row[2]:
                    canvasGain += 1
            print(tsvout)
        # MANTA
        # tally up Manta outputs
        elif "Manta" in row[2]:
            # find SV type
            varID = row[2].split(":")
            QUAL = "."
            vcffilter = "."
            end = "."
            # identify translocations
            if "MantaBND" in row[2]:
            # No Ends in BND cases as they are translocations - add N/A to end field for the TSV
                end = "."
                QUAL = "."
                vcffilter = row[6]
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaBND += 1

            # identify Dels - can be problematic
            elif "MantaDEL" in row[2]:
                # some Dels columns out of sync with no ref or QUAL field so need to find end value from appropated field
                if len(row) ==11:
                    end = row[7].split(';')[0]
                    QUAL = row[5]
                    vcffilter = row[6]

                elif len(row) < 11:
                    #check which column end and Filter criteria are in
                    if "END=" in row[5]:
                        end = row[5].split(';')[0]
                        QUAL = "."
                        vcffilter = row[4]
                    elif "END=" in row[6]:
                        end = row[6].split(';')[0]
                        QUAL = "."
                        vcffilter = row[4]
                    else:
                        pass
                else:
                    pass
                end = "".join(i for i in end if i.isdigit())
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaDEL += 1

            # Finds invs
            elif "MantaINV" in row[2]:
                end = row[7].split(';')[0]
                end = "".join(i for i in end if i.isdigit())
                QUAL = row[5]
                vcffilter = row[6]
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaINV += 1
            # find DUPs
            elif "MantaDUP" in row[2]:
                end = row[7].split(';')[0]
                end = "".join(i for i in end if i.isdigit())
                QUAL = row[5]
                vcffilter = row[6]
                # add to tallys
                if "PASS" in row:
                    MantaT += 1
                    MantaDUP += 1

            # write to tsv. If first record do not add the \n
            if len(tsvout) < 1:
                # sample chr start end QUAL Filter tool change
                tsvout += sample + "\t" + row[0] + "\t" + row[1] + "\t" + end + "\t" + QUAL + "\t" + vcffilter + \
                          "\tManta\t" + varID[0]
            else:
                tsvout += "\n" + sample + "\t" + row[0] + "\t" + row[1] + "\t" + end + "\t" + QUAL + "\t" + vcffilter + \
                          "\tManta\t" + varID[0]
            print(tsvout)
        else:
            pass


tally_out = open(output, 'w')
tally_out.write("Sample\tcanvas_Total\tcanvas_Gain\tcanvas_Loss\tManta_Total\tManta_BND\tManta_DEL\tManta_INV\tMantaDUP\n")
tally_out.write(str(sample) + "\t" + str(canvasT)+ "\t" + str(canvasGain) + "\t" + str(canvasLoss)+ "\t" +
                str(MantaT) + "\t" + str(MantaBND)+ "\t" + str(MantaDEL)+ "\t" +str(MantaINV) +"\t" + str(MantaDUP) + "\n")
tally_out.close()

"""


out.close()

exit(0)
