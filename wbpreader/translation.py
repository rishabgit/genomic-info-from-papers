#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
#import argparse
#import subprocess
#import requests
#import csv
import re
#import json
#import nltk.data
#import os
#from xml.dom import minidom
#from bs4 import BeautifulSoup
#import requests
import gffutils
from pybedtools import BedTool
from Bio.Seq import Seq



def translate_dict(mode):
    """Returns a dict with protein codes in mode IUPAC or all for easy lookup of translation FROM protein to seq"""

    prot={}
    table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T','ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R','AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P','CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R','CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A','GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E','GGA':'G', 'GGC':'G', 'GGG':'G','GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L','TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_','TGG':'W',
    }

    if mode=='all':
        prot['A']={}
        prot['C']={}
        prot['D']={}
        prot['E']={}
        prot['F']={}
        prot['G']={}
        prot['H']={}
        prot['I']={}
        prot['K']={}
        prot['L']={}
        prot['M']={}
        prot['N']={}
        prot['P']={}
        prot['Q']={}
        prot['R']={}
        prot['S']={}
        prot['T']={}
        prot['V']={}
        prot['W']={}
        prot['Y']={}
        prot['_']={}
        prot['ochre']={}
        prot['amber']={}
        prot['opal']={}
        prot['A']['GCA']=1
        prot['A']['GCC']=1
        prot['A']['GCG']=1
        prot['A']['GCT']=1
        prot['C']['TGC']=1
        prot['C']['TGT']=1
        prot['D']['GAC']=1
        prot['D']['GAT']=1
        prot['E']['GAA']=1
        prot['E']['GAG']=1
        prot['F']['TTC']=1
        prot['F']['TTT']=1
        prot['G']['GGA']=1
        prot['G']['GGC']=1
        prot['G']['GGG']=1
        prot['G']['GGT']=1
        prot['H']['CAC']=1
        prot['H']['CAT']=1
        prot['I']['ATA']=1
        prot['I']['ATC']=1
        prot['I']['ATT']=1
        prot['K']['AAA']=1
        prot['K']['AAG']=1
        prot['L']['CTA']=1
        prot['L']['CTC']=1
        prot['L']['CTG']=1
        prot['L']['CTT']=1
        prot['L']['TTA']=1
        prot['L']['TTG']=1
        prot['M']['ATG']=1
        prot['N']['AAC']=1
        prot['N']['AAT']=1
        prot['P']['CCA']=1
        prot['P']['CCC']=1
        prot['P']['CCG']=1
        prot['P']['CCT']=1
        prot['Q']['CAA']=1
        prot['Q']['CAG']=1
        prot['R']['AGA']=1
        prot['R']['AGG']=1
        prot['R']['CGA']=1
        prot['R']['CGC']=1
        prot['R']['CGG']=1
        prot['R']['CGT']=1
        prot['S']['AGC']=1
        prot['S']['AGT']=1
        prot['S']['TCA']=1
        prot['S']['TCC']=1
        prot['S']['TCG']=1
        prot['S']['TCT']=1
        prot['T']['ACA']=1
        prot['T']['ACC']=1
        prot['T']['ACG']=1
        prot['T']['ACT']=1
        prot['V']['GTA']=1
        prot['V']['GTC']=1
        prot['V']['GTG']=1
        prot['V']['GTT']=1
        prot['W']['TGG']=1
        prot['Y']['TAC']=1
        prot['Y']['TAT']=1
        prot['_']['TAA']=1
        prot['_']['TAG']=1
        prot['_']['TGA']=1
        prot['ochre']['TAA']=1
        prot['amber']['TAG']=1
        prot['opal']['TGA']=1

    elif mode=='IUPAC':
        prot['A']=['GCN'] # ACGT
        prot['C']=['TGY'] # CT
        prot['D']=['GAY'] # CT
        prot['E']=['GAR'] # AG
        prot['F']=['TTY'] # CT
        prot['G']=['GGN'] # ACGT
        prot['H']=['CAY'] # CT
        prot['I']=['ATH'] # ACT
        prot['K']=['AAR'] # AG
        prot['L']=['CTN','TTR'] # ACGT
        prot['M']=['ATG'] # only
        prot['N']=['AAY'] # CT
        prot['P']=['CCN'] # ACGT
        prot['Q']=['CAR'] # AG
        prot['R']=['AGR'] # AG
        prot['R']=['CGN'] # ACGT
        prot['S']=['AGY','TCN'] # CT
        prot['T']=['ACN'] # ACGT
        prot['V']=['GTN'] # ACGT
        prot['W']=['TGG'] # only
        prot['Y']=['TAY'] # CT
        prot['_']=['TAR'] # AG
        prot['_']=['TGA']
        prot['ochre']=['TAA']
        prot['amber']=['TAG']
        prot['opal']=['TGA']


    else :
        print ("Call subroutine translate_dict in mode IUPAC or all ")

    return prot

def sub_is(ref,varp):
    ''' This subroutine takes the referece seq and compares it to the mutated protein to figure out the substitution made '''
    refn=list(ref)
    
    prots=translate_dict('all')
    alt=prots[varp]
    mutdist={}

    # For each position in ref, calculate if there is a sub or not
    i=0
    

    for a in alt:
        mutdist[a]=0
        for p in refn:
            #i=i+i
            if p==a[i]:
                #print ("SAME",p, a[i], refn, a,i)
                i=i+1
            else:
                #print ("DIFF",p, a[i], refn, a,i)
                mutdist[a]=mutdist[a]+1
                i=i+1
        i=0    
    #print (mutdist)
    #mindist=min(mutdist, key=mutdist.get)
    temp = min(mutdist.values())
    mindist = [key for key in mutdist if mutdist[key] == temp]
    #print(mindist)
    alt=','.join(mindist)
    #ref=ref.lower()
    #alt=alt.lower()
    return (ref,alt)


def extract_data_from_snp(snp):

        try:
            #print(snp['variation'][0], snp.seqid, snp.start, snp['substitution'][0],snp['aachange'][0] ,snp['aa_position'][0] ,snp['hgvsc'][0] )
            tsc=snp['hgvsc'][0].split(':')
            mut=snp['aachange'][0]
            mut=str(mut)
            mutss=mut.split('/')
            muts= mutss[0] + snp['aa_position'][0] + mutss[1]
            val= str(snp['variation'][0]) + ':' + str(snp.seqid) + ':' + str(snp.start)  + ':' + str(snp.end)
            return (tsc[0],muts, val)

        except:
            pass



def vars_from_GFF(gff):
    '''Returns a dict with protein coding changes, by gene and AA'''

    gvar={}
    
    fn = gffutils.example_filename(gff)     
    db = gffutils.create_db(fn, ':memory:')
   
    for snp in db.features_of_type('SNP', order_by='start'):

        a=extract_data_from_snp(snp)
        if a is not None:
            if a[0] in gvar:
                gvar[a[0]][a[1]]=a[2]
            else:
                gvar[a[0]]={}
                gvar[a[0]][a[1]]=a[2]

    for snp in db.features_of_type('point_mutation', order_by='start'):
        a=extract_data_from_snp(snp)
        if a is not None:
            if a[0] in gvar:
                gvar[a[0]][a[1]]=a[2]
            else:
                gvar[a[0]]={}
                gvar[a[0]][a[1]]=a[2]
    
    return gvar


def clone_from_ts(gff,db):
    ''' Extract which clone overlaps with which gene '''

    for i in db.children(gene, featuretype='assembly_component', order_by='start'):
        print (i)

    genes = BedTool(gff)
    def clone_filter(feature):
        return feature[7] == 'INS'


    
    # Grab all the clones, and make a tmp file

    # Grab all the genes and make a tmp file

    # bedtools



#18221 SNP
#76286 point_mutation



def var_abspos(gene,ts,varelstart,db):
    '''Returns chromosome, abspos'''

    varabstart=0
    geneabstart=0
    generelpos=0
    geneori=gene.strand
    #geneori=re.split('(\t)',geneori)
    #varelstart=varelstart
    chro=0
    abspos=0

    #print ("GENE",gene,geneori)
    ts='Transcript:'+ts

    # If gene is on positive strand

    if geneori=="+":
        #print ("Gene pos")

        for i in db.children(gene, featuretype='CDS', order_by='start'):
            
            
            # Now step through the exons until you have reached further than the var

            if re.search(str(ts), str(i['Parent'])):
                geneabstart=i.start
                # Make length of exon
                exlen=int(i.stop)-int(i.start)
                #exlen=exlen * -1
                #print(str(ts), i['Parent'], i.start, i.stop, exlen, sep='\t')

                # If the varabstart is calculated
                if int(varabstart)>0:
                    #relgenepos=relgenepos+exlen
                    print ("DONE",  i.start, i.stop, exlen, generelpos, varabstart , varelstart , sep='\t' )
                    pass

                # If the length of the exon is more than the var = var is in this exon
                elif generelpos+exlen > varelstart:
                    # Calculate the abs position of the var
                    varabstart=int(i.start)+(varelstart-generelpos)
                    #varabstart=100
                    print ("HERE",  i.start, i.stop, exlen, generelpos,  varabstart, varelstart , sep='\t' )
                    #print ("HERE", varelstart,  generelpos, geneabstart, int(i.start), varelstart, varabstart,   sep='\t' )
                    #relgenepos=relgenepos+exlen
                    abspos=varabstart-3
                    chro=i.seqid

                # elif the length of the exon is less than the var = var is in a subsequent exon
                elif generelpos+exlen < varelstart:
                    print ("NEXT", i.start, i.stop, exlen, generelpos, varabstart, varelstart , sep='\t' )
                    #generelpos=generelpos+exlen
                    #varelstart=varelstart-exlen
                    #print (generelpos)
                    pass

                # Other cases
                else:
                    print("Should not happen")
                    # blah

                #Update position in gene
                generelpos=generelpos+exlen+1

            else:
                # other transcripts
                pass
                #print (str(ts),i['Parent'])

    elif  geneori=="-":
        #print ("Gene neg")
        varelstart=varelstart * -1
            # Now step through the exons until you have reached further than the var
        for i in db.children(gene, featuretype='CDS', order_by='start',reverse=True):
            
            if re.search(str(ts), str(i['Parent'])):
                geneabstart=i.start
                # Make length of exon
                exlen=int(i.stop)-int(i.start)
                #print(str(ts), i['Parent'], i.start, i.stop, exlen, sep='\t')

                # If the varabstart is calculated
                if int(varabstart)>0:
                    #relgenepos=relgenepos+exlen
                    #print ("DONE",  i.start, i.stop, exlen, generelpos, varabstart , varelstart , sep='\t' )
                    pass

                # If the length of the exon is more than the var = var is in this exon
                elif (generelpos-exlen) < varelstart:
                    # Calculate the abs position of the var
                    varabstart=int(i.stop)+(varelstart-generelpos)
                    #varabstart=100
                    #print ("HERE",  i.start, i.stop, exlen, generelpos,  varabstart, varelstart , sep='\t' )
                    abspos=varabstart
                    chro=i.seqid

                # elif the length of the exon is less than the var = var is in a subsequent exon
                elif (generelpos-exlen) > varelstart:
                    #print ("NEXT", i.start, i.stop, exlen, generelpos, varabstart, varelstart , sep='\t' )
                    pass

                # Other cases
                else:
                    print("Should not happen")
                    # blah

                #Update position in gene
                generelpos=generelpos-exlen-1



            else:
                # other transcripts
                pass
                #print (str(ts),i['Parent'])



    else:
        print("WARNING: Houston we have a problem - this gene does not have an orientation, try another gff file")
        exit(1)

    #print("RES",chro,int(abspos),"HGVSg CHROMOSOME_II:g.7861183C>T", '\n', '\n' )
    return(chro,int(abspos),geneori)




def calculate_flanks(varstv,fas,prots,oriprot):
    '''Returns chromosome, abspos'''
    varseq='n'
    left='n'
    right='n'

    #seq = Seq("TCGGGCCC")
    #seq=seq.reverse_complement()


    if varstv['ORI']=='+':
        #print("ORI PLUS")
        #for s in vars[t][v]:
        #print (t, v, s)
        if ('CHR' in varstv) and 'ABSPOS' in varstv :
                #print (t, v , vars[t][v]['CHR'], vars[t][v]['ABSPOS'] )
                varseq=fas[  varstv['CHR']  ][ varstv['ABSPOS']-1  : varstv['ABSPOS']+2 ].seq
                left_flank=fas[  varstv['CHR']  ][ varstv['ABSPOS']-31 : varstv['ABSPOS']-1 ].seq
                right_flank=fas[  varstv['CHR']  ][ varstv['ABSPOS']+2 : varstv['ABSPOS']+32 ].seq
                protein = prots[oriprot]

                if varseq in prots[oriprot]:
                    print ('MATCH',left_flank, varseq, right_flank, protein ,oriprot, varstv['ABSPOS'], sep='\t')
                    left=left_flank.lower()
                    right=right_flank.lower()
                else:
                    print ('NO MATCH', left_flank, varseq, right_flank, protein ,oriprot, varstv['ABSPOS'], sep='\t')
                    pass

        #print (v , vars[t][v['CHR'])
        #if v['CHR'] in fas:
        #    print ( v,  v['CHR'], fas[ v['CHR']][200:230].seq)
        else:
            print ("MISS MATCH",oriprot, varstv['ABSPOS'])
    elif varstv['ORI']=='-':

        #print("ORI MINUS")

        if ('CHR' in varstv) and 'ABSPOS' in varstv :
                #print (t, v , vars[t][v]['CHR'], vars[t][v]['ABSPOS'] )
                varseq=fas[  varstv['CHR']  ][ varstv['ABSPOS']  : varstv['ABSPOS']+3 ].seq
                left_flank=fas[  varstv['CHR']  ][ varstv['ABSPOS']-30 : varstv['ABSPOS'] ].seq
                right_flank=fas[  varstv['CHR']  ][ varstv['ABSPOS']+3 : varstv['ABSPOS']+33 ].seq
                protein = prots[oriprot]
                # ATTGGAGT
                # Revtrans
                varseqr=Seq(varseq)
                varseqr=varseqr.reverse_complement()
                #left=Seq(left)
                #left=left.reverse_complement()
                #right=Seq(right)
                #right=right.reverse_complement()

                if varseqr in prots[oriprot]:
                    print ('MATCH',left_flank, varseq, right_flank, protein ,oriprot, sep='\t')
                    left=left_flank.lower()
                    right=right_flank.lower()
                else:
                    print ('NO MATCH', left_flank, varseq, right_flank, protein ,oriprot, sep='\t')
                    pass
        



        pass
    else:
        print("ERROR: A gene is missing orientation in the GFF file")

    return (varseq,left,right)


