import os
from Bio import Entrez,SeqIO,Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pylab
import math
from Bio.SeqUtils import GC
from cLev import dist

from Primerselectiontools_py3 import *

def printPrimerSequences(primers):
    for primer in primers:
        for prim in primer:
            print(str(prim))
        print("--------------------")

def printNumOfPrimers(primers):
    numoftotprimers = 0
    for primer in primers:
        numoftotprimers += len(primer)
    print("Total 12mers (unique): " + str(len(primers)) + " (" + str(numoftotprimers) + ")")

def parsePrimerFile(filename,output_primer_length,inioligolength):
    #enter all possible 25mers in a 2d array
    records = SeqIO.parse(open(filename), "fasta")
    primers = []
    for record in records:
        recseq = record.seq
        possibleprimers = []
        #choose all possible 12mers that can be made from each 20mer
        for i in range(inioligolength-output_primer_length):
            tempseq = recseq[i:output_primer_length+i]
            possibleprimers.append(tempseq)
        primers.append(possibleprimers)
    return primers

def restrictionEnzymeFilter(filename, primers):
    records = SeqIO.parse(open(filename), "fasta")
    seqtoremovelist = []
    for record in records:
        seqtoremovelist.append(record.seq)
        
    newprimers = []
    for primer in primers:
        newprimer = []
        for prim in primer:
            keep = 1
            for seqtoremove in seqtoremovelist:
                if not prim.find(seqtoremove)==-1:
                    keep = -1
            if keep==1:
                newprimer.append(prim)        
                
        if len(newprimer)>0:
            newprimers.append(newprimer)
    
    return newprimers

def restrictionEnzymeFilterExt(filename, primers):
    records = SeqIO.parse(open(filename), "fasta")
    seqtoremovelist = []
    for record in records:
        seqtoremovelist.append(record.seq)     
    newprimers = []
    for primer in primers:
        newprimer = []
        for prim in primer:
            newprim = "GCTCTTCG" + prim + "CGAAGAGC"
            keep = 1
            for seqtoremove in seqtoremovelist:
                if not newprim.find(seqtoremove)==-1:
                    keep = -1
            if keep==1:
                newprimer.append(prim)        
                
        if len(newprimer)>0:
            newprimers.append(newprimer)
    
    return newprimers

def restrictionEnzymeFilterBspQI(filename, primers):
    records = SeqIO.parse(open(filename), "fasta")
    seqtoremovelist = []
    for record in records:
        seqtoremovelist.append(record.seq)     
    newprimers = []
    for primer in primers:
        newprimer = []
        for prim in primer:
            newprim = "GCTCTTCG" + prim + "CGAAGAGC"
            keep = 1
            for seqtoremove in seqtoremovelist:
                if not newprim.count_overlap(seqtoremove)<=1:
                    keep = -1
            if keep==1:
                newprimer.append(prim)        
                
        if len(newprimer)>0:
            newprimers.append(newprimer)
    
    return newprimers

def temperatureFilter(primers, lowtemp, hightemp):
    newprimers = []
    for primer in primers:
        newprimer = []
        for prim in primer:
            Tm = oligoTm(prim)
            if Tm >= lowtemp and Tm <= hightemp:
                newprimer.append(prim)
        if len(newprimer)>0:
            newprimers.append(newprimer)
    
    return newprimers

def selfDimers(primers, cutoff):
    newprimers=[]
    for primer in primers:
        newprimer = []
        for prim in primer:
            score = primerdimers(prim,prim)
            if not score>cutoff:
                newprimer.append(prim)
        if len(newprimer)>0:
            newprimers.append(newprimer)
    return newprimers

def secStructure(primers, cutoff):
    newprimers=[]
    for primer in primers:
        newprimer = []
        for prim in primer:
            score = calcSecondaryStructure(prim)
            if score>cutoff:
                newprimer.append(prim)
        if len(newprimer)>0:
            newprimers.append(newprimer)
    return newprimers

def GCcutoff(primers, cutoffmin, cutoffmax):
	newprimers=[]
	GCall = []
	for primer in primers:
		newprimer = []
		for prim in primer:
			score = GC(prim)
			if score>cutoffmin and score<cutoffmax:
				newprimer.append(prim)
			GCall.append(score)
		if len(newprimer)>0:
			newprimers.append(newprimer)

	generate_histogram = False
	fileGCout = "GCbarcode"
	if generate_histogram:
		print("min GC: " +str(min(GCall)))
		print("max GC: " +str(max(GCall)))
		pylab.hist(GCall, bins=60)
		pylab.title("%i oligo GC percentage\nfrom %i to %i" \
					% (len(GCall),min(GCall),max(GCall)))
		pylab.xlabel("GC content (%)")
		pylab.ylabel("Count")
		pylab.savefig(fileGCout+'_hist.png', bbox_inches='tight')
		pylab.savefig(fileGCout+'_hist.pdf', bbox_inches='tight')
		#pylab.show()
		pylab.close()
	return newprimers

def temperatureMcalc(primers, fileTMout):
	primersTM = []
	TMall = []
	for primer in primers:
		newprimerTM = []
		for prim in primer:
			Tm = oligoTm(prim)
			TMall.append(Tm)
	
# 	listfile = open(fileTMout+'.txt', 'w')
# 	for item in TMall:
# 		listfile.write("%s\n" % item)
# 	listfile.close()
	
	generate_histogram = True
	if generate_histogram:
		print("min TM: " +str(min(TMall)))
		print("max TM: " +str(max(TMall)))
		pylab.hist(TMall, bins=40)
		pylab.title("%i oligo melting temps\nfrom %i to %i" \
					% (len(TMall),min(TMall),max(TMall)))
		pylab.xlabel("Melting temp (C)")
		pylab.ylabel("Count")
		pylab.savefig(fileTMout+'hist.png', bbox_inches='tight')
		pylab.savefig(fileTMout+'hist.pdf', bbox_inches='tight')
		#pylab.show()
		pylab.close()

def modLev_cutoff(primers, cutoffmin, generate_histogram):
	newprimers=[]
	Lev_dist_all = []
	
	count_primers = 0
	del_list = []
	#loop over all primers
	for primer in primers:
		newprimer = []
		#loop over each sub-primer
		for prim in primer:
			count_primers2 = 0
			Lev_dist_temp = []
			#loop over each primer
			for primer2 in primers:
				#dont compare the same primers
				if count_primers != count_primers2:
					#loop over sub-primers
					for prim2 in primer2:
						if not (prim2 in del_list):
							#determine the Lev distance
							#score = dist(str(str(prim), "utf-8"), str(str(prim2), "utf-8"))
							score = dist(str(prim), str(prim2))
							Lev_dist_all.append(score)
							Lev_dist_temp.append(score)
				count_primers2 += 1
			#for each pair. make sure this primer is greater distance than cutoff
			if min(Lev_dist_temp)>=cutoffmin:
				newprimer.append(prim)
			else:#otherwise don't use it in future
				del_list.append(prim)
		if len(newprimer)>0:
			newprimers.append(newprimer)
		count_primers += 1
		#print(str(count_primers) + " checked, "+str(len(primers)-count_primers)+" remain")
	
	generate_histogram = False
	fileGCout = "Lev_dist"
	if generate_histogram:
		print("min Lev dist: " +str(min(Lev_dist_all)))
		print("max Lev dist: " +str(max(Lev_dist_all)))
		pylab.hist(Lev_dist_all, bins=20)
		pylab.title("%i oligo Lev_dist\nfrom %i to %i" \
					% (len(Lev_dist_all),min(Lev_dist_all),max(Lev_dist_all)))
		pylab.xlabel("Lev dist")
		pylab.ylabel("Count")
		pylab.savefig(fileGCout+'_hist.png', bbox_inches='tight')
		pylab.savefig(fileGCout+'_hist.pdf', bbox_inches='tight')
		#pylab.show()
		pylab.close()
	return newprimers

def modLev_plot(primers,fileLevout):
	newprimers=[]
	Lev_dist_all = []
	
	count_primers = 0
	del_list = []
	for primer in primers:
		count_primers2 = 0
		for primer2 in primers:
			if (count_primers != count_primers2):
				score = dist(str(primer[0]), str(primer2[0]))
				Lev_dist_all.append(score)
			count_primers2 += 1
		count_primers += 1

	generate_histogram = True
	if generate_histogram:
		print("min Lev dist: " +str(min(Lev_dist_all)))
		print("max Lev dist: " +str(max(Lev_dist_all)))
		pylab.hist(Lev_dist_all, bins=20)
		pylab.title("%i oligo Lev_dist\nfrom %i to %i" \
					% (len(Lev_dist_all),min(Lev_dist_all),max(Lev_dist_all)))
		pylab.xlabel("Lev dist")
		pylab.ylabel("Count")
		pylab.savefig(fileLevout+'_hist_after_filt.png', bbox_inches='tight')
		pylab.savefig(fileLevout+'_hist_after_filt.pdf', bbox_inches='tight')
		#pylab.show()
		pylab.close()
	return newprimers

##########################################
#this is the list of primers to scan through
input_sequences = "unused_elledge_25mers.fasta"
#unused_elledge_25mers_first_50k.fasta
#unused_elledge_25mers_first_24k.fasta
#unused_elledge_25mers
#this input files is unused elledge primers (no skpp20, or skpp15, or, em284, or em1536)

#enter all possible 25mers in a 2d array
# the above file is from Elledge's paper "Design of 240,000 orthogonal 25mer DNA barcode probes": http://www.pnas.org/content/106/7/2289.full
# see http://elledgelab.bwh.harvard.edu/Barcode/

#this file contains restriction site sequences
restriction_site_file = "restriction_minimal_Fall2021.fasta"

#same as above but without BspQI, BspQI (rc)
restriction_site_file2 = "restriction_minimal_Fall2021_noBspQI.fasta"

#output length
output_primer_length = 12
#input length
inioligolength = 25
cutoff_GC_min = 45
cutoff_GC_max = 55

lowtemp = 38
hightemp = 44
selfdimercutoff = 4

Lev_dist_check = False

#######################
#load the primers in:
primers = parsePrimerFile(input_sequences,output_primer_length,inioligolength)
print("Initial Setup")
printNumOfPrimers(primers)

#filter out primers with restriction sites
primers = restrictionEnzymeFilter(restriction_site_file, primers)

#after adding flanking Nt.BspQI sites, new restriction sites may be introduced
#filter out primers + Nt.BspQI sites that now have restriction sites
primers = restrictionEnzymeFilterExt(restriction_site_file2, primers)

#primers should only have 1 BspQI, 1 BspQI (rc) - filter out those that have more
# this file contains only BspQI, BspQI (rc)
filename3 = "BspQI.fasta"
#filter out primers + Nt.BspQI sites that now have extra BspQI sites
primers = restrictionEnzymeFilterBspQI(filename3, primers)

#how many primers are left?
print("Restriction Enzyme Filter")
printNumOfPrimers(primers)

#check GC content and filter out primers outside
primers = GCcutoff(primers, cutoff_GC_min, cutoff_GC_max)
print("GC content Filter")
printNumOfPrimers(primers)

#plot the Tm distribution before filter
fileTMout = "oligo_Tm_list"
temperatureMcalc(primers, fileTMout)

#filter primers based on Tm
primers = temperatureFilter(primers, lowtemp, hightemp)
print("Melting Temp Filter")
printNumOfPrimers(primers)

#plot the Tm distribution after filter
fileTMout = "oligo_Tm_list_filtered"
temperatureMcalc(primers, fileTMout)

#filter for self-dimers
primers = selfDimers(primers,selfdimercutoff)
print("Self Dimerization Cutoff")
printNumOfPrimers(primers)

# # this was commented
# secstructcutoff = -2
# primers = secStructure(primers,secstructcutoff)
# print "Secondary Structure Cutoff"
# printNumOfPrimers(primers)
if Lev_dist_check == True:
	# check Lev. distance
	#dist_min = 2
	#primers = modLev_cutoff(primers, dist_min,True)
	#modLev_plot(primers,"Lev_dist_after_2")
	#print("Mod Levenshtein Filter")
	#printNumOfPrimers(primers)

	dist_min = 3
	primers = modLev_cutoff(primers, dist_min,False)
	modLev_plot(primers,"Lev_dist_after_3")
	print("Mod Levenshtein Filter")
	printNumOfPrimers(primers)

#Output Primers
if Lev_dist_check == True:
	filename = "filt_prim_"+str(output_primer_length)+"nt_Lev_"+str(dist_min)+"_Tm_"+str(lowtemp)+"_"+str(hightemp)+"_GC_"+str(cutoff_GC_min)+"_"+str(cutoff_GC_max)+"_SD_"+str(selfdimercutoff)+".fasta"
else:
	filename = "filt_prim_"+str(output_primer_length)+"nt_Tm_"+str(lowtemp)+"_"+str(hightemp)+"_GC_"+str(cutoff_GC_min)+"_"+str(cutoff_GC_max)+"_SD_"+str(selfdimercutoff)+".fasta"

fileout = open(filename,'w')
for i in range(len(primers)):
    fileout.write(">cp12mer-"+str(i+1)+"\n")
    fileout.write(str(primers[i][0])+"\n")

fileout.close()
