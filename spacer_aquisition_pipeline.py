#!/usr/bin/env python3

"""Extract real spacers from 50 bp sequencing"""

#Imported modules
from subprocess import call
import os
from Bio import SeqIO
from Bio import Seq
import sys
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
import math
import random


#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2018, Laboratory of Dr. Kimberley Seed, UC Berkeley"
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Angus Angermeyer"
__email__ = "angermeyer@berkeley.edu"
__status__ = "Prototype"


#Input path of working directory
path = ""

#Make sure the reference file is in a folder called reference in the path
#This file is the reference for bowtie mapping. You can include any regions of interest that you would like to map to.
#The included reference file is chrI, chrII, PLE1 (modified to mutate the original lead spacer), and the ICP1 background genome.
index = path+"reference/reference.fasta"

os.chdir(path)

for filename in os.listdir(path):
	if filename.endswith("R1.fastq"):

		trim_log = path+filename[:filename.rindex(".fastq")]+"_trimlog.txt"
		output = path+filename[:filename.rindex(".fastq")]+"_trimmed.fastq"
		conditions = "MINLEN:50 AVGQUAL:20"
		
		trim = f"trimmomatic SE -quiet -threads 4 -phred33 -trimlog {trim_log} {filename} {output} {conditions}"
		print(f"Trimming:\t{filename}\t{conditions}")
		call(trim, shell=True)


for filename in os.listdir(path):
	if filename.endswith("_trimmed.fastq"):
		log = path+filename[:filename.rindex("_trimmed.fastq")]+"_log.log"
		spacers = []
		spacers34 = []
		
		pre_filtered = 0
		bp31_spacers = 0
		bp32_spacers = 0
		bp33_spacers = 0
		bp34_spacers = 0
		
		motif_rejects = 0
		tract_rejects = 0
		tail_rejects = 0
		N_rejects = 0
		
		print(f"Parsing spacers: {filename}")
		
		for record in SeqIO.parse(filename, "fastq"):
			pre_filtered += 1
			sequence = str(record.seq)

			if sequence[0:16] == "ATAGGCTGCTTAAAGA":
				
				if "CCCCCCCC" in sequence or "GGGGGGGG" in sequence:
					tract_rejects += 1
					continue
				
				else:		
					if "N" in sequence:
						N_rejects += 1
						continue
					
					else:
						#rename header to only the non-redundant parts (can still be back-tracked later if needed)
						header = record.id.split()
						unique_id =  header[0]+str(random.randint(1,10000))
						record.id = unique_id
						record.description = unique_id
					
						#31bp spacers
						if sequence[-4:] == "GTTA":
							spacers.append(sequence[16:sequence.rindex("GT")])

							letter_annotations = record.letter_annotations["phred_quality"][16:sequence.rindex("GT")]
							record.letter_annotations = {}
							record.seq = Seq.Seq(sequence[16:sequence.rindex("GT")])
							record.letter_annotations = {'phred_quality':letter_annotations}
							bp31_spacers += 1

							with open(filename[:filename.rindex(".fastq")]+"_spacers.fastq", "a") as output:
								SeqIO.write(record, output, "fastq")
					
						#32bp spacers
						elif sequence[-3:] == "GTT":
							spacers.append(sequence[16:sequence.rindex("GT")])

							letter_annotations = record.letter_annotations["phred_quality"][16:sequence.rindex("GT")]
							record.letter_annotations = {}
							record.seq = Seq.Seq(sequence[16:sequence.rindex("GT")])
							record.letter_annotations = {'phred_quality':letter_annotations}
							bp32_spacers += 1
						
							with open(filename[:filename.rindex(".fastq")]+"_spacers.fastq", "a") as output:
								SeqIO.write(record, output, "fastq")						
					
						#33bp spacers
						elif sequence[-2:] == "GT":
							spacers.append(sequence[16:sequence.rindex("GT")])

							letter_annotations = record.letter_annotations["phred_quality"][16:sequence.rindex("GT")]
							record.letter_annotations = {}
							record.seq = Seq.Seq(sequence[16:sequence.rindex("GT")])
							record.letter_annotations = {'phred_quality':letter_annotations}
							bp33_spacers += 1
						
							with open(filename[:filename.rindex(".fastq")]+"_spacers.fastq", "a") as output:
								SeqIO.write(record, output, "fastq")						
							
						
						elif sequence[-1:] == "G":
							spacers34.append(sequence[16:-1])
						
							letter_annotations = record.letter_annotations["phred_quality"][16:-1]
							record.letter_annotations = {}
							record.seq = Seq.Seq(sequence[16:-1])
							record.letter_annotations = {'phred_quality':letter_annotations}
							bp34_spacers += 1
						
							with open(filename[:filename.rindex(".fastq")]+"_spacers34.fastq", "a") as output:
								SeqIO.write(record, output, "fastq")
							
						else:
							tail_rejects += 1
				
			else:
				motif_rejects += 1
				continue
			
		with open(log, "a") as log:
			log.write("Filtering results for: "+filename+"\n")
			log.write("Pre_filtered reads:\t"+str(pre_filtered)+"\n")
			log.write("motif_rejects:\t"+str(motif_rejects)+"\n")
			log.write("tract_rejects:\t"+str(tract_rejects)+"\n")
			log.write("N_rejects:\t"+str(N_rejects)+"\n")
			log.write("tail_rejects:\t"+str(tail_rejects)+"\n\n")

			log.write("34bp_spacers:\t"+str(bp34_spacers)+"\n")
			log.write("31bp_spacers:\t"+str(bp31_spacers)+"\n")
			log.write("32bp_spacers:\t"+str(bp32_spacers)+"\n")
			log.write("33bp_spacers:\t"+str(bp33_spacers)+"\n")
			log.write("total(31-33)_spacers:\t"+str(len(spacers))+"\n\n\n")
			
		



#####Mapping with bowtie
#Before this step you need to create a bowtie index:
#bowtie-build reference.fasta reference_bowtie_index
for filename in os.listdir(path):
	if filename.endswith("spacers.fastq"): #Ignoring 34bp spacers currently
		output = filename[:filename.rindex(".fastq")]+"_bowtie.txt"
		
		index_file = "reference/CRISPR_combined"
		
		seed_length = 31
		max_total_mismatches = 0
		max_alignments = 1
				
		log = path+filename[:6]+"_log.log"
		
		#Choose one of the following mapping approaches:
		
		#Bowtie mapping limiting to max_alignments number of hits per spacer
# 		bowtie_call = f"bowtie -a -l {seed_length} -v {max_total_mismatches} -m {max_alignments} {index_file} -q {filename} {output} >> {log} 2>&1"
		
		#Bowtie mapping allowing all spacer matches
		bowtie_call = f"bowtie -a -l {seed_length} -v {max_total_mismatches} {index_file} -q {filename} {output} >> {log} 2>&1"
		
		with open(log, "a") as log:
			log.write("Mapping for file:\t"+filename+"\n")
			log.write(bowtie_call+"\n")
		
		print(f"bowtie mapping vs reference: {index}")
		call(bowtie_call, shell=True)
		



###Generate weblogos and flanking and spacer fasta files#####
for filename in os.listdir(path):
	if filename.endswith("bowtie.txt"):
		print("Parsing bowtie mappings:\t"+filename)

		mapped_spacer_count=0
		
		chromosome_1 = []
		chromosome_2 = []
		PLE1_mutated = []
		ICP1_2011_A  = []
		ICP1_2011_A_wo_array = []

		with open(filename, "r") as maps:
			for line in maps:
				mapped_spacer_count+=1		
				line = line.split()
				
				spacer_ident    = line[0]
				spacer_strand   = line[1]
				spacer_location = int(line[3])
				spacer_sequence = line[4]
				spacer_length   = len(spacer_sequence)
				spacer_count    = 1
				
				spacer_in_list  = False
				
				spacer_info = [spacer_sequence, spacer_location, spacer_strand, spacer_length, spacer_count]
				
				
				if 0 <= spacer_location <= 3070910:
					
					index_pos = 0
					for spacer in chromosome_1:
						if spacer[0] == spacer_sequence and spacer[1] == spacer_location:
							chromosome_1[index_pos][4] +=1
							spacer_in_list = True
							break
						index_pos += 1
					if spacer_in_list == False:
						chromosome_1.append(spacer_info)


				elif 3070920 <= spacer_location <= 4162285:

					index_pos = 0
					for spacer in chromosome_2:
						
						if spacer[0] == spacer_sequence and spacer[1] == spacer_location:
							chromosome_2[index_pos][4] +=1
							spacer_in_list = True
							break
						index_pos += 1
					if spacer_in_list == False:
						chromosome_2.append(spacer_info)


				elif 4162320 <= spacer_location <= 4180365:

					index_pos = 0
					for spacer in PLE1_mutated:
						
						if spacer[0] == spacer_sequence and spacer[1] == spacer_location:
							PLE1_mutated[index_pos][4] +=1
							spacer_in_list = True
							break
						index_pos += 1
					if spacer_in_list == False:
						PLE1_mutated.append(spacer_info)


				elif 4180380 <= spacer_location <= 4306730:

					index_pos = 0
					for spacer in ICP1_2011_A_wo_array:
						
						if spacer[0] == spacer_sequence and spacer[1] == spacer_location:
							ICP1_2011_A_wo_array[index_pos][4] +=1
							spacer_in_list = True
							break
						index_pos += 1
					if spacer_in_list == False:
						ICP1_2011_A_wo_array.append(spacer_info)


		reference_list = [["chromosome_1",chromosome_1, 3070910],["chromosome_2",chromosome_2, 1091365],["PLE1_mutated",PLE1_mutated, 18045], ["ICP1_2011_A_wo_array", ICP1_2011_A_wo_array, 126350]]


		record = SeqIO.read(index, "fasta")
		for reference in reference_list:
			y=1
			with open(filename[:filename.rindex("_trimmed_spacers")]+f"_{reference[0]}_flanking.fasta","w") as fasta_flanking:
				pam_list=[]
				pam_count=[]
				pam_cumulative_count=[]
				length_list = [0,0,0]
				plus_count = 0
				minus_count = 0
				plus_cumulative = 0
				minus_cumulative = 0
				
				for loc in reference[1]:
					strand   = loc[2]
					site     = loc[1]
					length   = loc[3]
					count    = loc[4]
					sequence = loc[0]
				
					length_list[length-31]+=1*count
						
					
					if strand == "+":
						plus_count += 1
						plus_cumulative += count
						pam = str(record.seq[site-2:site].reverse_complement())
						flanking = record.seq[site-10:site] +"NNN"+ record.seq[site+length:site+length+10]
						flanking = flanking.reverse_complement()

						fasta_flanking.write(f">{reference[0]}_{length}_{strand}_{site}_{count}_{y}\n{flanking}\n")

						if pam in pam_list:
							pam_count[pam_list.index(pam)]+=1
							pam_cumulative_count[pam_list.index(pam)]+=count
						else:
							pam_list.append(pam)
							pam_count.append(1)
							pam_cumulative_count.append(count)



					elif strand == "-":
						minus_count += 1
						minus_cumulative += count
						pam = str(record.seq[site+length:site+length+2])
						flanking = record.seq[site-10:site] +"NNN"+ record.seq[site+length:site+length+10]

						fasta_flanking.write(f">{reference[0]}_{length}_{strand}_{site}_{count}_{y}\n{flanking}\n")
					
						if pam in pam_list:
							pam_count[pam_list.index(pam)]+=1
							pam_cumulative_count[pam_list.index(pam)]+=count
						else:
							pam_list.append(pam)
							pam_count.append(1)
							pam_cumulative_count.append(count)
					


					
					with open(filename[:filename.rindex("_trimmed_spacers")]+f"_{reference[0]}_spacers.fasta","a") as fa_out:
						fa_out.write(f">{record.id}_{count}_{length}_{strand}_{site}_{y}\n{sequence}\n")
				
					y+=1
			
			
			log = path+filename[:6]+"_log.log"
			with open(log, "a") as log:
				log.write("\n\n"+reference[0]+"\n")
				log.write("PAM\tRep\tAll\n")
				for i in range(len(pam_list)):
					log.write(str(pam_list[i])+"\t"+str(pam_count[i])+"\t"+str(pam_cumulative_count[i])+"\n")
				log.write("\n")
				log.write("+\t"+str(plus_count)+"\t"+str(plus_cumulative)+"\n")
				log.write("-\t"+str(minus_count)+"\t"+str(minus_cumulative)+"\n")
				log.write("Tot\t"+str(sum(pam_count))+"\t"+str(sum(pam_cumulative_count))+"\n")
				z=31
				for count in length_list:
					log.write(str(z)+"bp spacer count:\t"+ str(count)+"\n")
					z+=1



			#Generate flanking weblogos
			file = filename[:filename.rindex("_trimmed_spacers")]+f"_{reference[0]}_flanking.fasta"
			title = "Sample:"+filename[:filename.rindex(".R1")]+"\ Target:" + reference[0]
			file_out = file+".eps"
			
			annotation = "10,,,,,5,,,,1,[,SPACER,],1,,,,,5,,,,10"
			
			weblogo = f"weblogo -A dna -D fasta -E YES -P '' --errorbars NO -U bits -t {title} --title-fontsize 10 --number-fontsize 4 --annotate '{annotation}' -f {file} -o {file_out}"
			call(weblogo, shell=True)





#####Output *.graph files for brig#####
index_records = list(SeqIO.parse(index, "fasta"))

	
for filename in os.listdir(path):
	log = path+filename[:6]+"_log.log"
	if filename.endswith("spacers.fasta"):
		ref = filename[filename.index("R1")+3:filename.rindex("_spacers.fasta")]
		for ref_record in reference_list:
			if ref_record[0] == ref:
				template_len = ref_record[2]
				pos_list = []
				pos_list2 = []
				
				for i in range(template_len):
					pos_list.append(0.0)
					pos_list2.append(0.0)
				
				with open(filename, "r") as spacer_file:
					count1 = 0
					count2 = 0
					for spacer_record in SeqIO.parse(spacer_file, "fasta"):
						spacer_length = int(spacer_record.id.split("_")[-4])
						strand = spacer_record.id.split("_")[-3]
						count = int(spacer_record.id.split("_")[-5])
						
						site =  int(spacer_record.id.split("_")[-2])
						if ref_record[0] == "chromosome_1":
							site = site-0
							
						elif ref_record[0] == "chromosome_2":
							site = site-3070920
							
						elif ref_record[0] == "PLE1_mutated":
							site = site-4162320
							
						elif ref_record[0] == "ICP1_2011_A_wo_array":
							site = site-4180380
							

						if strand == "+":
							for j in range(spacer_length):
								pos_list[site+j]+=1*count
							count1+=1
						
						elif strand == "-":
							for j in range(spacer_length):
								pos_list2[site-j]+=1*count
							
							count2+=1	
								
					with open(log, "a") as log:
						log.write("\n\n"+ref_record[0]+"\n")
						log.write("Max + strand spacers:\t"+str(max(pos_list))+"\n")
						log.write("Max - strand spacers:\t"+str(max(pos_list2))+"\n")

				
				with open(filename+"_for.graph", "w") as graph_out:
					x=0
					for pos in pos_list:
						graph_out.write(f"{str(x)}\t{str(x+1)}\t{pos}\n")
						x+=1
						
				with open(filename+"_rev.graph", "w") as graph_out:
					x=0
					for pos2 in pos_list2:
						graph_out.write(f"{str(x)}\t{str(x+1)}\t{pos2}\n")
						x+=1


######Create brigD3 images#####
for filename in os.listdir(path+"reference/"):
	if filename.endswith(".gbk"):
		record_name = path+"reference/"+filename

		ref_name = filename[:filename.index(".gbk")]

		
		for filename in os.listdir(path):
			if filename.endswith(".graph") and ref_name in filename and "_for" in filename:
				record = SeqIO.read(record_name, "genbank")
				record_ori = SeqIO.read(record_name, "genbank")
				
				scale_height = 0
				
				if "ICP1" in filename:
					scale_height = 10
				elif "chromosome" in filename:
					scale_height = 1000
				else:
					scale_height = 50000
				
				for_data = []
				with open(filename, "r") as for_graph_data:
					for line in for_graph_data:
						line=line.split()
						for_data.append((int(line[1]),float(line[2])))
						if float(line[2]) > scale_height:
							scale_height=float(line[2])
					for_data.append((len(record.seq)+5,scale_height))		
						
					

				rev_data = []
				with open(filename.replace("_for.graph", "_rev.graph"), "r") as rev_graph_data:
					for line in rev_graph_data:
						line=line.split()
						rev_data.append((int(line[1]),float(line[2])))
						if float(line[2]) > scale_height:
							scale_height=float(line[2])
					rev_data.append((len(record.seq)+5,scale_height))
				print(filename)
				print("Max height: "+ str(scale_height))
				
				if "ICP1" in filename:
					scale_height = 10
				elif "chromosome" in filename:
					scale_height = 1000
				else:
					scale_height = 50000
					
					
					
					
				
				gd_diagram = GenomeDiagram.Diagram("temp")
				gd_track_for_features = gd_diagram.new_track(5, name="Annotated Features", height=1, scale_ticks=0, greytrack=1, greytrack_labels=0, scale=0)
				gd_feature_set = gd_track_for_features.new_set(type='feature')

				gd_track_for_feature_names = gd_diagram.new_track(2, name="Annotated names", height=1, scale_ticks=1, greytrack=0, greytrack_labels=0, scale=1, hide=1)
				gd_feature_names_set = gd_track_for_feature_names.new_set(type='feature')

				gd_track_for_feature_scale = gd_diagram.new_track(1, name="Annotated scale", hide=1, height=1, scale_ticks=1, greytrack=0, greytrack_labels=0, scale=1, scale_color=(0,0,0), scale_largetick_interval=20000, scale_largeticks=0.5, scale_smalltick_interval=100000)
				gd_feature_scale_set = gd_track_for_feature_scale.new_set(type='feature')		

				#ORF name lines
				for feature in record.features:
					if feature.type != "CDS":
						continue
					feature.strand=-1
					
					if "chromosome" in filename:
						gd_feature_set.add_feature(feature, color=(245,245,245), label=False, sigil="BIGARROW", label_position="middle", label_size=10) #Color change from 254 (white) to ring grey
					elif "ICP1" in filename:
						gd_feature_set.add_feature(feature, color=(245,245,245), label=True, sigil="BIGARROW", label_position="middle", label_size=6)
					else:
						gd_feature_set.add_feature(feature, color=(245,245,245), label=True, sigil="BIGARROW", label_position="middle", label_size=17) #Color change from 254 (white) to ring grey


				#ORF arrows
				record = record_ori
	
				for feature in record.features:
					if feature.type != "CDS":
						continue
					if feature.strand == 1: 
						color = (255,255,255)
					else:
						color = (255,255,255)
					
					if "chromosome" in filename:
						gd_feature_set.add_feature(feature, color=color, label=True, sigil="BIGARROW", arrowshaft_height=1.0, label_position="middle", label_size=10, name="", arrowhead_length=0)
					else:
						gd_feature_set.add_feature(feature, color=color, label=True, sigil="BIGARROW", arrowshaft_height=1.0, label_position="middle", label_size=10, name="")
					
					feature.strand=1


				gd_track_for_graph = gd_diagram.new_track(6, name="graph", height=5, scale_ticks=0, greytrack=0, greytrack_labels=0, scale=0)
				gd_graph_set = gd_track_for_graph.new_set(type='graph')
				if "chromosome" in filename: 
					gd_graph_set.new_graph(for_data, name="graph test", style="line", color=colors.red, altcolor=colors.white, linewidth=.2)
				else:
					gd_graph_set.new_graph(for_data, name="graph test", style="line", color=colors.red, altcolor=colors.white, linewidth=1)

				gd_graph_set = gd_track_for_graph.new_set(type='graph')
				
				if "chromosome" in filename:
					gd_graph_set.new_graph(rev_data, name="graph test2", style="line", color=colors.blue, altcolor=colors.white, linewidth=.2)
				else:
					gd_graph_set.new_graph(rev_data, name="graph test2", style="line", color=colors.blue, altcolor=colors.white, linewidth=1)


				scale_bar_data = []

				for num in range(len(record)):
					if num == 0:
						scale_bar_data.append((num,scale_height))
					elif num == round(len(record)/4):
						scale_bar_data.append((num,scale_height))
					elif num == round(len(record)/2):
						scale_bar_data.append((num,scale_height))
					elif num == round(len(record)*3/4):
						scale_bar_data.append((num,scale_height))
					else:
						scale_bar_data.append((num,0))

				gd_graph_set = gd_track_for_graph.new_set(type='graph')
				gd_graph_set.new_graph(scale_bar_data, name="scale bar", style="line", color=colors.black, altcolor=colors.white)
				
				try:
					gd_diagram.draw(format="circular", circular=True, pagesize=(50*cm,50*cm), start=0, end=len(record), circle_core=.3, fragments=1)
					gd_diagram.write(filename+"_"+ref_name+"_spacer_map.pdf", "pdf")

				except:
					print(sys.exc_info())
					pass

				name = filename+"_"+ref_name+"_spacer_map.pdf"
				convert_call = f"convert -density 1200 {name} {name[:name.index('.pdf')]}.jpg"
				call(convert_call, shell=True)
				
				

weblogo_group = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -r300 -dDEVICEWIDTHPOINTS=300 -dDEVICEHEIGHTPOINTS=120 -sOutputFile=weblogos_motifs.pdf *.eps"
call(weblogo_group, shell=True)

tar_call = "tar -c -f spacer_and_flanking_sequences.tar *.fasta"
call(tar_call, shell=True)

tar_call = "tar -czvf graph_files_for_brig.tar *.graph"
call(tar_call, shell=True)

for filename in os.listdir(path):
	if ".eps" in filename or filename.endswith(".graph") or filename.endswith("flanking.fasta") or filename.endswith("spacers.fasta"):
		os.remove(filename)


