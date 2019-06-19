__author__ = 'Colin Wu'
# This file creates ribosome density data structures.
import GFF, argparse, pysam
from Bio import Seq

class densebuilder(object):
	def __init__(self,GFFgen,bamfile1,bamfile2,densitypathandfilestring,wigpathandfile,totreads,assignment,riboshiftdict,bamfileout,softclipped):
		self.GFFgen= GFFgen
		self.bamfile1= bamfile1
		self.bamfile2= bamfile2
		self.densitypathandfilestring= densitypathandfilestring
		self.wigpathandfile= wigpathandfile
		self.totreads= totreads
		self.assignment= assignment
		self.riboshiftdict= riboshiftdict
		self.bamfileout= bamfileout
		self.softclipped= softclipped

	# For alignments using bowtie2 and output bam files 
	def setdense(self):
		outputdata1, outputdata2= {}, {} # +1 strand and -1 strand
		
		 # Create BioSeq files from the GFF file for modification of the letter annotations field.
		for sequence in self.GFFgen:
			outputdata1[sequence.id]= [0 for x in range(len(sequence))]
			outputdata2[sequence.id]= [0 for x in range(len(sequence))]

		if self.bamfile2!= '-1':	bamfiles= [self.bamfile1, self.bamfile2]
		else:	bamfiles= [self.bamfile1]

		mapped= self.setdense_5or3assignment(outputdata1,outputdata2,bamfiles,self.riboshiftdict,self.assignment,self.bamfileout,self.softclipped)	

		if self.totreads== '-1':	
			self.totreads= mapped[0]
			print "Total reads mapped: "+ str(self.totreads)
			print "Using sum of chromosome and splice mappings for normalization: "+ str(mapped[0])
		else:
			print str(self.totreads)+ " total reads entered. Using this for normalization."
		
		#Now we normalize data.
		self.norm_m(outputdata1,(self.totreads))
		self.norm_m(outputdata2,(self.totreads))
		
		# Write to file as float or pickle them...
		self.writecountsf(outputdata1, self.densitypathandfilestring+"_plus_")
		self.writecountsf(outputdata2, self.densitypathandfilestring+"_minus_")

		if self.wigpathandfile != '-1': 
			self.countstowig(outputdata1, self.wigpathandfile+"_plus")
			self.countstowig(outputdata2, self.wigpathandfile+"_minus")


		# Write output file of comments.
		fc= open(self.densitypathandfilestring+"_output.txt","w")
		fc.write("Density was built with parameters:\n")
		fc.write("riboshiftdict="+str(self.riboshiftdict)+"\n")
		fc.write("assignment="+str(self.assignment)+"\n")
		fc.write(str(mapped[0])+" total reads mapped to the genome and splice junctions."+ "\n")
		fc.write(str(mapped[1])+" reads of interested length mapped to the genome and splice junctions."+ "\n")


	###### Setdense for 5'end alignment.
	def setdense_5or3assignment(self,outputdata1,outputdata2,bamfiles,riboshiftdict,assignment,bamfileout,softclipped):
		mappedreads= 0
		mappedreadslocal= 0
		tooshortlongreads= 0
		dumppedreads= 0
		junctreads= 0
		ntadd= 0
		for bam in bamfiles:
			print "bam file is "+ bam
			bamgen= pysam.AlignmentFile(bam, "rb")

			for read in bamgen.fetch():	# Iterate through all reads
				if read.flag== 4 or read.flag== 256 or read.flag== 272:	continue # A bowtie mismatch. These are the same as those previously saved in the files bowtie output as reads without matches and reads exceeding the -m limit. No need to write as fastq.
				
				chrom= read.reference_name
				startp= int(read.pos)   # in bam files, first position is ZERO. STAR outputs soft clipped startpos

				totaljunctlen= 0
				quitsignal= 0
				mappedreads+= 1

				readlen= 0
				for element in read.cigar:
					if element[0]== 0:	# (M) alignemnt match
						readlen+= element[1]
					elif element[0]== 1: # (I) insertion to reference 
						#readlen-= element[1] 
						quitsignal= 1
					elif element[0]== 2: # (D) deletion from the reference
						#readlen+= element[1]
						quitsignal= 1
					elif element[0]== 3: # (N) skipped region from the reference, --> junction(s) in mapped reads... should be okay for ~30nt 
						totaljunctlen+= element[1]
						junctreads+= 1
					elif element[0]== 4:	# (S) soft clipping 
						if element[1] <= 2:	readlen-= element[1] 	# For nontemplated addition...
						if element[1] >2 and element[1]<= int(softclipped):	readlen+= element[1]  #pass     # Don't count clipped As #
						if element[1]> int(softclipped):	quitsignal= 1
					else:	quitsignal= 1

				if quitsignal== 1:
					dumppedreads+= 1
					continue

				if self.riboshiftdict.has_key(readlen):
					riboshift= self.riboshiftdict[readlen][0]
				else:
					tooshortlongreads+= 1
					continue

				if read.flag== 0 or read.flag== 16:	mappedreadslocal+= 1
 
				if self.assignment== '5':
					if (read.flag== 0):             # Align to the forward strand
						ribojunctlen= self.junctlen_for_riboshift(read.cigar, riboshift)  
						outputdata1[chrom][startp+ ribojunctlen+ riboshift]+= 1 

					if (read.flag== 16):            # Align to the reverse strand
						ribojunctlen= self.junctlen_for_riboshift(read.cigar[::-1], riboshift) 
						outputdata2[chrom][startp+ totaljunctlen+ (readlen- 1)- ribojunctlen- riboshift]+= 1 

				if self.assignment== '3':
					if (read.flag== 0):
						ribojunctlen= self.junctlen_for_riboshift(read.cigar[::-1], abs(riboshift))  
						outputdata1[chrom][startp+ totaljunctlen+ (readlen- 1)- ribojunctlen- abs(riboshift)]+= 1  
					
					if (read.flag== 16):
						ribojunctlen= self.junctlen_for_riboshift(read.cigar, abs(riboshift)) 
						outputdata2[chrom][startp+ ribojunctlen+ abs(riboshift)]+= 1  

				self.bamfileout.write(read)

		print str(mappedreadslocal)+" reads mapped to the genome."
		print str(junctreads)+" reads mapped to splice junctions."
		print str(tooshortlongreads)+" reads are too long/short."
		return [mappedreads,mappedreadslocal+junctreads]


	# Function to determine the combined junction length for riboshift. 
	def junctlen_for_riboshift(self,readcigar,riboshift):
		junctnum= 0
		mappedexonlen= 0
		insert= 0
		for x in range(len(readcigar)):
			if readcigar[x][0]== 0:	# Match
				mappedexonlen+= readcigar[x][1]
			if mappedexonlen-1 < riboshift:
				junctnum+= 1

		for y in range(junctnum):
			if readcigar[y][0]== 3:	# Insert
				insert+= readcigar[y][1]
		return insert


	# Divide every read value by number of mapped reads. Values then in rpm.
	# No need to return readcounts since it is a mutable type being a list.
	def norm_m(self,readcounts,reads):
	    for chrom in readcounts.keys():
	        for position in range(len(readcounts[chrom])):
	            readcounts[chrom][position]/=float(reads)
	            readcounts[chrom][position]*= 1E6     


	# Function to write counts files to disk - Use for float.
	def writecountsf(self,resultlist,filestring):  #Resultlist it is the list of counts, filestring is the file prefix for each chr to be written.
	    import struct
	    f2= open(filestring+ "keys","w")
	    for chrom in resultlist.keys():
	        f= open(filestring+ chrom, "wb")
	        for position in resultlist[chrom]:
	            f.write(struct.pack("f",position))
	        f.close()  
	        f2.write(chrom+"\n")
	    f2.close()    


if __name__== '__main__':
	parser= argparse.ArgumentParser()
	parser.add_argument('--GFFfile', help= 'coding GFF or GTF file')
	parser.add_argument('--bamfile1', help= 'bamfile from first genome mapping')
	parser.add_argument('--bamfile2', help= 'bamfile from polyA removal')
	parser.add_argument('--densitypathandfilestring', help= 'density file output path')
	parser.add_argument('--wigpathandfile', help= 'wig file output path')
	parser.add_argument('--totreads', default= -1, help= 'total reads for normalization')
	parser.add_argument('--assignment', help= '5 or 3 end', required= True)
	parser.add_argument('--riboshiftdict', help= 'dictionary of riboshifts')
	parser.add_argument('--bamfileoutput', help= 'output bam file')
	parser.add_argument('--softclipped', help= 'number of soft clipped allowed')
	args = parser.parse_args()

	import ast
	riboshiftdict= ast.literal_eval(args.riboshiftdict) #convert string into dictionary
	GFFgen= GFF.parse(args.GFFfile)
	bamgen0= pysam.AlignmentFile(args.bamfile1, "rb")
	bamfileout= pysam.AlignmentFile(args.bamfileoutput, "wb", template= bamgen0)
	rfpdense= densebuilder(GFFgen,args.bamfile1,args.bamfile2,args.densitypathandfilestring,args.wigpathandfile,args.totreads,args.assignment,riboshiftdict,bamfileout,args.softclipped)
	rfpdense.setdense()
	    
