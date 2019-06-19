import GFF, csv
from Bio.Seq import Seq

def writegene_wf2(shift5,shift3,riboshift,densityfile,feature,gfffile,utrgfffilename,outfile):
	GFFgen= GFF.parse(gfffile)
	counts1= readcountsf(densityfile+"_plus_")
	counts2= readcountsf(densityfile+"_minus_")
	counts= [counts1,counts2]
	idtable= makeidtable2(GFFgen)
	GFFgen= GFF.parse(gfffile)
	GFFlist= makeGFFlist(GFFgen)
	goodgenes= 2
	print feature
	if utrgfffilename== "-1":	utrtable= {}
	else:
		utrtable= utrgffgen=GFF.parse(utrgfffilename)
		utrtable= genometools.makeutrtable(utrgffgen)
	bp1= shift5
	chromosome= idtable[feature][2]
	featurenum= idtable[feature][1]
	longfeature= idtable[feature][0].id

	if utrtable.has_key(longfeature):	bp2= utrtable[longfeature][1]-utrtable[longfeature][0]+ shift3
	else:	bp2= 0+shift3

	bp= [bp1,bp2,riboshift]
	retval= givegene(chromosome,featurenum,GFFlist,counts,bp,goodgenes)
	if retval[0]== -1:	print "Not a good gene..."

	t= []
	t.append(["pos", "rpm"])
	i= -shift5
	while i < len(retval[0])-shift5:
		newline= [i,retval[0][i+shift5]]
		t.append(newline)
		i+=1

	fcsv= open(outfile+"_"+feature+".csv", "w")
	writer = csv.writer(fcsv)
	writer.writerows(t)
	fcsv.close()

def readcountsf(filestring):
    import struct
    keys= []
    resultlist= {}
    f2= open(filestring+"keys","r")
    for line in f2:	keys.append(line.rstrip('\n')) 
    for chrom in keys:
        resultlist[chrom]=[]
        with open(filestring+chrom,"rb") as f:
            nextval = f.read(4) 
            while nextval != "": 
                resultlist[chrom].append(struct.unpack("f",nextval)[0])
                nextval = f.read(4)
	f2.close()
    return resultlist


def makeGFFlist(GFFgen):
	GFFlist={}
	for chr in GFFgen:	GFFlist[chr.id]=chr
	return GFFlist


def makeidtable2(GFFgen):
	idtable= {}
	for chrom in GFFgen:
		if chrom.id== "gi|49175990|ref|NC_000913.2|":	keyname="Name"
		else:	keyname= "Alias"
		num_feat= 0
		for feature in chrom.features:
			name_feat= feature.id
			for item in feature.qualifiers:
				if keyname in feature.qualifiers:
					alias_feat= feature.qualifiers[keyname][0]
					idtable[alias_feat]= [feature,num_feat,chrom.id]
			num_feat+=1
	return idtable


def givegene(chromosome,feature,GFFs,counts,bp,goodgenes):
	if type(GFFs)==list:
		GFFlist=GFFs[0]
		utrtable5=GFFs[1]
		utrtable3=GFFs[2]
	else:
		GFFlist=GFFs
		utrtable5={}
		utrtable3={}		
	if type(bp)!=int:
		if len(bp)==2:
			bp.append(0)
	else:
		bp=[bp,bp,0]
	
	if bp[0]<0 or bp[1]<0:
		print "Error - bp is negative!"
		return [-1,-1]
	
	if(GFFlist[chromosome].features[feature].strand== 1):	strand=0
	else:	strand=1
	
	if goodgenes > 0:
		#if (GFFlist[chromosome].id == 'chrMito' or GFFlist[chromosome].id == '2-micron'):
		if GFFlist[chromosome].id == '2-micron':	return [-1,-1]
		if (GFFlist[chromosome].features[feature].type!='gene'):	return [-1,-1]
		if (GFFlist[chromosome].features[feature].qualifiers.has_key("orf_classification")):
			if (GFFlist[chromosome].features[feature].qualifiers["orf_classification"][0]=="Dubious"):	return [-1,-1]

	start= GFFlist[chromosome].features[feature].location.start.position
	end= GFFlist[chromosome].features[feature].location.end.position
	
	if goodgenes == 1:
		feats=neighbors(GFFlist,chromosome,feature)
		prevfeat=feats[0]
		nextfeat=feats[1]
		prevend=GFFlist[chromosome].features[prevfeat].location.end.position
		nextstart=GFFlist[chromosome].features[nextfeat].location.start.position
		prevdir=GFFlist[chromosome].features[prevfeat].strand==GFFlist[chromosome].features[feature].strand
		nextdir=GFFlist[chromosome].features[nextfeat].strand==GFFlist[chromosome].features[feature].strand

		prevfeatid= GFFlist[chromosome].features[prevfeat].id
		nextfeatid= GFFlist[chromosome].features[nextfeat].id

		if utrtable5.has_key(prevfeatid):	prev5=utrtable5[prevfeatid][1]-utrtable5[prevfeatid][0]
		else:	prev5=0

		if utrtable3.has_key(prevfeatid):	prev3=utrtable3[prevfeatid][1]-utrtable3[prevfeatid][0]
		else:	prev3=0

		if utrtable5.has_key(nextfeatid):	next5=utrtable5[nextfeatid][1]-utrtable5[nextfeatid][0]
		else:	next5=0

		if utrtable3.has_key(nextfeatid):	next3=utrtable3[nextfeatid][1]-utrtable3[nextfeatid][0]
		else:	next3=0

		if strand==0:
			if prevdir==True:	prevend+=prev3
			if nextdir==True:	nextstart-=next5
		else:
			if prevdir==True:	prevend+=prev5
			if nextdir==True:	nextstart-=next3
			
		if strand==0:
			if (start-bp[0])<prevend and prevdir==True:	return[-2,-2]
			if (end+bp[1])>nextstart and nextdir==True:	return[-2,-2]
		else:
			if (start-bp[1])<prevend and prevdir==True:	return[-2,-2]
			if (end+bp[0])>nextstart and nextdir==True:	return[-2,-2]
	
	splicedseq= Seq('')
	splicedcounts= []
	for item in GFFlist[chromosome].features[feature].sub_features:
		if item.type == 'CDS': 
			start_feat = int(item.location.start.position)
			end_feat = int(item.location.end.position)
			splicedseq+=(GFFlist[chromosome][start_feat:end_feat]).seq
			splicedcounts+=counts[strand][chromosome][start_feat:end_feat]
			
	if splicedcounts==[] or str(splicedseq)=='':	return [-1,-1]
	
	if strand==1:	bp=[bp[1],bp[0],bp[2]]
		
	splicedcounts=counts[strand][chromosome][start-bp[0]:start]+splicedcounts
	splicedseq=GFFlist[chromosome][start-bp[0]:start].seq+splicedseq
	splicedcounts+=counts[strand][chromosome][end:end+bp[1]]
	splicedseq+=GFFlist[chromosome][end:end+bp[1]].seq

	if bp[2]!=0:
		if strand==0:
			if bp[2]>0:	splicedcounts=counts[strand][chromosome][start-bp[0]-bp[2]:start-bp[0]]+splicedcounts[:-bp[2]]
			else:	splicedcounts=splicedcounts[-bp[2]:]+counts[strand][chromosome][end+bp[1]:end+bp[1]-bp[2]]
		else:
			if bp[2]>0:	splicedcounts=splicedcounts[bp[2]:]+counts[strand][chromosome][end+bp[1]:end+bp[1]+bp[2]]
			else:	splicedcounts=counts[strand][chromosome][start-bp[0]+bp[2]:start-bp[0]]+splicedcounts[:bp[2]]

	if strand==1:
		splicedcounts.reverse()
		splicedseq=splicedseq.reverse_complement()
	
	if len(splicedcounts)!=len(splicedseq):
		return [-1,-1]
	return [splicedcounts,splicedseq]


def neighbors(GFFlist,chromosome,feature):
	if feature<0 or feature >=len(GFFlist[chromosome].features):
		print "Illegal feature given to neigbors program."
		return[-1,-1]
	
	i=1
	checkvar1=True
	checkvar2=True
	while(checkvar1==True or checkvar2==True):
		if((feature-i)<=0 and checkvar1==True):
			checkvar1=False
			prevfeat=0
		if((feature+i)>=len(GFFlist[chromosome].features) and checkvar2==True):
			checkvar2=False
			nextfeat=0
		if(checkvar1==True):
			if GFFlist[chromosome].features[feature-i].qualifiers.has_key("orf_classification"):
				if (GFFlist[chromosome].features[feature-i].qualifiers["orf_classification"][0]!='Dubious'):	dub=False
				else:	dub=True
			else:	dub=False
			if (GFFlist[chromosome].features[feature-i].type=='gene' and dub==False):
				prevfeat=feature-i
				checkvar1=False
		
		if(checkvar2==True):
			if GFFlist[chromosome].features[feature+i].qualifiers.has_key("orf_classification"):
				if (GFFlist[chromosome].features[feature+i].qualifiers["orf_classification"][0]!='Dubious'):	dub=False
				else:	dub=True
			else:	dub=False			
			if (GFFlist[chromosome].features[feature+i].type=='gene' and dub==False):
				nextfeat=feature+i
				checkvar2=False
		i+=1
		
	return [prevfeat,nextfeat]	   


def chrpostomrnapos(chrpos,chrom,featurenum,GFFlist):
	sublist=[]
	for subfeature in GFFlist[chrom].features[featurenum].sub_features:
		if subfeature.type== 'CDS':
			start=subfeature.location.start.position
			end=subfeature.location.end.position
			sublist.append([start,end])
	prevlength=0
	
	if(GFFlist[chrom].features[featurenum].strand==-1):	sublist.reverse()
	if len(sublist)==0:	return -1
	for item in sublist:
		start=item[0]
		end=item[1]
		length=end-start
		
		if(GFFlist[chrom].features[featurenum].strand==1):
			if chrpos >= start and chrpos < end:		
				mrnapos=prevlength+chrpos-start
				return mrnapos
			else:	prevlength+=length
		elif(GFFlist[chrom].features[featurenum].strand==-1):
			if chrpos < end and chrpos >= start:
				mrnapos=prevlength+end-chrpos-1
				return mrnapos
			else:	prevlength+=length
		else:
			print "Should not see this."

	if(GFFlist[chrom].features[featurenum].strand==1):	mrnapos=prevlength+chrpos-start-length	
	else:	mrnapos=prevlength+end-chrpos-1-length	
	return mrnapos




