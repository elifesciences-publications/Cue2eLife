import writegene_main, sys, multiprocessing

exp,file= sys.argv[1],sys.argv[2]

gfffile= '/home/xxx/yyy/yeast_KD2.gff'

utrgfffile= '-1'

footprints= [['27to32',5,'A'],['20to22',5,'A'],['15to17',3,'']]

shift5= 0    # Positive value gives extra on end.

shift3= 50    # Positive value puts extra on end.

gene= 'GFP'

def outputgene(footprints):
    sizerange,assignment,posshift= footprints[0],footprints[1],footprints[2]
    densityfile= '/home/xxx/Root/%s/%s/density%s_%s%sM01/%s_%sf%s/%s_%sf%s' %(exp,file,assignment,sizerange,posshift,file,sizerange,posshift,file,sizerange,posshift)
    outfilepath= '/home/xxx/data/genemodel/%s_%s%s' %(file,sizerange,posshift)
    print densityfile
    print outfilepath
    commandstring= writegene_main.writegene_wf2(shift5,shift3,0,densityfile,gene,gfffile,utrgfffile,outfilepath)

pool= multiprocessing.Pool(len(footprints))
Genewriter = pool.map(outputgene, footprints)
pool.close()
pool.join()

