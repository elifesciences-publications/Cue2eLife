import sys, os, multiprocessing

exp,file= sys.argv[1], sys.argv[2]

mismatch= '01'

GFFfile= '/home/xxx/yyy/yeast_KD2.gff'

rootpath= '/home/xxx/Root/%s/%s' %(exp,file)

bamfile1= '%s/%s_starM%s/%s_match.bam' %(rootpath,file,mismatch,file)

bamfile2= '-1'

footprints= [[27,32,5,'{27:[16],28:[16],29:[17],30:[17],31:[17],32:[17]}'],[20,22,5,'{20:[16],21:[17],22:[17]}'],[15,17,3,'{15:[0],16:[0],17:[0]}']]

totalreads= '-1'

softclipped= '2'

def densebuilder(footprints):
    minlen,maxlen,assignment,riboshiftdict= footprints[0],footprints[1],footprints[2],str(footprints[3])
    bamfileoutput= '%s/%s_starM%s/%s_%sto%s_match.bam' %(rootpath,file,mismatch,file,minlen,maxlen)
    densityfilepath= '%s/density%s_%sto%sM%s/%s_%sto%sf' %(rootpath,assignment,minlen,maxlen,mismatch,file,minlen,maxlen)
    densitypathandfilestring= '%s/%s_%sto%sf' %(densityfilepath,file,minlen,maxlen)
    wigpath= '%s/wigfile%s_%sto%sM%s' %(rootpath,assignment,minlen,maxlen,mismatch)
    wigpathandfile= '%s/%s_%sto%sM%s' %(wigpath,file,minlen,maxlen,mismatch)

    if os.path.exists(bamfile1):
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        if not os.path.exists(wigpath):    os.makedirs(wigpath)

    commandstring= 'python yeastdensebuilder_main.py --GFFfile %s --bamfile1 %s --bamfile2 %s --densitypathandfilestring %s --wigpathandfile %s --totreads %s --assignment %s --riboshiftdict %s --bamfileoutput %s --softclipped %s' %(GFFfile,bamfile1,bamfile2,densitypathandfilestring,wigpathandfile,totalreads,assignment,riboshiftdict,bamfileoutput,softclipped)

    print commandstring
    os.system(commandstring)

pool= multiprocessing.Pool(len(footprints))
Dense = pool.map(densebuilder, footprints)
pool.close()
pool.join()

