import sys
import json
sys.path.append('cadnano2')
import cadnano
import main
from model.document import Document
from model.io.decoder import decode
from model.io.encoder import encode
import util
import copy
from model.strand import Strand
import numpy as np

class analysisHelpers(object):
    def __init__(self, fileName):
        self._fN=fileName
        self._doc=Document()
        self._numStaples=0
        self._scaffoldLength=0
        decode(self._doc,file(self._fN).read())
        self.getModelInfo()
        self.printModelInfo()
        self.vhLoop()
        self.numberScaffold()
	
    ##Initialization Functions
    def vhLoop(self):
    	for vh in self._vhList:
    	    vh._renumber=-1
    	    vh._neighbors=self._part.getVirtualHelixNeighbors(vh)
    
    def getModelInfo(self):
		self._part=self._doc.selectedPart()
		self._vhList=self._part.getVirtualHelices()
		self._oligos=list(self._part.oligos())
	
    def printModelInfo(self):
		self._numVH=len(self._vhList)
		self._numOligos=len(self._oligos)
		print "Design with %s virtual helices and %s oligos loaded" %(self._numVH,self._numOligos)
	
    #Sequence Application Functions	
    def readSeqFile(self,sequenceFile):
		sqFile=file(sequenceFile)
		self._fullSeq=sqFile.read().strip()
		self._fullSeqLength=len(self._fullSeq)
		print "Sequence length %s" %(self._fullSeqLength)
	
    def applySequenceToScaffold(self, sequenceFile):
		self.readSeqFile(sequenceFile)
		#Find scaffold
		self.findScaffold()
		self._scaffold.applySequence(self._fullSeq)
		print "Applied sequence with %s length\n" %(len(self._scaffold.sequence()))
		self._scaffoldLength=len(self._scaffold.sequence())
        
    def findScaffold(self):
		for x in self._oligos:
			scaf=x.isStaple()
			if scaf==False:
				self._scaffold=x
                #print "scaffold found"
    
    def getCircleMapData(self):
        self.numberScaffold()
        self.numberStaples()
        self.printCircleData()

    
    def numberScaffold(self):
        #Get start info
        self.findScaffold()
        count=1
        first=0
        for strand in self._scaffold.strand5p().generator3pStrand():
            sLowIdx, sHighIdx = strand._baseIdxLow, strand._baseIdxHigh
            subn=strand.insertionLengthBetweenIdxs(sLowIdx, sHighIdx)
            difference=strand.length()+subn
            vh=strand.virtualHelix()
            if vh._renumber==-1:
                #print "Renumbering helix %s as %s\n" %(vh.number(), first)
            	vh._renumber=first
            	first=first+1
            #print "length is %s" %difference
            count+=(difference)
            idx5p=strand.idx5Prime()
            idx3p=strand.idx3Prime()
            #if strand.isDrawn5to3
            
            strand.num5Prime=count-difference
            strand.num3Prime=count-1
            if not self._scaffoldLength:
            	self._scaffoldLength=strand.num3Prime
        #scaffvh=self._scaffold._strand5p.virtualHelix()
        
    def numberStaples(self):
        for x in self._oligos:
            if x.isStaple()==True:            	
                self._numStaples+=1
                for strand in x.strand5p().generator3pStrand():
                    strand.num5Prime=[]
                    strand.num3Prime=[]
                    #print "Staple part\n"
                    scafComp=strand.getComplementStrands()
                    if len(scafComp)>0:
                    	for p in range(0,len(scafComp)):
                        	dir1=scafComp[p].isDrawn5to3()
                            	diff53p=np.abs(scafComp[p].idx5Prime()-strand.idx3Prime())
                                diff35p=np.abs(scafComp[p].idx3Prime()-strand.idx5Prime())
                                strand.num5Prime.append(scafComp[p].num3Prime-diff35p)
                            	strand.num3Prime.append(scafComp[p].num5Prime+diff53p)


    def printCircleData(self):
        fileString="CD_"
        fileString+=self._fN
        fileString+=".csv"
        self._circleFileName=fileString
        self._circleFile=open(fileString,'w')
        self._circleFile.write("Staple locations on scaffold 5'->3'\n")
        for x in self._oligos:
            if x.isStaple()==True:
                stapleString=""
                count=0
                for strand in x.strand5p().generator3pStrand():
                    if count>0:
                        stapleString+=","
                    scafComp=strand.getComplementStrands()
                    if len(scafComp)>0:
                   	n=len(strand.num5Prime)
                   	for p in range(0,n):
                   		if p!=0:
                   			stapleString+=","
                   		stapleString+=("%s, %s") %(strand.num5Prime[p],strand.num3Prime[p])

                    count+=1
                if stapleString:
                    stapleString+="\n"
                    #print stapleString
                    self._circleFile.write(stapleString)
        
        self._circleFile.close()
        
    def getVirtualHelixSequences(self):      
        for vh in self._vhList:
            seq=''
            first=1
            minLoc=0
            maxLoc=0
            for scaffStr in vh.scaffoldStrandSet():
            	scaffStr._minMax=[]
            	seq+=scaffStr.sequence()
            	if minLoc==0 and maxLoc==0:
            	   minLoc=scaffStr.idx5Prime()
            	   maxLoc=scaffStr.idx3Prime()
            	if scaffStr.isDrawn5to3():
            	  if scaffStr.idx5Prime()<minLoc:
            	     minLoc=scaffStr.idx5Prime()
            	  if scaffStr.idx3Prime()>maxLoc:
            	     maxLoc=scaffStr.idx3Prime()
            	else:
            	  if scaffStr.idx3Prime()<minLoc:
            	     minLoc=scaffStr.idx3Prime()
            	  if scaffStr.idx5Prime()>maxLoc:
            	     maxLoc=scaffStr.idx5Prime()     
                scaffStr._minMax.append([minLoc,maxLoc])
            vh._minMax=[minLoc,maxLoc]
            vh._sequence=seq
        
    
    def printScaffoldLength(self):
        return self._scaffoldLength
    
    def printCircleFileName(self):
        return self._circleFileName
    
    def setScaffoldLength(self,sfL):
    	self._scaffoldLength=sfL
   
    def getVirtualHelixList(self):
   	return self._vhList
   
    def getScaffoldStrand(self):
   	return self._scaffold
    


                        
                
            
