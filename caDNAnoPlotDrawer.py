##caDNAno Plot Drawer##
# This file contains definition of two classes.
# 1. caDNAnoCirclePlotDrawer draws circle plots
# 2. caDNAnoBasePlotDrawer

from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
import plotly.offline as offline
from plotly.graph_objs import *
from analysisHelpers import analysisHelpers
from scipy import signal

import pylab as pl
import numpy as np
import math
import csv
import matplotlib.cm as cm

class caDNAnoCirclePlotDrawer(object):
    def __init__(self, arg, *kwarg):
        if isinstance(arg,analysisHelpers):
            self._fN=arg.printCircleFileName()
            self._scaffoldLength=arg.printScaffoldLength()

        else:
            print "is not instance"
            self._fN=arg
            self._scaffoldLength=kwarg[0]
            
        self._lengthList=[]
        self._rad=self._scaffoldLength/(np.pi*2)
        self.readFileCSV()
        
        self._layout={
            'xaxis': {
            'range': [-self._rad-100, self._rad+100],
            'zeroline': False,
            'showgrid': False,
        },
        'yaxis': {
            'range': [-self._rad-100, self._rad+100],
            'showgrid': False,
            'zeroline': False,
        },
        'width': 600,
        'height': 600,
        'hovermode' : 'closest',
        'showlegend' : False,
        'shapes': [
            {
            'type': 'circle',
            'xref': 'x',
            'yref': 'y',
            'x0': -self._rad,
            'y0': -self._rad,
            'x1': self._rad,
            'y1': self._rad,
            'line': {
                'color': 'rgba(0,0,0,0.5)',
            }
            }
            ]
    
        }
        
        self.makeFigData()
        self._fig = dict( data=self._data , layout=self._layout)
        
    
    ##Functions for Bezier Lines
    def bezierPoint(self,b,t):
        Btx=(1-t)*((1-t)*b[0]+t*b[2])+t*((1-t)*b[2]+t*b[4])
        Bty=(1-t)*((1-t)*b[1]+t*b[3])+t*((1-t)*b[3]+t*b[5])
        return [Btx,Bty]         

    def computeBezier(self,pointsx, pointsy, Btype):
        midPointX=np.mean(pointsx)
        midPointY=np.mean(pointsy)
        distanceXY=((pointsx[0]-pointsx[1])**2+(pointsy[0]-pointsy[1])**2)**0.5
        #print distanceXY
        ref=0
        if distanceXY<(self._rad/2):
            ref=self._rad-distanceXY
        else:
            ref=0
        if (Btype==1):
            ref=self._rad

        angleAT=np.arctan2(midPointY,midPointX)
        refX=ref*np.cos(angleAT)
        refY=ref*np.sin(angleAT)
        #print "%s %s" %(refX,refY)
        if (True):
            ##Points for Bezier Curve - If making arc
            x0=np.cos(np.pi/2)
            y0=np.sin(np.pi/2)
            x1=(4-x0)/3
            y1=((1-x0)*(3-x0))/(3*y0)
            angleAT=np.arctan2(midPointY,midPointX)
            ref1x=ref*np.cos(angleAT)
            ref1y=ref*np.sin(angleAT)
            x2=x1
            y2=-y1
            angleAT=np.arctan2(midPointY,midPointX)
            ref2x=self._rad*np.cos(angleAT)
            ref2y=self._rad*np.sin(angleAT)
            x3=x1
            y3=-y0
            bezierPoints=[pointsx[0],pointsy[0],refX,refY,pointsx[1],pointsy[1]]
            return bezierPoints

    ##Read csv file
    def readFileCSV(self):
        f=open(self._fN)
        #circThreshold=(numThreshold*360.0)/self._scaffoldLength
        rowCount=sum(1 for row in csv.reader(f))
        self._numStaples=rowCount-1
        secsPerStaple=0
        gradient=np.linspace(0.2, 1, 2*self._rad)
        multfactor=5
        
        
        if self._numStaples<5:
            multfactor=10

		#This creates arrays and then sets all values in the array to 2*rad.
        self._xPlot=np.zeros((self._numStaples+1,self._numStaples*multfactor))
        self._xPlot=[x+2*self._rad for x in self._xPlot]
        
        
        self._yPlot=np.zeros((self._numStaples+1,self._numStaples*multfactor))
        self._yPlot=[x+2*self._rad for x in self._yPlot]
        
        self._circlePoints=np.zeros((self._numStaples+1, self._numStaples*multfactor))
        self._circlePoints=[x+2*self._rad+1 for x in self._circlePoints]
        
        #colorList=np.zeros((self._numStaples+1,self._numStaples*multfactor))
        self._labels=np.zeros((self._numStaples+1,self._numStaples*multfactor))
        self._labels2=np.zeros((self._numStaples+1,self._numStaples*multfactor))

        count=0
        f=open(self._fN)

        #Read data from file
        for row in csv.reader(f):
            count=count+1
            secsPerStaple=len(row)/2	#Number of sections to staple.
            if (count>1):
                x1=[int(row[x]) for x in range(0,len(row),2)]
                #x1=[i-self._scaffoldLength if (i>self._scaffoldLength) else i for i in x1]
                x1theta=[x*360.0/self._scaffoldLength for x in x1]
                x2=[int(row[x]) for x in range(1,len(row),2)]
                #x2=[i-self._scaffoldLength if (i>self._scaffoldLength) else i for i in x2]
                x2theta=[x*360.0/self._scaffoldLength for x in x2]
                for i in range(0,len(x1)+1):
                    
                    if (i==0):
                        keepNext=x2theta[i]
                        keepNextLabel=x2[i]
                        self._xPlot[count-2][0]=np.cos(math.radians(x1theta[i]))*(self._rad)		
                        self._yPlot[count-2][0]=np.sin(math.radians(x1theta[i]))*(self._rad)
                        self._circlePoints[count-2][0]=x1[i]
                        self._labels[count-2][0]=x1[i]
                        c=0
                    
                    else:
                        c=c+2
                        self._xPlot[count-2][c-1]=np.cos(math.radians(keepNext))*(self._rad)		#Uses previous value for first point
                        self._yPlot[count-2][c-1]=np.sin(math.radians(keepNext))*(self._rad)
                        #colorList[count-2][c-1]=np.ceil(keepNext/circThreshold)
                        self._labels[count-2][c-1]=(keepNextLabel)
                        self._circlePoints[count-2][c-1]=keepNextLabel
                        
                        if (i!=secsPerStaple):
                            self._xPlot[count-2][c]=np.cos(math.radians(x1theta[i]))*(self._rad)
                            self._yPlot[count-2][c]=np.sin(math.radians(x1theta[i]))*(self._rad)
                            #colorList[count-2][c]=np.ceil(x1theta[i]/circThreshold)
                            self._circlePoints[count-2][c]=x1[i]
                            self._labels[count-2][c]=x1[i]
                            keepNext=x2theta[i]
                            keepNextLabel=x2[i]
						
						
                            

    def makeFigData(self):
        self._layout=dict(self._layout)
        start=0
        xP=[]
        yP=[]
        self._data=[]
        self._histData=[]
        self._connectionLength=[]
        pathString=""

        for i in range(0,self._numStaples):
          labelTmp=np.where(self._labels[i][:]!=0)
          labelT=self._labels[i][labelTmp]
          colorsP=[]
          tmpX=np.where(self._xPlot[i][:]!=2*self._rad)    #x-coordinates for staple points
          tmpY=np.where(self._yPlot[i][:]!=2*self._rad)    #y-coordinates for staple points
          xP=(self._xPlot[i][tmpX])
          yP=(self._yPlot[i][tmpY])
          
          for n in range(0,len(labelT)):
            #color1=np.ceil(float(labelT[n]/numThreshold))
            #colorsP.append(colorArray[6])   ##Get colors for scatter points
            if (n>0):
                BezType=n%2 #If one, draw Bezier Curve along scaffold. If 0, draw along inside
                passX=[xP[n],xP[n-1]]
                passY=[yP[n],yP[n-1]]
                labelXY=[labelT[n],labelT[n-1]]
                connectL=np.abs(labelT[n]-labelT[n-1])
                
                if (connectL>self._scaffoldLength/2):
                    connectL=connectL-self._scaffoldLength/2

                if BezType==0:
                    self._connectionLength.append(connectL)
                
                passString=self.computeBezier(passX,passY,BezType)

                pathString = 'M %s,%s Q %s,%s %s,%s' %(passString[0],passString[1],passString[2],passString[3],passString[4],passString[5])
                pathString=pathString.replace('"','')
                pathString=pathString.replace("'","")
                shapeString = { 'type': 'path', 'path' : pathString, 'line': {'color': 'rgb(0, 0, 0)'}}
                #print layout['shapes'].append
                self._layout['shapes'].append(shapeString)
                bPdraw=self.bezierPoint(passString,0.5)


                trace=Scatter(
                    x=[bPdraw[0]],
                    y=[bPdraw[1]],
                    mode = 'markers',
                    marker = dict(
                    size = 1,
                    color = 'rgba(0, 0, 0, .8)',
                    ),
                    text=str(labelXY[0])+' to '+str(labelXY[1]),
                    hoverinfo='text',
                )
                self._data.append(trace)


        #All points for one staple are stored in a single trace object
          trace1=Scatter(
            x=xP[:].tolist(),
            y=yP[:].tolist(),
            mode='markers',
            marker=Marker( size=10, color=colorsP ),
            text=[str(labelT[p]) for p in range(0,len(labelT))],
            hoverinfo='text'
            )
          self._data.append(trace1)
        
        trace2=Histogram(
            x=self._connectionLength,
            histnorm='count',
            autobinx=False,
            xbins=dict(
               start=0,
               end=self._scaffoldLength/2,
              size=2
            ),
            marker=dict(
            color='black',
            ),
            opacity=0.75
		    )
        self._histData.append(trace2)
        
    def writeHTMLfile(self):
		fstring="circle.html"
		f=open(fstring,'w')
		f.write('<!DOCTYPE html>\n<html>\n<head>\n<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>\n</head>\n\n<body>\n<div id="myDiv" style="width:800px;height:800px;">\n</div>\n<script>\n')
		f.write('var data = ')
		fstringData=str(self._data)
		fstringData=fstringData.replace('False','false')
		fstringData=fstringData.replace('True','true')
		fstringLayout=str(self._layout)
		fstringLayout=fstringLayout.replace('False','false')
		fstringLayout=fstringLayout.replace('True','true')
		f.write(fstringData)
		f.write('\n\n')
		f.write('var layout= ')
		f.write(fstringLayout)
		f.write('\n\nPlotly.newPlot("myDiv",data,layout);')
		f.write('</script>\n</body>\n')
		f.close()
		
		fstring="circleNanohub.html"
		f=open(fstring,'w')
		f.write('<!DOCTYPE html>\n<html>\n<head>\n<script src="cadnanovis/src/plotly-latest.min.js"></script>\n</head>\n\n<body>\n<div id="myDiv" style="width:500px;height:500px;">\n</div>\n<script>\n')
		f.write('var data = ')
		fstringData=str(self._data)
		fstringData=fstringData.replace('False','false')
		fstringData=fstringData.replace('True','true')
		fstringLayout=str(self._layout)
		fstringLayout=fstringLayout.replace('False','false')
		fstringLayout=fstringLayout.replace('True','true')
		f.write(fstringData)
		f.write('\n\n')
		f.write('var layout= ')
		f.write(fstringLayout)
		f.write('\n\nPlotly.newPlot("myDiv",data,layout);')
		f.write('</script>\n</body>\n')
		f.close()
		
		fstring="histogramNanohub.html"
		f=open(fstring,'w')
		f.write('<!DOCTYPE html>\n<html>\n<head>\n<script src="cadnanovis/src/plotly-latest.min.js"></script>\n</head>\n\n<body>\n<div id="myDiv" style="width:500px;height:500px;">\n</div>\n<script>\n')
		f.write('var data = ')
		fstringData=str(self._histData)
		fstringData=fstringData.replace('False','false')
		fstringData=fstringData.replace('True','true')
		f.write(fstringData)
		f.write('\n\n')
		f.write('\n\nPlotly.newPlot("myDiv",data);')
		f.write('</script>\n</body>\n')
		
		fstring="histogram.html"
		f=open(fstring,'w')
		f.write('<!DOCTYPE html>\n<html>\n<head>\n<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>\n</head>\n\n<body>\n<div id="myDiv" style="width:800px;height:800px;">\n</div>\n<script>\n')
		f.write('var data = ')
		fstringData=str(self._histData)
		fstringData=fstringData.replace('False','false')
		fstringData=fstringData.replace('True','true')
		f.write(fstringData)
		f.write('\n\n')
		f.write('\n\nPlotly.newPlot("myDiv",data);')
		f.write('</script>\n</body>\n')
		f.close()
          
class caDNAnoBaseMapDrawer(object):
    def __init__(self, arg,arg2):
        if isinstance(arg,analysisHelpers):
            arg.getVirtualHelixSequences()
            self._part=arg._part
            self._numVH=arg._numVH
            self._smooth=arg2
            self._vhList=arg.getVirtualHelixList()
            self._scaffold=arg.getScaffoldStrand()
            self._colorscale = [[0, '#0505FF'],[0.5 ,'#DCDCDC'], [1, '#FF0505']]
            self._vhNums=[]
            self._x=[]
            self._y=[]
            self._z=[]
            self._scafX=[]
            self._scafY=[]
            self.numberVH()            

            
        else:
            print "input is not not instance of analysisHelpers Class"
            #Write function here to read data from file.
    
    def makeScaffoldMap(self):
		self.drawScaffoldRoute()
		self.makeScaffoldFigData()
    
    def drawScaffoldRoute(self):
		lookupY=self._vhLookUp
		for strand in self._scaffold.strand5p().generator3pStrand():
			vh=strand.virtualHelix()
			coord=vh.coord()
			vhNum=vh._number
			ycoord=lookupY[vhNum]
			stapleStrands=vh.stapleStrandSet()
			low, high=strand.idxs()
			fiveToThree=strand.isDrawn5to3()

			if fiveToThree==0:
				self._scafX.append(high)
				self._scafY.append(ycoord)
				self._scafX.append(low)
				self._scafY.append(ycoord)
			else:
				#reverse order
				self._scafX.append(low)
				self._scafY.append(ycoord)
				self._scafX.append(high)
				self._scafY.append(ycoord)			
		
    def makeBaseMap(self):
		self.fillDataStructure2()
		self.makeFigData()
		return 0
	
	
		
    def initDataStructure(self):
    	self._map=[None]*(self._vhRange+1)
    	for i in range(0,self._vhRange+1):
    	    self._map[i]=[None]*(self._range+1)
    	
    	print "Created list of dimension %s by %s" %(self._range+1,self._vhRange+1)
    
    
    def numberVH(self):
        first=1
        ycoord1=0
        vhList=[]
        yList=[]
        lookupY=[None]*(self._numVH)
        coord=None
        
    	for strand in self._scaffold.strand5p().generator3pStrand():
    		vh=strand.virtualHelix()
    		vhNum=vh._number
    		coord=vh.coord()
    		stapleStrands=vh.stapleStrandSet()
    		if first!=1:
				diff=[coord[0]-coordP[0],coord[1]-coordP[1]]
				if abs(sum(diff))>1:
					ycoord1+=1

    		if vhNum not in vhList:
				vhList.append(vhNum)
				yList.append(ycoord1)
				lookupY[vhNum]=ycoord1
				ycoord1+=1
    		
    		coordP=coord
    		first=0
		self._vhLookUp=lookupY
			
    
    def fillDataStructure2(self):
        lookupY=self._vhLookUp
    	for strand in self._scaffold.strand5p().generator3pStrand():
    		vh=strand.virtualHelix()
    		coord=vh.coord()
    		vhNum=vh._number

    		ycoord=lookupY[vhNum]
    		stapleStrands=vh.stapleStrandSet()
    		low, high=strand.idxs()
    		
    		if strand.sequence():
    		   #print "Strand has sequence\n"
    		   sequence=strand.sequence()
    		   compS=strand.getComplementStrands()
    		   numericSequence=self.convertSeqtoNum(sequence)
    		   fiveToThree=strand.isDrawn5to3()
    		   #print "Start set to %s; End set to %s\n" %(start,end)

    		if fiveToThree==0:
    		   numericSequence=numericSequence[::-1]
    		bcount=0
    		for x in range(low,high+1):
				subn=strand.insertionLengthBetweenIdxs(x, x)
				if (stapleStrands.hasStrandAt(x,x) and subn!=-1):
					self._z.append(numericSequence[bcount])
					self._y.append(ycoord)
					self._x.append(x)
					bcount+=1
				elif subn==-1:
					self._z.append('nan')
					self._y.append(ycoord)
					self._x.append(x)
				else:
					self._z.append(None)
					self._y.append(ycoord)
					self._x.append(x)
		
   		
    	xmax=max(self._x)
    	xmin=min(self._x)
    	ymax=max(self._y)
    	self._map=self.initDataStructure2(xmax,ymax)
    	for n in range(len(self._x)):
			self._map[self._y[n]][self._x[n]]=self._z[n]

    	for y in range(len(self._map)):
			noneSt=-1
			addEnd=0
			for x in range(len(self._map[y])):
				if self._map[y][x]!=None and noneSt==-1:
					noneSt=x
				if self._map[y][x]=='nan':
					del(self._map[y][x])
					if addEnd%2==1:
						self._map[y].append(None)
					else:
						self._map[y].insert(0,None)
					addEnd+=1

    	for p in range(len(self._map)):
			mRow=signal.savgol_filter(self._map[p],self._smooth,2).tolist()
			#Sgolay filter will have values greater than 1 or less than 0. This shouldn't
			#be possible, so fix it here.
			mRow=[x if np.isnan(x)==True or (np.isnan(x)==False and x<1) else 1 for x in mRow]
			mRow=[x if np.isnan(x)==True or (np.isnan(x)==False and x>0) else 0 for x in mRow]
			self._map[p]=mRow
				
	
    def initDataStructure2(self,x,y):	
		xyMap = [ [ None for i in range(x+1) ] for j in range(y+1) ]
		return xyMap
    		       
    
    def convertSeqtoNum(self,sequence):
    	seqList=[]
        #print sequence
    	for base in sequence:
     	   if (base is 'A') or (base is 'T'):
     	      baseNum=0
     	   elif (base is 'G') or (base is 'C'):
     	      baseNum=1
     	   else:
     	      baseNum=None
     	   seqList.append(baseNum)
     	return seqList

              
    def makeScaffoldFigData(self):
		self._scafData=[]
		
		trace1=Scatter(
			x=self._scafX,
			y=self._scafY,
			name="scaffold",
			mode='lines',
			marker = dict(
                color = 'rgba(0, 0, 0, .9)',
              ),
		
		)
		
		trace2=Scatter(
			x=[self._scafX[0]],
			y=[self._scafY[0]],
			name="5' scaffold",
			mode = 'markers',
            marker = dict(
                size = 10,
                color ='rgba(0.6, 0, 0, .8)',
              ),
            
		)
		
		trace3=Scatter(
			x=[self._scafX[-1]],
			y=[self._scafY[-1]],
			name="3' scaffold",
			mode = 'markers',
            marker = dict(
                size = 10,
                color = 'rgba(255, 182, 193, .9)',
              ),
            
		
		)
		
		self._scafData.append(trace1)
		self._scafData.append(trace2)
		self._scafData.append(trace3)
		self._figScaf=dict(data=self._scafData)
		
    
    def makeFigData(self):
       self._data=[{
          'z': self._map,
          'type': 'heatmap',
          'showscale' : True,
          'connectgaps' : False,
          'colorscale': self._colorscale,
          'colorbar' : {
          	'title':'GC Content',
          	'titleside':'right',
          	'titlefont' : {
          		'size': 16,
          	}
          
          },
        },
        ]
        
             
       self._layout={
    	 	"autosize": True,
    	 	"margin" : {
    	 	   't':0,
    	 	   'b':0,
    	 	   'l':0,
    	 	   'r':0,
    	 	
    	 	} ,
    	  
    	} 
    	
       self._meshData=Mesh3d(x=self._x,y=self._y,z=self._z)
    	
       self._fig2 = [self._meshData]
       self._fig = dict( data=self._data , layout=self._layout)
       

    def writeHTMLfile(self):
		fstringData=str(self._data)
		fstringData=fstringData.replace('False','false')
		fstringData=fstringData.replace('True','true')
		fstringData=fstringData.replace('None','null')
		fstringData=fstringData.replace('nan','null')
		
		fstringLayout=str(self._layout)
		fstringLayout=fstringLayout.replace('False','false')
		fstringLayout=fstringLayout.replace('True','true')
		
		fstring="baseMap.html"
		f=open(fstring,'w')
		f.write('<!DOCTYPE html>\n<html>\n<head>\n<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>\n</head>\n\n<body>\n<div id="myDiv" style="width:800px;height:800px;">\n</div>\n<script>\n')
		f.write('var data = ')
		f.write(fstringData)
		f.write('\n\n')
		f.write('var layout= ')
		f.write(fstringLayout)
		f.write('\n\nPlotly.newPlot("myDiv",data,layout);')
		f.write('</script>\n</body>\n')
		f.close()
		
		fstring="baseMapNanohub.html"
		f=open(fstring,'w')
		f.write('<!DOCTYPE html>\n<html>\n<head>\n<script src="cadnanovis/src/plotly-latest.min.js"></script>\n</head>\n\n<body>\n<div id="myDiv" style="width:500px;height:500px;">\n</div>\n<script>\n')
		f.write('var data = ')
		f.write(fstringData)
		f.write('\n\n')
		f.write('var layout= ')
		f.write(fstringLayout)
		f.write('\n\nPlotly.newPlot("myDiv",data,layout);')
		f.write('</script>\n</body>\n')
		f.close()

            
       
        
        
            
        
    
