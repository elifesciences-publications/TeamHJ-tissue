#script to create an nxn grid of square cells with two wall compartments separating each cell.

import random
n=2  # To create an n by n grid



global nCells, nWALLS, nVERTEX
c=0
nCells= 0
nWALLS=0
nVERTEX=0
chan = open("gridN.init","w")
nCells= n**2
nWALLS=4+6*(n-1)+2*(n-1)**2
nVERTEX=(n+1)**2
ncount=n-1 
vertexone=0
vertextwo=1
ncount3=3*n
chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   #prints first line of init file
while c<=ncount: #first row of vertical walls
   chan.write(" "+str(c)+" "+str(c)+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
   c=c+1
   vertexone=vertexone+1
   vertextwo=vertextwo+1


cellone=-1
celltwo=0
vertexone=n+1
vertextwo=celltwo
wallcounter=1
celloneA=0
celltwoA=celltwo
vertexoneA=n+1
vertextwoA=n+2


while wallcounter<=n-1:
 ncount=ncount+n+1
 chan.write(" "+str(c)+" "+str(celltwo)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
 c=c+1
 vertexone=vertexone+1
 vertextwo=vertextwo+1
 celltwo=celltwo+1
 cellone=cellone+1
 while c<=ncount-1: #horizontal walls
   chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
   c=c+1
   vertexone=vertexone+1
   vertextwo=vertextwo+1
   celltwo=celltwo+1
   cellone=cellone+1
 chan.write(" "+str(c)+" "+str(cellone)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
 vertexone=vertexone+1
 vertextwo=vertextwo+1
 ncount=ncount+n   
 c=c+1
 celltwoA=celltwo

 while c<=ncount: #internal vertical walls
  chan.write(" "+str(c)+" "+str(celloneA)+" "+str(celltwoA)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
  c=c+1
  vertextwoA+=1
  vertexoneA+=1
  celltwoA+=1
  celloneA+=1
 vertextwoA+=1
 vertexoneA+=1
 wallcounter+=1

ncount=ncount+n+1
chan.write(" "+str(c)+" "+str(celltwo)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
c=c+1
vertexone=vertexone+1
vertextwo=vertextwo+1
celltwo=celltwo+1
cellone=cellone+1
while c<=ncount-1: #horizontal walls
   chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
   c=c+1
   vertexone=vertexone+1
   vertextwo=vertextwo+1
   celltwo=celltwo+1
   cellone=cellone+1
chan.write(" "+str(c)+" "+str(cellone)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
vertexone=vertexone+1
vertextwo=vertextwo+1
ncount=ncount+n   
c=c+1

while c<=nWALLS-1:
  chan.write(" "+str(c)+" "+str(celloneA)+" "+str(-1)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
  c=c+1
  vertextwoA=vertextwoA+1
  vertexoneA=vertexoneA+1
  celltwoA=celltwoA+1
  celloneA=celloneA+1

chan.write(" \n")
chan.write(" "+str(nVERTEX)+" 2 \n")
ycoord=0  
xcoord=0
while xcoord<=n:
  while ycoord<=n:
    chan.write(""+str(xcoord)+" "+str(ycoord)+"\n")
    ycoord+=1
  ycoord=0
  xcoord+=1

chan.write("\n"+str(nWALLS)+" 1  6 \n ")
wallcounter=1

while wallcounter<=nWALLS:
  chan.write("1 0 0 0.1 0.1 0.1 0.1\n")  # Prints wall variabales
  wallcounter+=1
chan.write("\n"+str(nCells)+" 9\n ")
wallcounter=1

while wallcounter<=nCells:
  auxin=random.uniform(1, 1.2)
  pin=auxin+random.uniform(0, 0.1)
  chan.write("1 0 1 0.9 "+str(auxin)+ " "+str(pin)+ " 1 1 1\n") # Prints cell variables
  wallcounter+=1
chan.close()



