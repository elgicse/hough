#!/usr/bin/env python
import ROOT as r
import numpy as np
import sys

inrootfile = "data/UTHits.root"
outrootfile = "data/out/out.root"

TREENAME = "UTHits/MatchedTracks"
XY,XZ,YZ = 0,1,2
X,Y,Z = 0,1,2

class mydict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

def rootPreProcessing():
    r.gROOT.ProcessLine(".x /home/elena/Documenti/University/PhD/UZH/LHCb/lhcbstyle.C")

def mgInit():
    return r.TMultiGraph()

def matrixInit():
    return np.matrix()

def isGood(track):
    return track.reconstructible_asLong and track.true_p > 3000 and track.true_pt > 500 and track.true_eta > 2 and track.true_eta < 5 and abs(track.pid) != 22 and track.fromB

def HGraph(n,x,y):
    gr = r.TGraph(n,x,y)
    gr.SetMarkerStyle(7)
    #gr.SetLineColor()
    #gr.SetMarkerColor()
    gr.SetName("tracks_graph") 
    gr.SetTitle("tracks_graph") 
    return gROOT    

def hitGraph(track,plane):
    if plane == XY:
        gr = HGraph(track.vplen,track.velo_x_hit,track.velo_y_hit)
    if plane == XZ:
        gr = HGraph(track.vplen,track.velo_z_hit,track.velo_x_hit)   
    if plane == YZ:
        gr = HGraph(track.vplen,track.velo_z_hit,track.velo_y_hit)   
    return gr


def createHitsGraph(mg,tracks):
    for t in filter(isGood,tracks):
            mg[XY].Add(hitGraph(t,XY))
            mg[XZ].Add(hitGraph(t,XZ))
            mg[YZ].Add(hitGraph(t,YZ))
    return True

def getXYZ(hit,track):
    index = track.index(hit)
    x = track.velo_x_hit[index]
    y = track.velo_y_hit[index]
    z = track.velo_z_hit[index]
    #x,y,z = 0,0,0 #calc x,y,z
    hit = x,y,z
    return x,y,z

def calcRho(plane,th,x,y,z):
    if plane == XY:
        rho = transform(x,y,th)
    if plane == XZ:
        rho = transform(z,x,th)
    if plane == YZ:
        rho = transform(z,y,th)
    return rho

def transform(x1,x2,theta):
    return x1*np.cos(theta) + x2*np.sin(theta)

def setDictionaries(dictionaries,tracks):
    for t in tracks:
        if isGood(t):
            #"GET HIT LIST"
            hits = []*track.vplen
            for h in hits:
                x,y,z = getXYZ(h,track)
                rho = [],[],[]  
                for th in theta:
                    rho[XY] = calcRho(XY,th,x,y,z)
                    dictionaries[XY].setdefault((rho[XY],th),dict()).setdefault(X,list()).append(x)
                    dictionaries[XY].setdefault((rho[XY],th),dict()).setdefault(Y,list()).append(y)
                    dictionaries[XY].setdefault((rho[XY],th),dict()).setdefault(Z,list()).append(z)

                    rho[XZ] = calcRho(XZ,th,x,y,z)
                    dictionaries[XZ].setdefault((rho[XZ],th),dict()).setdefault(X,list()).append(x)
                    dictionaries[XZ].setdefault((rho[XZ],th),dict()).setdefault(Y,list()).append(y)
                    dictionaries[XZ].setdefault((rho[XZ],th),dict()).setdefault(Z,list()).append(z)

                    rho[YZ] = calcRho(YZ,th,x,y,z)
                    dictionaries[YZ].setdefault((rho[YZ],th),dict()).setdefault(X,list()).append(x)
                    dictionaries[YZ].setdefault((rho[YZ],th),dict()).setdefault(Y,list()).append(y)
                    dictionaries[YZ].setdefault((rho[YZ],th),dict()).setdefault(Z,list()).append(z)
    return True

def dict2Matrix(dictionaries,matrices):
    pass

def searchPeaks(matrices,peaklists):
    pass

def calcParams(peaklists,paramlists):
    pass

def createHoughGraph(mg,tracks):
    setDictionaries(Dictionaries,tracks)
    dict2Matrix(Dictionaries,Matrices)
    searchPeaks(Matrices,PeakLists)
    calcParams(PeakLists,ParamLists)
    "create graphs"
    "add graph in multigraf"
    "save in root file"
    "create DTP Dictionary"
    return True


if __name__ == "__main__":
    print "HOUGH START..."       

    rootPreProcessing()

        #init data structures
    MultiGraphs = mgInit(),mgInit(),mgInit()
    PeakLists = [],[],[]
    ParamLists = [],[],[]
    Dictionaries = mydict(),mydict(),mydict()
    #Matrices = matrixInit(),matrixInit(),matrixInit()
    theta = range(180)




        #root file descriptor
    ifd = r.TFile(inrootfile) #input file descriptor 
    ofd = r.TFile(outrootfile,"recreate") #output file descriptor 
        #"SET WRITE FLAG AND CREATE OUT TREE"

        #get tree
    Tracks = ifd.Get(TREENAME)

    if not createHitsGraph(MultiGraphs,Tracks):
        print "ERROR: createHitsGraph"
        sys.exit(1)
    sys.exit(0)
    if not createHoughGraph(MultiGraphs,Tracks):
        print "ERROR: createHoughGraph"
        sys.exit(1)

    sys.exit(0)
    print "HOUGH END!"
