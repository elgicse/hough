#!/usr/bin/env python
import ROOT as r
import numpy as np
import itertools as itools
import sys

inrootfile = "data/UTHits.root"
outrootfile = "data/out/out.root"

TREENAME = "UTHits/MatchedTracks"
XY,XZ,YZ = 0,1,2
YX,ZX,ZY = 0,1,2
X,Y,Z,W = 0,1,2,3

class mydict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

def rootPreProcessing():
    r.gROOT.ProcessLine(".x /home/elena/Documenti/University/PhD/UZH/LHCb/lhcbstyle.C")

def mgInit():
    return r.TMultiGraph()

def matrixInit():
    return np.matrix()

def makeThetaArray(mint,maxt,binning):
    theta = np.arange(mint,maxt,binning)
    theta = theta * np.pi / 180
    return theta

def isGood(track):
    return (track.reconstructible_asLong 
        and track.true_p > 3000 
        and track.true_pt > 500 
        and track.true_eta > 2 
        and track.true_eta < 5 
        and abs(track.pid) != 22 
        and track.fromB)

def HGraph(n,x,y):
    gr = r.TGraph(n,x,y)
    gr.SetMarkerStyle(20)
    gr.SetName("tracks_graph") 
    gr.SetTitle("tracks_graph") 
    return gr

def hitGraph(track,plane):
    if plane == XY:
        gr = HGraph(track.vplen, track.velo_x_hit, track.velo_y_hit)
    if plane == XZ:
        gr = HGraph(track.vplen, track.velo_z_hit, track.velo_x_hit)   
    if plane == YZ:
        gr = HGraph(track.vplen, track.velo_z_hit, track.velo_y_hit)   
    return gr

def createHitsGraph(mg,tracks):
    for t in itools.ifilter(isGood,tracks):
        mg[XY].Add(hitGraph(t,XY))
        mg[XZ].Add(hitGraph(t,XZ))
        mg[YZ].Add(hitGraph(t,YZ))
    return True

def getXYZ(hit,track):
    x = track.velo_x_hit[hit]
    y = track.velo_y_hit[hit]
    z = track.velo_z_hit[hit]
    return x,y,z

def transform(x1,x2,theta):
    return x1*np.cos(theta) + x2*np.sin(theta)


class myRho():
    """raw and binned Rho values"""
    def __init__(self, raw, minRho = 0,maxRho = 1000,rhoRes = 1):
        self.raw = raw
        self.minRho = minRho
        self.maxRho = maxRho
        self.rhoRes = rhoRes
        #self.nrBins = int((self.maxRho - self.minRho) / self.rhoRes)
    def nrBins(self):
        return int((self.maxRho - self.minRho) / self.rhoRes)

    def binned(self):
        self.makeBinning()
        return self.digiRho

    def setMaxRho(self,val):
        self.maxRho = val

    def setRhoResolution(self,val):
        self.rhoRes = val

    def makeBinning(self):
        val = [self.raw]
        bins = np.linspace(self.minRho, self.maxRho, self.nrBins())
        digi = np.digitize(val,bins)
        out = bins[digi-1] + np.diff(bins)[digi-1]/2
        self.digiRho = round(out[0],2)
        return True

def calcRho(plane,th,x,y,z):
    if plane == XY:
        rho = transform(x,y,th)
    if plane == XZ:
        rho = transform(z,x,th)
    if plane == YZ:
        rho = transform(z,y,th)
    return myRho(rho)

def setDictionaries(dictionaries,tracks):
    for t in itools.ifilter(isGood,tracks):
        hits = range(t.vplen)
        for h in hits:
            x,y,z = getXYZ(h,t)
            rho = [None]*3  
            for th in theta:
                rho[XY] = calcRho(XY,th,x,y,z).binned()
                rho[XZ] = calcRho(XZ,th,x,y,z).binned()
                rho[YZ] = calcRho(YZ,th,x,y,z).binned()

                dictionaries[XY].setdefault((rho[XY],th), {X:[],Y:[],Z:[],W:0})
                dictionaries[XY][(rho[XY],th)][X].append(x)
                dictionaries[XY][(rho[XY],th)][Y].append(y)
                dictionaries[XY][(rho[XY],th)][Z].append(z)
                dictionaries[XY][(rho[XY],th)][W]+=1

                dictionaries[XZ].setdefault((rho[XZ],th), {X:[],Y:[],Z:[],W:0})
                dictionaries[XZ][(rho[XZ],th)][X].append(x)
                dictionaries[XZ][(rho[XZ],th)][Y].append(y)
                dictionaries[XZ][(rho[XZ],th)][Z].append(z)
                dictionaries[XZ][(rho[XZ],th)][W]+=1

                dictionaries[YZ].setdefault((rho[YZ],th), {X:[],Y:[],Z:[],W:0})
                dictionaries[YZ][(rho[YZ],th)][X].append(x)
                dictionaries[YZ][(rho[YZ],th)][Y].append(y)
                dictionaries[YZ][(rho[YZ],th)][Z].append(z)
                dictionaries[YZ][(rho[YZ],th)][W]+=1



                

    print "DICTS CREATED"
    return True

def dict2Matrix(dictionaries,matrices):
    pass

def searchPeaks(matrices,peaklists):
    pass

def calcParams(peaklists,paramlists):
    pass

def createHoughGraph(mg,tracks):
    return setDictionaries(Dictionaries,tracks)
    dict2Matrix(Dictionaries,Matrices)
    searchPeaks(Matrices,PeakLists)
    calcParams(PeakLists,ParamLists)
    #"create graphs"
    #"add graph in multigraph"
    #"save in root file"
    #"create DTP Dictionary"
    return True

def makeTracklets(params,trackletsContainer):
    pass

def matchTracklets(trackletsContainer,dictionaries,matchedTracklets):
    pass

#################### MAIN PROGRAM ##############################

if __name__ == "__main__":
    print "HOUGH START..."       

    rootPreProcessing()

        #init data structures
    MultiGraphs = mgInit(),mgInit(),mgInit()
    PeakLists = [],[],[]
    ParamLists = [],[],[]
    Tracklets = [],[],[] #tracks on every plane
    Matched = [],[],[] #tracks on every plane
    Dictionaries = mydict(),mydict(),mydict()
    #Matrices = matrixInit(),matrixInit(),matrixInit()
    minTheta, maxTheta, thetaBinning = 0.0, 180.0, 0.2
    theta = makeThetaArray(minTheta,maxTheta,thetaBinning)


        #root file descriptor
    ifd = r.TFile(inrootfile) #input file descriptor 
    ofd = r.TFile(outrootfile,"recreate") #output file descriptor 
        #"SET WRITE FLAG AND CREATE OUT TREE"

        #get tree
    Tracks = ifd.Get(TREENAME)

    if not createHitsGraph(MultiGraphs, Tracks):
        print "ERROR: createHitsGraph"
        sys.exit(1)
    
    #sys.exit(0) #debug
    
    if not createHoughGraph(MultiGraphs, Tracks):
        print "ERROR: createHoughGraph"
        sys.exit(1)
    
    if not makeTracklets(ParamLists, Tracklets):
        print "ERROR: makeTracklets"
        sys.exit(1)

    if not matchTracklets(Tracklets, Dictionaries, Matched):
        print "ERROR: matchTracklets"
        sys.exit(1)

    sys.exit(0)
    print "HOUGH END!"
