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
    index = track.index(hit)
    x = track.velo_x_hit[index]
    y = track.velo_y_hit[index]
    z = track.velo_z_hit[index]
    return x,y,z

def transform(x1,x2,theta):
    return x1*np.cos(theta) + x2*np.sin(theta)

def makeBinning(val,bot,top,nbins):
    val = [val]
    bins = np.linspace(bot,top,nbins)
    digi = np.digitize(val,bins)
    out = bins[digi-1] + np.diff(bins)[digi-1]/2
    return round(out[0],2)

class myRho():
    """raw and binned Rho values"""
    self.minRho = 0
    self.maxRho = 1000
    self.rhoRes = 1
    def __init__(self, raw):
        super(myRho, self).__init__()
        self.raw = raw
        #self.nrBins = int((self.maxRho - self.minRho) / self.rhoRes)
    def self.nrBins(self):
        return int((self.maxRho - self.minRho) / self.rhoRes)
    def self.binned(self):
        self.digiRho = makeBinning(self.raw, self.minRho, self.maxRho, self.nrBins(self))
        return self.digiRho
    def self.setMaxRho(self,val):
        self.maxRho = val
    def self.setRhoResolution(self,val):
        self.rhoRes = val

def calcRho(plane,th,x,y,z):
    if plane == XY:
        rho = transform(x,y,th)
    if plane == XZ:
        rho = transform(z,x,th)
    if plane == YZ:
        rho = transform(z,y,th)
    #return binRho(rho,maxRho,rhoBinning)
    rho = myRho(rho)
    return rho.binned()

def setDictionaries(dictionaries,tracks):
    for t in itools.ifilter(isGood,tracks):
        #"GET HIT LIST"
        hits = []*t.vplen
        for h in hits:
            x,y,z = getXYZ(h,t)
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
    print "Dicts created"
    sys.exit(0)
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
