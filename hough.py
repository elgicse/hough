#!/usr/bin/env python
import ROOT as r
import numpy as np
import itertools as itools
import sys
import array as arr
import time

inrootfile = "data/UTHits.root"
outrootfile = "data/out/out.root"

TREENAME = "UTHits/MatchedTracks"
XY,XZ,YZ = 0,1,2
YX,ZX,ZY = 0,1,2
X,Y,Z,W = 0,1,2,3


def new_numpy1d_with_pointer(size):
    np_a = np.zeros( size, dtype=np.float32 )
    pointer, read_only_flag = np_a.__array_interface__['data']
    return np_a, pointer

def key2index(key,rhomin,rhomax,thetamin,thetamax):
    #key[1] is theta, key[0] is rho
    #rows are thetas, columns are rhos

    #col = int(np.round(key[1]*(rhomax-rhomin)/binsy))
    #row = int(np.round(key[0]*(thetamax-thetamin)/binsx))

    rowIdx = int(np.round( (key[1]-thetamin)*(binsx-1)/(thetamax-thetamin) ))
    colIdx = int(np.round( (key[0]-rhomin)*(binsy-1)/(rhomax-rhomin) ))
    return rowIdx,colIdx


def getKeyList(dict):
    keylist = dict.keys()
    del keylist[keylist.index("source")]
    del keylist[keylist.index("dest")]
    del keylist[keylist.index("c_source")]
    del keylist[keylist.index("c_dest")]

def makeMatrices(dictionary):
    
    dictionary.update({ "source":{},
                        "c_source":arr.array( 'L', [0]*binsx ) ,
                        "dest":{},
                        "c_dest":arr.array( 'L', [0]*binsx )
                        })

    for i in xrange(binsx):
        dictionary.source[i], dictionary.c_source[i] = new_numpy1d_with_pointer( binsy )
        dictionary.dest[i], dictionary.c_dest[i]     = new_numpy1d_with_pointer( binsy )    

    keys = getKeyList(dictionary)

    global rhomax = max(keys,key=lambda k: k[0])[0]
    global rhomin = min(keys,key=lambda k: k[0])[0]
    global thetamax = max(keys,key=lambda k: k[1])[1]
    global thetamin = min(keys,key=lambda k: k[1])[1]

    for key in keys:
        i,j = key2index(key,rhomin,rhomax,thetamin,thetamax)
        dictionary.source[i][j] = dictionary[key][W]



    sourceHist = r.TH2F("sourceHist", "sourceHist", binsx, 0, binsx, binsy, 0, binsy)
    for j in range(binsx):
        for k in range(binsy):
            sourceHist.SetBinContent(j+1,k+1,dictionary.source[j][k])

    sourceHist.Draw("surf2")
    sourceHist.Write()

    while True:
        time.sleep(5)

    



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
    for g in xrange(2):
        mg[g].Draw("ap")
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
    def __init__(self, raw, minRho = -1000,maxRho = 1000,rhoRes = 1):
        self.raw = raw
        self.minRho = minRho
        self.maxRho = maxRho
        self.rhoRes = rhoRes
        self.recalcBins()

    def recalcBins(self):
        self.nrBins = int((self.maxRho - self.minRho) / self.rhoRes)
        self.bins = np.linspace(self.minRho, self.maxRho, self.nrBins)

    def binned(self):
        self.makeBinning()
        return self.digiRho

    def setMaxRho(self,val):
        self.maxRho = val
        self.recalcBins()

    def setRhoResolution(self,val):
        self.rhoRes = val
        self.recalcBins()

    def makeBinning(self):
        val = [self.raw]
        #self.bins = np.linspace(self.minRho, self.maxRho, self.nrBins)
        digi = np.digitize(val,self.bins)
        out = self.bins[digi-1] + np.diff(self.bins)[digi-1]/2
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

def setDictionaries(dictionaries,tracks,theta):
    for t in itools.ifilter(isGood,tracks):
        hits = xrange(t.vplen)
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

    global binsx, binsy, rhoSpace
    rhoSpace = calcRho(YZ,th,x,y,z).bins
    binsy = calcRho(YZ,th,x,y,z).nrBins
    binsx = len(theta)   

    print "DICTS CREATED"
    return True

def linspaceToAxis(linspace, aMin, aMax):
    """LIST OF VALUES TO LIST OF BIN INDICES"""
    aBins = len(linspace)
    axis = np.linspace
    return axis
    """RIFARE"""

def sumMatrices(plane):
    pass

def searchPeaks(dictionaries,peaklists):
    """PSEUDO - CODE """
    sigma = 1.3
    threshold = 22
    pList = [], [], []
    paramList = [], [], []
    iPlane = -1
    for dict in dictionaries:
        iPlane += 1

        # Apply TSpectrum class methods
        spectrum = r.TSpectrum()
        nPeaks = spectrum.SearchHighRes(dict.c_source,dict.c_dest,binsx,binsy,sigma,threshold,r.kTRUE,3,r.kFALSE,3)

        # Draw and save the smoothed source histogram
        smoothedHist = r.TH2F("smoothedHist"+str(iPlane), "smoothedHist"+str(iPlane), binsx, 0, binsx, binsy, 0, binsy)
        for j in range(binsx):
            for k in range(binsy):
                smoothedHist.SetBinContent(j+1,k+1,dict.dest[j][k])
        smoothedHist.Draw("surf2")
        smoothedHist.Write()

        # Extract and write peak parameters
        peaksFile = open("data/peaks" + str(iPlane) + "_th" + str(threshold) + "_sigma" + str(sigma)+ ".dat","w")
        s = ("Found " + str(npeaks) + " peaks\n")
        peaksFile.write(s)
        s = ("pos1 \t pos2 \t ampl \t theta \t rho \t slope \t intercept\n")
        peaksFile.write(s)
        for p in xrange(nPeaks):
            hitList = []
            xh = spectrum.GetPositionX()[p]
            yh = spectrum.GetPositionY()[p]
            ampl = dict.source[int(xh)][int(yh)]
            #pList[iPlane].append( (xh,yh,ampl) )
            theta,rho,slope,intercept,hitList = analyzePeak(xh,yh,dict)
            s = (str(xh) + "\t" + str(yh) + "\t" + str(ampl) + "\t"
                + str(theta) + "\t" + str(rho) + "\t" 
                + str(slope) + "\t" + str(intercept) + "\n")
            peaksFile.write(s)
            paramList[iPlane].append( (xh,yh,ampl,theta,rho,slope,intercept) )
        peaksFile.close()


def analyzePeak(xh,yh,dict):
    th = thetamin + xh * (thetamax-thetamin)/binsx
    rh = rhomin + yh * (rhomax-rhomin)/binsy
    slope = -1 * np.cos(th)/np.sin(th)
    intercept = rh / np.sin(th)
    hitList = [dict[key] for key in getKeyList(dict) if areSimilarKeys(key, (xh,yh))]
    return th, rh, slope, intercept, hitList


def calcParams(peaklists,paramlists):  #un ciclo di troppo?
    pass

def createHoughGraph(mg,tracks,theta):
    setDictionaries(Dictionaries,tracks,theta)
    makeMatrices(Dictionaries[XZ])
    return
    searchPeaks(dictionaries,PeakLists)
    #calcParams(PeakLists,ParamLists)
    #"create graphs"
    #"add graph in multigraph"
    #"save in root file"
    #"create DTP Dictionary"
    return True

def makeTracklets(mgContainer,params,trackletsContainer):
    verticalThreshold = 0.15
    iPlane = -1
    for parList in params: #for every plane
        iPlane += 1
        low, high = mgContainer[iPlane].GetXaxis().GetXmin(), mgContainer[iPlane].GetXaxis().GetXmax()
        p = 0
        for pars in parList: #pars[5] = slope, pars[6] = intercept
            if np.abs(pars[5]) > verticalThreshold:
                track = r.TF1("tracklet"+str(iPlane)+"_"+str(p),"pol1",low,high) #completare costruttore!!!!!!!!
                track.SetParNames("Intercept","Slope")
                track.FixParameter(0,pars[6])
                track.FixParameter(1,pars[5])
                trackletsContainer[iPlane].append(track)
                p += 1


def areSimilarKeys(k1,k2):
    """check if two keys are similar"""
    diff1 = np.abs(k1[0]-k2[0])
    diff2 = np.abs(k1[1]-k2[1])
    if (diff1<0.05*k1[0]) and (diff2<0.05*k1[1]):
        return 1
    else:
        return 0


def matchTracklets(trackletsContainer,dictionaries,matchedTracklets):
    pass

#################### MAIN PROGRAM ##############################

if __name__ == "__main__":
    print "HOUGH START..."       

    binsx, binsy, rhoSpace = 0, 0, np.linspace(0, 1000, 1000)

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
    
    if not createHoughGraph(MultiGraphs, Tracks, theta):
        print "ERROR: createHoughGraph"
        sys.exit(1)
    
    if not makeTracklets(MultiGraphs, ParamLists, Tracklets):
        print "ERROR: makeTracklets"
        sys.exit(1)

    if not matchTracklets(Tracklets, Dictionaries, Matched):
        print "ERROR: matchTracklets"
        sys.exit(1)

    sys.exit(0)
    print "HOUGH END!"
