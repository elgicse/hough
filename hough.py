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
    return keylist



def makeMatrices(Dictionaries):
    global rhomax, rhomin, thetamax, thetamin
    rhomax = []
    rhomin = []
    plane = 0
    for dictionary in Dictionaries:
        dictionary.update({ "source":{},
                            "c_source":arr.array( 'L', [0]*binsx ) ,
                            "dest":{},
                            "c_dest":arr.array( 'L', [0]*binsx )
                            })
        for i in xrange(binsx):
            dictionary.source[i], dictionary.c_source[i] = new_numpy1d_with_pointer( binsy )
            dictionary.dest[i], dictionary.c_dest[i]     = new_numpy1d_with_pointer( binsy )    
        keys = getKeyList(dictionary)
        rhomax.append( max(keys,key=lambda k: k[0])[0] )
        rhomin.append( min(keys,key=lambda k: k[0])[0] )
        if plane is 0: # theta is always the same array
            thetamax = max(keys,key=lambda k: k[1])[1]
            thetamin = min(keys,key=lambda k: k[1])[1]
        for key in keys:
            i,j = key2index(key,rhomin[plane],rhomax[plane],thetamin,thetamax)
            dictionary.source[i][j] = dictionary[key][W]
        sourceHist = r.TH2F("sourceHist", "sourceHist", binsx, 0, binsx, binsy, 0, binsy)
        for j in range(binsx):
            for k in range(binsy):
                sourceHist.SetBinContent(j+1,k+1,dictionary.source[j][k])
        sourceHist.Draw("surf2")
        sourceHist.Write()
        plane += 1
# To view the plot do:
#    while True:
#        time.sleep(5)



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
        mg[g].SetTitle("HitGraph_plane_"+str(g))
        mg[g].Draw("ap")
        mg[g].Write()
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

#                dictionaries[XY].setdefault((rho[XY],th), {X:[],Y:[],Z:[],W:0})
#                dictionaries[XY][(rho[XY],th)][X].append(x)
#                dictionaries[XY][(rho[XY],th)][Y].append(y)
#                dictionaries[XY][(rho[XY],th)][Z].append(z)
#                dictionaries[XY][(rho[XY],th)][W]+=1
#
#                dictionaries[XZ].setdefault((rho[XZ],th), {X:[],Y:[],Z:[],W:0})
#                dictionaries[XZ][(rho[XZ],th)][X].append(x)
#                dictionaries[XZ][(rho[XZ],th)][Y].append(y)
#                dictionaries[XZ][(rho[XZ],th)][Z].append(z)
#                dictionaries[XZ][(rho[XZ],th)][W]+=1
#
#                dictionaries[YZ].setdefault((rho[YZ],th), {X:[],Y:[],Z:[],W:0})
#                dictionaries[YZ][(rho[YZ],th)][X].append(x)
#                dictionaries[YZ][(rho[YZ],th)][Y].append(y)
#                dictionaries[YZ][(rho[YZ],th)][Z].append(z)
#                dictionaries[YZ][(rho[YZ],th)][W]+=1

                hit = [x, y, z]

                dictionaries[XY].setdefault((rho[XY],th), {'hitlist':[], W:0})
                dictionaries[XY][(rho[XY],th)]['hitlist'].append(hit)
                dictionaries[XY][(rho[XY],th)][W]+=1

                dictionaries[XZ].setdefault((rho[XZ],th), {'hitlist':[], W:0})
                dictionaries[XZ][(rho[XZ],th)]['hitlist'].append(hit)
                dictionaries[XZ][(rho[XZ],th)][W]+=1

                dictionaries[YZ].setdefault((rho[YZ],th), {'hitlist':[], W:0})
                dictionaries[YZ][(rho[YZ],th)]['hitlist'].append(hit)
                dictionaries[YZ][(rho[YZ],th)][W]+=1


    global binsx, binsy, rhoSpace
    rhoSpace = calcRho(YZ,th,x,y,z).bins
    binsy = calcRho(YZ,th,x,y,z).nrBins
    binsx = len(theta)   
    print "DICTS CREATED"
    return True


#def linspaceToAxis(linspace, aMin, aMax):
#    """LIST OF VALUES TO LIST OF BIN INDICES"""
#    aBins = len(linspace)
#    axis = np.linspace
#    return axis
#    """RIFARE"""



class myTrack():
    """Track candidates"""
    def __init__(self, peakParam):
        self.xPeak = peakParam[0]
        self.yPeak = peakParam[1]
        self.peakAmplitude = peakParam[2]
        self.hitList = []
    def getHitList(self):
        return self.hitList
    def makeTF1(self, plane, index, lowMargin, highMargin):
        self.rootTF1 = r.TF1("tracklet"+str(plane)+"_"+str(index),"pol1",lowMargin,highMargin)
        self.rootTF1.SetParNames("Intercept","Slope")
        self.rootTF1.FixParameter(0,self.intercept)
        self.rootTF1.FixParameter(1,self.slope)
        self.rootTF1.SetLineColor(r.kRed)
        return self.rootTF1
    def analyzePeak(self,plane,dict):
        self.theta = thetamin + self.xPeak * (thetamax-thetamin)/binsx
        self.rho = rhomin[plane] + self.yPeak * (rhomax[plane]-rhomin[plane])/binsy
        self.slope = -1 * np.cos(self.theta)/np.sin(self.theta)
        self.intercept = self.rho / np.sin(self.theta)
        for key in getKeyList(dict):
            if areSimilarKeys(plane, key, (self.rho,self.theta)):
                self.hitList += dict[key]['hitlist']
        #self.hitList = [(hit for hit in dict[key]['hitlist']) for key in getKeyList(dict) if areSimilarKeys(plane, key, (self.rho,self.theta))]
        #self.hitList = [dict[key] for key in getKeyList(dict) if (key is (self.rho,self.theta))]
        


def areSimilarKeys(plane,k1,k2):
    """check if two keys are similar"""
    rhoToll = 2*(rhomax[plane]-rhomin[plane])/binsy
    thetaToll = 2*(thetamax-thetamin)/binsx

    diff1 = np.abs(k1[0]-k2[0])
    diff2 = np.abs(k1[1]-k2[1])

    if (diff1<rhoToll*np.abs(k1[0])) and (diff2<thetaToll*np.abs(k1[1])):
        return True
    else:
        return False



def searchPeaks(dictionaries,Tracklets):
    #trackLists = [], [], []
    sigma = 1.3
    threshold = 22
    #pList = [], [], []
    #paramList = [], [], []
    iPlane = -1
    for dict in dictionaries:
        iPlane += 1
        # Apply TSpectrum class methods
        spectrum = r.TSpectrum2()
        nPeaks = spectrum.SearchHighRes(dict.c_source,dict.c_dest,binsx,binsy,sigma,threshold,r.kTRUE,3,r.kFALSE,3)
        # Draw and save the smoothed histogram
        smoothedHist = r.TH2F("smoothedHist"+str(iPlane), "smoothedHist"+str(iPlane), binsx, 0, binsx, binsy, 0, binsy)
        for j in range(binsx):
            for k in range(binsy):
                smoothedHist.SetBinContent(j+1,k+1,dict.dest[j][k])
        smoothedHist.Draw("surf2")
        smoothedHist.Write()
        # Extract and write peak parameters
        peaksFile = open("data/peaks" + str(iPlane) + "_th" + str(threshold) + "_sigma" + str(sigma)+ ".dat","w")
        s = ("Found " + str(nPeaks) + " peaks\n")
        peaksFile.write(s)
        s = ("pos1 \t pos2 \t ampl \t theta \t rho \t slope \t intercept\n")
        peaksFile.write(s)
        # Make tracks from peaks
        for p in xrange(nPeaks):
            xh, yh = spectrum.GetPositionX()[p], spectrum.GetPositionY()[p]
            ampl = dict.source[int(xh)][int(yh)]
            peakParams = (xh, yh, ampl)
            trackCandidate = myTrack(peakParams)
            trackCandidate.analyzePeak(iPlane,dict)
            s = (str(xh) + "\t" + str(yh) + "\t" + str(ampl) + "\t"
                + str(trackCandidate.theta) + "\t" + str(trackCandidate.rho) + "\t" 
                + str(trackCandidate.slope) + "\t" + str(trackCandidate.intercept) + "\n")
            peaksFile.write(s)
            #trackLists[iPlane].append(trackCandidate)
            Tracklets[iPlane].append(trackCandidate)
        peaksFile.close()


def makeTracklets(MultiGraphs,trackLists,Matched):
    # Match tracks in every plane
    Matched = matchTracklets(trackLists)
    # Make TF1s out of matched tracks
    verticalThreshold = 0.15
    iPlane = -1
    for matchedList in Matched: #for every plane
        iPlane += 1
        MultiGraphs[iPlane].Draw("ap")
        low, high = MultiGraphs[iPlane].GetXaxis().GetXmin(), MultiGraphs[iPlane].GetXaxis().GetXmax()
        p = 0
        for track in matchedList:
            if np.abs(track.slope) > verticalThreshold: # ma serve ancora dopo il matching?
                func2plot = track.makeTF1(iPlane, p, low, high)
                funcGraph = r.TGraph(func2plot)
                MultiGraphs[iPlane].Add(funcGraph)
                p += 1
    return 1



def matchTracklets(TrackLists):
    matched = [], [], []
    list1 = TrackLists[XY]
    list2 = TrackLists[XZ]
    list3 = TrackLists[YZ]
    for track1 in list1:
        for track2 in list2:
            for track3 in list3:
                if correspondingHits(track1, track2, X) and correspondingHits(track1, track3, Y):
                    #matched.append( (track1, track2, track3) )
                    matched[XY].append(track1)
                    matched[XZ].append(track2)
                    matched[YZ].append(track3)
    return matched


def correspondingHits(track1, track2, axis):
    hitList1 = track1.getHitList()
    hitList2 = track2.getHitList()
    avgNumberOfHits = 0.5 * ( len(hitList1) + len(hitList2) )
    minHitsInCommon = int(0.7 * avgNumberOfHits) - 1
    nHitsInCommon = 0
    #Questo non va bene: le liste vanno ribaltate!
    for hit1 in hitList1:
        for hit2 in hitList2:
            if hit1[axis] == hit2[axis]:
                nHitsInCommon += 1

    # Hit list structure (dictionaries):
    # list of x ; list of y ; list of z ; weight
    #for i in xrange( len(hitList1[axis]) ):
    #    for j in xrange( len(hitList2[axis]) ):
    #        if hitList1[axis][i] == hitList2[axis][j]:
    #            nHitsInCommon += 1

    if nHitsInCommon > minHitsInCommon:
        return True
    else:
        return False


def createHoughGraph(mg,tracks,theta,Tracklets):
    setDictionaries(Dictionaries,tracks,theta)
    makeMatrices(Dictionaries)
    searchPeaks(Dictionaries,Tracklets)
    #calcParams(PeakLists,ParamLists)
    #"create graphs"
    #"add graph in multigraph"
    #"save in root file"
    #"create DTP Dictionary"
    return True



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
    inputTracks = ifd.Get(TREENAME)

    if not createHitsGraph(MultiGraphs, inputTracks):
        print "ERROR: createHitsGraph"
        sys.exit(1)
    
    if not createHoughGraph(MultiGraphs, inputTracks, theta, Tracklets):
        print "ERROR: createHoughGraph"
        sys.exit(1)
    
    if not makeTracklets(MultiGraphs, Tracklets, Matched):
        print "ERROR: makeTracklets"
        sys.exit(1)

    for mg in MultiGraphs:
        mg.Draw("alp")
        mg.Write()

    ofd.Close()
    #if not matchTracklets(Tracklets, Dictionaries, Matched):
    #    print "ERROR: matchTracklets"
    #    sys.exit(1)

    sys.exit(0)
    print "HOUGH END!"
