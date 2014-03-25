#!/usr/bin/env python
# IMPORTANT:
# Documentation stored at
# /afs/cern.ch/user/e/egraveri/public/VeloHough.pdf
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


def key2index(key,rmin,rmax,thmin,thmax,plane):
    rowIdx = int(np.round( (key[1]-thmin)*(binsx[plane]-1)/(thmax-thmin) ))
    colIdx = int(np.round( (key[0]-rmin)*(binsy[plane]-1)/(rmax-rmin) ))
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
    thetamax = []
    thetamin = []
    plane = 0
    for dictionary in Dictionaries:
        dictionary.update({ "source":{},
                            "c_source":arr.array( 'L', [0]*binsx[plane] ) ,
                            "dest":{},
                            "c_dest":arr.array( 'L', [0]*binsx[plane] )
                            })
        for i in xrange(binsx[plane]):
            dictionary.source[i], dictionary.c_source[i] = new_numpy1d_with_pointer( binsy[plane] )
            dictionary.dest[i], dictionary.c_dest[i]     = new_numpy1d_with_pointer( binsy[plane] )    
        keys = getKeyList(dictionary)
        rhomax.append( max(keys,key=lambda k: k[0])[0] )
        rhomin.append( min(keys,key=lambda k: k[0])[0] )
        thetamax.append( max(keys,key=lambda k: k[1])[1] )
        thetamin.append( min(keys,key=lambda k: k[1])[1] )

        for key in keys:
            if dictionary[key][W] > minHitsPerTrack: #remove background prior to TSpectrum2 evaluation
                i,j = key2index(key,rhomin[plane],rhomax[plane],thetamin[plane],thetamax[plane],plane)
                dictionary.source[i][j] = dictionary[key][W]
        sourceHist = r.TH2F("sourceHist", "sourceHist", binsx[plane], 0, binsx[plane], binsy[plane], 0, binsy[plane])
        for j in range(binsx[plane]):
            for k in range(binsy[plane]):
                sourceHist.SetBinContent(j+1,k+1,dictionary.source[j][k])
        cSaver.append(r.TCanvas())
        sourceHist.Draw("surf2")
        sourceHist.Write()
        plane += 1


class mydict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


def rootPreProcessing():
    r.gROOT.ProcessLine(".x /home/elena/Documenti/University/PhD/UZH/LHCb/lhcbstyle.C")


def mgInit():
    return r.TMultiGraph()


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

def averageNumOfHits(tracks):
    hitsHist = r.TH1F("hitsHist","hitsHist",30,0,29)
    #hitsHist = r.TH1F()
    for t in itools.ifilter(isGood,tracks):
        #print str(t.vplen)
        hitsHist.Fill(t.vplen)
    cSaver.append(r.TCanvas())
    hitsHist.Draw()
    avg, rms = hitsHist.GetMean(), hitsHist.GetRMS()
    print "Average number of hits for selected tracks: %s (sigma: %s)"%(avg,rms)
    minHitsPerTrack = int(avg - 3*rms)
    print "Selected %s as threshold number of hits"%minHitsPerTrack
    return minHitsPerTrack

def createHitsGraph(mg,tracks):
    for t in itools.ifilter(isGood,tracks):
        mg[XY].Add(hitGraph(t,XY))
        mg[XZ].Add(hitGraph(t,XZ))
        mg[YZ].Add(hitGraph(t,YZ))
    for g in xrange(3):
        cSaver.append(r.TCanvas())
        mg[g].SetTitle("HitGraph_plane_"+str(g))
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
    def __init__(self, raw, plane, minRho = [-1000,-400,-400],maxRho = [1000,400,400],rhoRes = [1,0.5,0.5]):
        self.raw = raw
        self.minRho = minRho[plane]
        self.maxRho = maxRho[plane]
        self.rhoRes = rhoRes[plane]
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
        self.recalcBins()
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
    return myRho(rho,plane)



def setDictionaries(dictionaries,tracks,theta):
    global binsx, binsy
    binsx = [None]*3
    binsy = [None]*3
    print "X and Y slopes of selected source tracks:"
    for t in itools.ifilter(isGood,tracks):
        print str(t.tx_velo) + "\t" + str(t.ty_velo)
        hits = xrange(t.vplen)
        for h in hits:
            x,y,z = getXYZ(h,t)
            rho = [None]*3
            hit = [x, y, z]

            for plane in xrange(3):
                for th in theta[plane]:
                    rho[plane] = calcRho(plane,th,x,y,z).binned()
                    dictionaries[plane].setdefault( (rho[plane],th), {'hitlist':[], W:0} )
                    dictionaries[plane][(rho[plane],th)]['hitlist'].append(hit)
                    dictionaries[plane][(rho[plane],th)][W]+=1
                binsx[plane] = len(theta[plane])
                binsy[plane] = calcRho(plane,th,x,y,z).nrBins
    print "DICTS CREATED"
    return True


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
        self.rootTF1.SetMarkerStyle(1)
        return self.rootTF1
    def analyzePeak(self,plane,dict):
        self.theta = thetamin[plane] + self.xPeak * (thetamax[plane]-thetamin[plane])/(binsx[plane]-1)
        self.rho = rhomin[plane] + self.yPeak * (rhomax[plane]-rhomin[plane])/(binsy[plane]-1)
        self.slope = -1 * np.cos(self.theta)/np.sin(self.theta)
        self.intercept = self.rho / np.sin(self.theta)
        for key in getKeyList(dict):
            if areSimilarKeys(plane, key, (self.rho,self.theta)):
                self.hitList += dict[key]['hitlist']
        


def areSimilarKeys(plane,k1,k2):
    """check if two keys are similar"""
    tollBins = 0.75
    rhoToll = tollBins*(rhomax[plane]-rhomin[plane])/binsy[plane]
    thetaToll = tollBins*(thetamax[plane]-thetamin[plane])/binsx[plane]
    diff1 = np.abs(k1[0]-k2[0])
    diff2 = np.abs(k1[1]-k2[1])
    if (diff1<rhoToll) and (diff2<thetaToll):
        return True
    else:
        return False



def searchPeaks(dictionaries,Tracklets):
    plane = -1
    for dict in dictionaries:
        plane += 1
        print "Searching tracks for plane "+str(plane)+"..."
        spectrum = r.TSpectrum2()
        #sigma = 1.3
        #threshold = 23
        sigma = 1.3
        threshold = 10
        #bgRemove = r.kTRUE
        bgRemove = r.kFALSE
        markovReplace = r.kFALSE
        if plane is XY:
            #threshold = 50
            sigma = 4
        nPeaks = spectrum.SearchHighRes(dict.c_source,dict.c_dest,binsx[plane],binsy[plane],sigma,threshold,bgRemove,3,markovReplace,3)
        print "Plane " + str(plane) + ": " + str(nPeaks) + " peaks found in Hough accumulator."
        smoothedHist = r.TH2F("smoothedHist"+str(plane), "smoothedHist"+str(plane), binsx[plane], 0, binsx[plane], binsy[plane], 0, binsy[plane])
        for j in range(binsx[plane]):
            for k in range(binsy[plane]):
                smoothedHist.SetBinContent(j+1,k+1,dict.dest[j][k])
        cSaver.append(r.TCanvas())
        smoothedHist.Draw("surf2")
        smoothedHist.Write()
        peaksFile = open("data/peaks_" + str(plane) + "_th" + str(threshold) + "_sigma" + str(sigma)+ ".dat","w")
        s = ("Found " + str(nPeaks) + " peaks\n")
        peaksFile.write(s)
        s = ("pos1 \t pos2 \t ampl \t theta \t rho \t slope \t intercept\n")
        peaksFile.write(s)
        for p in xrange(nPeaks):
            xh, yh = spectrum.GetPositionX()[p], spectrum.GetPositionY()[p]
            ampl = dict.source[int(xh)][int(yh)]
            peakParams = (xh, yh, ampl)
            trackCandidate = myTrack(peakParams)
            trackCandidate.analyzePeak(plane,dict)
            s = (str(xh) + "\t" + str(yh) + "\t" + str(ampl) + "\t"
                + str(trackCandidate.theta) + "\t" + str(trackCandidate.rho) + "\t" 
                + str(trackCandidate.slope) + "\t" + str(trackCandidate.intercept) + "\n")
            peaksFile.write(s)
            Tracklets[plane].append(trackCandidate)
        peaksFile.close()


def makeTracklets(MultiGraphs,trackLists):
    Matched = matchTracklets(trackLists)
    verticalThreshold = 0.15
    low = [None]*3
    high = [None]*3
    for mtrack in Matched:
        if (np.abs(mtrack.slopes[XZ]) > verticalThreshold) or (np.abs(mtrack.slopes[YZ]) > verticalThreshold):
            Matched.remove(mtrack)
    nMatchedNonVerticalTracks = len(Matched)

    print "Found "+str(nMatchedNonVerticalTracks)+" non-vertical tracks."
    print "X and Y slopes of found tracks:"
    for mtrack in Matched:
        print str(mtrack.slopes[XZ]) + "\t" + str(mtrack.slopes[XY])

    outFileXY = open("data/tracks_XY.dat","w")
    outFileXZ = open("data/tracks_XZ.dat","w")
    outFileYZ = open("data/tracks_YZ.dat","w")
    outFiles = [outFileXY, outFileXZ, outFileYZ]
    s = "slope \t intercept \n"
    for plane in (XY,XZ,YZ):
        cSaver.append(r.TCanvas())
        MultiGraphs[plane].Draw("alp")
        low[plane], high[plane] = MultiGraphs[plane].GetXaxis().GetXmin(), MultiGraphs[plane].GetXaxis().GetXmax()
        outFiles[plane].write(s)
    p = 0
    for mtrack in Matched:
        for plane in (XY,XZ,YZ):
            if (np.abs(mtrack.slopes[plane]) < verticalThreshold) or (plane is XY):
                s = (str(mtrack.slopes[plane]) + "\t" + str(mtrack.intercepts[plane]) + "\n")
                outFiles[plane].write(s)
                func2plot = mtrack.projections[plane].makeTF1(plane, p, low[plane], high[plane])
                funcGraph = r.TGraph(func2plot)
                MultiGraphs[plane].Add(funcGraph)
                p += 1
    for plane in (XY,XZ,YZ):
        outFiles[plane].close()
    return Matched



def matchTracklets(TrackLists):
    matchedtrk = []
    list1 = TrackLists[XY]
    list2 = TrackLists[XZ]
    list3 = TrackLists[YZ]
    matched3DTracks = []
    minHitsForMatching = 4
    print "Now matching tracklets..."
    for track1 in list1:
        for track2 in list2:
            for track3 in list3:
                # Change minhitspertrack to 70% of the number of hits of the tracklet with less hits
                minHitsForMatching = min( minHitsPerTrack, int(0.7 * min(len(track1.hitList), len(track2.hitList), len(track3.hitList))) )
                print minHitsForMatching
                if (correspond(track1.hitList, track2.hitList, minHitsForMatching)
                    and correspond(track1.hitList, track3.hitList, minHitsForMatching)):
                    matchedtrk.append( (track1, track2, track3) )
                    t3d = my3DTrack()
                    t3d.setProjections(track1,track2,track3)
                    matched3DTracks.append(t3d)
    matchedtrk = list(set(matchedtrk))
    matched3DTracks = list(set(matched3DTracks))
    print "Found "+str(len(matched3DTracks))+" matching tracks"
    return matched3DTracks


def correspond(hlist1, hlist2, minHits):
    if (len(hlist1) > minHits) and (len(hlist2) > minHits):
        nCommon = 0
        for hit1 in hlist1:
            for hit2 in hlist2:
                if hit1 == hit2:
                    nCommon += 1
        if nCommon > minHits:
            return True
    else:
        return False


def createHoughGraph(mg,tracks,theta,Tracklets):
    setDictionaries(Dictionaries,tracks,theta)
    makeMatrices(Dictionaries)
    searchPeaks(Dictionaries,Tracklets)
    return True


class my3DTrack():
    """Track candidates"""
    def __init__(self):
        self.hitList = []
        self.slopes = [None]*3
        self.intercepts = [None]*3
        #self.projections = [None]*3
    def setProjection(self,plane,track2D):
        if plane is XY:
            self.trackXY = track2D
        if plane is XZ:
            self.trackXZ = track2D
        if plane is YZ:
            self.trackYZ = track2D
        self.hitList.append(track2D.hitList)
    def setProjections(self,t1,t2,t3):
        self.trackXY = t1
        self.trackXZ = t2
        self.trackYZ = t3
        self.hitList.extend(t1.hitList)
        self.hitList.extend(t2.hitList)
        self.hitList.extend(t3.hitList)
        self.hitList = list(set(map(tuple,self.hitList)))
        self.projections = [self.trackXY, self.trackXZ, self.trackYZ]
        self.computeParameters()
    def addHit(self,hit):
        self.hitList.append(hit)
    def getHitList(self):
        return self.hitList
    #    return list(set(self.hitList))
    def computeParameters(self):
        for plane in xrange(3):
            self.slopes[plane] = self.projections[plane].slope
            self.intercepts[plane] = self.projections[plane].intercept
        #self.slopes[XZ] = self.trackXZ.slope
        #self.slopes[YZ] = self.trackYZ.slope
        #self.slopes[XY] = self.trackXY.slope
        #self.intercepts[XZ] = self.trackXZ.intercept
        #self.intercepts[YZ] = self.trackYZ.intercept
        #self.intercepts[XY] = self.trackXY.intercept




#################### MAIN PROGRAM ##############################

if __name__ == "__main__":
    print "HOUGH START..."

    global cSaver, minHitsPerTrack
    cSaver = []
    #minHitsPerTrack = 4

    rootPreProcessing()

        #init data structures
    MultiGraphs = mgInit(),mgInit(),mgInit()
    PeakLists = [],[],[]
    ParamLists = [],[],[]
    Tracklets = [],[],[] #tracks on every plane
    Dictionaries = mydict(),mydict(),mydict()
    minTheta = [0.0, 70.0, 70.0]
    maxTheta = [180.0, 110.0, 110.0]
    thetaBinning = [0.5, 0.1, 0.1]
    theta = [None]*3
    theta[XY] = makeThetaArray(minTheta[XY],maxTheta[XY],thetaBinning[XY])
    theta[XZ] = makeThetaArray(minTheta[XZ],maxTheta[XZ],thetaBinning[XZ])
    theta[YZ] = makeThetaArray(minTheta[YZ],maxTheta[YZ],thetaBinning[YZ])

        #root file descriptor
    ifd = r.TFile(inrootfile) #input file descriptor 
    ofd = r.TFile(outrootfile,"recreate") #output file descriptor 

        #get tree
    inputTracks = ifd.Get(TREENAME)

    minHitsPerTrack = averageNumOfHits(inputTracks)

    if not minHitsPerTrack:
        print "ERROR: averageNumOfHits"
        sys.exit(1)

    if not createHitsGraph(MultiGraphs, inputTracks):
        print "ERROR: createHitsGraph"
        sys.exit(1)
    
    if not createHoughGraph(MultiGraphs, inputTracks, theta, Tracklets):
        print "ERROR: createHoughGraph"
        sys.exit(1)
    
    Matched = makeTracklets(MultiGraphs, Tracklets)
    if not Matched:
        print "ERROR: makeTracklets"
        sys.exit(1)

    print "Saving data..."
    for mg in MultiGraphs:
        cSaver.append(r.TCanvas())
        mg.Draw("alp")
        mg.Write()

    ofd.Close()

    print "HOUGH END."
    sys.exit(0)
