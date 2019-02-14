#functions for manipulating sequences and alignments, working with sliding windows, doing population genetics etx.

import numpy as np
import math
from copy import deepcopy


##################################################################################################################
#Bits for intyerpreting and manipulating sequence data

DIPLOTYPES = ['A', 'C', 'G', 'K', 'M', 'N', 'S', 'R', 'T', 'W', 'Y']
PAIRS = ['AA', 'CC', 'GG', 'GT', 'AC', 'NN', 'CG', 'AG', 'TT', 'AT', 'CT']
diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))

def haplo(diplo): return diploHaploDict[diplo]

def diplo(pair): return haploDiploDict[pair]


#convert one ambiguous sequence into two haploid pseudoPhased sequences

def pseudoPhase(sequence, seqType = "diplo"):
    if seqType == "diplo": pairs = [haplo(s) for s in sequence]
    else: pairs = sequence
    return [[p[0] for p in pairs], [p[1] for p in pairs]]


################################################################################################################

#modules for working with and analysing alignments

numSeqDict = {"A":0,"C":1,"G":2,"T":3,"N":np.NaN}

class Alignment:
    def __init__(self, sequences = None, names=None, groups = None, length = None, numArray = None):
        assert sequences is not None or length is not None, "Specify either sequences or length of empty sequence object."
        if sequences is not None:
            assert isinstance(sequences, (list,tuple,np.ndarray)), "Sequences must be a list, tuple or numpy array."
            if isinstance(sequences, np.ndarray): self.array = sequences
            else: self.array = np.array([list(seq) for seq in sequences])
        else:
            self.array = np.empty((0,length))
            self.numArray = np.empty((0,length))
        
        if numArray is not None:
            assert numArray.shape == sequences.shape, "Numeric array is different shape from sequence array."
            self.numArray = numArray
        else:
            self.numArray = np.array([[numSeqDict[b] for b in seq] for seq in sequences])
         
        self.nanMask = ~np.isnan(self.numArray)
        
        self.N,self.l = self.array.shape
        
        if not names: names = range(self.N)
        else: assert len(names) == self.N, "Incorrect number of names."
        self.names = names
        
        if not groups: groups = [None]*self.N
        else: assert len(groups) == self.N, "Incorrect number of groups."
        self.groups = groups
                
    def addSeq(self, sequence, name = None, group = None):
        self.array = np.vstack((self.array, np.array(list(sequence))))
        self.numArray = np.vstack((self.numArray, np.array([numSeqDict[b] for b in sequence])))
        self.N += 1
        if not name: name = self.N
        self.names.append(name)
        self.groups.append(group)
        self.nanMask = ~np.isnan(self.numArray)
    
    def subset(self, indices = [], names = [], groups = []):
        indices = [i for i in xrange(self.N) if i in indices or
                   (self.names[i] and self.names[i] in names) or
                   (self.groups[i] and self.groups[i] in groups)]
        
        new = Alignment(sequences = self.array[indices], numArray=self.numArray[indices],
                        names=[self.names[i] for i in indices], groups=[self.groups[i] for i in indices])
        return new
    
    def column(self,x): return self.array[:,x]
    
    def numColumn(self,x): return self.numArray[:,x]
    
    def distMatrix(self):
        distMat = np.zeros((self.N,self.N))
        for i in range(self.N - 1):
            for j in range(i + 1, self.N):
                distMat[i,j] = distMat[j,i] = numHamming(self.numArray[i,:], self.numArray[j,:])
        return distMat
    
    def varSites(self): return np.where([np.unique(self.numArray[:,x][self.nanMask[:,x]]) > 1 for x in xrange(self.l)])[0]
    
    def biSites(self): return np.where([len(np.unique(self.numArray[:,x][self.nanMask[:,x]]) == 2) for x in xrange(self.l)])[0]
        
    def siteFreqs(self, sites):
        if not sites: sites = xrange(self.l)
        if type(sites) is not list: sites = [sites]
        return [binFreqs(self.numArray[:,x][self.nanMask[:,x]].astype(int)) for x in sites]




def callsToAlignment(seqs, sampleData, seqType = "diplo"):
    seqNames = []
    groups = []
    pseudoPhasedSeqs = []
    #first pseudo phase all seqs
    for indName in seqs.keys():
        if sampleData.ploidy[indName] == 2:
            pseudoPhasedSeqs += pseudoPhase(seqs[indName], seqType)
            seqNames += [indName + "A", indName + "B"]
            groups += [sampleData.getPop(indName)]*2
        else:
            pseudoPhased[indName] = seqs[indName]
            seqNames.append(indName)
            groups.append(sampleData.getPop(indName))
    return Alignment(sequences=pseudoPhasedSeqs, names = seqNames, groups=groups)



def binFreqs(numArr):
    n = len(numArr)
    if n == 0: return np.array([np.NaN]*4)
    else: return 1.* np.bincount(numArr, minlength=4) / n



## a distance matrix method that uses numerical arrays 
## there is considerable overhead in making the arrays,
## so this isn't fast for pair-wise distance, but is good for sets,
## as you onlyhave to make the array once for each.

def numHamming(numArrayA, numArrayB):
    dif = numArrayA - numArrayB
    return np.nanmean(dif[~np.isnan(dif)] != 0)


def distMatrix(sequences):
    numSeqs = [[numSeqDict[b] for b in seq] for seq in seqs]
    DNAarray = np.array(numSeqs)
    N,ln = DNAarray.shape
    distMat = np.zeros((N,N))
    for i in range(N - 1):
        for j in range(i + 1, N):
            distMat[i,j] = distMat[j,i] = numHamming(DNAarray[i,:], DNAarray[j,:])
    return distMat


class SampleData:
    def __init__(self, indNames = [], popNames = None, popInds = [], popNumbers = None, ploidy = None):
        if not popNumbers:
            popNumbers = range(len(popInds))
        if not popNames:
            popNames = [str(x) for x in popNumbers]
        assert len(popNames) == len(popInds) == len(popNumbers), "Names, inds and numbers should be same length."
        self.popNames = popNames
        self.popNumbers = popNumbers
        self.popInds = {}
        for x in range(len(popInds)):
            for indName in popInds[x]:
                if indName not in indNames:
                    indNames.append(indName)
            self.popInds[popNames[x]] = popInds[x]
            self.popInds[popNumbers[x]] = popInds[x]
        self.indNames = indNames
        if not ploidy:
            ploidy = [2]*len(indNames)
        assert len(indNames) == len(ploidy), "ploidy and indNames should be same length."
        self.ploidy = dict(zip(indNames, ploidy))
    
    def getPop(self, indName):
        pop = [p for p in self.popNames if indName in self.popInds[p]]
        if len(pop) == 0: return None
        elif len(pop) == 1: return pop[0]
        else: return tuple(pop)
    
    def getPopNumber(self, popName):
        if popName in self.popNames:
            return self.popNumbers[self.popNames.index(popName)]


def popDiv(Aln):
    distMat = Aln.distMatrix()
    np.fill_diagonal(distMat, np.NaN) # set all same-with-same to Na
    
    pops,indices = np.unique(Aln.groups, return_inverse = True)
    nPops = len(pops)
    assert nPops > 1, "At least two populations required."
    
    #get population indices - which positions in the alignment correspond to each population
    # this will allow indexing specific pops from the matrix.
    popIndices = [list(np.where(indices==x)[0]) for x in range(nPops)]
    
    output = {}
    
    #pi for each pop
    for x in range(nPops):
        output["pi_" + pops[x]] = np.nanmean(distMat[np.ix_(popIndices[x],popIndices[x])])
    
    #pairs
    for x in range(nPops-1):
        for y in range(x+1, nPops):
            #dxy
            output["dxy_" + pops[x] + "_" + pops[y]] = output["dxy_" + pops[y] + "_" + pops[x]] = np.nanmean(distMat[np.ix_(popIndices[x],popIndices[y])])
            
            #fst
            n_x = len(popIndices[x])
            n_y = len(popIndices[y])
            w = 1.* n_x/(n_x + n_y)
            pi_s = w*(output["pi_" + pops[x]]) + (1-w)*(output["pi_" + pops[y]])
            pi_t = np.nanmean(distMat[np.ix_(popIndices[x]+popIndices[y],popIndices[x]+popIndices[y])])
            output["fst_" + pops[x] + "_" + pops[y]] = output["fst_" + pops[y] + "_" + pops[x]] = 1 - pi_s/pi_t
    
    return output


def ABBABABA(Aln, P1, P2, P3, P4):
    #subset by population
    P1Aln = Aln.subset(groups=[P1])
    P2Aln = Aln.subset(groups=[P2])
    P3Aln = Aln.subset(groups=[P3])
    P4Aln = Aln.subset(groups=[P4])
    P123Aln = Aln.subset(groups=[P1,P2,P3,P4])
    ABBAsum = BABAsum = maxABBAsum = maxBABAsum = 0.0
    #get derived frequencies for all biallelic siites
    for i in P123Aln.biSites():
        allFreqs = Aln.siteFreqs(i)[0] #an array with 4 values, the freq for A,C,G and T
        # get frequencies for wach pop
        P1Freqs,P2Freqs,P3Freqs,P4Freqs = [A.siteFreqs(i)[0] for A in (P1Aln, P2Aln, P3Aln, P4Aln)]
        #check for bad data
        if np.any(np.isnan(P1Freqs)) or np.any(np.isnan(P2Freqs)) or np.any(np.isnan(P3Freqs)) or np.any(np.isnan(P4Freqs)): continue
        #if the outgroup is fixed, then that is the ancestral state - otherwise the derived state is the most common allele overall
        if np.max(P4Freqs) == 1.:
            anc = np.where(P4Freqs == 1)[0][0] #ancetral allele is which is fixed (get the index)
            der = [i for i in np.where(allFreqs > 0)[0] if i != anc][0] # derived is the index that is > 0 but not anc
        else:der = np.argsort(allFreqs)[-2] # the less common base overall
        #derived allele frequencies
        P1derFreq = P1Freqs[der]
        P2derFreq = P2Freqs[der]
        P3derFreq = P3Freqs[der]
        P4derFreq = P4Freqs[der]
        PDderFreq = max(P2derFreq,P3derFreq)
        # get weigtings for ABBAs and BABAs
        ABBAsum += (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq)
        BABAsum += P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        maxABBAsum += (1 - P1derFreq) * PDderFreq * PDderFreq * (1 - P4derFreq)
        maxBABAsum += P1derFreq * (1 - PDderFreq) * PDderFreq * (1 - P4derFreq)
    #calculate D, fd
    output = {}
    try: output["D"] = (ABBAsum - BABAsum) / (ABBAsum + BABAsum)
    except: output["D"] = np.NaN
    try:
        if output["D"] >= 0: output["fd"] = (ABBAsum - BABAsum) / (maxABBAsum - maxBABAsum)
        else: output["fd"] = np.NaN
    except: output["fd"] = np.NaN
    output["ABBA"] = ABBAsum
    output["BABA"] = BABAsum
    
    return output


################################################################################################

#modules for working with sliding windows

#Window object class, stores names, sequences and window information

class SeqWindow: 
    def __init__(self, scaffold = None, start = None, end = None, seqs = None, names = None, positions = None, ID = None):
        if not names and not seqs:
            names = []
            seqs = []
        elif not names:
            names = [None]*len(seqs)
        elif not seqs:
            seqs = [[] for name in names]
        assert len(names) == len(seqs)
        if not positions:
            positions = []
        if len(seqs) > 0:
            #print len(seqs[0]), len(positions)
            assert len(seqs[0]) == len(positions) # ensure correct number of positions is given
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs:
                assert type(seq) is list # added sequences must each be a list
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.names = names
        self.positions = positions
        self.seqs = seqs
        self.n = len(self.names)
        self.ID = ID
    
    #method for adding
    def addBlock(self, seqs, positions):
        assert len(seqs) == self.n # ensure correct number seqs is added
        if len(seqs) > 0:
            assert len(seqs[0]) == len(positions) # ensure correct number of positions is given
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs:
                assert type(seq) is list # added sequences must each be a list
        for x in range(len(seqs)):
            self.seqs[x] += seqs[x]
        self.positions += positions
    
    def addSite(self, calls, position):
        assert len(calls) == self.n # ensure correct number seqs is added
        for x in range(self.n):
            self.seqs[x].append(calls[x])
        self.positions.append(position)
    
    def seqLen(self):
        return len(self.positions)
    
    def firstPos(self):
        return min(self.positions)
    
    def lastPos(self):
        return max(self.positions)
    
    def slide(self,step=None,newStart=None,newEnd=None):
        #function to slide window along scaffold
        assert step != None or newStart != None 
        if step:
            newStart = self.start + step
            newEnd = self.end + step            
        self.start = newStart
        if newEnd:
            self.end = newEnd
        #find first position beyon newStart
        i = 0
        while i < len(self.positions) and self.positions[i] < newStart:
            i += 1
        #slide positions
        self.positions = self.positions[i:]
        #slide seqs
        self.seqs = [seq[i:] for seq in self.seqs]
    
    def seqDict(self):
        return dict(zip(self.names,self.seqs))
    
    def midPos(self):
        try:
            return int(round(sum(self.positions)/len(self.positions)))
        except:
            pass


#site object class for storing the information about a single site
class Site:
    def __init__(self,scaffold=None, position=None, calls=[]):
        self.scaffold = scaffold
        self.position = position
        self.calls = calls

#function to parse a clls line into the Site class
def parseCallsLine(line):
    site = Site()
    objects = line.split()
    if len(objects) >= 3:
        site.scaffold = objects[0]
        site.position = int(objects[1])
        site.calls = objects[2:]
    return site

#sliding window generator function
def slidingWindows(callsFile, windSize, stepSize, names = None, include = None, exclude = None, skipDeepcopy = False):
    #get file headers
    headers = callsFile.readline().split()
    allNames = headers[2:]
    if not names:
        names = allNames
    columns = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #window counter
    windowsDone = 0
    #initialise an empty window
    window = SeqWindow(names = names, ID = 0)
    #read first line
    line = callsFile.readline()
    site = parseCallsLine(line)
    while line:
        #build window
        while site.scaffold == window.scaffold and site.position <= window.end:
            #add this site to the window
            window.addSite(calls=[site.calls[columns[name]] for name in names], position=site.position)
            #read next line
            line = callsFile.readline()
            site = parseCallsLine(line)
        
        '''if we get here, the line in hand is incompatible with the currrent window
            If the window is not empty, yeild it'''
        
        if window.scaffold is not None:
            windowsDone += 1
            
            if skipDeepcopy: yield window
            else: yield deepcopy(window)
        
        #now we need to make a new window
        #if on same scaffold, just slide along
        if site.scaffold == window.scaffold:
            window.slide(step = stepSize)
            window.ID = windowsDone + 1
        
        #otherwise we're on a new scaffold (or its the end of the file)
        else:
            #if its one we want to analyse, start new window
            if (not include and not exclude) or (include and site.scaffold in include) or (exclude and site.scaffold not in exclude):
                window = SeqWindow(scaffold = site.scaffold, start = 1, end = windSize, names = names, ID = windowsDone + 1)
            
            #if its a scaf we don't want, were going to read lines until we're on one we do want
            else:
                badScaf = site.scaffold
                while site.scaffold == badScaf or (include and site.scaffold not in include and site.scaffold is not None) or (exclude and site.scaffold in exclude and site.scaffold is not None):
                
                    line = callsFile.readline()
                    site = parseCallsLine(line)
            
        #if we've reached the end of the file, break
        if len(line) <= 1:
            break
    

