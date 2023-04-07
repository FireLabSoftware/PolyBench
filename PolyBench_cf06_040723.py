#!/usr/bin/env python -i

## PolyBench version cf05 03-31-23
## PolyBench graphs coverage and variants for a several kb reference
##    k-mers matching a reference sequence are used to count coverage
##    k-mers with a mismatch at the central base are used to count variants
##    k-mers with multiple variants are only counted as representing the central base
## PolyBench is intended to visualize
##    Quality of material from a defined RNA/DNA preparation
##    Quality of data from a defined sequence pipeline
## Optional filters (these can be applied improve fidelity)- 
##    Read-level filters: Filters that ignore potentially errant reads
##       AmbiguityFilter: Any read pair with an N in either read will be ignored
##       DuplicationFilter: Reads that have identical starts for both strands are ignored (deduplication)
##       SmarterStrandedFilter: A specific filter that removes reads that artefactually map to the wrong strand in SmarterStranded seq
##    Kmer-level filters
##       DualStrandRequire: Requires that a mutant k-mer be present in R1 and (in reverse complement) in R2
##       MinQScore: Requires a minimum quality score
##    Only k-mers observed as complements in both sequence reads (R1 and R2) are counted
## Inputs are as follows (command line, Key=Value syntax)
##     RefFile= <FastA file with list of sequences to match k-mers from>
##     Data= <List of .fastq files or .fasta files>
##        Single file, or a list of files that is comma delimited with no spaces
##          .fasta or .fastq files can be gzip compressed, albeit slowing the program 
##          .fasta and .fasta.gz files must have exactly one sequence per line (no multiline sequences)
##        * Wildcards are allowed here, or list of files in a file with the extension .files
##        Providing a directory here will search all files in this directory or subdirectory for fasta and fastq data files
##        Some sets of filenames may produce errors in "guessing" the correspondence of R1 and R2 files
##        If this happens, you can explicitly specify read files with R1File=xxx and R2File=yyy and avoid such errors
##
## Optional Parameters (will default to reasonable values if not set)
##   Input Data Handling
##     R1Only= <default false> Setting this to true tells Polybench to only look at R1 data
##     Trim5/Trim3= <Default is 0> No real need to trim since the sought variants are always in the middle of a k-mer
##         Setting this to a positive integer will trim that number of bases off each end of each read
##     AllowSecondaryMutation <default True> Setting this to true instructs Polybench to count variants where there is up to one additional snp in the k-mer relative to reference
##         (The central base in the k-mer is the counted position, this just allows counting of the k-mer if there is another error somewhere else [default=True])
##     klen = <How long are the k-mers used>.  Default klen = 21.  Recommend an odd number between 21 and 33 depending on complexity of reference  
##     Circular = <Set to True to force every Reference sequence to be treated as a circle>
##     MaxFileReads = <Default is 0 [no maximum]>  Setting this to a fixed value stops reading after a certain number of read pairs from each file
##   Output-
##     Output consists of
##       A graph with metadata and parameters added (generally an svg file [strongly recommended] but can be set to other modes, e.g. png with Graphmode=)
##       A text file (.tdv) with results from the run in a line-by-line table
##       A log file with details of the run
##       Optionally: a fasta file with reads that have unexpected strandedness properties
##     A few output parameters can be set
##       DisplaySmoothingWindow= <Default 100>  This sets the smoothing (averaging) window for the graph lines
##       DisplayGranularity= <Default 100>  This sets the window for the graph (the distance in base pairs between plotted points)
##       OutFileBase = <Default will be files names based on input files>  This can be an optionally user-set character String to Label Output Files with>
##       OutDir= <Default is current working directory>  This sets the output directory
##       GraphTitle= <Character String to Label graph with, will use a default based on input file names if left onset>
##       DisplayVariants= <default True> Setting this to true displays variants that meet the criteria in FractionThreshold and CountThreshold on the graph and reports them in the output text file
##         FractionThreshold= <default 0.0025> Fraction of total matches at a given base that need to be variant to display/report that variant
##            Note that indels in homopolymers have several possible positions; therefor hp length is used to discount frequency before deciding whether to display these variants
##         CountThreshold= <default 4>  Minimum number of variant instances per base required to display/report a variant>
##         LookoutVariants= <default none> Will instruct Polybench to look out for a set of additional variants that then will be displayed even if threshold is not met (standard format for input OriginalPositionMutant (e.g. G234A)
##       ReadRecapture= <default '' [none]>  Setting this captures reads with k-mers with a specified strandedness into a new fasta file, can be 'a' (capture antisense reads' or 's' (capture sense reads), 'b' capture reads with both sense and antisense k-mers or any combination
## Running the program:
##   python PolyBench##.py
##     RefFile=<MyRefFile>
##     Data=MyFastA1.fasta,MyFastA2.fasta.gz,MyFastQ*.fastq
##     <Other_Parameters> 
## ********
## Note about Complex reference sequences
##   This package is not designed to handle reference sequences with substantial internal repeats
##   Sequences with possible repetitive character, are, however flagged and dealt with in a consistent way
##   For any k-mer that is present in the reference as a perfect (unmutated) sequence, the first occurence is treated as the address to assign occurences
##         Note that indel variants in homopolymer run are intrinsically ambiguous in their position.  Polybench reports such variants at central consistent position in any homopolymer run
##         This reflects a need to make some choice to maintain integral counts; note, however, a consequent increase in indel frequencies at the
##           beginning of homopolymer runs due to the fact that these positions "steal" indels from later positions in the homopolymer

## End Help

## Here are default parameters
RefFile1 = '' ## (ref=) (reffile=) reference file for alignments
DataFile0 = [''] ## (data=) (datafiles=) (datafile=) high throughput sequence data (R1 file),  Can be file or list of files
GeneMap1 = 'default' ## (genemap=) A list of start and stop sites for genes within the reference: format [GeneName1, Start1, End1, GeneName2, Start2, End2, ... (0 based)] Default is start-to-end for whole gene fragment
OrfMap1 = 'default' ## (orfmap=) A list of start and stop sites for orfs within each gene [GeneName1,OrfStart1,OrfEnd1,GeneName2,OrfStart2,OrfEnd2, ...]
R1Only1 = False ## (r1only=) Set to true if only one read is provided
klen1 = 21  ## (klen=) length of k-mers used.  Odd numbers are preferred due to avoiding ambiguity of self-complementary k-mers
Trim5 = 0   ## (trim5=) whether to trim sequences from the beginning of each R1 read sequence (end of R2 read)
Trim3 = 0   ## (trim3=) whether to trim sequences from the end of each R1 read sequence (end of R1); Trim to sequence if this is a sequence (TGGAATTCT for TruSeq sRNA)
ReadCutoff1 = 0 ## (readcutoff=) A nonzero value n of ReadCutoff1 will truncate each read to n bases (done on both R1 and R2 and before any trimming)
AllowSecondaryMutation1 = True ## (allowsecondarymutation=) Setting this to true looks allow both primary(central) and secondary (anywhere else) mutations in each k-mer.  This also decreases memory footprint for ,medium/large reference sequences
RequireR1Match1 = False ## (requirer1match=) or (legacy=) Setting this to True will only tally R2 reads if there is a tallied match in R1
AssertSense1 = False ## (assertsense=) Setting this to true only searches for K-mers on the sense strand (ignores any matches to antisense strand)
SpanTolerance1 = 3 ## (spantolerance=) This is the maximum trim distance for any read to be counted in assembling a histogram of fragment sizes.  Set to 3 (default) means that there must be a unique k-mer within the first three bases of the sequence on each side to be counted.
MaxFileReads1 = 0 #100000 #0 ## (maxread=) Maximum number of reads per file (default is to analyze every read MaxFileReads1=0)
Circular1 = False ## (circ=) sets default of whether individual reference sequences are considered linear or circular
DisplaySmoothingWindow1 = 'default' ## (displaysmooth=) Smoothing window for fraction output
DisplayGranularity1 ='default' ## (displaygran=) Granularity for fraction output
CoverageSmoothingWindow1 = 'default' ## (coveragesmooth=) Smoothing window for coverage plot (default will be 'klen1//2')
CoverageGranularity1 ='default' ## (coveragegran=) Granularity for coverage plot (default will be 'klen1//2')
ReportGranularity1 = 200000 ## (reportgran=)(reportgranularity=) How often to report progress to the user (default is once per 200000 reads)
OutputFileBase0 = 'default'  ## (outfilebase=)(outputfilebase=) Setting this will instruct the program as to the basename for files that are produced as output
OutputDirectory1 = ''  ## (outdir=)(outdirectory=) Output will sent to the specified directory.   Default is current working directory
UnusualReadRecapture1 = '' ## (unusualreadrecapture=) (readrecapture=) (readcapture=) (unusualreadcapture=).  Setting this instructs PolyBench to create an output fasta file with reads that have unusual character (s=sense reads from refereence, a=antisense reads, b=both orientations in same read) 
## Setting equalassign=True places equal weight on each position where the read could be assigned.  False (default) assigns to the first possible position on the plus strand
GraphTitle1 = 'default' ## (graphtitle=) Setting this sets the title of the resulting graph (suggestion: avoid commas)
## Setting GraphTitle1 to "VaccineProject" imposes a title on each graph that is related to the 2021-2022 FireLab vacccine sequencing project (very specific to that project)
VScale1 = 600 ## (vscale=) Constant for vertical scaling of graph (pixels per log unit)
HScale1 = 'default' ## (hscale=) Constant for horizontal scaling of graph (pixels per base)
AutoScale1 = True ## (autoscale=) Setting this to False will use a fixed scaling of HScale=1.5, DisplayGranularity=100, DisplaySmoothing=100
Delimiter1 = '\r' ## (delimiter) Delimiter for text output (generally '\r' or '\n')
LookoutVariants0 = [] ## (lookoutvariants=) This is a list of variants to be on the lookout for.  These are inputted at the nucleic acid level as strings <refbase><position><variantbase>, e.g, G2445A is a G->A transition at base 2445
DisplayVariants1 = True ## (displayvariants=) Setting this to true displays variants on the graph and reports them in the output text file
FractionThreshold1 = 0.0025 ## (fthresh=) Fraction of total matches at a given base that need to be variant to display that variant as a box
CountThreshold1 = 4 ## (cthresh=) Minimum number of variant instances required to display a variant as a box
MinYAxisDisplay1 = 0.000001 ## (miny=) Minimum value for display on Y axis (zero will be displayed with this value)
CleanGraph1 = False ## (cleangraph=)  Setting this to true suppresses the output of numerical data on the graphic summary
DetailsOnGraph1 = False ## (cleangraph=)  Setting this to true displays detailed information on run on the graph (lots of extra lines of info)
DisplayGraphic1 = True ## (displaygraphic=) Setting this to True tries to display graphic immediately; False suppresses immediate display of graphic
BaseOffset1 = 0 ## (baseoffset=) Fixed distance to offset plots (setting to 100 will add 100 to all plots on X axis)
GraphStartBase1 = -1 ## (graphstartbase=) Default (-1) only includes positions where there are recovered k-mers.  Setting this to a specific value will expand plotting rightward 
GraphEndBase1 = -1 ## (graphendbase=)  Default (-1) only includes positions where there are recovered k-mers.  Setting this to a specific value will expand plotting leftward
AggregateStarts1 = False ## (aggregatestarts=) Aggregates R1 and R2 starts on the StartEnd plot (otherwise they are shown separately)
GraphicMode1 = 'svg' ## (graphicmode=) Sets mode for graphics file (.svg is recommended but can be png, jpg, tif, eps, or [on some systems] pdf
## Filtering or reads (ie conditions under which entire reads are skipped)
IndexFilter1 = False ## (indexfilter=) Filter by index-- find the most common index for each read file and insist that only these are counted
IndexList0 = 'default' ## (indexlist=) Implements filtering of reads by illumina index (last item on description line of each read.  0 or '' does no filtering, 'default' filters for just the most prevalent index for each file
## a list in the form 'AGAGA+GAGAG,GAGAG+CACAC;GAGAG+ACACA,GGAAG+TTTTT' filters the first file for only those reads with AGAGA+GAGAG or GAGAG+CACAC, and so on
DuplicationFilter1 = False   ## (duplicationfilter=) Setting this to "True" ignores any read pair for which the end nucleotides (collapseK in legth) were already encoutered
DuplicationK1 = 'default'  ## (duplicationk=) Length of each of a pair of k-mers on each side of a read pair to use to detect duplicates (default or no give a value of klen1)
AmbiguityFilter1 = True ## (ambiguityfilter=) Setting this to true skips any read where any base in R1 or R2 is ambiguous (ie "N" in the FastA/FastQ data file)
SmarterStrandedFilter1 = False ## (smarterfilter=) Setting this to true will implement a basic filtering routine for smarter stranded sequencing
## The requirements are
## 1.  R1 must start with three (C or G) bases,
## 2.  The k-mer starting with the fourth base must be in the variant dictionary for the template (at most one mismatch or indel
## 3.  When aligned to that k-mer, the first three bases of the read can have at most one match to the reference
## 4.  Insertion or deletion of 1 base can't align the upstream three bases.
## This ignores some reads that are potentially correct, but does so in a reasonably agnostic manner, while removing many
## wrong stranded reads from the Smarter Stranded protocol.
## ********
## Filtering of K-mers
MinQScore1 = 36 ## (minqscore=) Requires a minimum quality score at the individual base to call a variant (default =36).  Zero for no filter
DualStrandRequire1 = False ## (dual=) (dualstrandrequire=) Setting this to true will require that any k mer that is counted is observed in both R1 and (in cpmplement) in R2 
DissonanceFilter1 = True ## (dissonancefilter=) Setting this to true will ignore rare cases where a k-mer in the same position but different variants is observed in R1 and R2 
## Filtering of reads through insisting that the ends match cleanly (or mismatch)
R1MinTrim5 = -1 ## (r1mintrim5=) Minimum R1 trim 5' side allowed before throwing out read.  Default (-1) doesn't throw out reads. 
R2MinTrim5 = -1 ## (r2mintrim5=) Minimum R1 trim 5' side allowed before throwing out read.  Default (-1) doesn't throw out reads. 

R1MaxTrim5 = -1 ## (r1maxtrim5=) Maximum R1 trim 5' side allowed before throwing out read.  Default (-1) doesn't throw out reads. 
R2MaxTrim5 = -1 ## (r2maxtrim5=) Maximum R1 trim 5' side allowed before throwing out read.  Default (-1) doesn't throw out reads. 
MinSpan0 = 0 ## (minspan=) Minimum span between R1 start and R2 end.  Default (-1) doesn't throw out anything
MaxSpan0 = 0 ## (maxspan=) Maximum span between R1 start and R2 end.  Default (-1) doesn't throw out anything
NubbinCount0 = 3 ## (nubbincount=) How much of a nubbin is counted when cataloging matching (and unmatching) 5' end regions

## Start-End plot This is an optional plot that highlights features reflective of capturable or blocked termini in NGS sequence data (thanks to William Wang and Nimit Jain for recognition of these_
StartEndPlot1 = False ## (startendplot=) Construct a plot with strand bias and start points
spStranded1 = True    ## (spstranded=) True to calculate separate statistics for R1 and R2, consistent with strand-specific input material and library methodology
spCoverage1 = True    ## (spcoverage=) True to show a detailed + and - strand coverage plot below the start site graph
spGraphTitle1='default' ## (spgraphtitle=) Allows the user to specify the start+end graph title (default will construct a title from filenames)
TerminusStrokeWidth1 = 20 ## (terminusstrokewidth=) Width of strokes at 1x magnification on StartEnd plot for ratios
RatioStrokeWidth1 = 10 ## (ratiostrokewidth=) Width of strokes at 1x magnification on StartEnd plot for ratios
CoverageStrokeWidth1 = 3 ## (coveragestrokewidth=) Width of strokes at 1x magnification on StartEnd plot for coverage
## OtherRandomStuff
FastQDumpProgram1 = '' ## (fastqdump=) (fasterqdump=) (fastq-dump=) (fasterq-dump=) Full path of a program that will download files from SRA if needed, or "True" to search hard drive for program.
## An ability to set filenames (or lists) explicitly is added (version bq1).
## This is is only needed for nonstandard filenames.  If your filed are named xxx_R1 and xxx_R2, you shouldn't need to explicitly set F2 names
## But there are cases where this fails, and the abilityu to set R1 and R2 file names explicitly is needed
R1File1 = '' ## (r1file=)(r1files=) Allows user to specify R1 file(s) explicitly.   
R2File1 = '' ## (r2file=)(r2files=) Allows user to specify R2 file(s) explicitly
ReadsToProcess1 = -1 ## (readstoprocess=) Setting this to a positive integer limits the number of reads that are processed

from VSG_ModuleEV import *
from itertools import chain, zip_longest
import statistics, glob, subprocess
from collections import Counter
import gc
import sys
vcommand()

vclear()
vset(bg=white)
vset(svgstrokezoom=True)

LabelsD1 = {('w','f'):['black','black',24,12,'Reference_Base_Fraction'],
           ('t','rc'):['red','red',8,4,'Overall_Coverage'],
           ('c','rc'):['gray80','gray80',40,0,'Expected_Representation'],
           ('m','f'):['magenta','magenta',20,10,'Antisense_Fraction'],
           ('p','f'):['gray','gray',16,8,'Sense_Fraction'],
           ('s','f'):['green','lightgreen',18,8,'Alternate_Base_Fraction'],
           ('i','f'):['brown','yellow',16,8,'Inserted_Base_Fraction'],
           ('d','f'):['blue','cyan',14,8,'Deleted_Base_Fraction']}

## key is <variant type, 'c' for coverage, 'f' for fraction relative to wild type, 'a' is antisense fraction [wt+var],
## value is fill, stroke, radius, Label
           
## a fast genetic code parser.  Everything needs to be upper case and needs to be a real base
## can pre-filter also, or at least make sure of upper case with translate(s.upper()), etc
a1list=[]
for a1 in ['A','G','C','T']:
    for a2 in ['A','G','C','T']:
        for a3 in ['A','G','C','T']:
            a1list.append(a1+a2+a3)
aa='KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF'
gc1={a1list[inde]:aa[inde] for inde in range(64)}
def translate(s):
    return ''.join(gc1[s[i:i+3]] for i in range(0,len(s),3))

SmarterFilterD1 = set(['GGG','GGC','GCG','GCC','CGG','CGC','CCG','CCC'])

## Some tools to look for differential patterns in sense and antisense read sets
SenseTrimR1SeqCounter = Counter()
SenseTrimR2SeqCounter = Counter()
SenseTrimR1LenCounter = Counter()
SenseTrimR2LenCounter = Counter()
SenseStartR1SeqCounter = Counter()
SenseStartR2SeqCounter = Counter()
SenseSpanCounter0 = Counter()
AntiTrimR1SeqCounter = Counter()
AntiTrimR2SeqCounter = Counter()
AntiTrimR1LenCounter = Counter()
AntiTrimR2LenCounter = Counter()
AntiStartR1SeqCounter = Counter()
AntiStartR2SeqCounter = Counter()
AntiSpanCounter0 = Counter()
Starts1 = Counter()
Ends1 = Counter()
kPlus1 = Counter()
kMinus1 = Counter()

OutBase2 = "Default_NoInputFound" ## This will change as soon as Polybench finds an input file
class GraphPath1:
    def __init__(self,gcolor, glinewidth, gname, gpriority):
        self.color = gcolor
        self.linewidth = glinewidth
        self.name = gname
        self.points = []  ## individual x:y points for graphing
        self.data = {}    ## raw data values from which the xy points were generated.
        self.graphed = {}    ## raw data values from which the xy points were generated.
        self.priority = gpriority
    def add(self,p,q=(None,None)): ## p = tuple (xvlaue,yvalue), q can be a raw data which is the basis for this
        self.points.append(p)
        self.data[q[0]]= q[1]
        self.graphed[q[0]] = p[1]
    def draw(self,VSGCanvas):
        if self.points:
            VSGCanvas.vpolyline(points=self.points,strokewidth=self.linewidth,fill=self.color,stroke=self.color,colorkey=self.name,priority=self.priority,svgstrokezoom=True)
        
Counters12 = ["SenseTrimR1SeqCounter","SenseTrimR2SeqCounter","SenseTrimR1LenCounter","SenseTrimR2LenCounter","SenseStartR1SeqCounter","SenseStartR2SeqCounter","SenseSpanCounter0","AntiTrimR1SeqCounter","AntiTrimR2SeqCounter","AntiTrimR1LenCounter","AntiTrimR2LenCounter","AntiStartR1SeqCounter","AntiStartR2SeqCounter","AntiSpanCounter0"]
def find1(n):
    ''' find a file on this machine'''
    o = []
    for r,d,f in os.walk('/'): 
        if n in f:
            o.append(os.path.join(r, n))
    return o

def findfastqdump(candidate):
    ver = 0
    try:
        ver = subprocess.check_output([candidate,'-V'])
        return candidate
    except:
        pass
    if os.path.isfile('~/FastQDumpLocation.txt'):
        candidate = open('~/FastQDumpLocation.txt', mode='rt').read()
        try:
            ver = subprocess.check_output([candidate,'-V'])
            return candidate
        except:
            pass
    vLog('Looking for fasterq-dump-- if this fails, provide a location in the command line (fastqdump=<path>)')
    vLog('or reinstall and allow execution (chmod +X <path>)')
    vLog('Note this can take some real time, get a cup of coffee or tea')
    NewCands = find1('fasterq-dump')
    NewItems = []
    for candidate in NewCands:
        try:
            ver = subprocess.check_output([candidate,'-V'])
            NewItems.append([versioner1(ver),candidate])
        except:
            pass
    if NewItems:
        candidate = sorted(NewItems, reverse=True)[0][1]
        try:
            open('~/FastQDumpLocation.txt', mode='w').write(candidate)
        except:
            pass
        return candidate
    vLog('Unable to find fast-q dump.  Recommend that you reinstall this from ncbi or download fastq/fasta files directly')
    return '' 

kRangeS1 = tuple((2*(klen1-1-x) for x in range(klen1)))
kRangeA1 = tuple((2*x for x in range(klen1)))
kMask1 = 4**(klen1-1)-1
kSam1 = 4**(klen1-1)
kSam3 = 3*kSam1
BaseCS1 = Counter({'G':0,'A':1,'T':2,'C':3,'g':0,'a':1,'t':2,'c':3})
BaseCS2 = Counter({'G':0,'A':kSam1,'T':2*kSam1,'C':3*kSam1,'g':0,'a':kSam1,'t':2*kSam1,'c':3*kSam1})
BaseCA1 = Counter({'C':0,'T':1,'A':2,'G':3,'c':0,'t':1,'a':2,'g':3})
BaseCA2 = Counter({'C':0,'T':kSam1,'A':2*kSam1,'G':3*kSam1,'c':0,'t':kSam1,'a':2*kSam1,'g':3*kSam1})

def sEncode1(sss):  ##last base is in least significant position in numerical encoding
    return sum([BaseCS1[css]<<iss for css,iss in zip(sss,kRangeS1)])
def aEncode1(ass):
    return sum([BaseCA1[css]<<iss for css,iss in zip(ass,kRangeA1)])
    

def CounterReport12(CounterNameList):
    cR = ''
    for c1 in CounterNameList:
        PName = c1.split('Counter')[0]
        PValue = globals()[c1]
        Count1 = sum(PValue.values())
        cR += PName + ' = ' + str(globals()[c1]) + Delimiter1
        if not('seq' in PName.lower()):
            PList = list(chain.from_iterable([[p]*PValue[p] for p in PValue]))
            if Count1>0:
                cR += PName + '_median: ' + '%.1f'%statistics.median(PList) + Delimiter1
                cR += PName + '_mean: ' + '%.4f'%statistics.mean(PList) + Delimiter1
            if Count1>1:
                cR += PName + '_stdev: ' + '%.4f'%statistics.stdev(PList) + Delimiter1
                cR += PName + '_sterr: ' + '%.4f'%(statistics.stdev(PList)/(Count1**0.5)) + Delimiter1
    return cR
        
## CounterReport12(Counters12) is the full report of all the various counters    
def smoothD1(dn, dd, w, g):
    ## dn = dictionary where indices are keys, numerators are values, dd=dict where indices are keys, denominators are values
    ## w is the window and g is the granularity of the output
    ## (list, window [odd]) Smoothed ratio betweeen lists with a window size of w and granularity of g, returns a dictionary where keys are average index, values are running average
    ## ignores any value for which counts of denominator are zero
    if type(dd) in (int,float):
        dd = {x:dd for x in dn}
    d2 = {}
    plist = list(dd.keys())
    pmin = min(plist)
    pmax = max(plist)
    for i in range(pmin,pmax+1,g):
        accu = 0
        num = 0
        for j in range(i,i+w):
            if (j in dd) and dd[j]>0:
                accu += dn.get(j,0)/dd[j]
                num += 1
        if num>0:
            d2[i+w//2] = accu/num
        else:
            d2[i+w//2] = -1
    return d2

def ReportValue1(myvalue):
    if type(myvalue)==int:
        return str(myvalue)+'): '
    if type(myvalue)==bool:
        if myvalue==False:
            return 'Not_Applied): '
        else:
            return 'Applied): '
    else:
        return  str(myvalue)+'): '

       

IndexReport2 = [] ## list of all indices as enforced by the program
    
class variant():
    def __init__(self, sE, reference,sequence,position, original, mutant):
        ## position is the position in the sequence of the base mutation (0 based) in the reference
        ## position of mutated base in the reference is always klen2 (0 based)
        ## There would be klen2 upstream bases (all matching)
        if AllowSecondaryMutation1:
            self.s = sequence
        self.sE = sE ## sEncode value of sequence of k-mer with variant position in 'center' (center-0.5 if klen is even)
        self.p = position ## position of center base (1 based)
        self.o = original ## reference base at center position
        self.m = mutant   ## variant base or bases at center position
        hpU = 0  ## homopolymer length upstream
        hpD = 0  ## homopolymer length downstream
        self.h = 1.0 ## homopolymer number for insertion/deletion variants... 1 for no homopolymer
        if mutant==original:
            self.t = 'w' ## 'wild type'
        elif len(mutant)==1:
            self.t = 's' ## 'snp'
        if mutant=='':
            for temppos in range(klen2-1,0,-1):
                if reference[temppos]!=original: break
                hpU += 1
            for temppos in range(klen2+1,klen1):
                if reference[temppos]!=original: break
                hpD += 1
            if (hpU==hpD) or ((hpU==hpD+1) and (self.o in 'AG')) or ((hpU==hpD-1) and (self.o in 'CT')):
                self.h += hpU+hpD
                self.t = 'd' ## self.t is variant type, d is 'deletion' where the central base is
            else:
                self.t = 'o' ## k-mer has 1b deletion in homopolymer run but not central (so deletion is not assigned here)    
        if len(mutant)==2:
            for temppos in range(klen2,0,-1):
                if reference[temppos]!=mutant[1]: break
                hpU += 1
            for temppos in range(klen2+1,klen1):
                if reference[temppos]!=mutant[1]: break
                hpD += 1
            if (hpU==hpD) or ((hpU==hpD+1) and (self.o in 'AG')) or ((hpU==hpD-1) and (self.o in 'CT')):
                self.h += hpU+hpD
                self.t = 'i' ## self.t is variant type, d is 'deletion' where the central base is
            else:
                self.t = 'e' ## k-mer has 1b insertion in homopolymer run but not central (so insertion is not assigned here)    
        self.d0 = False ## True if there are other (perfect) match "duplicated" positions (in the form (position,strand) where position is zero based and strand is 'p' or 'm'
        self.d1 = False ## True if there are (one-base-off) match "duplicated" positions (in the form (position,strand) where position is zero based and strand is 'p' or 'm'
        self.f = 0 ## count of "forward" instances
        self.b = 0 ## count of "backwards" instances
        self.q = 0 ## a sum of quality scores


if OutputDirectory1 == '':
    OutputDirectory1 = os.path.dirname(OutputFileBase0)
if OutputFileBase0=='default':
    OutputFileBase1=os.path.basename(RefFile1).split('.')[0].split('_')[0]
else:
    OutputFileBase1=os.path.basename(OutputFileBase0)
    

if not(os.path.isfile(RefFile1)):
    if RefFile1=='':
        vLog('\n***Error: No RefFile file specified.  Please run again with a valid FastA reference file specified in command line, e.g. \n\n     RefFile="<YourReferenceFile.fasta>" \n\n')
    else:
        vLog('\n***Error: Unable to locate reference file: '+RefFile1+'\n Please check the file name/path, enclose in quotes (RefFile="<YourReferenceFile.fasta>), and try again')
    exit()

RefName1, S1 = list(vFastAToDict(RefFile1, upper=True).items())[0] ## Uses VSG tool to read reference file to dictionary
vLog('Reading sequence from file '+RefFile1)
vLog('Keeping First Reference Sequence: '+ RefName1)
vLog('Reference Length: '+str(len(S1)))
vLog('Indicated Circularity: '+str(Circular1))
RealBases1 = S1.count('G')+S1.count('A')+S1.count('T')+S1.count('C')
if RealBases1!=len(S1):
    vLog('Nonstandard bases encountered--(GATC:'+str(RealBases1)+ '  U:'+S1.count('U')+ '  N:'+S1.count('N')+ '  OtherBases:'+str(len(S1)-S1.count('N')-S1.count('U')-RealBases1))
    vLog('Converting "U" to "T" and all other nonstandard bases to G to run Polybench')
    S1x = ''
    for cx1 in S1:
        if cx in 'TU':
            s1x += 'T'
        if cx1 in 'AC':
            S1x += cx1
        else:
            S1x  += 'G'
    S1= S1x
    


geneD0 = {}  ## Keys are gene names, values are gene objects
geneD1 = {}  ## Keys are position in S1, values are gene id (or null value if not in a gene)
class gene():
    def __init__(self,name,start,end):
        self.name = name
        self.start = start
        self.end = end
        self.dict = {}  ## keys are positions in S1, values are relative positions in the gene
        if end>start:
            self.seq = S1[start:end+1]
            for i in range(self.start,self.end+1):
                self.dict[i] = i-start
                geneD1[i] = self
        else:
            self.seq = vantisense(S1[end:start+1])
            for i in range(start,end-1,-1):
                self.dict[i] = start-i
                geneD1[i] = self
        self.firstorf() ## default... can be also be set to specific value in the OrfMap variable, which would override this
    def aachange(self,p,newbase):
    ##input: relative position [relative to gene], new base); output amino acid change
        rp = self.dict[p]
        if rp<self.orfstart:
            return("5'UTR")
        if rp>self.orfend:
            return("3'UTR")
        ProteinIndex = (rp-self.orfstart)//3
        WildtypeTriplet = self.seq[self.orfstart+3*ProteinIndex:self.orfstart+3*ProteinIndex+3]
        OldAA = gc1[WildtypeTriplet]
        if len(newbase)==1:
            VariantTriplet = list(WildtypeTriplet)
            VariantTriplet[(rp-3*ProteinIndex)%3]=newbase
            NewAA = gc1[''.join(VariantTriplet)]
        else:
            NewAA = "fs"
        return OldAA+str(ProteinIndex+1)+NewAA
    def firstorf(self):
        stops = set(('TAA','TGA','TAG'))
        seqU = self.seq.upper()
        if not('ATG') in seqU:
            self.orfstart,self.orfend = 0,len(seqU)
        p = self.seq.find('ATG')
        q = p
        while q<=len(self.seq)-3 and not(self.seq[q:q+3] in stops):
            q+=3
        self.orfstart,self.orfend = p,q
if not(GeneMap1) or GeneMap1=='default':
    GeneMap1 = (RefName1,0,len(S1))
for i in range(0,len(GeneMap1),3):
    GeneName1 = GeneMap1[i]
    GeneStart1 = int(GeneMap1[i+1])
    GeneEnd1 = int(GeneMap1[i+2])
    geneD0[GeneName1] = gene(GeneName1,GeneStart1,GeneEnd1)
if OrfMap1 and not(OrfMap1=='default'):
    for i in range(0,len(OrfMap1),3):
        GeneName1 = OrfMap1[i]
        OrfStart1 = int(OrfMap1[i+1])
        OrfEnd1= int(OrfMap1[i+2])
        geneD0[GeneName1].orfstart = OrfStart1
        geneD0[GeneName1].orfend = OrfEnd1  
    
if Circular1: ## Handle circular sequences so that k-mers that span the circle junction are included in analysis
    S2 = S1+S1[:klen1-1]
    S3 = S1*3 ## buffered by bases on each side for finding flanking upstream and downstream segments
else:
    S2 = S1
    S3 = ('N'*len(S1))+S1+('N'*len(S1))
klen2 = klen1//2

LookoutVariants1 = set()  ## not debugged yet [06-19-22]
for lv1 in LookoutVariants0:
    lb1 = '' ## original base
    lp1 = '' ## position of change
    lb2 = '' ## new base
    for lp0,lc1 in enumerate(lv1):
        if lc1 in 'GATC':
            lb1+=lc1
        else:
            break
    for lp0,lc1 in enumerate(lv1[lp0-1:]):
        if lc1.isnum():
            lp1 += lc1
        else:
            break
    for lc1 in enumerate(lv1[lp0:]):
        if lc1 in 'GATC':
            lb2 += lc1
        else:
            break
    lp11 = int(lp1)-1
    varKseq1 = S3[len(S1)+lp11-klen2:len(S1)+lp11]+lb1+S3[len(S1)+lp11+1:len(S1)+lp11+klen1*2][:klen1]
    LookoutVariants1.add(varKseq1)
        
    
    
VariantD1 = {}   ## keys are k-mers, values are 2-ple: (o,v), o=orientation ('p' or 'm'), v=variant objects corresponding to that k-mer.
DMut1 = {'G':['A','T','C','GG','GA','GT','GC',''],
        'A':['G','T','C','AG','AA','AT','AC',''],
        'T':['G','A','C','TG','TA','TT','TC',''],
        'C':['G','A','T','CG','CA','CT','CC','']}
SwitchD1 = {'p':'m', 'm':'p'} ## switch from minus to plus strand and vice versa
SenseR1 = 0 ## Total number of sense reads encountered (sense reads have at least one sense k-mer and no antisense k-mer)
AntiR1 = 0  ## Total number of antisense reads encountered (sense reads have at least one antisense k-mer and no sense k-mer)
BothR1 = 0  ## Total number of ambiguous polarity reads encountered (both sense and antisense k-mers)
CountD1 = {x:Counter() for x in 'wsdipmtc'} ## Keys to each are positions, values are counts in each category
CountD1['o'] = CountD1['w']
CountD1['e'] = CountD1['w']
## Categories here and elsewhere are
     ##  'w': wild type (unmutant)
     ##  's': single nucleotide polymorphism (no indel)
     ##  'd': single nucleotide deletion ('o' are homopolymer positions that are accounted for as wild type since the deletion is considered at the earliest possible position for counting purposes) 
     ##  'i': single nucleotide insertion ('e' are homopolymer positions that are accounted for as wild type since the insertion is considered at the earliest possible position for counting purposes) 
     ##  'p': plus strand
     ##  'm': minus strand
     ##  't': total
     ##  'c': copies

for i1 in range(len(S2)-klen1+1):
    s1 = S2[i1:i1+klen1]
    b0 = S2[i1+klen2]
    sE1 = sEncode1(s1)
    for t1 in CountD1:
        CountD1[t1][(i1+klen2)%len(S1)] = 0
    if not(sE1 in VariantD1):
        CountD1['c'][(i1+klen2)%len(S1)] += 1
        VariantD1[sE1] = ['p',variant(sE1,s1, s1, (i1+klen2)%len(S1), s1[klen2], s1[klen2])]
        if not(AssertSense1):
            aE1 = aEncode1(s1)
            if not(aE1 in VariantD1):
                VariantD1[aE1] = ['m',VariantD1[sE1][1]]
    else:
        CountD1['c'][VariantD1[sE1][1].p] += 1
        VariantD1[sE1][1].d0 = True
for i1 in range(len(S2)-klen1+1):
    if CountD1['c'][(i1+klen2)%len(S1)]==0: continue
    b0 = S2[i1+klen2]
    s1 = S2[i1:i1+klen1]
    sE1 = sEncode1(s1)
    for b1 in DMut1[b0]:
        if len(b1)==1:
            si1 = S2[i1:i1+klen2]+b1+S2[i1+klen2+1:i1+klen1]
        elif b1=='':
            if i1+klen1==len(S2): continue
            si1 = S2[i1:i1+klen2]+S2[i1+klen2+1:i1+klen1+1]
        else:
            si1 = S2[i1:i1+klen2]+b1+S2[i1+klen2+1:i1+klen1-1]
        siE1 = sEncode1(si1)
        if not(siE1 in VariantD1):
            VariantD1[siE1] = ['p',variant(siE1,s1, si1, (i1+klen2)%len(S1), s1[klen2], b1)]
            if not(AssertSense1):
                aiE1 = aEncode1(si1)
                if not(aiE1 in VariantD1):
                    VariantD1[aiE1] = ['m',VariantD1[siE1][1]]
        else:
            VariantD1[siE1][1].d0 = True
VariantD2 = {}
if AllowSecondaryMutation1:    
    for siE1 in VariantD1:
        o1,V1 = VariantD1[siE1]
        if o1=='m': continue
        si1 = V1.s
        del(V1.s)
        b1 = si1[klen2]
        for j1 in range(klen1):
            if j1==klen2: continue
            if j1==klen2+1 and V1.t =='i': continue
            b2 = si1[j1]
            for b3 in 'GATC':
                if b3==b2:continue
                si2 = si1[:j1]+b3+si1[j1+1:]
                siE2 = sEncode1(si2)
                if siE2 in VariantD1:
                    V1.d1 = True
                elif siE2 in VariantD2:
                    VariantD2[siE2][1].d1 = True
                else:
                    VariantD2[siE2] = ['p',V1]
                    if not(AssertSense1):
                        aiE2 = aEncode1(si2)
                        if not(aiE2 in VariantD2):
                            VariantD2[aiE2] = ['m',V1]

## Count unique reference and variant k-mers to give the user a sense of whether their k-mers are long enough
UniqueReferenceCount1 = 0; NonUniqueReferenceCount1 = 0; UniqueVariantCount1 = 0; NonUniqueVariantCount1 = 0
OneBaseOffReferenceCount1 = 0; OneBaseOffVariantCount1 = 0
for s1 in VariantD1:
    o1,V1 = VariantD1[s1]
    if o1=='m': continue
    if V1.t=='w':
        if not(V1.d0):
            UniqueReferenceCount1 += 1
        else:
            NonUniqueReferenceCount1 += 1
        if V1.d1:
            OneBaseOffReferenceCount1 += 1
    else:
        if not(V1.d0):
            UniqueVariantCount1 += 1
        else:
            NonUniqueVariantCount1 += 1
        if V1.d1:
            OneBaseOffVariantCount1 += 1
gc.collect()
vLog('Total Reference Positions: '+str(len(S2)-klen1+1))
vLog('Total Reference k-mers: '+str(UniqueReferenceCount1+NonUniqueReferenceCount1))
vLog('Unique Reference k-mers: '+str(UniqueReferenceCount1))
vLog('Non-Unique Reference k-mers: '+str(NonUniqueReferenceCount1))
vLog('One-Base-Off Reference k-mers: '+str(OneBaseOffReferenceCount1))
vLog('Unique Variant k-mers: '+str(UniqueVariantCount1))
vLog('Non-Unique Variant k-mers: '+str(NonUniqueVariantCount1))
vLog('One-Base-Off Variant k-mers: '+str(OneBaseOffVariantCount1))
AllReads1 = 0
PassIndexFilter1 = 0  ## Read pairs that have the indicated index (all if no index filter is applied)
PassDuplicationFilter1 = 0  ## Read pairs that don't have ends matching a previous read (ie reads remaining after deduplication)
PassAmbiguityFilter1 = 0  ## Read pairs with no N in either read 
PassSmarterFilter1 = 0  ## Read pairs that match the Smarter Stranded strand accuracy filter
ReadsWithSomeAlignment1 = 0  ## Read pairs with at least one k-mer marching the reference
## Note that filters are applied serially, so reads filtered out by earlier filters won't be counted
AllMatchedKmer1 = 0  ## Count of k-mers in all reads
PassQualityFilterKmer1 =0   ## Count of K-mers matching the reference and meeting the minimal quality standard in the central base
PassDualReadFilterKmer1 =0  ## Count of K-mers matching the reference that are detected on both strands for a given read pair (limits detection to the overlap between strands)

DataFile1 = []
if type(DataFile0)==str:
    if '*' in DataFile0:        
        ggdf1 = glob.glob(DataFile0)
        DataFile1.extend(ggdf1)
    else:
        DataFile1.extend( DataFile0.strip().strip('[').strip(']').strip('(').strip(')').split(','))
else:
    DataFile1 = DataFile0
if R1File1 or R2File1:
    R1F1 = []
    for r1f1 in R1File1.split(','):
        R1F1.extend(sorted(glob.glob(r1f1)))
    R2F1 = []
    for r2f1 in R2File1.split(','):
        R2F1.extend(sorted(glob.glob(r2f1)))
    for (r1f1,r2f1) in zip_longest(R1F1,R2F1):
        if not(r1f1): r1f1=''
        if not(r2f1): r2f1=''
        DataFile1.append((r1f1,r2f1))
    if R1File1 and not(R2File1):
        R1Only1 = True

for fn1 in DataFile1:
    if not(os.path.isfile(fn1)) and not(os.path.isfile(fn1+'_1.fastq')) and not('.' in fn1) and FastQDumpProgram1:
        FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
        break
    
if not(DataFile1):
    vLog('******* Error: No Data File Found, check your data filename and put in quotes if needed '+str(DataFile0))
    exit()

def TurnItIntoAFileName1(s):
    s1 = ''
    for s0 in s:
        if s0.isalnum():
            s1 += s0
        elif s0=='*':
            s1 += 'star'
        else:
            if s1 and s1[-1]!='_':
                s1 += '_'
    return s1[:240]

    

for f1, Fn1 in enumerate(DataFile1):
    if not(Fn1): continue
    if type(Fn1)==tuple:
        Fn1,Fn2 = Fn1
        Fn11 = os.path.basename(Fn1)
    else:
        if not(os.path.isfile(Fn1)) and not('.' in Fn1):
            if FastQDumpProgram1 and not(os.path.isfile(Fn1+'_1.fastq')) and not('.' in Fn1):
                vLog(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
                vLog("If this isn't what you intended, a likely possibility is that the filename is not provided or not embedded in quotes-- ")
                vLog("Suggested remedy is to check the file name and rerun with quotations around it (e.g., R1File='MyData_R1.fastq'")
                vLog("But Polybench will try to download from NCBI [type ctrl-c to stop this if this is not what you intended]")
                vLog("Preparing to download sequence read set "+Fn1+" from NCBI")
                TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                    Fn1])
                vLog("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
            if os.path.isfile(Fn1+'.fastq'):
                Fn1 = Fn1+'.fastq'
                R1Only1 = True
            else:
                Fn1 = Fn1+'_1.fastq'
        if not(os.path.isfile(Fn1)):
            vLog('DataFile not found error.  Unable to find a file for for: '+Fn1)
            continue
        Fn11 = os.path.basename(Fn1)
        if ('.2R_' in Fn11[::-1]) and (Fn11[::-1].replace('.2R_','.1R_',1)[::-1] in DataFile1) and not(R1Only1):
            vLog("Redundant R2 Entry in DataFileList [notice only, will process with R1 file]: "+Fn11)
            continue
        if ('.2_' in Fn11[::-1]) and (Fn11[::-1].replace('.2_','.1_',1)[::-1] in DataFile1) and not(R1Only1):
            vLog("Redundant R2 Entry in DataFileList [notice only, will process with R1 file]: "+Fn11)
            continue
        if ('_2R_' in Fn11[::-1]) and (Fn11[::-1].replace('_2R_','_1R_',1)[::-1] in DataFile1) and not(R1Only1):
            vLog("Redundant R2 Entry in DataFileList [notice only, will process with R1 file]: "+Fn11)
            continue
        if ('_2_' in Fn11[::-1]) and (Fn11[::-1].replace('_2_','_1_',1)[::-1] in DataFile1) and not(R1Only1):
            vLog("Redundant R2 Entry in DataFileList [notice only, will process with R1 file]: "+Fn11)
            continue

        Fn2 = Fn1[::-1].replace('.1R_','.2R_',1)[::-1]
        if Fn2==Fn1 or not(os.path.isfile(Fn2)):
            Fn2 = Fn1[::-1].replace('.1_','.2_',1)[::-1]
        if Fn2==Fn1 or not(os.path.isfile(Fn2)):
            Fn2 = Fn1[::-1].replace('_1R_','_2R_',1)[::-1]
        if Fn2==Fn1 or not(os.path.isfile(Fn2)):
            Fn2 = Fn1[::-1].replace('_1_','_2_',1)[::-1]
        if Fn2==Fn1 or not(os.path.isfile(Fn2)):
            Fn2 = Fn1[::-1].replace('_1R','_2R',1)[::-1]
    if Fn2==Fn1 and not(R1Only1):
        vLog('**Warning***: Unable to parse filename and find file for R2 corresponding to your R1 file '+Fn1)
        vLog('**Warning***: Continuing with R1 file only but please note that you can specify R1File=xxx R2File=yyy on command line if needed'+Fn1)
        Fn2 = ''
    if Fn1.lower().endswith('fasta') or Fn1.lower().endswith('fasta.gz') or Fn1.lower().endswith('fa') or Fn1.lower().endswith('fa.gz'):
        LineGran1 = 2
    else:
        LineGran1 = 4
    F2 = iter(())
    if not(os.path.isfile(Fn2)):
        Fn2=''
        R1Only1 = True
    if Fn1.lower().endswith('gz'):       
        F1 = gzip.open(Fn1, mode='rt')
        if Fn2 and not(R1Only1):
            F2 = gzip.open(Fn2, mode='rt')       
    else:
        F1 = open(Fn1, mode='rt')
        if Fn2 and not(R1Only1):
            F2 = open(Fn2, mode='rt')
    ## This is a very specific routine for setting Graph names based on parameters from vaccine work (June 2022).  Could be modified for any specific
    ## project application, or names can be set from the command line
    if OutputFileBase0 == 'default':
        if '_S' in Fn11:
            Mnemonic1 = Fn11[::-1].split('S_',1)[-1][::-1]
        elif ('_' in Fn11) and len(Fn11.split('_')[0])>2:
            Mnemonic1=Fn11.split('_')[0]
        else:
            Mnemonic1=Fn11.split('R1')[0] ## best guess for a non-standard data file name
        if Mnemonic1.endswith('_T') and len(Mnemonic1)>2:
            Mnemonic1 = Mnemonic1[:-2]
        OutputFileBase1 = Mnemonic1+'_kmerCoverage_'+OutputFileBase1+'__'+vnow.replace('_','')
    else:
        Mnemonic1 = OutputFileBase1
    GraphTitle1 = str(GraphTitle1)
    if GraphTitle1.lower()=='vaccineproject':
        GraphTitle1 = 'PolyBench CoverageAndVariants:'
        if "moderna" in (os.path.basename(RefFile1)+RefName1).lower():
            if "dna" in (os.path.basename(RefFile1)+RefName1).lower():
                GraphTitle1 += "_ModernaNIAID_DNA\r"
            else:
                GraphTitle1 += "_ModernaNIAID_RNA\r"
        elif "pfizer" in (os.path.basename(RefFile1)+RefName1).lower():
            if "dna" in (os.path.basename(RefFile1)+RefName1).lower():
                GraphTitle1 += "_BNTPfizer_DNA\r"
            else:
                GraphTitle1 += "_BNTPfizer_RNA\r"
        elif "ms2" in (os.path.basename(RefFile1)+RefName1).lower():
            GraphTitle1 += "_MS2RNA\r"
        elif "phi8" in (os.path.basename(RefFile1)+RefName1).lower().replace('-',''):
            GraphTitle1 += "_phi8kanRNA\r"
        elif "lambda" in (os.path.basename(RefFile1)+RefName1).lower().replace('-',''):
            GraphTitle1 += "_lambdaDNA\r"
        if '*' in str(DataFile0):
            GraphTitle1 += str(DataFile0)
        else:
            GraphTitle1 += ' + '.join([os.path.basename(x).split('_S')[0] for x in DataFile1])
        if R1Only1:
            GraphTitle1 += "  R1Only"
        GraphTitle1 += '\r'
        if "NJA" in str(DataFile1):
            GraphTitle1 += "LigationCapture"
        elif R1Only1:
            GraphTitle1 += "TruSeq"
        else:
            GraphTitle1 += "SmarterStranded"
        if SmarterStrandedFilter1:
            GraphTitle1+= " + SmarterFilter"
        if DualStrandRequire1:
            GraphTitle1+= " + DualFilter"
        if MinQScore1:
            GraphTitle1+= " + Q>"+str(MinQScore1)
    elif GraphTitle1.lower()=='default':
        GraphTitle1 = OutputFileBase1
    OutBase2 = TurnItIntoAFileName1(GraphTitle1)
    GraphicMode1 = GraphicMode1.lower().strip().strip('.')
    UnusualReadRecaptureA1 = 'a' in UnusualReadRecapture1.lower()
    UnusualReadRecaptureS1 = 's' in UnusualReadRecapture1.lower()
    UnusualReadRecaptureB1 = 'b' in UnusualReadRecapture1.lower()
    if UnusualReadRecapture1:
        UnusualReadFileName1 = 'Captured_'+'Sense_'*UnusualReadRecaptureS1+'Anti_'*UnusualReadRecaptureA1+'Both_'*UnusualReadRecaptureB1+OutBase2+'.fasta'
        UnusualReadFile1 = open(os.path.join(OutputDirectory1,UnusualReadFileName1), mode='w')
    if DuplicationFilter1:
        if DuplicationK1 in (0, 'default'):
            DuplicationK1 = klen1
        DuplicationD1 = {}
    if IndexFilter1:
        if IndexList0 == 'default':
            IndexCounter1 = Counter()
            for i0,L000 in enumerate(F1):
                if i0%LineGran1==0:
                    IndexCounter1[L000.strip().split(':')[-1]] += 1
            F1.seek(0)
            IndexReport1 = IndexCounter1.most_common(5)
            vLog('Index Report '+Fn11+': '+str(IndexReport1))
            IndexSet1 = set([IndexReport1[0][0],])
        elif IndexList0:
            IndexSet1 = set(IndexList0.strip().split(';')[f1].split(','))
        else:
            IndexSet1 = set()
        IndexReport2.append(str(IndexReport1))
    RegisterWarned1 = False
    QualityWarned1 = False
    for i1,(L00,M00) in enumerate(zip_longest(F1,F2,fillvalue='')):
        if (i1%(ReportGranularity1*LineGran1)==0):
            vLog('Processing dataset '+Mnemonic1+', Line '+str(i1//LineGran1))
        if (ReadsToProcess1>=0) and (i1>ReadsToProcess1*LineGran1): break
        if MaxFileReads1 and i1>LineGran1*MaxFileReads1: break
        if IndexFilter1:
            if i1%LineGran1==0:
                if L00.strip().split(':')[-1] in IndexSet1:
                    CorrectIndex1 = True
                else:
                    CorrectIndex1 = False
        else:
            CorrectIndex1 = True
        if i1%LineGran1==0 and not(RegisterWarned1):
            if (not(L00.startswith('>')) and not(L00.startswith('@'))) or (not(R1Only1) and not(M00.startswith('>')) and not(M00.startswith('@'))):
                RegisterWarned1 = True
                print('*****')
                print('***** Warning: Nonstandard fasta or fastq data format detected.')
                print('*****   File1 being analyzed=: '+Fn1)
                if not(R1Only1):
                    print('*****   File2 being analyzed=: '+Fn2)                      
                print('*****   Expected ">" or "@" at beginnings of declarations')
                print('*****   Expected a repeating file structure with '+str(LineGran1)+' lines per entry')
                print('*****   This is the entry observed at expected declaration lines # '+str(i1))
                print('*****     R1: '+L00.strip())
                if not R1Only1:
                    print('*****     R2: '+M00.strip())
                print('*****       Only the first break in the expected output format is reported')
                print('*****       so some or all additional entries may also be out of register')
                print('*****')
        if i1%LineGran1==1:
            L0 = L00
            M0 = M00
            if ReadCutoff1:
                L0 = L0[:ReadCutoff1]
                M0 = M0[:ReadCutoff1]
        if i1%LineGran1==LineGran1-1:
            AllReads1 += 1
            if not(CorrectIndex1): continue
            PassIndexFilter1 += 1
            if DuplicationFilter1:
                MyBarcode1 = L00[:DuplicationK1]+M00[:DuplicationK1]
                if MyBarcode1 in DuplicationD1:
                    DuplicationD1[MyBarcode1] += 1
                    continue
                else:
                    DuplicationD1[MyBarcode1] = 1
            PassDuplicationFilter1 += 1
            if AmbiguityFilter1:
                if (('N' in L0) or ('N' in M0)): continue
            PassAmbiguityFilter1 += 1
            SenseK1 = 0  ## perfect match k-mers on sense strand
            AntiK1 = 0    ## perfect match k-mers on antisense strand
            SenseK0 = 0  ## any match k-mers on sense strand
            AntiK0 = 0    ## any match k-mers on antisense strand
            ## Comments 060722  This filters on a variety of indicators that a read is inappropriately assigned to the wrong strand
            ## Correctly strand-assigned reads from smarter strand have 3-4 bases of nontemplated bases (almost always G/C at the 5' end
            ## Inspection of a sample dataset indicated that a 3 base extension was much more common (>97%)
            ## Artefactual reads apparently have some mismatches in the first 3 bases (a fraction have one mismatch in position 1 or 3, but almost none have two mismatces
            ## This filter removes reads on either apparent strand that lack a non-templated 5' extension (criterion >2/3 mismatches in 5' 3 bases)
            ## The next k-mer (starting at the fourth base is required to match with the possible exception of the central base
            ## This routine treats the two strand equally but may remove some correct reads, particularly in highly gc rich regions.  Coverage is not likely to be
            ## Strongly affected, but for prcise coverage estimates it may be beneficial to turn off the filter (which is less relevant for converage estimates                
            if SmarterStrandedFilter1:
                L3 = L0[:3]
                if not(L3 in SmarterFilterD1): continue
                L4 = L0[3:3+klen1]
                if not(len(L4)) == klen1: continue
                sE1 = sEncode1(L4)
                if sE1 in VariantD1:
                    o4,V4 = VariantD1[sE1]
                elif AllowSecondaryMutation1 and (sE1 in VariantD2):
                    o4,V4 = VariantD2[sE1]
                else:
                    continue
                if o4=='p':
                    s3start = len(S1)+V4.p-klen2-3
                elif o4=='m':
                    s3start = len(S1)+V4.p-klen2+klen1
                    L3 = vantisense(L3)
                if (L3==S3[s3start+1:s3start+4]) or (L3==S3[s3start-1:s3start+2]) or (sum([L3[j]==S3[s3start+j] for j in (0,1,2)]))>1: continue
            PassSmarterFilter1 += 1
            Ls0 = L0.strip()
            if Trim5: Ls0 = Ls0[Trim5:]
            if Trim3:
                if type(Trim3)==int:
                    Ls0 = Ls0[:-Trim3]
                else:
                    rpos3 = Ls0.rfind(Trim3)
                    if rpos3>0:
                        Ls0 = Ls0[:rpos3]
                    else:
                        continue
            RDict1 = {}  ## Hits from current read.  Keys are position, values are R1-relative orientation, Variant Object, qualityR1,qualityR2 
            FirstMatchR1 = -1; FirstMatchR2 = -1; LastMatchR1 = -1
            FirstPosR1 = -1; FirstPosR2 = -1; LastPosR1 = -1
            Span0 = None
            sE1 = 0
            for i0,c0 in enumerate(Ls0):
                sE1 = ((sE1&kMask1)<<2) + BaseCS1[c0]
                i = i0-klen1+1
                if i<0: continue
                if LineGran1==4 and len(L00)>Trim5+i+klen2:
                    q1 = ord(L00[Trim5+i+klen2])-33  #L00 is the quality string, Q values are L00 entry-33 as per Illumina standards [error rate=1/(10**(Q/10))]
                else:
                    q1 = 1
                    if not(QualityWarned1):
                        print("####")
                        print("#### Warning: No Quality Scores or Length Mismatch between BaseCall and Quality Strings")
                        print("####   This may reflect a fastA input or compromised or nonstandard fastq input file")
                        print("####   Line Number = "+str(i1))
                        print("####   Input File = "+Fn1)
                        print("####   BaseCalls = "+L0.strip())
                        print("####   Qualities = "+L00.strip())
                        print("####")
                        QualityWarned1 = True
                o1 = ''
                if (sE1 in VariantD1):
                    o1,V1 = VariantD1[sE1]
                    RDict1[V1.p] = [V1,o1,q1,0]
                    if V1.t == 'w' and o1 == 'p' and not(V1.d0):
                        SenseK1 += 1
                    elif V1.t == 'w' and o1 == 'm' and not(V1.d0):
                        AntiK1 += 1
                elif AllowSecondaryMutation1 and (sE1 in VariantD2):
                    o1,V1 = VariantD2[sE1]
                    RDict1[V1.p] = [V1,o1,q1,0]
                elif FirstMatchR1==-1 and i>R1MaxTrim5>=0:
                    break
                else:
                    continue
                if i==0:
                    if o1=='p':
                        Starts1[((V1.p-klen2)%len(S1),o1)] += 1
                    elif o1=='m':
                        Starts1[((V1.p-klen2+klen1-1)%len(S1),o1)] += 1
                if not(V1.d0) and (i+klen1+SpanTolerance1<=len(Ls0)):
                    LastPosR1 = V1.p
                    LastMatchR1 = i+klen1 ## once the loop has completed, this will be the the first mismatched position (zero based)
                if o1 == 'p':
                    SenseK0 += 1
                    efb1 = V1.sE>>(2*klen1-2) ## GATC=3210 for first base of sense strand
                else:
                    AntiK0 += 1
                    efb1 = 3-(V1.sE & 3) ## GATC=3210 for first base of antisense strand
                if (FirstMatchR1==-1) and (sE1>>(2*klen1-2))==efb1:
                    FirstMatchR1 = i
                if (FirstPosR1==-1) and ((sE1>>(2*klen1-2))==efb1) and not(V1.d0) and (i<=SpanTolerance1): ## efb= "Expected first base"
                    FirstPosR1 = V1.p
            if RequireR1Match1 and ((FirstMatchR1==-1) or (FirstMatchR1<R1MinTrim5) or (FirstMatchR1>R1MaxTrim5>=0)): continue
            if (R1Only1 or not(Fn2)):
                if o1=='p':
                    Ends1[((V1.p-klen2+klen1-1)%len(S1),o1)] += 1
                elif o1=='m':
                    Ends1[((V1.p-klen2)%len(S1),o1)] += 1
                if (FirstPosR1>=0) and (LastPosR1>=0):
                    Span0 = abs(FirstPosR1-LastPosR1)+klen1
            else:
                Ms0 = M0.strip()
                if Trim5: Ms0 = Ms0[:-Trim5]
                if Trim3: Ms0 = Ms0[Trim3:]
                aE1 = 0
                for i0,c0 in enumerate(Ms0):
                    aE1 = (aE1>>2) + BaseCA2[c0]
                    i = i0-klen1+1
                    if i<0: continue
                    if LineGran1==4 and len(M00)>Trim3+i+klen2:
                        q2 = ord(M00[Trim3+i+klen2])-33  #M00 is the quality string, Q values are L00 entry-33 as per Illumina standards [error rate=1/(10**(Q/10))]
                    else:
                        q2=1
                    o2 = ''
                    if (aE1 in VariantD1):
                        o2,V2 = VariantD1[aE1]
                        if V2.t == 'w' and o2 == 'p' and not(V2.d0):
                            SenseK1 += 1
                        elif V2.t == 'w' and o2 == 'm' and not(V2.d0):
                            AntiK1 += 1
                    elif AllowSecondaryMutation1 and (aE1 in VariantD2):
                        o2,V2 = VariantD2[aE1]
                    elif FirstMatchR2==-1 and i>R2MaxTrim5>=0:
                        break
                    else:
                        continue
                    if i==0:
                        if o2=='m':
                            Ends1[((V2.p-klen2)%len(S1),'p')] += 1
                        elif o2=='p':
                            Ends1[((V2.p-klen2+klen1-1)%len(S1),'m')] += 1
                    if V2.p in RDict1:
                        if V2==RDict1[V2.p][0]:
                            RDict1[V2.p][3] = q2
                        elif DissonanceFilter1:                            
                            del(RDict1[V2.p])  ## ignore k-mers that are present on both strands but different at the mutational level
                    else:
                        RDict1[V2.p] = [V2,o2,0,q2]
                    if o2 == 'p':
                        SenseK0 += 1
                        efb2 =  3-(V2.sE&3)
                    else:
                        AntiK0 += 1
                        efb2 = V2.sE>>(2*klen1-2)
                    if FirstMatchR2==-1 and 3-(aE1&3)==efb2:
                        FirstMatchR2 = i
                    if FirstPosR2==-1 and 3-(aE1&3)==efb2 and not(V2.d0) and (i<=SpanTolerance1):
                        FirstPosR2 = V2.p
                if (FirstMatchR2==-1) or (FirstMatchR2<R2MinTrim5) or (FirstMatchR2>R2MaxTrim5>=0): continue
                if (FirstPosR1>=0 and FirstPosR2>=0):
                    Span0 = abs(FirstPosR2-FirstPosR1)+klen1
            if (Span0!=None) and ((MaxSpan0 and Span0>MaxSpan0) or (Span0<MinSpan0)):
                continue
            for p1 in RDict1:
                [V1,o1,q1,q2] = RDict1[p1]
                if StartEndPlot1:
                    if not(q1 and q2):
                        if (o1=='p') == bool(q1):
                            kPlus1[V1.p] +=1
                        else:
                            kMinus1[V1.p] +=1
                AllMatchedKmer1 += 1
                if MinQScore1 and (q1+q2<MinQScore1): continue
                PassQualityFilterKmer1 += 1
                if DualStrandRequire1 and (q1*q2==0): continue
                PassDualReadFilterKmer1 += 1
                if o1 == 'p':
                    V1.f += 1
                    V1.q += q1+q2
                elif o1 == 'm':
                    V1.b += 1
                    V1.q += q1+q2
                CountD1[V1.t][V1.p] += 1
                CountD1['t'][V1.p] += 1
                CountD1[o1][V1.p] += 1
            if SenseK0>0 or AntiK0>0:
                ReadsWithSomeAlignment1 += 1
            if SenseK1>0 and AntiK0==0:
                SenseR1+=1
                SenseTrimR1SeqCounter[Ls0[:FirstMatchR1][:NubbinCount0]] += 1
                SenseTrimR1LenCounter[FirstMatchR1] += 1
                SenseStartR1SeqCounter[Ls0[FirstMatchR1:FirstMatchR1+NubbinCount0]] += 1
                if Span0!=None: SenseSpanCounter0[Span0] += 1
                if R1Only1 or not(Fn2):
                    SenseTrimR2SeqCounter[Ls0[LastMatchR1+1:][:NubbinCount0]] += 1
                    SenseTrimR2LenCounter[len(Ls0)-LastMatchR1] += 1
                    SenseStartR2SeqCounter[Ls0[LastMatchR1-NubbinCount0:LastMatchR1]] += 1
                else:
                    SenseTrimR2SeqCounter[Ms0[:FirstMatchR2][:NubbinCount0]] += 1
                    SenseTrimR2LenCounter[FirstMatchR2] += 1
                    SenseStartR2SeqCounter[Ms0[FirstMatchR2:FirstMatchR2+NubbinCount0]] += 1
                if UnusualReadRecaptureS1:
                    UnusualReadFile1.write('>SenseMatch_'+str(SenseR1)+'_'+Fn11+Delimiter1)
                    UnusualReadFile1.write(L0.strip()+Delimiter1)
                    if not(R1Only1) and Fn2:
                        UnusualReadFile1.write(M0.strip()+Delimiter1)
            elif SenseK0==0 and AntiK1>0:
                AntiR1+=1
                AntiTrimR1SeqCounter[Ls0[:FirstMatchR1][:NubbinCount0]] += 1
                AntiTrimR1LenCounter[FirstMatchR1] += 1
                AntiStartR1SeqCounter[Ls0[FirstMatchR1:FirstMatchR1+NubbinCount0]] += 1
                if Span0!=None: AntiSpanCounter0[Span0] += 1
                if R1Only1 or not(Fn2):
                    AntiTrimR2SeqCounter[Ls0[LastMatchR1+1:][:NubbinCount0]] += 1
                    AntiTrimR2LenCounter[len(Ls0)-LastMatchR1] += 1
                    AntiStartR2SeqCounter[Ls0[LastMatchR1-NubbinCount0:LastMatchR1]] += 1
                else:
                    AntiTrimR2SeqCounter[Ms0[:FirstMatchR2][:NubbinCount0]] += 1
                    AntiTrimR2LenCounter[FirstMatchR2] += 1
                    AntiStartR2SeqCounter[Ms0[FirstMatchR2:FirstMatchR2+NubbinCount0]] += 1
                if UnusualReadRecaptureA1:
                    UnusualReadFile1.write('>AntisenseMatch_'+str(AntiR1)+'_'+Fn11+Delimiter1)
                    UnusualReadFile1.write(L0.strip()+Delimiter1)
                    if not(R1Only1) and Fn2:
                        UnusualReadFile1.write(M0.strip()+Delimiter1)
            elif SenseK0>0 and AntiK0>0:
                BothR1+=1
                if UnusualReadRecaptureB1:
                    UnusualReadFile1.write('>BothStrandMatch_'+str(BothR1)+'_'+Fn11+Delimiter1)
                    UnusualReadFile1.write(L0.strip()+Delimiter1)
                    if not(R1Only1) and Fn2:
                        UnusualReadFile1.write(M0.strip()+Delimiter1)
    F1.close()
    if Fn2 and not(R1Only1):
        F2.close()
    if UnusualReadRecapture1:
        UnusualReadFile1.close()
    vLog('Finishing file '+Fn11+'  High Confidence Sense Reads: '+str(SenseR1)+'  High Confidence Antisense Reads: '+str(AntiR1)+'  Dual Strand:'+str(BothR1))

TotalsD1 = {x:sum(CountD1[x].values()) for x in CountD1} 

GSB1 = max(0,GraphStartBase1)
GEB1 = len(S1)
if GraphEndBase1>0: GEB1 = GraphEndBase1

if HScale1 == 'default':
    if not(AutoScale1):
        HScale1 = 1.5
    else:
        HScale1 = min(15,15000/(GEB1-GSB1))
if DisplayGranularity1 == 'default':
    if not(AutoScale1):
        DisplayGranularity1 = 100
    else:
        DisplayGranularity1 = max(21,(GEB1-GSB1)//50)
if DisplaySmoothingWindow1 == 'default':
    if not(AutoScale1):
        DisplaySmoothingWindow1 = 100
    else:
        DisplaySmoothingWindow1 = DisplayGranularity1

VariantList1 = []
if DisplayVariants1:
    for s1 in VariantD1:
        o1,V1 = VariantD1[s1]
        if o1=='m' or V1.t in 'woe': continue
        if (s1 in LookoutVariants1) or (CountD1['t'][V1.p] and ((V1.f+V1.b)/(CountD1['t'][V1.p])>(V1.h*FractionThreshold1)) and (V1.f+V1.b)>CountThreshold1):
            v1stroke, v1fill, v1rad, v1width, v1TypeLabel = LabelsD1[(V1.t,'f')]
            v1SpotLabel1 = V1.o+str(V1.p+1)
            if V1.t == 's' or V1.t=='i':
                v1SpotLabel1 += V1.m
            elif V1.t == 'd':
                v1SpotLabel1 += 'del'
            if V1.p in geneD1:
                v1SpotLabel1 += '('+geneD1[V1.p].aachange(V1.p,V1.m)+')'
            meanq12 = 1.0
            if V1.f+V1.b>0:
                meanq12 = V1.q/(V1.f+V1.b)
            v1SpotLabel1 += "Q%.0f"%meanq12
            VariantList1.append(v1SpotLabel1+'@'+'%.0f'%(V1.f+V1.b)+'/'+'%.0f'%CountD1['t'][V1.p])
            vsquare(xc=V1.p*HScale1,yc=log10((V1.f+V1.b)/(CountD1['t'][V1.p]))*VScale1,r=v1rad, fill=v1fill,stroke=v1stroke,strokewidth=v1width,label=v1SpotLabel1, labelfont="DejaVuSerif 36 Bold", priority=2)
if GraphStartBase1>=0:
    vcircle(xc=(GraphStartBase1-BaseOffset1)*HScale1,yc=0,r=0,xg=GraphStartBase1)
if GraphEndBase1>=0:
    vcircle(xc=(GraphEndBase1-BaseOffset1)*HScale1,yc=0,r=0,xg=GraphEndBase1)

TotalsD1['t'] = max(TotalsD1['t'],MinYAxisDisplay1)  # Avoid later divide by zero

if CoverageSmoothingWindow1=='default':
    CoverageSmoothingWindow1 = 1 ## klen2
if CoverageGranularity1=='default':
    CoverageGranularity1 = 1 ## klen2


for t1 in ('w','m','p','s','i','d'):
    if t1=='a' and AssertSense1: continue
    v1stroke, v1fill, v1rad, v1width, v1TypeLabel = LabelsD1[(t1,'f')]
    smoothV1 = smoothD1(CountD1[t1],CountD1['t'],DisplaySmoothingWindow1, DisplayGranularity1)
    lastp1 = -1
    lastv1 = -1
    for p1 in smoothV1:
        v1 = smoothV1[p1]
        if v1<0:
            lastv1 = -1
        else:
            v1 = max(MinYAxisDisplay1,v1)
            vcircle(xc=p1*HScale1,yc=log10(v1)*VScale1,r=v1rad, fill=v1fill,stroke=v1stroke,strokewidth=v1width,xg=p1+BaseOffset1+1,yg=v1, colorkey=v1TypeLabel, priority=2)
            if (lastp1>=0 and lastv1>=0):
                vline(x1=lastp1*HScale1, x2=p1*HScale1, y1=log10(lastv1)*VScale1, y2=log10(v1)*VScale1, color=v1stroke,strokewidth=v1width, priority=1)
            lastp1 = p1
            lastv1 = v1


v1stroke, v1fill, v1rad, v1width, v1TypeLabel = LabelsD1[('c','rc')]
lastv1 = -1
lastp1 = -1
p1max = max(CountD1['c'].keys())
p1min = min(CountD1['c'].keys())
for p1 in CountD1['c']:
    v1 = CountD1['c'][p1]
    if v1!=lastv1 or p1==p1max:
        if p1>p1min:
            lastv2 = max(MinYAxisDisplay1,lastv1)
            vrect(x1=lastp1*HScale1, x2=(p1-1)*HScale1, yc=log10(lastv2)*VScale1, yr=v1rad, fill=v1fill, stroke=v1stroke, colorkey=v1TypeLabel, priority=0)
        lastp1 = p1
        lastv1 = v1

v1stroke, v1fill, v1rad, v1width, v1TypeLabel = LabelsD1[('t','rc')]
smoothV1 = smoothD1(CountD1['t'],TotalsD1['t']/len(CountD1[t1]),CoverageSmoothingWindow1, CoverageSmoothingWindow1)
smoothC1 = smoothD1(CountD1['c'],1.0,CoverageSmoothingWindow1, CoverageSmoothingWindow1)
for p1 in smoothV1:
    v1 = smoothV1[p1]
    c1 = smoothC1[p1]
    v1 = max(MinYAxisDisplay1,v1)
    if c1>0:
        vcircle(xc=p1*HScale1,yc=log10(v1)*VScale1,r=v1rad, fill=v1fill,stroke=v1stroke,strokewidth=v1width,xg=p1+BaseOffset1+1,yg=v1, colorkey=v1TypeLabel, priority=2)

##v1stroke, v1fill, v1rad, v1width, v1TypeLabel = LabelsD1[('t','rc')]
##for p1 in CountD1['t']:    
##    if TotalsD1['t']:
##        v1 = CountD1['t'][p1]*(UniqueReferenceCount1+NonUniqueReferenceCount1)/TotalsD1['t']
##    else:
##        v1 = 0
##    v1 = max(v1,MinYAxisDisplay1)
##    vcircle(xc=p1*HScale1,yc=log10(v1)*VScale1,r=v1rad, fill=v1fill,stroke=v1stroke,strokewidth=v1width,xg=p1,yg=v1,colorkey=v1TypeLabel, priority=2)
##v1stroke, v1fill, v1rad, v1width, v1TypeLabel = LabelsD1[('c','rc')]
##for p1 in CountD1['c']:    
##    v1 = max(CountD1['c'][p1],MinYAxisDisplay1/1.2)
##    vcircle(xc=p1*HScale1,yc=log10(v1)*VScale1,r=v1rad, fill=v1fill,stroke=v1stroke,strokewidth=v1width,xg=p1,yg=v1,colorkey=v1TypeLabel, priority=4)
   


vgrid(gxlabel="Position", gmajorcolor=cyan, gminorcolor=cyan, gylog=True, gylabel="Relative Coverage", gtitle=GraphTitle1, gtitlefontsize=270)

## Draw a very simple map
GeneWidth1 = 40
yGene1 = VSG.ymin ## Vertical position of map
MapFontSize1 = 192
for g0 in geneD1:
    g1 = geneD1[g0]
    myName = g1.name.replace('_',' ')
    if g1.start<g1.end:
        vtext(text=myName+'->', x1=(g1.orfstart+g1.start)*HScale1,y1=yGene1-MapFontSize1//2, font="DejaVuSerif "+str(MapFontSize1)+" Bold", color='gray40')
        vrect(x1=g1.start*HScale1, x2=g1.end*HScale1, y1=yGene1-MapFontSize1//2-2*GeneWidth1, y2=yGene1-MapFontSize1//2-4*GeneWidth1, fill='gray40', stroke=none, strokewidth=0)
        if g1.orfstart>-1:
            vrect(x1=(g1.orfstart+g1.start)*HScale1, x2=(g1.orfend+g1.start)*HScale1, y1=yGene1-MapFontSize1//2-GeneWidth1, y2=yGene1-MapFontSize1//2-5*GeneWidth1, fill='gray40', stroke=none, strokewidth=0)            
    else:
        vtext(text='<-'+myName, x2=(g1.start-g1.orfstart)*HScale1,y2=yGene1-6*GeneWidth1, font="DejaVuSerif "+str(MapFontSize1)+" Bold", color='gray60')
        vrect(x1=g1.end*HScale1, x2=g1.start*HScale1, y1=yGene1-1*GeneWidth1, y2=yGene1-3*GeneWidth1, fill='gray60', stroke=none, strokewidth=0)
        if g1.orfstart>-1:
            vrect(x1=(g1.start-g1.orfend)*HScale1, x2=(g1.start-g1.orfstart)*HScale1, y1=yGene1, y2=yGene1-4*GeneWidth1, fill='gray60', stroke=none, strokewidth=0)            
            

TotalR1 = SenseR1+AntiR1 ##+BothR1 (For now we report fraction of cleearly callable reads)
vcolorkey()
tD1 = max(1,TotalsD1['t'])
tR1 = max(1,TotalR1)
if not(CleanGraph1):
    if AssertSense1:
        vlegend(text="\rSense reads: "+str(SenseR1)+  "  (%0.4f%%)"%(SenseR1*100.0/tR1)+'\r'+
        "(AssertSense is ON so no Antisense k-mers searched)", font="DejaVu Serif 192 Bold")
    else:
        vlegend(text="\rSense reads: "+str(SenseR1)+  "  (%0.4f%%)"%(SenseR1*100.0/tR1)+'\r'+
        "Antisense reads: "+str(AntiR1)+  "  (%0.4f%%)"%(AntiR1*100.0/tR1)+
            "\r (Strand-ambiguous reads: "+str(BothR1)+')',
            font="DejaVu Serif 192 Bold")
    vlegend(text="\rReference-bases: %0.0f"%TotalsD1['w'] + "  (%0.4f%%)"%(TotalsD1['w']*100.0/tD1)+'\r'
        + "Substitution-base-variants: %0.0f"%TotalsD1['s'] + "  (%0.4f%%)"%(TotalsD1['s']*100.0/tD1)+'\r'
        + "Deleted-base-variants: %0.0f"%TotalsD1['d'] + "  (%0.4f%%)"%(TotalsD1['d']*100.0/tD1)+'\r'
        + "Inserted-base-variants: %0.0f"%TotalsD1['i'] + "  (%0.4f%%)"%(TotalsD1['i']*100.0/tD1), font="DejaVu Serif 192 Bold")
if DetailsOnGraph1:
    vlegend(text=vSysLogInfo1, font="DejaVuSerif 48 Bold")

details1 = vSysLogInfo1.split("Full External Variable List (after substitutions from vCommand):")

DataFile0 = []
for df11 in DataFile1:
    if type(df11)==str:
        DataFile0.append(df11)
    else:
        DataFile0.extend(df11)
details2 = ''
for dL in details1[0].splitlines():
    details2 += dL.strip().replace('\t',' ').replace('\r',' ').replace('\n',' ')+';'
details2+= Delimiter1+'Variable_List'+Delimiter1
for dL in details1[-1].strip().splitlines():
    vname1 = dL.split('=')[0].strip()
    details2 += vname1+' = '+str(globals()[vname1]).replace('\r','\\r')+Delimiter1
details2 += "AllFilesRead = "+';'.join(DataFile0)+Delimiter1
details2 += "Indices_Used_For_PolyBench_Precise_Demultiplexing = "+';'.join(IndexReport2)
report0 = "##PolyBench_Output_at_TimeStamp: "+vnow+Delimiter1
report0 += "Program_Version: "+os.path.basename(sys.argv[0])+Delimiter1
report0 += "Reference_Sequence: "+os.path.basename(RefFile1)+Delimiter1
report0 += "Data_Files: "+str([os.path.basename(DataFile11) for DataFile11 in DataFile0])+Delimiter1
report0 += "Command_Line_Settings: "+' '.join(sys.argv[1:])+Delimiter1
report0 += "Total_Reads_Processed: "+str(AllReads1)+Delimiter1
report0 += "Reads_Passing_Index_Filter ("+ReportValue1(IndexFilter1)+str(PassIndexFilter1)+Delimiter1
report0 += "Reads_Passing_Duplication_Filter ("+ReportValue1(DuplicationFilter1)+str(PassDuplicationFilter1)+Delimiter1
report0 += "Reads_Passing_Ambiguity_Filter ("+ReportValue1(AmbiguityFilter1)+str(PassAmbiguityFilter1)+Delimiter1
report0 += "Reads_Passing_Smarter_Stranded_Stringent_Filters ("+ReportValue1(SmarterStrandedFilter1)+str(PassSmarterFilter1)+Delimiter1
report0 += "Reads_With_Some_Alignment: "+str(ReadsWithSomeAlignment1)+Delimiter1
report0 += "Matched_Kmers: "+str(AllMatchedKmer1)+Delimiter1
report0 += "Kmers_Passing_Quality_Filter (Q>"+ReportValue1(MinQScore1)+str(PassQualityFilterKmer1)+Delimiter1
report0 += "Kmers_Passing_Dual_Read_Filter ("+ReportValue1(DualStrandRequire1)+str(PassDualReadFilterKmer1)+Delimiter1

report1 = Delimiter1.join(["Sense_Reads: "+str(SenseR1),
          "Sense_Read_Percentage: "+"%0.4f"%(SenseR1*100.0/tR1),
          "Antisense_Reads: "+str(AntiR1),
          "AntiSense_Read_Percentage: "+'%0.4f'%(AntiR1*100.0/tR1),
          "Mixed_Strand_Reads: "+str(BothR1),
          'Reference_Bases: '+"%0.0f"%TotalsD1['w'],
          'Reference_Base_Percentage: '+'%0.4f'%(TotalsD1['w']*100.0/tD1),
           'Substitution_Base_Variants: '+"%0.0f"%TotalsD1['s'],
           'Substitution_Base_Percentage: '+"%0.4f"%(TotalsD1['s']*100.0/tD1),
           "Deleted-base-variants: "+"%0.0f"%TotalsD1['d'],
           'Deleted_Base_Percentage: '+"%0.4f"%(TotalsD1['d']*100.0/tD1),
           "Inserted_base_variants: "+"%0.0f"%TotalsD1['i'],
           'Inserted_Base_Percentage: '+"%0.4f"%(TotalsD1['i']*100.0/tD1),
           'Prevalent_Variant_List: '+','.join(VariantList1)])+Delimiter1
ReportFile1 = open(OutBase2+"_"+vnow+"pyb.tdv", mode='w')
ReportFile1.write(report0+report1+CounterReport12(Counters12)+details2)
ReportFile1.close()
gc.collect()
if StartEndPlot1:
    from BayesRatioAH import ConservativeFoldDifference  ## (PExpt,PRef,TExpt,TRef,FDR=0.05,Sp=2)

    SEPlot1 = VSGcanvas()
    SEPlot2 = VSGcanvas()
    SEPlot1.vset(bg=white)
    SEPlot2.vset(bg=white)
    spHScale1 = HScale1
    spVScale1 = VScale1
    spVScale3 = spVScale1/2 ## Vertical scale for coverage plots
    spYAxisStartToBiasConversion1 = 50 ## This sets the graphing equivalence between Y axis units of log2 (for differential strand presence) and enrichment (for start sites)
    spVScale2 = VScale1/spYAxisStartToBiasConversion1
    spLabelFontSize1 = 120
    VAdjust1 = spLabelFontSize1/(spVScale2*1.5)
    MinLog2VScaleSpan1 = 6 ##Sets a minimum span on the Y axis
    MinLog2VScaleSpan2 = 4 ##Sets a minimum span on the Y axis (CoveragePlots)
    spFlags1 = {} ## keys are position (zero based) values are  ##pStart (R1 start on positive strand, pEnd (R2 start on positive strand), mStart (R1 start on negative strand), mEnd (R2 start on negative strand)
    spLabelFont1='DejaVuSansBold '+str(spLabelFontSize1)
    spSequenceFont1='DejaVuSansMono Bold '+str(spLabelFontSize1//2)
    kT1 = max(0.5,sum(kPlus1.values())+sum(kMinus1.values()))    
    MinFoldEnrichmentForEndLabel1 = 10
    spN1 = 2*(UniqueReferenceCount1+NonUniqueReferenceCount1)
    StartsT1 = max(1,sum(Starts1.values()))
    EndsT1 = max(1,sum(Ends1.values()))
    if spGraphTitle1=='default':
        spGraphTitle1 = 'Starts-and-Ends:  '+'ref='+os.path.basename(RefFile1).split('.')[0].split('_')[0]+'  data='+Mnemonic1
    GraphPathD1 = {'+ Termini':GraphPath1('magenta1',TerminusStrokeWidth1,'+ Termini',3),
                   '- Termini':GraphPath1('purple1',TerminusStrokeWidth1,'- Termini',3),
                   '+ Termini-R1':GraphPath1('magenta1',TerminusStrokeWidth1,'+ Termini-R1',3),
                   '- Termini-R1':GraphPath1('purple1',TerminusStrokeWidth1,'- Termini-R1',3),
                   '+ Termini-R2':GraphPath1('magenta4',TerminusStrokeWidth1,'+ Termini-R2',2),
                   '- Termini-R2':GraphPath1('purple4',TerminusStrokeWidth1,'- Termini-R2',2),                   
                   '+/- Coverage Ratio':GraphPath1('red',RatioStrokeWidth1,'+/- Coverage Ratio',1),
                   '(R1s+R2a)/(R1a+R2s) Coverage Ratio':GraphPath1('orange',RatioStrokeWidth1,'(R1s+R2a)/(R1a+R2s) Coverage Ratio',1),
                   '+ Coverage':GraphPath1('blue',CoverageStrokeWidth1,'+ Coverage',0),
                   '- Coverage':GraphPath1('green',CoverageStrokeWidth1,'- Coverage',0)}
                   
    for spi in range(len(S1)):
        if spStranded1 and (CountD1['p'][spi]+CountD1['m'][spi]):
            d1 = ConservativeFoldDifference(CountD1['p'][spi],
                                CountD1['m'][spi],
                                CountD1['p'][spi]+CountD1['m'][spi],
                                CountD1['p'][spi]+CountD1['m'][spi],
                                1/(len(S1)+1))
            d1 = min(max(-MinLog2VScaleSpan1,log(d1,2)),MinLog2VScaleSpan1)
            GraphPathD1['(R1s+R2a)/(R1a+R2s) Coverage Ratio'].add((spi*spHScale1,d1*spVScale1),(spi,(CountD1['p'][spi],CountD1['m'][spi])))

        if kPlus1[spi]+kMinus1[spi]:
            d2 = ConservativeFoldDifference(kPlus1[spi],
                                            kMinus1[spi],
                                            kPlus1[spi]+kMinus1[spi],
                                            kPlus1[spi]+kMinus1[spi],
                                            1/(len(S1)+1))
            d2 = min(max(-MinLog2VScaleSpan1,log(d2,2)),MinLog2VScaleSpan1)
            GraphPathD1['+/- Coverage Ratio'].add((spi*spHScale1,d2*spVScale1),(spi,(kPlus1[spi],kMinus1[spi])))
        if spCoverage1:
            if kMinus1[spi]:
                covm1 = min(max(-MinLog2VScaleSpan2,log((spN1*kMinus1[spi])/kT1,2)),MinLog2VScaleSpan2)
                GraphPathD1['- Coverage'].add((spi*spHScale1,covm1*spVScale3),(spi,kMinus1[spi]))
            if kPlus1[spi]:
                covp1 = min(max(-MinLog2VScaleSpan2,log((spN1*kPlus1[spi])/kT1,2)),MinLog2VScaleSpan2)
                GraphPathD1['+ Coverage'].add((spi*spHScale1,covp1*spVScale3),(spi,kPlus1[spi]))                  
        spLabel1 = str(spi+BaseOffset1+1)
        if spStranded1:
            pStart1 = spN1*Starts1[(spi,'p')]/StartsT1
            mStart1 = spN1*Starts1[(spi,'m')]/StartsT1
            pEnd1 = spN1*Ends1[(spi,'p')]/EndsT1
            mEnd1 = spN1*Ends1[(spi,'m')]/EndsT1
            pStart2= min(pStart1,spYAxisStartToBiasConversion1*MinLog2VScaleSpan1)
            mStart2= min(mStart1,spYAxisStartToBiasConversion1*MinLog2VScaleSpan1)
            pEnd2= min(pEnd1,spYAxisStartToBiasConversion1*MinLog2VScaleSpan1)
            mEnd2= min(mEnd1,spYAxisStartToBiasConversion1*MinLog2VScaleSpan1)
            GraphPathD1['+ Termini-R1'].add((spi*spHScale1,pStart2*spVScale2),(spi,Starts1[(spi,'p')]))     
            GraphPathD1['- Termini-R1'].add((spi*spHScale1,-mStart2*spVScale2),(spi,Starts1[(spi,'m')]))     
            GraphPathD1['+ Termini-R2'].add((spi*spHScale1,pEnd2*spVScale2),(spi,Ends1[(spi,'p')]))     
            GraphPathD1['- Termini-R2'].add((spi*spHScale1,-mEnd2*spVScale2),(spi,Ends1[(spi,'m')]))     
        else:
            pStart1 = spN1*(Starts1[(spi,'p')]+Ends1[(spi,'p')])/(StartsT1+EndsT1)
            mStart1 = spN1*(Starts1[(spi,'m')]+Ends1[(spi,'m')])/(StartsT1+EndsT1)
            pStart2= min(pStart1,-3+spYAxisStartToBiasConversion1*MinLog2VScaleSpan1)
            mStart2= min(mStart1,-3+spYAxisStartToBiasConversion1*MinLog2VScaleSpan1)
            mAdjust1 = 0.0; pAdjust1 = 0.0
            GraphPathD1['+ Termini'].add((spi*spHScale1,pStart2*spVScale2),(spi,Starts1[(spi,'p')]+Ends1[(spi,'p')]))     
            GraphPathD1['- Termini'].add((spi*spHScale1,-mStart2*spVScale2),(spi,Starts1[(spi,'m')]+Ends1[(spi,'m')]))     
        if pStart1>MinFoldEnrichmentForEndLabel1:            
            if not spi in spFlags1:
                spFlags1[spi] = []
            spFlags1[spi].append('pStart_Flag')
            chaser1 = ''+(' [%.1fx]'%pStart1)*(pStart1>pStart2)
            if spStranded1:
                chaser1='(r1)'+chaser1
                if 0<=pStart2-pEnd2<spLabelFontSize1/spVScale2/1.5:
                    pStart2+=VAdjust1-(pStart2-pEnd2)
            SEPlot1.vtext(x1=spi*spHScale1,
                      yc=pStart2*spVScale2,
                      font=spLabelFont1,
                      color='magenta1',
                      text='->'+spLabel1+chaser1,
                      priority=4)
        if mStart1>MinFoldEnrichmentForEndLabel1:
            if not spi in spFlags1:
                spFlags1[spi] = []
            spFlags1[spi].append('mStart_Flag')
            chaser1 = ''+(' [%.1fx]'%mStart1)*(mStart1>mStart2)
            if spStranded1:
                chaser1='(r1)'+chaser1
                if 0<=mStart2-mEnd2<spLabelFontSize1/spVScale2/1.5:
                    mStart2+=VAdjust1-(mStart2-mEnd2)
            SEPlot1.vtext(x2=spi*spHScale1,
                      yc=-mStart2*spVScale2,
                      font=spLabelFont1,
                      color='purple1',
                      text=spLabel1+chaser1+'<-',
                      priority=4)
        if spStranded1 and pEnd1>MinFoldEnrichmentForEndLabel1:
            if not spi in spFlags1:
                spFlags1[spi] = []
            spFlags1[spi].append('pEnd_Flag')
            chaser1 = ''+(' [%.1fx]'%pEnd1)*(pEnd1>pEnd2)
            if 0<pEnd2-pStart2<spLabelFontSize1/spVScale2/1.5:
                pEnd2+=VAdjust1-(pEnd2-pStart2)
            SEPlot1.vtext(x1=spi*spHScale1,
                          yc=pEnd2*spVScale2,
                          font=spLabelFont1,
                          color='magenta4',
                          text='->'+spLabel1+'(r2)'+chaser1,
                          priority=4)
        if spStranded1 and mEnd1>MinFoldEnrichmentForEndLabel1:
            if not spi in spFlags1:
                spFlags1[spi] = []
            spFlags1[spi].append('mEnd_Flag')
            chaser1 = ''+(' [%.1fx]'%mEnd1)*(mEnd1>mEnd2)
            if 0<mEnd2-mStart2<spLabelFontSize1/spVScale2/1.5:
                mEnd2+=VAdjust1-(mEnd2-mStart2)
            SEPlot1.vtext(x2=spi*spHScale1,
                          yc=-mEnd2*spVScale2,
                          font=spLabelFont1,
                          color='purple4',
                          text=spLabel1+'(r2)'+chaser1+'<-',
                          priority=4)
    if GraphEndBase1<spi+BaseOffset1:
        GraphEndBase1 = spi+BaseOffset1
    for spg1 in GraphPathD1:
        if spg1 in ('+ Coverage','- Coverage'):
            if spCoverage1:
                GraphPathD1[spg1].draw(SEPlot2)
        else:
            GraphPathD1[spg1].draw(SEPlot1)
    ## Draw a grid for log values (automatic from VSG
    SEPlot1.vgrid(gtitle=spGraphTitle1,gtitlefontsize=400,gxlabel='Reference Position'*(not(spCoverage1)), gylabel='MinEst-log2(coverage-ratio)',gmajorwidth=10,gmajorcolor='cyan',gylabelcolor='red',
                                  gminorwidth=5,gminorcolor='cyan',gxminor=True,gyminor=True,
                                  gyaxisfontsize=300,gylabelfontsize=300,gxaxisfontsize=300,gxlabelfontsize=300,
                                  gxmin=(GraphStartBase1,(GraphStartBase1-BaseOffset1)*spHScale1),
                                  gxmax=(GraphEndBase1-GraphStartBase1,(GraphEndBase1-BaseOffset1)*spHScale1),
                                  gymin=(MinLog2VScaleSpan1+0.25,(MinLog2VScaleSpan1+0.25)*spVScale1),
                                  gymax=(-MinLog2VScaleSpan1-0.25,(-MinLog2VScaleSpan1-0.25)*spVScale1))
    SEPlot2.vgrid(gtitle=' ',gxlabel='Reference Position', gylabel='log2(coverage)',gmajorwidth=10,gmajorcolor='cyan',gylabelcolor='bluegreen',
                                  gminorwidth=1,gminorcolor='cyan',gxminor=True,gyminor=True,
                                  gyaxisfontsize=300,gylabelfontsize=300,gxaxisfontsize=300,gxlabelfontsize=300,
                                  gxmin=(GraphStartBase1,(GraphStartBase1-BaseOffset1)*spHScale1),
                                  gxmax=(GraphEndBase1,(GraphEndBase1-BaseOffset1)*spHScale1),
                                  gymax=(MinLog2VScaleSpan2+0.25,(MinLog2VScaleSpan2+0.25)*spVScale3),
                                  gymin=(-MinLog2VScaleSpan2-0.25,(-MinLog2VScaleSpan2-0.25)*spVScale3))
    ## , and a manually set up grid for linear values [right side])
    LabelStart1 = max(spi,GraphEndBase1-BaseOffset1)*spHScale1        
    SEPlot1.vtext(text=str(0)+'x',x1=LabelStart1,yc=0, font='DejaVuSansBold 300', color='magenta2')
    for spl1 in range(1,MinLog2VScaleSpan1+1):
        SEPlot1.vtext(text=str(spYAxisStartToBiasConversion1*spl1),x1=LabelStart1,yc=spl1*spVScale1, font='DejaVuSansBold 300', color='magenta')
        SEPlot1.vtext(text=str(spYAxisStartToBiasConversion1*spl1),x1=LabelStart1,yc=-spl1*spVScale1, font='DejaVuSansBold 300', color='purple')
    SEPlot1.vtext(text='(+) Start enrichment',xc=LabelStart1+1000,yc=MinLog2VScaleSpan1*VScale1/2,rotate=-90, font='DejaVuSansBold 300', color='magenta')
    SEPlot1.vtext(text='(-) Start enrichment',xc=LabelStart1+1000,yc=-MinLog2VScaleSpan1*VScale1/2,rotate=-90, font='DejaVuSansBold 300', color='purple')
    SEPlot1.ymin = min(SEPlot1.ymin,SEPlot2.ymin-spVScale1*MinLog2VScaleSpan1-spVScale3*MinLog2VScaleSpan2-600)
vdisplay(OutBase2+"."+GraphicMode1,display=DisplayGraphic1)
if StartEndPlot1:
    if spCoverage1:
        SEPlot1 += SEPlot2.transform(yd=-spVScale1*MinLog2VScaleSpan1-spVScale3*MinLog2VScaleSpan2-600)
        SEPlot1.vcircle(xc=0,yc=0,r=0,fill='blue',stroke='blue',colorkey='+ Coverage')
        SEPlot1.vcircle(xc=0,yc=0,r=0,fill='green',stroke='green',colorkey='- Coverage')
    SEPlot1.vcolorkey(ckfont='DejaVuSansBold 300')

    SEPlot1.vdisplay('StartEndPlot_'+OutBase2+"."+GraphicMode1,display=DisplayGraphic1)
    SEData1 = open('StartEndData_'+OutBase2+".tdv",mode='w')
    SEData1.write('Position\tBase\t')
    for featurename1 in GraphPathD1:
        if GraphPathD1[featurename1].points:
            SEData1.write("'"+featurename1+'_GraphedScaled'+"'"+'\t'+"'"+featurename1+'_RawData'+"'"+'\t')
    SEData1.write('Sense_Coverage\tAntiSense_Coverage\t')
    SEData1.write('Flags\r')
    for spi in range(len(S1)):
        SEData1.write(str(spi+BaseOffset1+1)+'\t')
        SEData1.write(S1[spi]+'\t')
        for fn1 in GraphPathD1:
            if GraphPathD1[fn1].points:                
                fnD1 = GraphPathD1[fn1].graphed
                fnD2 = GraphPathD1[fn1].data
                if spi in fnD1:
                    SEData1.write(str(fnD1[spi])+'\t')
                    SEData1.write(str(fnD2[spi])+'\t')
                else:
                    SEData1.write('\t\t')
        SEData1.write(str(CountD1['p'][spi])+'\t')
        SEData1.write(str(CountD1['m'][spi])+'\t')        
        if spi in spFlags1:
            SEData1.write(str(spFlags1[spi])) 
        SEData1.write('\r')
    SEData1.close()
gc.collect()
testMod1 = str(sys.modules).lower()
if not('idle' in testMod1) and not('thonny' in testMod1):
    sys.exit(0)
##
##  Copyright 2021-2023, Andrew Fire and Stanford University
##  Revision List (partial)
##  Version BQ improves the detection of R1 and R2 filenames and adds an ability to explicitly specify file names if imputation is resulting
##     in "file not found" or other errors.
##  Version BR adds an option RequireR1Match, which was implicit in previous versions.  This requires that there be an identified match for Read1
##     (in the interval before MaxR1MaxTrim5) or Read2 will be ignored.  Previously this was defaulted to "True".   Current default is to "False",
##     So now R2 reads will be counted even if R1 matches nothing.  This can also be set to True with the keyword Legacy "Legacy=True")
##  Version CD03 uses numerical representation of sequences to limit memory footprint, and fixes the AssertSense option (thanks to Drew Galls)
##     Adds a failsafe to convert nonstandard bases to "G" (thanks to William Wang)
##     And includes Start/End plots, and a somewhat more modest memory footprint
##  Important but slight change as of version CD05, "N" bases in the reference sequence will be converted to Gs 
