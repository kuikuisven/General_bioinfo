"""
    Nicolas Rapin
    
    MINiML import funciton helper
    does a lot of things related to bioinfrmatics anlayses.
    
    
"""
#import MacOS
from lxml import etree
#from BeautifulSoup import UnicodeDammit
import h5py
import numpy
import os
import shutil
import ftplib
import time
import gzip
import gc
import sys

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

from scipy.stats import stats
import hcluster
from rpy2.robjects.packages import importr
import Bio.Cluster
import KEGG
import pandas as pa
from Queue import Queue
from threading import Thread
import traceback
import pandas as pa
#from SC import  *

__author__ = "Nicolas Rapin (nrapin@binf.ku.dk)"
__version__ = "$0.1$ "
__date__ = "$Date: 22/2/10 $"
__copyright__ = "Copyright ?"
__license__ = "Python"

#Some global variables
home = os.environ['HOME']
MINiML_header='{http://www.ncbi.nlm.nih.gov/projects/geo/info/MINiML}'
#GEO_data_dl_dir='/Volumes/Donnees/GEO/tmp/'
#tmp_file_path='/Users/nrapin/iris/GEO/_GEO_db.hdf5'
#db_file_path='/Users/nrapin/iris/GEO/GEO_db.hdf5'
GEO_data_dl_dir=home+'/GEO/tmp/'
tmp_file_path=home+'/GEO/'
db_file_path=home+'/GEO/GEO_db.hdf5'
refseq2gene=0
gene2refseq=0
platform_id=''

def smart_zscore_p(x):
    from  numpy import  argwhere, max ,min
    from scipy import  nanmean, nanstd
    not_null = argwhere(x != 0).ravel()
    if len(not_null) > 0:
        mean_gene = nanmean( x[not_null])
        std_gene  = nanstd( x[not_null])
        if std_gene != 0:
            x -= mean_gene
            x /= std_gene
            mini = min(x)
            x -= mini
        maxi = max(x)
        x /= maxi
    return x
    pass

def smart_zscore(X,parallel = True):
    """docstring for smart_zscore
        does zscoring of data, without including zeros in the computation.
        also rescale between [0-1]
    """
    from  numpy import  argwhere, max ,min
    from scipy import  nanmean, nanstd
    from  numpy import  arange, array
    import scipy
    import pandas as pa
    from joblib import Parallel, delayed
    import multiprocessing
    if parallel:
        X = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(smart_zscore_p)(X[gene_i,:]) for gene_i in arange(X.shape[0]))  # n_jobs = number of pro
    else :
        for gene_i in arange(X.shape[0]):
            #print gene_i
            not_null = argwhere(X[gene_i,:] != 0).ravel()
            if len(not_null) > 0:
                mean_gene = nanmean( X[gene_i,:][not_null])
                std_gene  = nanstd( X[gene_i,:][not_null])
                if std_gene != 0:
                    X[gene_i,:] -= mean_gene
                    X[gene_i,:] /= std_gene
                    mini = min(X[gene_i,:])
                    X[gene_i,:] -= mini
                maxi = max(X[gene_i,:])
                X[gene_i,:] /= maxi
                if False in  pa.DataFrame(X[gene_i,:]).notnull().values.ravel():
                    break
    return array(X)
    pass
def smart_ranks(X):
    """docstring for smart_zscore
        does zscoring of data, without including zeros in the computation.
        also rescale between [0-1]
    """
    from  numpy import  argwhere, max ,min
    from scipy import  nanmean, nanstd
    from scipy.stats import rankdata
    from  numpy import  arange, array
    import scipy
    import pandas as pa
    from joblib import Parallel, delayed
    import multiprocessing
    if 0:
        X = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(smart_zscore_p)(X[gene_i,:]) for gene_i in arange(X.shape[0]))  # n_jobs = number of pro
    else :
        for gene_i in arange(X.shape[0]):
            #print gene_i
            not_null = argwhere(X[gene_i,:] != 0).ravel()
            if len(not_null) > 0:
                ranks = rankdata(X[gene_i,not_null] ,method='dense')
                ranks = max(ranks) - ranks  + 1.
                ranks /= max(ranks)
                X[gene_i,not_null] = ranks
    return array(X)



def et_oui():
    print 'Et oui!'
    pass

def import_kallisto_outpout(data_file =  'kallisto_table.csv',sep=' ',convert_gene_names=True):
    """
    Function to import outpout from sleuth, kallisto analysis lib in R.
    default file is:
    kallisto_table.csv
    """
    import mygene #pip install mygene
    g=pa.read_csv(data_file,sep=sep)
    samples = g['sample'].unique()
    dat=g[g['sample'] == samples[0]][['target_id','est_counts']]
    for s in samples[1:]:
        dat = pa.merge(dat,g[g['sample'] == s][['target_id','est_counts']],on = 'target_id')
    dat.columns = numpy.concatenate([['target_id'], samples])
    if convert_gene_names:
        mg = mygene.MyGeneInfo()
        out = mg.querymany(dat.target_id.tolist(), scopes='symbol', fields='ensembltranscript', species='human')
        c=pa.DataFrame([[b['query'],b['entrezgene']] for b in out if b.has_key('entrezgene') ], columns=['id','target_id'])
        return pa.merge(c,dat,on='target_id')
    else:
        return dat
    pass

def make_colormap(seq,ncolors=10):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).

    example:
    c = mcolors.ColorConverter().to_rgb
    rvb = make_colormap([c('red'), c('violet'), 0.33, c('violet'), c('blue'), 0.66, c('blue')])
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict, ncolors)


def import_bcbio_output_Deseq2(exprs='annotated_combined.counts',summary = 'project-summary.yaml',return_raw_counts=False):
    """ imports bcbio_nextgen's output file (just raw expression)
        this is the file where the expression is stored: annotated_combined.counts
        process the file with DESeq2
    """
    import yaml
    import re
    from numpy import array, argwhere, mean
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    DESeq2 = import_rlib(['DESeq2'])[0]

    if not os.path.isfile(exprs):
        print "expression file not present... check path."
        return



    data2 = pa.read_csv(exprs,sep='\t')

    names =  array(data2.columns[1:-1])
    try:
        names = array([re.search('(^[A-Za-z]*[_-][A-Za-z]*)', i ).group(0)  for i in names ])
    except Exception,e:
        print e



    design = pa.DataFrame(array(names),columns=['type'])
    design.index = array(data2.columns[1:-1])
    d2 = data2[data2.columns[1:-1]]
    #d2.columns=names
    d2.index=array(data2.id)
    d2r = pandas2ri.py2ri(d2)



    if return_raw_counts:
        vstMat = numpy.array(d2)
        counts_data = array(vstMat)
        print vstMat.shape
    else:
        GenomicRanges = import_rlib(['GenomicRanges'])[0]
        BiocGenerics = import_rlib(['BiocGenerics'])[0]
        dds = DESeq2.DESeqDataSetFromMatrix(countData = d2r,colData=design,design = robjects.r(" ~ type"))
        #rld = DESeq2.rlog(dds)
        vsd = DESeq2.varianceStabilizingTransformation(dds)
        #rlogMat = GenomicRanges.assay(rld)
        vstMat = GenomicRanges.assay(vsd)
        counts_data = array(vstMat)
    non_empty =  argwhere(mean(counts_data , axis=1) != min(mean(array(vstMat), axis=1)) ).ravel()
    dat  = {}
    dat['data'] = counts_data[non_empty,:]
    dat['colnames_orig'] = array(data2.columns[1:-1])
    dat['rownames'] = array(data2.symbol)[non_empty]
    # a bit of precessing...
    dat['colnames'] = dat['colnames_orig'] # array([re.search('(^[A-Za-z]*[_-][A-Za-z]*)', i ).group(0)  for i in dat['colnames_orig'] ])
    dat['stemness'] = stemness_by_class(dat)
    pickle_save(dat,'annotated_combined.counts.DESeq2.pkl')
    return dat



def import_bcbio_output(exprs='annotated_combined.counts',summary = 'project-summary.yaml'):
    """ imports bcbio_nextgen's output file (just raw expression)
        this is the file where the expression is stored: annotated_combined.counts
        this is the file for the summary: project-summary.yaml
    """
    import yaml
    import re
    from numpy import array, argwhere, mean
    limma = import_rlib(['limma'])[0]

    if not os.path.isfile(exprs):
        print "expression file not present... check path."
        return
    if not os.path.isfile(summary):
        print "summary file not present... check path."
        return

    data = pa.read_csv(exprs, sep = '\t')
    summary = yaml.load(open(summary,'r'))

    counts_data = array(data[data.columns[1:-1]])
    non_empty =  argwhere(mean(array(data[data.columns[1:-1]]) , axis=1) !=0).ravel()
    dat  = {}
    dat['data'] = counts_data[non_empty,:]
    dat['colnames_orig'] = array(data.columns[1:-1])
    dat['rownames'] = array(data.symbol)[non_empty]
    # a bit of precessing...
    dat['data'] = array(limma.voom(dat['data'])[0])
    dat['colnames'] = array([re.search('(^[A-Za-z]*[_-][A-Za-z]*)', i ).group(0)  for i in dat['colnames_orig'] ])
    dat['stemness'] = stemness_by_class(dat)
    pickle_save(dat,'annotated_combined.counts.pkl')
    return dat
    pass


def import_bcbio_output_isoforms(exprs='combined.isoform.fpkm'):
    """ imports bcbio_nextgen's output file (just raw expression)
        this is the file where the expression is stored: annotated_combined.counts
        this is the file for the summary: project-summary.yaml
    """
    import mygene #pip install mygene
    import yaml
    import re
    from numpy import array, argwhere, mean
    limma = import_rlib(['limma'])[0]

    if not os.path.isfile(exprs):
        print "expression file not present... check path."
        return

    data = pa.read_csv(exprs, sep = '\t')
    counts_data = array(data[data.columns[1:]])
    non_empty =  argwhere(mean(array(data[data.columns[1:]]) , axis=1) !=0).ravel()
    dat  = {}
    dat['data'] = counts_data[non_empty,:]
    dat['colnames_orig'] = array(data.columns[1:])
    dat['rownames'] = array(data.id)[non_empty]
    # a bit of precessing...
    dat['data'] = array(limma.voom(dat['data'])[0])
    dat['colnames'] = array([re.search('(^[A-Za-z]*[_-][A-Za-z]*)', i ).group(0)  for i in dat['colnames_orig'] ])
    dat['stemness'] = stemness_by_class(dat)
    mg = mygene.MyGeneInfo()
    out = mg.querymany(array(ds.From).tolist(), scopes='ensembltranscript', fields='symbol', species='mouse')
    c=pa.DataFrame([[b['query'],b['symbol']] for b in out if b.has_key('symbol')])
    dd = {i[0]:i[1] for i in array(c)}
    symbols = [dd[i] if i in dd.keys() else i for i in a.id]
    dat['rownames'] = array(symbols)[non_empty]
    pickle_save(dat,'annotated_combined.counts.pkl')
    return dat
    pass





def probes_to_genes_hugene_10st(p):
    """docstring for probes_to_genes"""


    robjects.r('''
  	hugene1_probe2symb <- function(probes){
  	if (!"hugene10sttranscriptcluster.db" %in% installed.packages()) {source("http://bioconductor.org/biocLite.R"); biocLite("hugene10sttranscriptcluster.db")}
  	library(hugene10sttranscriptcluster.db);
  	keys=keys(hugene10sttranscriptcluster.db) # get it all
  	symb=select(hugene10sttranscriptcluster.db, keys, "SYMBOL", "PROBEID" ) # 3468 probes have multiple targets or annotations...
  	symb=symb[!duplicated(symb[,1]),]
  	return(symb[match(probes, symb[,1]),][,2])
  }
    ''')
    r_human_probe2symb =robjects.globalenv['hugene1_probe2symb']

    genes = r_human_probe2symb(p)
    return genes
    pass


def probes_to_genes_mogene_10st(p):
    """docstring for probes_to_genes"""


    robjects.r('''
  	mogene1_probe2symb <- function(probes){
  	if (!"mogene10sttranscriptcluster.db" %in% installed.packages()) {source("http://bioconductor.org/biocLite.R"); biocLite("mogene10sttranscriptcluster.db")}
  	library(mogene10sttranscriptcluster.db);
  	keys=keys(mogene10sttranscriptcluster.db) # get it all
  	symb=select(mogene10sttranscriptcluster.db, keys, "SYMBOL", "PROBEID" ) # 3468 probes have multiple targets or annotations...
  	symb=symb[!duplicated(symb[,1]),]
  	return(symb[match(probes, symb[,1]),][,2])
  }
    ''')
    r_human_probe2symb =robjects.globalenv['mogene1_probe2symb']

    genes = r_human_probe2symb(p)
    return genes
    pass


def probes_to_genes_430(p):
    """docstring for probes_to_genes"""
    robjects.r('''m430_probe2symb <- function(probes){
  	if (!"mouse4302.db" %in% installed.packages()) {source("http://bioconductor.org/biocLite.R"); biocLite("mouse4302.db")}
  	library(mouse4302.db);
  	keys=keys(mouse4302.db) # get it all
  	symb=select(mouse4302.db, keys, "SYMBOL", "PROBEID" ) # 3468 probes have multiple targets or annotations...
  	symb=symb[!duplicated(symb[,1]),]
  	return(symb[match(probes, symb[,1]),][,2])
}''')
    r_human_probe2symb =robjects.globalenv['m430_probe2symb']
    genes = r_human_probe2symb(p)
    return genes
    pass

def probes_to_genes_plus_2(p):
    """docstring for probes_to_genes"""
    import collections
    robjects.r('''
  	p2_human_probe2symb <- function(probes){
  	if (!"hgu133plus2.db" %in% installed.packages()) {source("http://bioconductor.org/biocLite.R"); biocLite("hgu133plus2.db")}
  	library(hgu133plus2.db);
  	keys=keys(hgu133plus2.db) # get it all
  	symb=select(hgu133plus2.db, keys, "SYMBOL", "PROBEID" ) # 3468 probes have multiple targets or annotations...
  	symb=symb[!duplicated(symb[,1]),]
  	return(symb[match(probes, symb[,1]),][,2])
  }
    ''')
    r_human_probe2symb =robjects.globalenv['p2_human_probe2symb']

    genes = r_human_probe2symb(p)
    gene_to_probes = collections.defaultdict(list)
    return genes
    pass

def genes_to_probes_plus_2(p):
    import collections
    ref=pickle_load('/Volumes/Hdd4/BRIC/HemaExplorer/datasets/AML_MILE_VERHAAK_vs_nl.pkl')['rownames']
    genes = probes_to_genes_plus_2(ref)
    gene_to_probes = collections.defaultdict(list)
    for i,names in enumerate(genes):
            try:
                gene_to_probes[names].append(ref[i])
            except Exception,e:
                pass
    return [gene_to_probes[i] for i in p]
    pass

def genes_to_probes_430(p):
    import collections
    ref= pickle_load('/Volumes/Hdd4/BRIC/HemaExplorer/datasets/nl_mouse_data.pkl')['rownames']
    genes = probes_to_genes_430(ref)
    gene_to_probes = collections.defaultdict(list)
    for i,names in enumerate(genes):
            try:
                gene_to_probes[names].append(ref[i])
            except Exception,e:
                pass

    return [gene_to_probes[i] for i in p]
    pass
def genes_to_probes_mogene_10st(p):
    import collections
    ref= pickle_load('/Volumes/Hdd4/BRIC/HemaExplorer/datasets/immgen_Stromal_cells.pkl')['rownames']
    genes = probes_to_genes_mogene_10st(ref)
    gene_to_probes = collections.defaultdict(list)
    for i,names in enumerate(genes):
            try:
                gene_to_probes[names].append(ref[i])
            except Exception,e:
                pass

    return [gene_to_probes[i] for i in p]
    pass

def genes_to_probes_hugene_10st(p):
    import collections
    ref= pickle_load('/Users/nrapin/Dropbox/BRIC/People/Kristina/Mette_Levinsen_CEL-files/normalized_data.pkl')['rownames']
    genes = probes_to_genes_hugene_10st(ref)
    gene_to_probes = collections.defaultdict(list)
    for i,names in enumerate(genes):
            try:
                gene_to_probes[names].append(ref[i])
            except Exception,e:
                pass

    return [gene_to_probes[i] for i in p]
    pass



class Worker(Thread):
    """Thread executing tasks from a given tasks queue"""
    def __init__(self, tasks):
        Thread.__init__(self)
        self.ncalls  = 0
        self.nscalls = 0
        self.nfcalls = 0
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.tasks.get()
            if '_threadname' in kargs:
                self.name = kargs.pop('_threadname')
            self.ncalls += 1
            try:
                func(*args, **kargs)
                self.nscalls += 1
            except Exception, e:
                self.nfcalls += 1
                traceback.print_exception(sys.exc_info()[0], e, sys.exc_info()[2])
            self.tasks.task_done()

class ThreadPool:
    """Pool of threads consuming tasks from a queue"""
    def __init__(self, num_threads):
        self.tasks = Queue(num_threads)
        for _ in range(num_threads): Worker(self.tasks)

    def add_task(self, func, *args, **kargs):
        """Add a task to the queue"""
        self.tasks.put((func, args, kargs))

    def wait_completion(self):
        """Wait for completion of all the tasks in the queue"""
        self.tasks.join()




def kdtree( data, leafsize=10 ):
    """
    build a kd-tree for O(n log n) nearest neighbour search

    input:
        data:      2D ndarray, shape =(ndim,ndata), preferentially C order
        leafsize:   max. number of data points to leave in a leaf

    output:
        kd-tree:    list of tuples
    """

    ndim = data.shape[0]
    ndata = data.shape[1]

    # find bounding hyper-rectangle
    hrect = numpy.zeros((2,data.shape[0]))
    hrect[0,:] = data.min(axis=1)
    hrect[1,:] = data.max(axis=1)

    # create root of kd-tree
    idx = numpy.argsort(data[0,:], kind='mergesort')
    data[:,:] = data[:,idx]
    splitval = data[0,ndata/2]

    left_hrect = hrect.copy()
    right_hrect = hrect.copy()
    left_hrect[1, 0] = splitval
    right_hrect[0, 0] = splitval

    tree = [(None, None, left_hrect, right_hrect, None, None)]

    stack = [(data[:,:ndata/2], idx[:ndata/2], 1, 0, True),
             (data[:,ndata/2:], idx[ndata/2:], 1, 0, False)]

    # recursively split data in halves using hyper-rectangles:
    while stack:

        # pop data off stack
        data, didx, depth, parent, leftbranch = stack.pop()
        ndata = data.shape[1]
        nodeptr = len(tree)

        # update parent node

        _didx, _data, _left_hrect, _right_hrect, left, right = tree[parent]

        tree[parent] = (_didx, _data, _left_hrect, _right_hrect, nodeptr, right) if leftbranch \
            else (_didx, _data, _left_hrect, _right_hrect, left, nodeptr)

        # insert node in kd-tree

        # leaf node?
        if ndata <= leafsize:
            _didx = didx.copy()
            _data = data.copy()
            leaf = (_didx, _data, None, None, 0, 0)
            tree.append(leaf)
        # not a leaf, split the data in two
        else:
            splitdim = depth % ndim
            idx = numpy.argsort(data[splitdim,:], kind='mergesort')
            data[:,:] = data[:,idx]
            didx = didx[idx]
            nodeptr = len(tree)
            stack.append((data[:,:ndata/2], didx[:ndata/2], depth+1, nodeptr, True))
            stack.append((data[:,ndata/2:], didx[ndata/2:], depth+1, nodeptr, False))
            splitval = data[splitdim,ndata/2]
            if leftbranch:
                left_hrect = _left_hrect.copy()
                right_hrect = _left_hrect.copy()
            else:
                left_hrect = _right_hrect.copy()
                right_hrect = _right_hrect.copy()
            left_hrect[1, splitdim] = splitval
            right_hrect[0, splitdim] = splitval
            # append node to tree
            tree.append((None, None, left_hrect, right_hrect, None, None))

    return tree


class ProgressBar:
    def __init__(self, min_value = 0, max_value = 100, width=77,**kwargs):
        self.char = kwargs.get('char', '#')
        self.mode = kwargs.get('mode', 'dynamic') # fixed or dynamic
        if not self.mode in ['fixed', 'dynamic']:
            self.mode = 'fixed'

        self.bar = ''
        self.min = min_value
        self.max = max_value
        self.span = max_value - min_value
        self.width = width
        self.amount = 0    # When amount == max, we are 100% done
        self.update_amount(0)


    def increment_amount(self, add_amount = 1):
        """
        Increment self.amount by 'add_ammount' or default to incrementing
        by 1, and then rebuild the bar string.
        """
        new_amount = self.amount + add_amount
        if new_amount < self.min: new_amount = self.min
        if new_amount > self.max: new_amount = self.max
        self.amount = new_amount
        self.build_bar()


    def update_amount(self, new_amount = None):
        """
        Update self.amount with 'new_amount', and then rebuild the bar
        string.
        """
        if not new_amount: new_amount = self.amount
        if new_amount < self.min: new_amount = self.min
        if new_amount > self.max: new_amount = self.max
        self.amount = new_amount
        self.build_bar()


    def build_bar(self):
        """
        Figure new percent complete, and rebuild the bar string base on
        self.amount.
        """
        diff = float(self.amount - self.min)
        percent_done = int(round((diff / float(self.span)) * 100.0))

        # figure the proper number of 'character' make up the bar
        all_full = self.width - 2
        num_hashes = int(round((percent_done * all_full) / 100))

        if self.mode == 'dynamic':
            # build a progress bar with self.char (to create a dynamic bar
            # where the percent string moves along with the bar progress.
            self.bar = self.char * num_hashes
        else:
            # build a progress bar with self.char and spaces (to create a
            # fixe bar (the percent string doesn't move)
            self.bar = self.char * num_hashes + ' ' * (all_full-num_hashes)

        percent_str = str(percent_done) + "%"
        self.bar = '[ ' + self.bar + ' ] ' + percent_str


    def __str__(self):
        return str(self.bar)




def make_D3_network(data,out='HemaTree.json',prune_percent=10):
	"""docstring for make_network
		take the usual data structure with labels and stuffs ,and produce a networkX network out of that, based on correlation/clustering
	"""
	import pandas as pa
	from numpy import *
	import orange, orngClustering, Orange
	from Orange import orange
	dat_ave =   average_cells_by_stemness(data)
	#fast_cor(dat_ave, probes_to_use=p2u,plot=True,figure ='correlation.pdf')

	#and do some orange clustering:
	p2u=select_probes_by_variance(dat_ave,var_thr=3)
	root = orngClustering.hierarchicalClustering(orange.ExampleTable(dat_ave['data'][p2u,:].T),\
                    distance_constructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    #linkage=orange.HierarchicalClustering.Complete,\
                    linkage=orange.HierarchicalClustering.Ward)
	prune(root,root.height/prune_percent)
	f=open(out,'w')
	f.write(print_clustering3(root,dat_ave))
	f.close()
	export_HTML_tree(out)
	pass
def prune(cluster, h):
	if cluster.branches:
		if cluster.height < h:
			cluster.branches = None
		else:
			for branch in cluster.branches:
				prune(branch, h)
	return cluster

def print_clustering2(cluster):
	if cluster.branches:
		return "(%s %s)" % (print_clustering2(cluster.left), print_clustering2(cluster.right))
	else:
		return str(tuple(cluster))


def print_clustering3(cluster,data,level=0):
	'''
	model : { "name": "HSC", "children": [ { "name": "MPP", "children": [ {"name":"CMP"} ]} ] }

	'''
	if cluster.branches:
		return ' { "name":" " , "size":2 , "children":[%s,%s]\n }' % (print_clustering3(cluster.left,data,level=level+1), print_clustering3(cluster.right,data,level=level+1)) + ' '*level
	else:
		return ','.join([ '{"name":"' +  str(data['colnames'][i]) +'","size":10}\n'+ ' '*level for i in tuple(cluster)])

def export_HTML_tree(json):
	"""docstring for export_HTML_tree"""
	page = '''
	<!DOCTYPE html>
	<html>
	  <head>
	    <meta http-equiv="Content-Type" content="text/html;charset=utf-8"/>
	    <link type="text/css" rel="stylesheet" href="style.css"/>
	    <script type="text/javascript" src="d3/d3.js"></script>
	    <script type="text/javascript" src="d3/d3.layout.js"></script>
	    <style type="text/css">


	.node circle {
	  cursor: pointer;
	  fill: #fff;
	  stroke-opacity:.51;
	  stroke: steelblue;
	  stroke-width: 1.5px;
	}

	.node text {
	  font-size: 11px;
	}

	path.link {
	  fill: none;
	  stroke: #ccc;
	  stroke-width: 1.5px;
	}
	    </style>
	  </head>
	  <body>
	    <div id="body">
	      <div id="footer">
	        d3.layout.tree.test
	        <div class="hint">click or option-click to expand or collapse</div>
	      </div>
	    </div>
	    <script type="text/javascript">

	var m = [20, 120, 20, 120],
	    w = 900 - m[1] - m[3],
	    h = 500 - m[0] - m[2],
	    i = 0,
	    root;

	var tree = d3.layout.tree()
	    .size([h, w]);

	var diagonal = d3.svg.diagonal()
	    .projection(function(d) { return [d.y, d.x]; });

	var vis = d3.select("#body").append("svg:svg")
	    .attr("width", w + m[1] + m[3])
	    .attr("height", h + m[0] + m[2])
	  .append("svg:g")
	    .attr("transform", "translate(" + m[3] + "," + m[0] + ")");

	d3.json("%s", function(json) {//"flare.json", function(json) {
	  root = json;
	  root.x0 = h / 2;
	  root.y0 = 0;


		// Initialize the display to show a few nodes.
		//root.children.forEach(toggleAll);


	  // Initialize the display to show a few nodes.
//	  root.children.forEach(toggleAll);
	  //toggle(root.children[0]);
	//  toggle(root.children[1]);
	//  toggle(root.children[0].children[0]);
	//  toggle(root.children[0].children[0].children[0]);
	//  toggle(root.children[0].children[0].children[0].children[1]);
	// toggle(root.children[1].children[1]);
	//  toggle(root.children[9]);
	// toggle(root.children[9].children[0]);

	  update(root);
	});

	function update(source) {
	  var duration = d3.event && d3.event.altKey ? 5000 : 500;

	  // Compute the new tree layout.
	  var nodes = tree.nodes(root).reverse();

	  // Normalize for fixed-depth.
	  nodes.forEach(function(d) { d.y = d.depth * 60; });

	  // Update the nodes...
	  var node = vis.selectAll("g.node")
	      .data(nodes, function(d) { return d.id || (d.id = ++i); });

	  // Enter any new nodes at the parent's previous position.
	  var nodeEnter = node.enter().append("svg:g")
	      .attr("class", "node")
	      .attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
	      .on("click", function(d) { toggle(d); update(d); });

	  nodeEnter.append("svg:circle")
	      .attr("r", 1e-6)
	      .style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

	  nodeEnter.append("svg:text")
	      .attr("x", function(d) { return d.children || d._children ? -10 : 10; })
	      .attr("dy", ".35em")
	      .attr("text-anchor", function(d) { return d.children || d._children ? "end" : "start"; })
	      .attr("transform", function(d) { return (d.children || d._children )&& d.name.length > 10 ?  "rotate(-15)":"rotate(0)";})
	      .text(function(d) { return d.name; })
	      .style("fill-opacity", 1e-6);

	  // Transition nodes to their new position.
	  var nodeUpdate = node.transition()
	      .duration(duration)
	      .attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; });

	  nodeUpdate.select("circle")
	      .attr("r", function(d) {return d.size;} )//14.5)
	      .style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

	  nodeUpdate.select("text")
	      .style("fill-opacity", 1);

	  // Transition exiting nodes to the parent's new position.
	  var nodeExit = node.exit().transition()
	      .duration(duration)
	      .attr("transform", function(d) { return "translate(" + source.y + "," + source.x + ")"; })
	      .remove();

	  nodeExit.select("circle")
	      .attr("r", 1e-6);

	  nodeExit.select("text")
	      .style("fill-opacity", 1e-6);

	  // Update the links...
	  var link = vis.selectAll("path.link")
	      .data(tree.links(nodes), function(d) { return d.target.id; });

	  // Enter any new links at the parent's previous position.
	  link.enter().insert("svg:path", "g")
	      .attr("class", "link")
	      .attr("d", function(d) {
	        var o = {x: source.x0, y: source.y0};
	        return diagonal({source: o, target: o});
	      })
	    .transition()
	      .duration(duration)
	      .attr("d", diagonal);

	  // Transition links to their new position.
	  link.transition()
	      .duration(duration)
	      .attr("d", diagonal);

	  // Transition exiting nodes to the parent's new position.
	  link.exit().transition()
	      .duration(duration)
	      .attr("d", function(d) {
	        var o = {x: source.x, y: source.y};
	        return diagonal({source: o, target: o});
	      })
	      .remove();

	  // Stash the old positions for transition.
	  nodes.forEach(function(d) {
	    d.x0 = d.x;
	    d.y0 = d.y;
	  });
	}

	// Toggle children.
	function toggle(d) {
	  if (d.children) {
	    d._children = d.children;
	    d.children = null;
	  } else {
	    d.children = d._children;
	    d._children = null;
	  }
	}

	    </script>
	  </body>
	</html>

	'''%json
	f=open(json.replace('json','html'),'w')
	f.write(page)
	f.close()
	pass



def generate_GSEA_one_vs_one_data(data,data_to_use=None,out='GSEA_individual'):
    if data_to_use==None:
        data_to_use=numpy.arange(data['data'].shape[1])
    if 'GSEA_individual' in os.listdir('.'):
        import shutil
        shutil.rmtree('GSEA_individual')
        print 'removing GSEA_individual dir...'
    os.mkdir(out)
    os.chdir(out)
    tile=numpy.tile
    for i in data_to_use:
        name =data['cel_file'][i].replace('.CEL.gz','')
        print name
        sample={}
        sample['data'] =numpy.concatenate([tile(data['data'][:,i],(3,1)).T,tile(data['chimera_sample'].T[:,i],(3,1)).T],axis=1)
        sample['colnames'] = numpy.array([name,name,name,'nl','nl','nl'])
        sample['rownames'] = data['rownames']
        sample['stemness'] = stemness_by_class(sample)
        export_values_for_gsea(sample,filename=name+'.txt',log2=False)

    os.chdir('..')
    pass

def generate_GSEA_aml_vs_nl_data(data,data_to_use=None,out='GSEA_aml_vs_nl'):
    if data_to_use==None:
        data_to_use=numpy.arange(data['data'].shape[1])
    tile=numpy.tile
    name =out
    print name
    sample={}
    sample['data'] =numpy.concatenate([data['data'][:,data_to_use] , data['chimera_sample'].T[:,data_to_use]], axis=1)
    sample['colnames'] = numpy.concatenate([['aml']*len(data_to_use) ,['nl']*len(data_to_use) ])
    sample['rownames'] = data['rownames']
    sample['stemness'] = stemness_by_class(sample)
    for i in limma_multiclass(sample,p=1e-5,limit=500)['indexes']:
        print  sample['rownames'][i]
    export_values_for_gsea(sample,filename=name+'.txt',log2=True)
    return sample
    pass

def generate_GSEA_aml_vs_relapse_data(data,data_to_use=None,out='GSEA_aml_vs_relapse'):
    from numpy import array
    if data_to_use==None:
        print 'provide list of paired samples...'
    tile=numpy.tile
    name =out
    print 'Generating ' + name
    denovo=array(data_to_use.values())[:,0]
    relapse=array(data_to_use.values())[:,1]
    sample={}
    sample['data'] =numpy.concatenate([data['data'][:,denovo] , data['data'][:,relapse]], axis=1)
    sample['colnames'] = numpy.concatenate([['de_novo']*len(denovo) ,['relapse']*len(relapse) ])
    sample['rownames'] = data['rownames']
    sample['stemness'] = stemness_by_class(sample)
    for i in limma_multiclass(sample,p=1e-5,limit=500)['indexes']:
        print  sample['rownames'][i]
    #f=open()
    print 'colname'
    print sample['colnames'].shape
    print 'rewname'
    print sample['rownames'].shape
    print 'stemness'
    print sample['stemness'].shape
    print 'data', str(sample['data'].shape)
    export_values_for_gsea(sample,filename=name+'.txt',log2=True)

    pass


def Compute_pvals_survival_data(amls,genes_all,no_plot=True, plot_anyways=False):
    """docstring for Compute_pvals_survival_data"""
    import numpy, os
    #amls['fc'] = amls['raw_fold_change'].T

    data_to_use =  range(len(amls['colnames']))
    pvals=[]
    for index, genes in enumerate(genes_all):
        print 1.*index/len(genes_all)*100.
        p2u_bo =get_probe_index([genes],amls['rownames'])[0]
        #amls = pickle_load('aml_nk_clusters.pkl')
        #amls['data']=amls['fc']

        #x= amls['data'][p2u_bo, :][ :,data_to_use]
        x= amls['data'][p2u_bo, :][ :,data_to_use]
        order=get_sorted_index(x)

        sl = len(order) / 4
        #quartiles
        low = range(sl*0,sl*1)
        high = range(sl*3,sl*4)
        #half
        low = range(0,len(order) / 2)
        high = range(len(order)/2,len(order))


        data_to_use = numpy.array(data_to_use)
        amls['stemness'] = numpy.array(amls['stemness'] ,dtype='i')
        amls['stemness'][...] = 1
        amls['stemness'][data_to_use[order[high]]] = 21
        amls['stemness'][data_to_use[order[low]]] = 22

        #amls['EventOS'][numpy.argwhere(amls['OS'] >= 365).reshape(-1) ] = 0
        #amls['OS'] [numpy.argwhere(amls['OS'] >= 365).reshape(-1) ] = 365

        ##print genes, amls['OS'][data_to_use[order[high]]].mean() - amls['OS'][data_to_use[order[low]]].mean()
        day_diff = amls['OS'][data_to_use[order[high]]].mean() - amls['OS'][data_to_use[order[low]]].mean()
        #sd[genes] = [amls['OS'][data_to_use[order[high]]].mean() , amls['OS'][data_to_use[order[low]]].mean()]

        if numpy.abs(day_diff ) > 50:
            print 'possible target ..'
            try:
                p = plot_survival(amls, numpy.concatenate([data_to_use[order[high]],data_to_use[order[low]]]) , [21,22], no_plot=no_plot)
                print 'Ok. one Pvalue'
            except Exception, e :
                p=1.
                print e

        else:
            p=1.
        if plot_anyways:
        	plot_survival(amls, numpy.concatenate([data_to_use[order[high]],data_to_use[order[low]]]) , [21,22], no_plot=False)
        pvals.append([p,genes])
        if (no_plot == False):
            try:
                os.rename('OS_fc_color_many_probes.pdf','OS_fc_color_many_probes_%s.pdf'%genes)
                print 'this is a test'
            except Exception, e :
                print e
    return pvals
    pass


def expand_labels(sample,label=None,n=2):
    label_i = numpy.argwhere(sample['colnames']== label).ravel()
    indexes = []
    for i in range(sample['data'].shape[1]):
        if i in label_i:
            for j in range(n):
                indexes.append(i)
        else:
            indexes.append(i)
    sample['data'] = sample['data'][:,indexes]
    sample['colnames'] = sample['colnames'][indexes]
    return sample
    pass

def mannwhitneyu(x, y, sd, use_continuity=False):
    """Computes the Mann-Whitney rank test on samples x and y.


    Parameters
    ----------
        x : array_like 1d
        y : array_like 1d
        use_continuity : {True, False} optional, default True
            Whether a continuity correction (1/2.) should be taken into account.

    Returns
    -------
        u : float
            The Mann-Whitney statistics
        prob : float
            one-sided p-value assuming a asymptotic normal distribution.

    Notes
    -----
    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is
    significant if the u-obtained is LESS THAN or equal to the critical
    value of U.

    This test corrects for ties and by default uses a continuity correction.
    The reported p-value is for a one-sided hypothesis, to get the two-sided
    p-value multiply the returned p-value by 2.

    """
    x = stats.asarray(x)
    y = stats.asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = stats.rankdata(numpy.concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    #ranky = ranked[n1:]        # the rest are y-ranks

    u1 = n1*n2 + (n1*(n1+1))/2.0 - numpy.sum(rankx,axis=0)  # calc U for x
    u2 = n1*n2 - u1                            # remainder is U for y
    bigu = max(u1,u2)
    smallu = min(u1,u2)
#   #T = np.sqrt(tiecorrect(ranked))  # correction factor for tied scores
#   T = stats.tiecorrect(ranked)
#   if T == 0:
#       raise ValueError, 'All numbers are identical in amannwhitneyu'
    T=1.0
#   sd = numpy.sqrt(T*n1*n2*(n1+n2+1)/12.0)

    if use_continuity:
        # normal approximation for prob calc with continuity correction
        z = (bigu-0.5-n1*n2/2.0) / sd
    else:
        z = (bigu-n1*n2/2.0) / sd  # normal approximation for prob calc
    return smallu, z
    #, stats.distributions.norm.sf(abs((z))) # (1.0 - stats.zprob(z))
    pass


def mannwhitneyu_1(x, y, sd):
    """Another implementation of man withney...
    """
    x = stats.asarray(x)
    y = stats.asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = stats.rankdata(numpy.concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    #ranky = ranked[n1:]        # the rest are y-ranks
    u1  = numpy.sum(rankx,axis=0) - n1*(n1+1)/2
    m_U = n1*n2/2
#   sd = numpy.sqrt((n1*n2 * (n1 + n2 + 1))/12)
    z  = (u1 - m_U)/sd
    return u1, z
    pass



def read_signatures(file_path):
    """docstring for read_signatures:
    import signatures (gene sets) from text file.
    structure is: Name \t na \t gene1 \t ...\t gene_n
    returns a dictionary.
    """
    dic={}
#   t=open(file_path).read()

    for line in open(file_path).readlines():
        line_split=line.split('\t')
#       print line_split
        name=line_split[0]
        #striped_name=name.strip('_UP').strip('_DN')
        gene_set=line_split[2:-1]
        dic[name] = numpy.array([i for i in gene_set if i != ''])

        #if not(striped_name in dic):
        #   dic[striped_name]={}
        #if name.find('_UP')!=-1: #there is _up somewhere
        #   dic[striped_name]['up']=gene_set
        #elif name.find('_DN')!=-1:
        #   dic[striped_name]['dn']=gene_set
        #else:
        #   dic[striped_name]['both']=gene_set

    return dic
    pass



def quantile(x, q, qtype = 7, issorted = False):
    """ Args:
     x - input data
     q - quantile
     qtype - algorithm
    issorted- True if x already sorted.
    Compute quantiles from input array x given q.For median, specify q=0.5.
     References: http://reference.wolfram.com/mathematica/ref/Quantile.html http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile
    Author: Ernesto P.Adorio Ph.D. UP
    Extension Program in Pampanga, Clark Field.
    """
    from math import modf, floor

    if not issorted:
        y = sorted(x)
    else: y = x

    if not (1 <= qtype <= 9):
        return None # error!
        # Parameters for the Hyndman and Fan algorithm
    abcd = [(0, 0, 1, 0), # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
            (0, 0, 0, 1), # California linear interpolation, R type 4
            (0.5, 0, 0, 1), # hydrologists method, R type 5
            (0, 1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1, -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
            (1.0/3, 1.0/3, 0, 1), # median-unbiased , R type 8
            (3/8.0, 0.25, 0, 1) # normal-unbiased, R type 9.
            ]

    a, b, c, d = abcd[qtype-1]
    n = len(x)
    g, j = modf( a + (n+b) * q -1)
    if j < 0:
        return y[0]
    elif j > n:
        return y[n]

    j = int(floor(j))
    if g == 0:
        return y[j]
    else:
        return y[j] + (y[j+1]- y[j])* (c + d * g)






def AGC(x , copy=0, c=0.2, mode='log_transform'):
    """
    docstring for AGC
    Performs Array generation based gene centering, based on the paper from Jaakko Astola.

    Input:

    - c parameter from the AGC algorithm. (automatic value is 1/5)
    - a dictionarry where keys are data, row, names, colnames, platforms.
        row names are probes labels [1xm]
        col names are samples labels [1xn]
        data is the data [m x n]
        platform is the platform [1xn]

    Returns:
     The centered data. (you did guess that huh?)
    """
    #find all platforms.
    platform_list=[]
    for p in x['platforms']:
        if p not in platform_list:
            platform_list.append(p)
    #verify that there are no negative values in the expression data:
    if (numpy.min(x['data']) <0 ):
        print 'Bad. Really bad.. I think you have centered your data or done something..'
        print 'Found negative values.'
        print 'Carrying on with analysis anyways.'

    #log2 transform the data.
    if mode=='log_transform':
        log2x=numpy.log2(x['data'])
    elif mode=='exp_transform':
        log2x=numpy.exp2(x['data'])
    else:
        log2x=x['data']

    #get median
    mu_k=numpy.median(log2x, axis=1)


    #make matrices, median a_k and b_k for samples from same platform:
    sep_data={}
    for platform in platform_list:
        sep_data[platform] = {}
        data_to_use=[]
        print 'Extracting %s:\n' % platform
        print
        for i,pl in enumerate( x['platforms'] ):
#           print pl , platform
            if pl == platform :
                data_to_use.append(i)
        l=len(data_to_use)
        sep_data[platform]['colnames'] = x['colnames'][data_to_use]                 #samples names.
        sep_data[platform]['cel_file'] = x['cel_file'][data_to_use]
        sep_data[platform]['stemness'] = x['stemness'][data_to_use]

        sep_data[platform]['mu_i_k'] = numpy.median(log2x[:,data_to_use], axis=1)   #probe mean across platfrom i.
        mu_i_k=sep_data[platform]['mu_i_k']
        mu_i_k_mat=numpy.array(l*[mu_i_k]).transpose()
        mu_k_mat= numpy.array(l*[mu_k]).transpose()
        sep_data[platform]['data'] = log2x[:,data_to_use] - (mu_i_k_mat -mu_k_mat  ) #data, centered, but not yet ajusted.
        #a_k and b_k require more than one line to be calculated...
        #sort the gene expression data gene wise. Used to get 2% high and low in a supid way.. a smart way would be to fit the data to a norm distribution. and get the tail.
        print 'It would be very smart to fit a normal distribution here in the code.'
        slog2x=numpy.sort(log2x[:,data_to_use], axis=1)
        if l < 100: # use the mean of the lowest and highest two values, there is not a lot of data...
            print 'working with few data... be careful.'
            print 'estimating 2% quantile, done, with MINiML.quantile()'
#           sep_data[platform]['a_k'] = numpy.mean( slog2x[:,[0,1,2] ] ,axis=1)
#           sep_data[platform]['b_k'] = numpy.mean( slog2x[:,[-1,-2,-3] ] ,axis=1)

            sep_data[platform]['a_k'] = numpy.array([quantile(slog2x[i,:], 0.02 ) for i in range(len(slog2x[:,0]))])

            print 'estimating 98% quantile, done, with MINiML.quantile()'
            sep_data[platform]['b_k'] = numpy.array([quantile(slog2x[i,:], 0.98 ) for i in range(len(slog2x[:,0]))])

        else:
            print 'To be implemented!' #but still we want the function to work.
            sep_data[platform]['a_k'] = numpy.mean( slog2x[:,[0] ] ,axis=1)
            sep_data[platform]['b_k'] = numpy.mean( slog2x[:,[-1] ] ,axis=1)

        # adjust Values
        prog = ProgressBar(0, l , 20, mode='fixed')
        for i in range(l):
            for j in range(len(x['rownames'])):
                if sep_data[platform]['data'][j,i] > sep_data[platform]['b_k'][j]:
                    sep_data[platform]['data'][j,i] = c*(sep_data[platform]['data'][j,i] - sep_data[platform]['b_k'][j]) + sep_data[platform]['b_k'][j]
                if sep_data[platform]['data'][j,i] < sep_data[platform]['a_k'][j]:
                    sep_data[platform]['data'][j,i] =  sep_data[platform]['a_k'][j] - c*(sep_data[platform]['a_k'][j] - sep_data[platform]['data'][j,i] )
            oldprog = str(prog)
            prog.update_amount(i+1)
            if oldprog != str(prog):
                print prog, "\r",
                sys.stdout.flush()
                oldprog=str(prog)

        #convert back values:
        #log2 transform the data.
        if mode=='log_transform':
            sep_data[platform]['data']=numpy.exp2(sep_data[platform]['data'])
        if mode=='exp_transform':
            sep_data[platform]['data']=numpy.log2(sep_data[platform]['data'])

    #merge back for shipping!
    y={}
    y['data']=numpy.concatenate(([sep_data[i]['data'] for i in  sep_data.keys()]),axis=1)
    y['colnames']=numpy.concatenate([sep_data[i]['colnames'] for i in  sep_data.keys()])
    y['platforms']=numpy.concatenate([ [i]*len(sep_data[i]['colnames']) for i in  sep_data.keys()])


    #Sort order of cel files according to stemness: (this is already done in the excel file)
#   [(data.index(i),i[5],i[11],i[1]) for i in sorted(data, key=lambda i: int(i[11] )  )]
# this is why people should program in C!
    final_colnames=numpy.concatenate([sep_data[i]['colnames'] for i in  sep_data.keys()])
    final_stemness=numpy.concatenate([sep_data[i]['stemness'] for i in  sep_data.keys()])
    final_celfile=numpy.concatenate([sep_data[i]['cel_file'] for i in  sep_data.keys()])
    xx=numpy.array([final_colnames,final_stemness,final_celfile]).transpose().tolist()
#   [(x.index(i),i[0],i[1]) for i in sorted( x , key=lambda i:( int(i[1]),i[0],i[2] ) )]
    order=[xx.index(i) for i in sorted( xx , key=lambda i:( int(i[1]),i[0],i[2] ) )]

    y['cel_file']=numpy.array(final_celfile)[order]
    y['stemness']=numpy.array(final_stemness)[order]
    y['platforms']=y['platforms'][order]
    y['colnames']=numpy.array(final_colnames)[order]
    y['data']=numpy.copy(y['data'][:,order])
    y['rownames']=x['rownames']


    return y
    pass


def ComBat(data):
    """docstring for Combat
        Uses Johnson, WE, ComBat method to normalize batch effect.

        Reference: Johnson, WE, Rabinovic, A, and Li, C (2007).
        Adjusting batch effects in microarray expression data using Empirical Bayes methods.
        Biostatistics 8(1):118-127.

    """

    import pandas as pa
    import patsy
    import combat
    import numpy as np

    x1cat = []
    for i in data['X1']:
        if i not in x1cat:
            x1cat.append(i)
    x1=[]
    for i in data['X1']:
        x1.append(x1cat.index(i))

    data['X1'] = x1#data['X1'].tolist()
    #data['X1'].pop(-1)
    #data['X1'].append('19')
    data['X1'] = numpy.array(data['X1'], dtype='i')

    data['X2'] = x1#data['X2'].tolist()
    #data['X2'].pop(-1)
    #data['X2'].append('19')
    data['X2'] = numpy.array(data['X2'], dtype='i')

    person_lab=[]
    for i in data['cel_file']:
        if i[0] == 'H':
            #print 'new'
            person_lab.append(1)
        elif i[0] == 'J':
            person_lab.append(2)
            #print 'new2'
        else:
            person_lab.append(0)
            #print 'old'
    dat={}
    dat['cel'] = [ i.replace('.CEL.gz','') for i in data['cel_file']]
    #dat['Sample'] = data['colnames']
    dat['stemness'] =   numpy.array(data['stemness'], dtype='i')
    #print dat['stemness']
    dat['date'] = numpy.array(data['X1'], dtype='i')
    dat['person'] = numpy.array(person_lab, dtype='i')
    #dat['batch'] = numpy.array(dat['date'],dtype ='i')

    df_h = pa.DataFrame(dat)
    df_d=pa.DataFrame(data['data'])
    #df_d.columns = dat['cel']
    #df_h = df_h.T
    #df_h.rows = dat['cel']


    #mod_ = patsy.dmatrix("~ stemness ", df_h, return_type="dataframe")
    mod_ = patsy.dmatrix(" ~ stemness", df_h, return_type="dataframe")
    import time
    t = time.time()
    #ebat = combat.combat(df_d,dat['date'], mod_ ,)
    #print np.array(ebat)
    ebat = combat.combat(df_d,dat['person'], mod_ ,)

    #merged_data_before_AGC['data'] = array(ebat)
    #T,p,e=MINiML.run_PCA(merged_data_before_AGC, arange(len(merged_data_before_AGC['colnames']))[:40], select_probes_by_variance(merged_data_before_AGC, var_thr=2, data_to_use=arange(len(merged_data_before_AGC['colnames']))[:40]), fit_ellipse=0,Text_Labels=0)
    #dat['batch'] = dat['platform']
    #df_h = pa.DataFrame(dat)
    #mod_ = patsy.dmatrix("~ stemness ", df_h, return_type="dataframe")
    #ebat = combat.combat(ebat, df_h.batch, mod_,["stemness" ])

    print "%.2f seconds" % (time.time() - t)
    return numpy.array(ebat)

    pass
#

def ComBatR(data):
    """docstring for Combat
        Uses Johnson, WE, ComBat method to normalize batch effect.


        Reference: Johnson, WE, Rabinovic, A, and Li, C (2007).
        Adjusting batch effects in microarray expression data using Empirical Bayes methods.
        Biostatistics 8(1):118-127.

    """
    data['X1'] = data['X1'].tolist()
    data['X1'].pop(-1)
    data['X1'].append('19')
    data['X1'] = numpy.array(data['X1'], dtype='i')

    data['X2'] = data['X2'].tolist()
    data['X2'].pop(-1)
    data['X2'].append('19')
    data['X2'] = numpy.array(data['X2'], dtype='i')


    os.mkdir('_combat_%s'%os.getpid())
    os.chdir('_combat_%s'%os.getpid())

    sep='\t'
    row,col=data['data'].shape
    header_file='h%s.txt'%os.getpid()
    data_file = 'd%s.txt'%os.getpid()

    #Writes header file:
    h=open(header_file,'w')
#   h.write('Array Name\tSample Name\tBatch\tstemness\n')
#   for i in range(col):
#       h.write(str(i)+sep+data['colnames'][i]+sep+data['platforms'][i]+data['stemness'][i]+sep+data['stemness'][i]+'\n')
    if  data.has_key('Batch') :
#       h.write('Array Name\tSample_Name\tBatch\tstemness\n')
#       for i in range(col):
#           h.write(str(i)+sep+data['colnames'][i]+sep+str(data['Batch'][i])+sep+str(data['stemness'][i])+'\n')
        h.write('Array Name\tSample_Name\tBatch\n')
        for i in range(col):
            h.write(str(i)+sep+data['colnames'][i]+sep+str(data['Batch'][i])+'\n')


    else:
#       h.write('Array Name\tSample_Name\tBatch\tStemness\tX1\n')
        h.write('Array Name\tSample_Name\tBatch\tStemness\n')
        for i in range(col):
#           h.write(str(i)+sep+ str(data['colnames'][i])+sep+ str(data['stemness'][i])+sep+str(data['platforms'][i])+sep+str(data['X1'][i])+'\n')
            h.write(str(i)+sep+ str(data['colnames'][i])+sep+ str(data['platforms'][i])+sep+str(data['stemness'][i])+'\n')


    #Writes data file
    d=open(data_file,'w')
    for i in range(col):
        d.write(str(i)+sep)
    d.write('\n')
    for i in range(row):
        for j in range(col):
            d.write(str(data['data'][i][j]) + sep)
#           print data['data'][i][j]
        d.write('\n')
    h.close()
    d.close()

    ComBat = importr('ComBat')

    m = numpy.array (  ComBat.Combat('d%s.txt'%os.getpid(),'h%s.txt'%os.getpid()) )
#   produces: Adjusted_d.txt_.xls

### f=open('Adjusted_'+data_file+'_.xls')
### t=f.read()
### f.close()
### data_=[]
### i=0
### for sample in t.splitlines():
###     #skip first line
###     if i==0:
###         i+=1
###     else:
###         data_.append( sample.split('\t') )
###         for j,field in enumerate(data_[-1]):
###             data_[-1][j] = field.strip() #clears white spaces before and after..
###         i+=1
### print numpy.array(data_)
### m=numpy.zeros([row,col])
### for i in range(row):
###     for j in range(col):
###         m[i][j] = data_[i][j]
###
###
    os.chdir('..')
    import shutil
    shutil.rmtree('_combat_%s'%os.getpid())
    return m.T

    pass

#



def extract_data(data,index):
    """docstring for extract_data"""
    new = {}
    no_chamge=['chimera_sample','rownames', 'low_n', 'raw_fold_change', 'top_n','voisins']
    for key in data.keys():
        print key
        try:
            ndim=data[key].ndim
            if key not in no_chamge :
                if ndim==1:
                    new[key] = data[key][index]
                else:
                    new[key] = data[key][:,index]
            else:
                new[key] = data[key]
        except Exception, e:
            print e
    return new
    pass

def export_values_MIC(data,probes_to_use=None,data_to_use=None,out='mic_dat.csv' ):
    """docstring for export_values_MIC
        this generate the text file needed by MINE.
    """
    if data_to_use == None:
        data_to_use = range(data['data'].shape[1])
    if probes_to_use == None:
        probes_to_use = range(data['data'].shape[0])

    f=open(out,'w')
    for probe in probes_to_use:
        f.write('%s,'%data['rownames'][probe])
        f.write('%s\n'%  ','.join( numpy.array(data['data'][probe,:][:,data_to_use], dtype='str')))
    f.close()
    pass

def run_MINE(master_var_index):
    """docstring for run_MINE
        this runs the java program with the correct paremeters.
        uses popen. works in parallel.
    """
    import subprocess as sp
    p1 = sp.Popen(["java", "-Xmx25g", "-jar", "MINE.jar", "mic_dat.csv", "-masterVariable", str(master_var_index), "c=25"], stdout = sp.PIPE)
    results=p1.stdout.read()
    print results
    return master_var_index
    pass

def MINE(data,probes_to_use=None,data_to_use=None,ncpu=4, done_file = None):
    """docstring for MINE
        this is an attempt to run the MINE software, on many cores, in masterVariable mode so that
        we get a lot of files, and it hopefully runs.
    """
    from multiprocessing import Pool
    import shutil, platform
    if data_to_use == None:
        data_to_use = range(data['data'].shape[1])
    if probes_to_use == None:
        probes_to_use = range(data['data'].shape[0])

    os.mkdir('_MINE_%s'%os.getpid())
    os.chdir('_MINE_%s'%os.getpid())
    if platform.uname()[1] == "SGI (Porse)":
        shutil.copy('/home/nrapin/Dropbox/Python-stuff/MINE.jar', '.')
    else:
        shutil.copy('/Users/nrapin/Dropbox/Python-stuff/MINE.jar', '.')
    export_values_MIC(data,probes_to_use=probes_to_use,data_to_use=data_to_use,out='mic_dat.csv')

    pool = Pool(processes=ncpu)
    jobs_list=None
    pre_jobs_list = range(len( probes_to_use ))
    if done_file != None:
        if os.path.exists(done_file):
            done = numpy.array(numpy.array(import_sample_data(done_file, line_to_skip=1)).reshape(-1), dtype = int)
            jobs_list=[]
            for i in pre_jobs_list:
                if i not in done:
                    jobs_list.append(i)
    else:
        jobs_list=range(len( probes_to_use ))
    results=[]
    results=pool.map(run_MINE, jobs_list)


    os.chdir('..')

    pass

def export_values_MIC_cor(data,probes_to_use=None,data_to_use=None,out='mic_dat.csv' ):
    """docstring for export_values_MIC
        this generate the text file needed by MINE.
    """
    if data_to_use == None:
        data_to_use = range(data['data'].shape[1])
    if probes_to_use == None:
        probes_to_use = range(data['data'].shape[0])

    f=open(out,'w')
    for sample in data_to_use:
        f.write('%s,'%data['colnames'][sample])
        f.write('%s\n'%  ','.join( numpy.array(data['data'][probes_to_use,:][:,sample], dtype='str')))
    f.close()
    pass

def run_MINE_cor(master_var_index):
    """docstring for run_MINE
        this runs the java program with the correct paremeters.
        uses popen. works in parallel.
    """
    import subprocess as sp
    p1 = sp.Popen(["java", "-Xmx25g", "-jar", "MINE.jar", "mic_dat.csv", "-masterVariable", str(master_var_index), "c=25"], stdout = sp.PIPE)
    results=p1.stdout.read()
    print results
    return master_var_index
    pass


def MINE_cor(data,probes_to_use=None,data_to_use=None,ncpu=4, done_file = None):
    from multiprocessing import Pool
    import shutil, platform
    if data_to_use == None:
        data_to_use = range(data['data'].shape[1])
    if probes_to_use == None:
        probes_to_use = range(data['data'].shape[0])

    os.mkdir('_MINE_%s'%os.getpid())
    os.chdir('_MINE_%s'%os.getpid())
    if platform.uname()[1] == "SGI (Porse)":
        shutil.copy('/home/nrapin/Dropbox/Python-stuff/MINE.jar', '.')
    else:
        shutil.copy('/Users/nrapin/Dropbox/Python-stuff/MINE.jar', '.')
    export_values_MIC_cor(data,probes_to_use=probes_to_use,data_to_use=data_to_use,out='mic_dat.csv')

    pool = Pool(processes=ncpu)
    jobs_list=None
    pre_jobs_list = range(len( data_to_use ))
    if done_file != None:
        if os.path.exists(done_file):
            done = numpy.array(numpy.array(import_sample_data(done_file, line_to_skip=1)).reshape(-1), dtype = int)
            jobs_list=[]
            for i in pre_jobs_list:
                if i not in done:
                    jobs_list.append(i)
    else:
        jobs_list=range(len( data_to_use ))
    results=[]
    results=pool.map(run_MINE, jobs_list)


    os.chdir('..')

    pass


def export_genes_for_bonjo(all_data,probes_to_use,data_to_use=None,out='bongo_dat.txt', discretize=True):
    """docstring for export_genes_for_bonjo
        generates the data file for dynamic network generation in bonjo.


    """
    import collections

    if data_to_use == None:
        data_to_use = range(all_data['data'].shape[1])

    #first we get the probe names as gnenes.
    AnnotationDbi=importr('AnnotationDbi')
    hg133=importr('hgu133plus2.db')

    p2go=AnnotationDbi.as_list(hg133.hgu133plus2SYMBOL)

    gene_to_go = collections.defaultdict(list)
    go_to_gene = collections.defaultdict(list)
    p2go_names=numpy.array(p2go.names)
    for i,names in enumerate(p2go_names):
            try:
                for go in  p2go[i]:
                    go_term=str(go)
                    gene_to_go[names].append(go_term)
                    go_to_gene[go_term].append(names)

            except Exception,e:
                pass
    f=open(out,'w')
    for probe in probes_to_use:
        in_gene=gene_to_go[all_data['rownames'][probe]]
        if in_gene[0]=='NA':
            in_gene=[all_data['rownames'][probe]]
        print in_gene[0]

        f.write(in_gene[0]+'\t')
    f.write('\n')

    if discretize:
        matrix=numpy.zeros(all_data['data'].shape)
        for sample in data_to_use:
            for probe in probes_to_use:
                val=all_data['data'][probe,sample]
                if val < -4:
                    matrix[probe,sample] = int(0)
                elif val < 0:
                    matrix[probe,sample] = int(1)
                elif val < 4:
                    matrix[probe,sample] = int(2)
                else :
                    matrix[probe,sample] = int(3)

    for sample in data_to_use:
        for probe in probes_to_use:
            if discretize:
                f.write('%d\t'  % matrix[probe,sample])
            else:
                f.write('%5e\t'  % all_data['data'][probe,sample])

        f.write('\n')
    f.close()
    pass


def merge_datasets_on_highest(data1,data2,do_batch_correct=True,plots=False):
    from numpy import array, concatenate, argwhere , argmax
    overlap =  len(set(data1['rownames']).intersection(data2['rownames']))
    if overlap > 0 :
        print 'Found %d overlaping genes.'%(overlap)
    else:
        print 'Found no ovelap...'
        return

    all_genes  = array(list(set(data1['rownames']).intersection(data2['rownames'])))
    all_genes.sort()
    merged_data = {}
    merged_data['rownames'] = all_genes
    merged_data['colnames'] = concatenate([data1['colnames'],data2['colnames']])
    merged_data['data'] = []
    # find max value of genes in both datasets:

    prog = ProgressBar(0, len(all_genes) , 50, mode='fixed')
    for index,g in enumerate(all_genes):
        #print g
        data_gene_1 = data1['data'][argwhere(data1['rownames'] == g).ravel(),:]
        maxi_1 = argmax(abs(array(data_gene_1.max(axis=1))))
        data_gene_2 = data2['data'][argwhere(data2['rownames'] == g).ravel(),:]
        maxi_2 = argmax(abs(array(data_gene_2.max(axis=1))))
        merged_data['data'].append(concatenate([  data_gene_1[maxi_1], data_gene_2[maxi_2]  ]))

        oldprog = str(prog)
        prog.update_amount(index)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)

    merged_data['data'] =  array(merged_data['data'],dtype='f' )
    merged_data['stemness'] = stemness_by_class(merged_data)
    merged_data['stemness'][...] =0
    merged_data['stemness'][0:len(data1['colnames'])] =1
    p2u = select_probes_by_variance(merged_data,var_thr = 3  )
    if p2u == []:
        p2u = range(20)
    if plots:
        tsne_paper(merged_data,probes=p2u, out='no_batch correction.pdf',perplexity=3)
        plot_exprs_histo(merged_data)
        plot_cluster(merged_data,probes_to_use=p2u)
    if do_batch_correct:
        merged_data['data'] = batch_correct(merged_data, batch_key='stemness')
        p2u = select_probes_by_variance(merged_data,var_thr = 2  )
        if p2u == []:
            p2u = range(20)
        merged_data['stemness'] = stemness_by_class(merged_data)
        PCA_paper(merged_data,probes=p2u, out='batch correction.pdf')
        return merged_data
    else:
        return merged_data
    pass
def probes_to_something(p,platform,something):

    """docstring for probes_to_genes

    platform can be :
        mogene10sttranscriptcluster.db
        hugene10sttranscriptcluster.db
        ...

    something is :
        SYMBOL
        REFSEQ
        ...
    """

    robjects.r("""
  human_probe2symb <- function(probes){
    if (!"%s" %%in%% installed.packages()) {source("http://bioconductor.org/biocLite.R"); biocLite("%s")}
    library(%s);
    keys=keys(%s) # get it all
    symb=select(%s, keys, "%s", "PROBEID" )
    symb=symb[!duplicated(symb[,1]),]
    return(symb[match(probes, symb[,1]),][,2])
  }
    """%(platform,platform,platform,platform,platform,something)
    )
    r_human_probe2symb =robjects.globalenv['human_probe2symb']

    genes = r_human_probe2symb(p)
    return genes
    pass


def make_hemaexplorer_app(data,annot,name='BloodSpot'):
    """
    create a hemaexplorer java app from a dataset.

    annot is a dictionnary .... of genes:probes1 probes 2 etc..
    """
    # make a copy of the template program.
    #(located in /Users/nrapin/Dropbox/BRIC/other_projects/HemaExplorer/HemaExplorer-src/HemaExplorer_js_with_text/application.macosx/)
    import shutil,os
    initial_dir = os.path.abspath('.')
    shutil.copytree('/Users/nrapin/Dropbox/BRIC/other_projects/HemaExplorer/HemaExplorer-src/HemaExplorer_js_with_text/application.macosx/',\
                       './'+name )
    os.chdir('./'+name + '/HemaExplorer_js_with_text.app/Contents/Java/data')

    # then make the classes and annotation files

    f = open('classes_h.txt','w')
    f.write('\t'.join(data['colnames']))
    f.close()

    f = open('list_h.txt','w')
    for p in annot:
        f.write(p+' '+' '.join(annot[p])+'\n')
    f.write('start start\n')
    f.close()

    #finally export the data in the right folder.
    os.chdir('exp_data')
    prog = ProgressBar(0, len(data['rownames']) , 50, mode='fixed')
    for i, p in enumerate(data['rownames']):
        f=open(p.replace('/','_').upper()+'.txt','w')
        f.write( ' '.join(numpy.array(data['data'][i,:],dtype='str').tolist()) )
        f.close()
        oldprog = str(prog)
        prog.update_amount(i)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)
    f=open('start.txt','w')
    f.write( ' '.join(numpy.array(data['data'][i,:],dtype='str').tolist()) )
    f.close()

    #launch program
    os.chdir(initial_dir)

    pass



def batch_correct(data, batch_key = 'stemness',covariate=None):
    import combat
    import pandas as pa

    dat1 = pa.DataFrame({i:data['data'][:,i] for i in range(len(data['colnames'])) } , index=data['rownames'])
    if covariate== None:
        pheno1 = pa.DataFrame({'batch':data[batch_key]})
        return numpy.array(combat.combat(dat1,pheno1.batch,False))
    else:
        import patsy
        pheno1 = pa.DataFrame({'batch':data[batch_key],'covariate':covariate})
        mod = patsy.dmatrix(" ~ covariate", pheno1, return_type="dataframe")
        print mod
        ebat = numpy.array(combat.combat(dat1,pheno1.batch,mod))
        return ebat
    pass


def export_as_txt(data,filename='SM.txt',sep='\t'):
    """docstring for export_values_as_txt
        utils=importr('utils')
        affy2gene=gene2affy={}
    #   affyId2geneSymbol=robjects.r('toTable(hgu133plus2SYMBOL)')
        affyId2geneSymbol=numpy.array(robjects.r.toTable(hgu133plus2.hgu133plus2SYMBOL))
        for i in range(affyId2geneSymbol.shape[1]):
            affy2gene[affyId2geneSymbol[1,i]]=affyId2geneSymbol[0,i]
            gene2affy[affyId2geneSymbol[0,i]]=affyId2geneSymbol[1,i]
        converted_rownames=[]
        for i in merged_data['rownames']:
            try:
                converted_rownames.append(affy2gene[i])
            except Exception,e:
                converted_rownames.append(i)
    """
    row,col=data['data'].shape

    keys=['cel_file']

    f=open(filename,'w')
    f.write('\t'.join(data['cel_file']))
    f.write('\n')

    labels = data['cel_file'].tolist()
    labels.append("")
    labels=numpy.array(labels)
    numpy.c_[data['rownames'],data['data']]
    mat_tmp =  numpy.r_[ labels,numpy.c_[data['rownames'],data['data']]]

#    for i in range(row):
#        f.write(data['rownames'][i] + sep)
#        f.write('\t'.join(numpy.array(data['data'][i],dtype='str')))
#        f.write('\n')

    f.close()
    pass


def export_values_as_txt(data,filename='SM.txt',sep='\t',probes=None, log2=True):
    """docstring for export_values_as_txt
        utils=importr('utils')
        affy2gene=gene2affy={}
    #   affyId2geneSymbol=robjects.r('toTable(hgu133plus2SYMBOL)')
        affyId2geneSymbol=numpy.array(robjects.r.toTable(hgu133plus2.hgu133plus2SYMBOL))
        for i in range(affyId2geneSymbol.shape[1]):
            affy2gene[affyId2geneSymbol[1,i]]=affyId2geneSymbol[0,i]
            gene2affy[affyId2geneSymbol[0,i]]=affyId2geneSymbol[1,i]
        converted_rownames=[]
        for i in merged_data['rownames']:
            try:
                converted_rownames.append(affy2gene[i])
            except Exception,e:
                converted_rownames.append(i)
    """
    if probes!=None:
        data['data']= data['data'][probes,:]

    keys=['colnames','platforms','stemness','cel_file']
    keys=['colnames']
    keys=['colnames','cel_file']

    if data.has_key('OS'):
        keys.append('OS')
        keys.append('EventOS')
    f=open(filename,'w')
    row,col=data['data'].shape



    for key in keys:
        if data.has_key(key):
            f.write(str(key))
            f.write(sep)
            for i in range(col):
                f.write('%s'%str(data[key][i]))
                f.write(sep)
            f.write('\n')

    for i in range(row):
        f.write(data['rownames'][i] + sep)
        for j in range(col):
            if log2:
                f.write(str(data['data'][i][j]) + sep)
            else:
                f.write(str( numpy.exp2(data['data'][i][j] )) + sep)
        f.write('\n')

    f.close()
    pass
#
def export_values_for_gsea(data,filename='GSEA.txt',sep='\t',probes=None, log2=True):
    """export for GSEA -
    """

    if probes!=None:
        data['data']= data['data'][probes,:]

    keys=['colnames','platforms','stemness','cel_file']
#   keys=['colnames']
    f=open(filename,'w')
    row,col=data['data'].shape

    for key in keys:
        if data.has_key(key):
            f.write('Name'+sep+'Description')
            f.write(sep)
            for i in range(col):
                f.write('%s%d'% (str(data[key][i]) , i) )
                f.write(sep)
            f.write('\n')

    print row
    for i in range(row):
        f.write(data['rownames'][i] + sep + 'NA' + sep)
        for j in range(col):
            if log2:
                f.write(str(data['data'][i][j]) + sep)
            else:
                f.write(str( numpy.exp2(data['data'][i][j] )) + sep)
        f.write('\n')

    f.close()

    data_colnames=data['colnames']
    data_colnames = [i.replace('(', '').replace(')', '').replace(';', '') for i in data_colnames]
    cat=[]
    for i in data_colnames:
     if i not in cat:
        cat.append(i)
    f=open(filename+'.cls','w')
    f.write('%d %d 1\n'%(len(data_colnames),len(cat)))
    f.write('# ')
    for i in cat:
        f.write('%s '%i)
    f.write('\n')
    for i in data_colnames:
        f.write('%s '%i)
    f.write('\n')


    pass

def split_array(x,indexes):
    """docstring for split_array
        x=[n1, n2, n3, ..., nn]
        split_array(x, [n1,n3,nj, .. ,nk]) -> [ [n1,n3,nj, .. ,nk] , [n2,n4...] ]
    """
#   x=x.tolist()
    y=[]
    z=[]
    for i,j in enumerate(x):
        if i in indexes:
            y.append(j)
        else:
            z.append(j)

    return [y,z]
    pass

def make_clusters(linkage ,leaves_labels=None,  depth=2, cur_depth=0 , return_list=None):
    """docstring for make_clusters
    get a linkage and iterativelly separate the cluster by two, returning a list of leaves in each cluster, and their color

    the idea:
    make the dedrograme, and iterativelly find a threshold where you get only two colors.
    this gives you the first separation.
    rerun on each group, and repeat procedure until you reach depth.

    then, you get groups, and each group has a color.

    """
    import hcluster
    topGO=import_rlib(['topGO'])
    if leaves_labels==None:
        print 'Please provide Leaves labels.'
        return
    #create list at the begining.
    if return_list==None:
        return_list={}
    #get 2 clusters:
    while depth > 0 :
        thr=0
        c=hcluster.fcluster(b,2,criterion='maxclust')
        c=split_array(c, [i for i,j in enumerate(c) if j==1])
        depth-=1


    pass

def find_probes_DBI(gene,organism):
    """docstring for find_probes
        loads the annotation DBI R package converted by the generate_annotation_dics.py script.
        loops through all possible annotations, and the fuction hopefully returns the correct probes ;)
    """
    path_to_use = '/Users/nrapin/BRIC/databases'
    gene = gene.lower()
    if 0: #slow method.
        p=None
        if organism == 'human':
            annot_to_probe = pickle_load(path_to_use + '/h_annot_to_probe.pkl')
        else:
            annot_to_probe = pickle_load(path_to_use + '/m_annot_to_probe.pkl')


        for k in annot_to_probe:
            p = annot_to_probe[k].get(gene)
            if p != None:
                print 'found in '+k
                break

    else: #fast method
        if organism == 'human':
            r=h5py.File(path_to_use + '/h_annot_to_probe.hdf5','r')
        else:
            r=h5py.File(path_to_use + '/m_annot_to_probe.hdf5','r')

        root= r['root']
        p=None
        for k in root:
            if gene in root[k]:
                p = root[k].get(gene)[...]
            if p != None:
                break
        r.close()
    return p
    pass



def convert_affy_2_genes_in_data(data,array_type=None,compute_fc=True):
    """docstring for convert_affy_2_genes_in_data
        changes data['rownaes'] form affy to genenames so that they can be used in the msig database.
        leaves unknown genes from the data as affy probes.
        also, array_type can be used for other than humans arrays.
        use the R annotation dbi package name.
    """
    import copy
    import collections

    copydata={}

    if not data.has_key('rownames'):
        print 'missing rownames...'
        return

#   now reconstruct the matrix
#   And rownames , with gene names.
    try:
    	0 - 'a'
#   	    AnnotationDbi=importr('AnnotationDbi')
#   	    if array_type == None:
#   	        hg133=importr('hgu133plus2.db')
#   	        p2go=AnnotationDbi.as_list(hg133.hgu133plus2SYMBOL)
#
#   	    else:
#   	        hg133=importr(array_type)
#   	        name = array_type.split('.db')[0]
#   	        if name[-1] == '2':
#   	            db_name = 'hg133'+'.'+name[0:-1]+'2SYMBOL'
#   	        else:
#   	            try:
#   	                db_name = 'hg133'+'.'+name+'2SYMBOL'
#   	                p2go=eval('AnnotationDbi.as_list(%s)'%db_name)
#   	            except :
#   	                db_name = 'hg133'+'.'+name+'SYMBOL'
#
#   	        print db_name
#   	        p2go=eval('AnnotationDbi.as_list(%s)'%db_name)
#
#   	    gene_to_go = collections.defaultdict(list)
#   	    go_to_gene = collections.defaultdict(list)
#   	    p2go_names=numpy.array(p2go.names)
#   	    for i,names in enumerate(p2go_names):
#   	            try:
#   	                for go in  p2go[i]:
#   	                    go_term=str(go)
#   	                    gene_to_go[names].append(go_term)
#   	                    go_to_gene[go_term].append(names)
#
#   	            except Exception,e:
#   	                pass
    except Exception,e:
    	gene_to_go  = {data['rownames'][index]:[i] for index,i in enumerate(probes_to_genes_plus_2(data['rownames']))}

    data['rownames_converted']=[]
    rownames=[]
    for index,probe in enumerate(data['rownames']):
#       print index, probe
#       print gene_to_go[probe]
        if len(gene_to_go[data['rownames'][index]])>1:
            print data['rownames'][index]  + ' is ' + str(gene_to_go[data['rownames'][index]])
        try:

            if gene_to_go[probe][0] =='NA':
                rownames.append(probe)
            else:
                rownames.append(gene_to_go[probe][0])
        except Exception , e:
            print 'probe not on chip... ' + probe
            rownames.append(probe)

    data['rownames_converted']=numpy.array(rownames)
    probes_to_keep=[i for i,j in enumerate(data['rownames_converted']) if j!='NA']


    copydata['rownames']=data['rownames_converted'][probes_to_keep]
    copydata['data']=data['data'][probes_to_keep,:]
    copydata['colnames']=data['colnames']

    uniques_row_names=[]
    for j,g in enumerate(copydata['rownames']):
        if g not in uniques_row_names :
            uniques_row_names.append(g)

    copydata['uniques_row_names'] = uniques_row_names


    if data.has_key('raw_fold_change') and (compute_fc):
        data['fc']=data['raw_fold_change'].T
        matrix =  numpy.zeros([len(uniques_row_names),data['raw_fold_change'].shape[0]])
        for j,g in enumerate(uniques_row_names):
            print j
            if go_to_gene.has_key(g):
#               print 'gene:',str(g)
#               print 'probes',str(go_to_gene[g])
                probes_i = get_probe_index(go_to_gene[g], data['rownames'])
                small_data = data['fc'][ probes_i , : ]
                maxes = small_data.max(axis=0)
                t = small_data == maxes
                nb_probes = len(probes_i)
                nb_samples = data['raw_fold_change'].shape[0]
                #t can contain some genes that are maxed at several probes. rare, but it happens. to we just select one by random-
                T=t.T
                for index , c in enumerate(T):
                    if c.sum != 1:
                        several_maxes_index = numpy.arange(nb_probes)[c]
                        kept_index = several_maxes_index[0]
                        T[index][:]=False
                        T[index][kept_index]= True
                #print g
                #print [g]
                #print t
                #print t.shape
                #print probes_i
                #print nb_samples,nb_probes
                max_probes_index = numpy.array(probes_i.tolist()*nb_samples).reshape(nb_samples,nb_probes)[t.T]
                #print max_probes_index
                #print len(max_probes_index)
            else:
                max_probes_index = get_probe_index( [g], data['rownames']).tolist()*data['raw_fold_change'].shape[0]# len(data['colnames'])
            #print matrix.shape
            #print data['fc'].shape
            #print data['raw_fold_change'].shape[0]
            for i,k in enumerate(max_probes_index):
                matrix[j,i]= data['fc'][k,i]
        copydata['fc_genes']=matrix
    return copydata
    pass

def AML_NK_criteria(data, key, value, above=1):
    """docstring for AML_NK_criteria
        select patient base on some criteria.
    """
    if above==1:
        patients = (data[key] >= value)
        data_to_use = numpy.array(range(len(data[key])))[patients]

    else:
        patients = (data[key] <= value)
        data_to_use = numpy.array(range(len(data[key])))[patients]

    return data_to_use
    pass

def map_platforms_probes(l1,p1,p2):
    """docstring for map_platforms_probes
        the idea is to give a list of probes in one platform, and get the probes of on the second platform.
        Can be used to map mouse to human gene sets for example
        Arguments
        - l1 is the list of query probes
        - p1 is the source platform (as string with R Annotation DBI name)
        - p2 is the query platform (as string with R Annotation DBI name)

        example:
        l2 = map_platforms_probes(['1424032_at', '1438176_x_at', '1455040_s_at'], 'mouse4302.db' ,'hgu133plus2.db')
    """
    import copy
    import collections

    copydata={}


#   now reconstruct the matrix
#   And rownames , with gene names.

    AnnotationDbi=importr('AnnotationDbi')
    platform_1=importr(p1)
    platform_2=importr(p2)
    p1_2_SYMBOL=AnnotationDbi.as_list(eval('platform_1.'+p1[0:p1.find('.db')-1]+'2SYMBOL'))
    p2_2_SYMBOL=AnnotationDbi.as_list(eval('platform_2.'+p2[0:p2.find('.db')-1]+'2SYMBOL'))
    probes_to_genes_1 = collections.defaultdict(list)
    genes_to_probes_1 = collections.defaultdict(list)
    probes_to_genes_2 = collections.defaultdict(list)
    genes_to_probes_2 = collections.defaultdict(list)
    p1_2_SYMBOL_names=numpy.array(p1_2_SYMBOL.names)
    p2_2_SYMBOL_names=numpy.array(p2_2_SYMBOL.names)

    for i,names in enumerate(p1_2_SYMBOL_names):
            try:
                for go in  p1_2_SYMBOL[i]:
                    go_term=str(go).upper()
                    probes_to_genes_1[names].append(go_term)
                    genes_to_probes_1[go_term].append(names)
            except Exception,e:
                pass

    for i,names in enumerate(p2_2_SYMBOL_names):
            try:
                for go in  p2_2_SYMBOL[i]:
                    go_term=str(go).upper()
                    probes_to_genes_2[names].append(go_term)
                    genes_to_probes_2[go_term].append(names)
            except Exception,e:
                pass
    l2=[]
    for i in l1:
        if probes_to_genes_1[i][0] != 'NA':
            the_probes=genes_to_probes_2[probes_to_genes_1[i][0]]
            for p in the_probes:
                if p not in l2:
                    l2.append(p)

    return l2
    pass



######################################################################
class fit_ellipse_2d(object):
    """docstring for fit_ellipse_2d
    input:
    2d array.

    output:
    2x2matrix + 2x1 vector that describes the linear transformation of a circle to get the ellipsoid.

    """
    def __init__(self, XX, color='red'):
        import numpy
        import pylab
        from math import log, pi
        from cvxopt import blas, lapack, solvers, matrix, sqrt, mul, cos, sin
        solvers.options['show_progress'] = False

        super(fit_ellipse_2d, self).__init__()
        self.color = color
        self.XX = matrix(XX)
        self.XX = matrix([self.XX,self.XX[0,:]])
        self.X = matrix(self.XX)
#       self.X=self.X.T
        self.m = self.X.size[0] - 1

        # Inequality description G*x <= h with h = 1
        self.G, self.h = matrix(0.0, (self.m,2)), matrix(0.0, (self.m,1))
        self.G = (self.X[:self.m,:] - self.X[1:,:]) * matrix([0., -1., 1., 0.], (2,2))
        self.h = (self.G * self.X.T)[::self.m+1] # get first element in rows
        self.G = mul(self.h[:,[0,0]]**-1, self.G)
        self.h = matrix(1.0, (self.m,1))

        print self.X
        # Loewner-John ellipsoid
        #
        # minimize   log det A^-1
        # subject to   xk'*A*xk - 2*xk'*b + b'*A^1*b <= 1,  k=1,...,m
        #
        # 5 variables x = (A[0,0], A[1,0], A[1,1], b[0], b[1])

        sol = solvers.cp(self.F)
        A = matrix( sol['x'][[0, 1, 1, 2]], (2,2))
        b = sol['x'][3:]


        L = +A
        lapack.potrf(L)
        c = +b
        lapack.potrs(L, c)

        # 1000 points on the unit circle
        nopts = 1000
        angles = matrix( [ a*2.0*pi/nopts for a in xrange(nopts) ], (1,nopts) )
        circle = matrix(0.0, (2,nopts))
        circle[0,:], circle[1,:] = cos(angles), sin(angles)

        # ellipse = L^-T * circle + c
        blas.trsm(L, circle, transA='T')
        ellipse = circle + c[:, nopts*[0]]
        ellipse2 = 0.5 * circle + c[:, nopts*[0]]

        #pylab.plot(ellipse[0,:].T, ellipse[1,:].T, 'k-')
        pylab.fill(ellipse[0,:].T, ellipse[1,:].T, facecolor=self.color,alpha=.4)
#       pylab.legend('aaa')
#       pylab.plot(self.X[:,0], self.X[:,1], 'ko')


    def F(self,x=None, z=None):
        from math import log, pi
        from cvxopt import blas, lapack, solvers, matrix, sqrt, mul, cos, sin
        solvers.options['show_progress'] = False
        if x is None:
            return self.m, matrix([ 1.0, 0.0, 1.0, 0.0, 0.0 ])
        # Factor A as A = L*L'.  Compute inverse B = A^-1.
        A = matrix( [x[0], x[1], x[1], x[2]], (2,2))
        L = +A
        try: lapack.potrf(L)
        except: return None
        B = +L
        lapack.potri(B)
        B[0,1] = B[1,0]
        # f0 = -log det A
        f = matrix(0.0, (self.m+1,1))
        f[0] = -2.0 * (log(L[0,0]) + log(L[1,1]))
        # fk = xk'*A*xk - 2*xk'*b + b*A^-1*b - 1
        #   = (xk - c)' * A * (xk - c) - 1  where c = A^-1*b
        c = x[3:]
        lapack.potrs(L, c)
        for k in xrange(self.m):
            f[k+1] = (self.X[k,:].T - c).T * A * (self.X[k,:].T - c) - 1.0
        # gradf0 = (-A^-1, 0) = (-B, 0)
        Df = matrix(0.0, (self.m+1,5))
        Df[0,0], Df[0,1], Df[0,2] = -B[0,0], -2.0*B[1,0], -B[1,1]
        # gradfk = (xk*xk' - A^-1*b*b'*A^-1,  2*(-xk + A^-1*b))
        #       = (xk*xk' - c*c', 2*(-xk+c))
        Df[1:,0] = self.X[:self.m,0]**2 - c[0]**2
        Df[1:,1] = 2.0 * (mul(self.X[:self.m,0], self.X[:self.m,1]) - c[0]*c[1])
        Df[1:,2] = self.X[:self.m,1]**2 - c[1]**2
        Df[1:,3] = 2.0 * (-1.*self.X[:self.m,0] + c[0])
        Df[1:,4] = 2.0 * (-1.*self.X[:self.m,1] + c[1])
        if z is None: return f, Df
        # hessf0(Y, y) = (A^-1*Y*A^-1, 0) = (B*YB, 0)
        H0 = matrix(0.0, (5,5))
        H0[0,0] = B[0,0]**2
        H0[1,0] = 2.0 * B[0,0] * B[1,0]
        H0[2,0] = B[1,0]**2
        H0[1,1] = 2.0 * ( B[0,0] * B[1,1] + B[1,0]**2 )
        H0[2,1] = 2.0 * B[1,0] * B[1,1]
        H0[2,2] = B[1,1]**2
        # hessfi(Y, y)
        #    = ( A^-1*Y*A^-1*b*b'*A^-1 + A^-1*b*b'*A^-1*Y*A^-1
        #            - A^-1*y*b'*A^-1 - A^-1*b*y'*A^-1,
        #        -2*A^-1*Y*A^-1*b + 2*A^-1*y )
        #    = ( B*Y*c*c' + c*c'*Y*B - B*y*c' - c*y'*B,  -2*B*Y*c + 2*B*y )
        #    = ( B*(Y*c-y)*c' + c*(Y*c-y)'*B, -2*B*(Y*c - y) )
        H1 = matrix(0.0, (5,5))
        H1[0,0] = 2.0 * c[0]**2 * B[0,0]
        H1[1,0] = 2.0 * ( c[0] * c[1] * B[0,0] + c[0]**2 * B[1,0] )
        H1[2,0] = 2.0 * c[0] * c[1] * B[1,0]
        H1[3:,0] = -2.0 * c[0] * B[:,0]
        H1[1,1] = 2.0 * c[0]**2 * B[1,1] + 4.0 * c[0]*c[1]*B[1,0]  + \
                  2.0 * c[1]**2 + B[0,0]
        H1[2,1] = 2.0 * (c[1]**2 * B[1,0] + c[0]*c[1]*B[1,1])
        H1[3:,1] = -2.0 * B * c[[1,0]]
        H1[2,2] = 2.0 * c[1]**2 * B[1,1]
        H1[3:,2] = -2.0 * c[1] * B[:,1]
        H1[3:,3:] = 2*B

        return f, Df, z[0]*H0 + sum(z[1:])*H1



#####################################################################


def draw_pca_ellipse(test, color=None):
    """docstring for pca_ellipse
        test is some data.
        draws an allipse around this data.
    """

    if color==None:
        color='red'
    import pca_module
    from cvxopt import blas, lapack, solvers, matrix, sqrt, mul, cos, sin
    from math import log, pi
    import pylab
    #o be really sure...
    test=numpy.array(test)
#   print test
#   pylab.scatter( test.T[0],test.T[1])
    mtest=test-test.mean(axis=0)
    stdtest= matrix([0.,0.,0.,0.],(2,2))
    stdtest[0,0]= test.std(axis=0)[0]
    stdtest[1,1]= test.std(axis=0)[1]

    print stdtest
#   t,p,e=pca_module.PCA_nipals2(mtest,PCs=2)
    t,p,e=pca_module.PCA_nipals2(mtest,PCs=2)
    t,e,p=numpy.linalg.svd(test)
    c=+matrix(test.mean(axis=0))
    # 100 points on the unit circle
    nopts = 100
    angles = matrix( [ a*2.0*pi/nopts for a in xrange(nopts) ], (1,nopts) )
    circle = matrix(0.0, (2,nopts))
    circle[0,:], circle[1,:] = cos(angles), sin(angles)
    print stdtest.size
    print circle.size
    print matrix(p).size
    print (stdtest*circle).size
    ellipse = matrix(p)*(stdtest*circle) + c[:, nopts*[0]]
    pylab.fill(ellipse[0,:].T, ellipse[1,:].T, facecolor=color,alpha=.4)
    return t,e,p
    pass


def project_PCA(data,query,samples=None,probes=None,query_samples=None,out='ICA.pdf',interactive=False,Text_Labels=True,PCs=10,flow=None, query_labels=False):
    """docstring for project_PCA
        runs an project_PCA, ie. make a normal PCA space on data, and project query onto that.
    """
    import copy
    import mdp
    import pylab
    from scikits.learn import fastica

    if probes == None:
        probes = range(data['data'].shape[0])
        print probes
    if samples == None:
        samples = range(data['data'].shape[1])
        print samples
    if query_samples == None:
        query_samples = range(query['data'].shape[1])
        print query_samples

    ica = fastica.fastica

    markers= [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
    colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    X = numpy.matrix(data['data'][probes,:][:,samples])
    X = X - X.mean(0).repeat(X.shape[0], 0)
    X = numpy.array(X).T
    #center X
#   X = (X - numpy.mean(X, 0))


    #Rescale x to have zero mean and unit variance
    X = (X - numpy.mean(X, 0))/numpy.std(X, axis=0, ddof=0)
#   pylab.boxplot(X.T)
#   pylab.show()

# generate the PCA space.
    if flow==None:
#       flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.FastICANode(whitened = True)])
#       flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.FastICANode()

        flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs)])
        flow = mdp.nodes.PCANode(output_dim=PCs)

        flow.train(X)
        0
    #use first two components of ICA,

#   data_2d=numpy.array([flow.execute(X)[:,0], flow.execute(X)[:,1] ])
#   [K,W,S] = ica(X)

#   data_2d=numpy.array([K[0,:], K[2,:] ])
#   maxi=data_2d.max()

    maxi=1
#   t1, t2 = K[:,0]/maxi, K[:,1]/maxi

    t1, t2 = flow.execute(X)[:,0]/maxi, flow.execute(X)[:,1]/maxi

    for i in range(len(data['colnames'][samples])):
#       print data['colnames'][samples][i] , data['stemness'][samples][i]

        marker=markers[ int(data['stemness'][samples][i]  ) ]
        color  = colors[ int(data['stemness'][samples][i]  ) ]

        if Text_Labels:
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=6)
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=9, label=data['colnames'][samples][i], alpha=.5)
        else:
            if data.has_key('types')!=True:
                data['types']=data['colnames']
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=9, label=data['types'][samples][i], alpha=.5 )
#           pylab.legend([data['types'][samples][i]])
#           print data['colnames'][samples][i]
#           pylab.legend(data['colnames'][samples][i])
    if not Text_Labels: #then we need some sort of class detection an labelling.
        types = {} # determine types, based on names.
        if data.has_key('types'):
            for i,s in enumerate(data['types'][samples]):
                if s not in types.keys():
                    types[s]=int(data['stemness'][samples][i])
        else:
            types['data']=0
        #now we know which cathegories, and which number they have.
        used_types=[0]*len(types)
        types=types.keys() # temp hack.
        for i in range(len(data['colnames'][samples])):
            #find first of each type, and
            print data['types'][samples][i]
            print used_types
            print
            if used_types[ types.index(data['types'][samples][i]) ]==0 :#  and data['stemness'][samples][i] == '22' :
                pylab.text(t1[i],t2[i],data['types'][samples][i], fontsize=6)
                used_types[ types.index(data['types'][samples][i]) ]=1

# now project the query data into that space
    colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]



    if query.has_key('stemness'):
        cathegories=[]
        for s in query['stemness'][query_samples]:
            if s not in cathegories:
                cathegories.append(s)
    else:
        cathegories=['same']

    train_data=numpy.matrix(data['data'][probes,:][:,samples])
    Y = numpy.matrix(query['data'][probes,:][:,query_samples])
    print Y.shape
    Y = Y - train_data.mean(1)
    Y = numpy.array(Y).T
    #center X
    #   X = (X - numpy.mean(X, 0))
    #Rescale x to have zero mean and unit variance
    Y = (Y - numpy.mean(Y, 0))/numpy.std(Y, axis=0, ddof=0)


    for cindex , cat in enumerate(cathegories):
        indexes=[j for j,i in  enumerate(query['stemness'][query_samples]) if i == cat]

        print 'query %s as %d elements' %(cat,len(indexes))

#       s=numpy.array([Y[sample,:]])
#       print  Y[sample,:].shape
        t1, t2 = flow.execute(Y)[:,0]/maxi, flow.execute(Y)[:,1]/maxi
        print 't1 %d'%(t1.shape)
        for i,sample in enumerate(indexes):
#           pylab.scatter([t1],[t2],color=colors[cindex],alpha=.6,label=query['colnames'][sample])
            print sample, i
            pylab.plot(t1[sample],t2[sample], marker='s', color=colors[cindex+9+cindex*2], markersize=5,  alpha=.5 )
#           pylab.plot(t1[sample],t2[sample], marker='*', color=colors[cindex], markersize=5,  alpha=.5 )
            if query_labels:
                pylab.text(t1[sample],t2[sample],query['colnames'][query_samples][sample], fontsize=6)



    pylab.savefig(out,format='PDF')
    if interactive:
        pylab.show()
    else:
        pylab.close()

    return flow, flow.execute(Y)
    pass

def Classify(data , labels, folds=5):
    """docstring for Classify
    run classification hehe. with cross validation
    """
    from time import time
    import logging
    import pylab as pl

    from sklearn.cross_validation import StratifiedKFold
    from sklearn.cross_validation import KFold
#   from sklearn.datasets import fetch_lfw_people
    from sklearn.grid_search import GridSearchCV
    from sklearn.metrics import classification_report
    from sklearn.metrics import confusion_matrix
    from sklearn.decomposition import RandomizedPCA
    from sklearn.svm import SVC



    # Display progress logs on stdout
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


    ################################################################################
    # The data should be as numpy arrays

    X = data.T
    n_features = X.shape[1]
    n_samples = X.shape[0]
    cat=[]
    for i in labels:
        if i not in cat:
            cat.append(i)
    y = numpy.array([cat.index(i) for i in labels])
    # the label to predict is the id of the person
#   y = lfw_people.target
    target_names = numpy.array(cat)
    n_classes = target_names.shape[0]

    print "Total dataset size:"
    print "n_samples: %d" % n_samples
    print "n_features: %d" % n_features
    print "n_classes: %d" % n_classes


    ################################################################################
    # Split into a training set and a test set using a stratified k fold

    # split into a training and testing set
#   train, test = iter(StratifiedKFold(y, k=4)).next() #StratifiedKfold is of course better.

#   for train, test in iter(KFold(n=n_samples, k=2)):
    for train, test in iter(StratifiedKFold(y, n_folds=folds)):
#   train, test = iter(KFold(n=n_samples, k=10)).next()
        if folds > 1:
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]
        else:
            X_train, X_test = X, X
            y_train, y_test = y, y

        ################################################################################
        # Compute a PCA on the  dataset (treated as unlabeled
        # dataset): unsupervised feature extraction / dimensionality reduction
        n_components = 150

        print "Extracting the top %d eigensamples from %d samples" % (
        n_components, X_train.shape[0])
        t0 = time()
        pca = RandomizedPCA(n_components=n_components, whiten=True).fit(X_train)
        print "done in %0.3fs" % (time() - t0)

        eigenfaces = pca.components_#.reshape((n_components, h, w))

        print "Projecting the input data on the eigenfaces orthonormal basis"
        t0 = time()
        X_train_pca = pca.transform(X_train)
        X_test_pca = pca.transform(X_test)
        print "done in %0.3fs" % (time() - t0)


        ################################################################################
        # Train a SVM classification model

        print "Fitting the classifier to the training set"
        t0 = time()
        param_grid = {
         'C': [1, 5, 10, 50, 100],
         'gamma': [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1],
        }
        clf = GridSearchCV(SVC(kernel='rbf'), param_grid)#,
                          # fit_params={'class_weight': 'auto'})
        clf = clf.fit(X_train_pca, y_train)
        print "done in %0.3fs" % (time() - t0)
        print "Best estimator found by grid search:"
        print clf.best_estimator_


        ################################################################################
        # Quantitative evaluation of the model quality on the test set

        print "Predicting the AML classes in testing set"
        t0 = time()
        y_pred = clf.predict(X_test_pca)
        print "done in %0.3fs" % (time() - t0)

        print classification_report(y_test, y_pred, target_names=target_names)
        print confusion_matrix(y_test, y_pred, labels=range(n_classes))
    return X_test_pca

def make_mst(data,probes_to_use=None,data_to_use=None, max_node=5 , min_corr=.9,out='tree.pdf'):
    """docstring for make_mst
        compute correlation matrix for samples, and make a network,
        which is then purned as an minimum_spanning_tree.
    """
    import Bio
    import networkx as nx
    import pylab


    if probes_to_use == None:
        print 'addidng probes'
        probes_to_use = numpy.arange(data['data'].shape[0])
    if data_to_use == None:
        print 'addidng sample'
        data_to_use = numpy.arange(data['data'].shape[1])

    print 'Using %d probes'%len(probes_to_use)
    print 'Using %d samples'%len(data_to_use)

    names= data['colnames'][data_to_use]
    cats={}
    new_names=[]
    for i in names:
        if i not in cats:
            cats[i]=1
            new_names.append(i)
        else:
            cats[i]+=1
            new_names.append(i+'_' +str(cats[i]))


    fc=data['data'][probes_to_use,:][:,data_to_use]

    distance_matrix = Bio.Cluster.distancematrix(fc.T , dist="c")
    nb_s = len(data_to_use)
    print 'making contact matrix'
    prog = ProgressBar(0, nb_s , 50, mode='fixed')
    similarity_matrix=numpy.zeros( [nb_s, nb_s] , dtype=float)
    for i,vec in enumerate(distance_matrix):
    #   similarity_matrix[i,i] = False
    #   similarity_matrix[i,range(i)] = similarity_matrix[range(i),i] = numpy.abs(1. - vec) > cor_thr
        similarity_matrix[i,i] = 0
        similarity_matrix[i,range(i)] = similarity_matrix[range(i),i] = numpy.abs(1. - vec)

        #progress bar
        oldprog = str(prog)
        prog.update_amount(i)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)

    APL_G=nx.Graph()
    max_node = max_node
    thr =  min_corr
    c=similarity_matrix
    indexes = numpy.arange(len(data_to_use))
    for j,l in enumerate(c):
    #   print j
        b=numpy.sort(l)
        for i in b[-max_node:]:
            if i > thr and i <1:
                k = indexes[l == i][0]
                print new_names[j], '->' ,new_names[k], str(i)
                g1 = new_names[j]
                g2 = new_names[k]
                APL_G.add_node(g1)
                APL_G.add_node(g2)
                APL_G.add_edge(g1,g2)
                APL_G.edge[g1][g2]['weight'] = i

    #remove clusters of 2
    clustering = nx.clustering(APL_G)
    for i in clustering:
        if clustering[i] == 0:
            0
#            APL_G.remove_node(i)

    #analysis..
    mst = nx.minimum_spanning_tree(APL_G)
  #  mst = APL_G
    colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']


    #drawing can be made better....
    #nx.draw_graphviz(APL_G)


#    return APL_G

    try:
        pos =nx.graphviz_layout(mst)
        print 'graphviz-'
    except:
        pos=nx.spring_layout(mst,iterations=400)

    color ={}
    ii =0
    for i,n in enumerate(APL_G.nodes()):
        if n.split('_')[0] not in color:
            color[n.split('_')[0]] = ii
            ii += 1
    pylab.figure(figsize=(20,12))
    nx.draw_networkx_edges(mst,pos,alpha=0.3, edge_color='m')
#   nodesize=[wins[v]*50 for v in H]
    node_color = [ color[n.split('_')[0]] for n in APL_G.nodes()]#[ APL_G.node[i]['color'] for i in APL_G.nodes()]
#   node_size = [ numpy.abs(APL_G.node[i]['size']) for i in APL_G.nodes()]
    nx.draw_networkx_nodes(APL_G,pos,node_size=500,node_color=node_color,alpha=0.4)
#    edge_width = [ (APL_G.edge[i[0]][i[1]]['weight']-thr)*100.  for i in APL_G.edges()]
    nx.draw_networkx_edges(mst,pos,alpha=0.4,node_size=0,edge_color='k')
    nx.draw_networkx_labels(mst,pos,fontsize=6)


    pylab.savefig(out)
    pylab.close()
    return mst,APL_G
    pass

def propagation_similarities(data, probes_to_use=None, data_to_use=None ,s=None ):
    """docstring for propagation_similarities

        data: is the data , i.e. dataset['data']

        s: is the function to calculate the similarity.
            None uses the basic similarity found in scikit page
            http://scikit-learn.sourceforge.net/auto_examples/cluster/plot_affinity_propagation.html
    """
    from sklearn.cluster import AffinityPropagation
    from sklearn import metrics


    if probes_to_use == None:
        probes_to_use = range(data.shape[0])
    if data_to_use == None:
        data_to_use = range(data.shape[1])

    X= data[probes_to_use,:][:,data_to_use].T

    ################################################################################
    # Compute similarities
    ################################################################################
    if s == None:
        X_norms = numpy.sum(X*X, axis=1)
        S = - X_norms[:,numpy.newaxis] - X_norms[numpy.newaxis,:] + 2 * numpy.dot(X, X.T)
    else:
        S = 0
    p = .5*numpy.median(S)

    ################################################################################
    # Compute Affinity Propagation
    ################################################################################
    af = AffinityPropagation(preference=-50).fit(X)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_

    n_clusters_ = len(cluster_centers_indices)

    print('Estimated number of clusters: %d' % n_clusters_)
    #print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    #print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    #print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    #print("Adjusted Rand Index: %0.3f" % metrics.adjusted_rand_score(labels_true, labels))
    #print("Adjusted Mutual Information: %0.3f"% metrics.adjusted_mutual_info_score(labels_true, labels))
    #print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels, metric='sqeuclidean'))


    return labels
    pass
#


def compare_genes_aml(data, gene,data_to_use):
    '''
    compare expression of gene x in amls and in chimera_sample normal.
    '''
    from numpy import concatenate
    g = get_probe_index(gene,data['rownames'])
    g_aml = data['data'][g,data_to_use]
    g_nl = data['chimera_sample'].T[g,data_to_use]
    l=len(data_to_use)
    beeswarm(concatenate([g_aml,g_nl]),concatenate([['AML']*l,['NL']*l]) , title =str(gene),plot_median =1 , out=str(gene[0])+'.pdf' )

    pass

def pretty_print_significance(x,log_neg=False):
	if log_neg==True:
		x=pow(10,-x)
	if x > .05:
		print 'ns'
	elif x <=.05 and x >.001:
		print '*'
	elif x <=.001 and x >.0001:
		print '**'
	elif x <= pow(10,-5):
		print '***'

	pass

def get_sample_from_hdf5(name):
    """docstring for get_sample_from_hdf5
        provide a name (can be cel file, or GSM number and return the sample.)
    """
    import copy
    import gc
    import platform;
    print name
    name= name.replace('.CEL','').replace('.gz','')
    hdf5_files=[\
#               '/Users/nrapin/BRIC/AML_NK_study/amls_NK_w_new_voisins_bis.hdf5',\
                '/Users/nrapin/BRIC/MILE_study/amls_MILE_w_new_voisins_bis.hdf5',\
                '/Users/nrapin/BRIC/MILE_study/merged_Mile_w_p2.hdf5',\
                '/Users/nrapin/BRIC/BRIC_AML_Study/merged_bric.hdf5',\
                '/Users/nrapin/BRIC/BRIC_AML_Study/Kasumi/merged_Kasumi.hdf5',\
#               '/Users/nrapin/BRIC/OS_study/merged_os_plus2.hdf5',\
                '/Users/nrapin/BRIC/AML_NK_study/merged_AML_NK.hdf5',\
                '/Users/nrapin/BRIC/GSE6891/merged_GSE6891.hdf5'\
                ]
    hdf5_files=['/Volumes/hdd/BRIC/Processed_data/all_merged.hdf5']
    if platform.uname()[1] == "SGI (Porse)":
            hdf5_files=['/home/nrapin/BRIC/MILE_study/merged_Mile_w_p2.hdf5',\
                        '/home/nrapin/BRIC/BRIC_AML_Study/merged_bric.hdf5',\
                        '/home/nrapin/BRIC/BRIC_AML_Study/Kasumi/merged_Kasumi.hdf5',\
#                       '/home/nrapin/BRIC/OS_study/merged_os_plus2.hdf5',\
                        '/home/nrapin/BRIC/AML_NK_study/merged_AML_NK.hdf5',\
                        '/home/nrapin/BRIC/GSE6891/merged_GSE6891.hdf5'\
                        ]


    #create a dic with all the merged_data for individual sample

    all_dic={}
    i=0
    file_pointers=[]
    for db_f in list(hdf5_files):
        file_pointers.append(h5py.File(db_f,'r'))
    for db_f in file_pointers:
        print db_f
#       db_f = '/Users/nrapin/BRIC/BRIC_AML_Study/merged_bric.hdf5'
        i+=1
        print 'Looking in' + str(db_f)
        root=db_f['root']
#       for sample, l in entries:
        if name in root:
            print name + ' found in ' + db_f.filename.split('/')[-1]
            a= hdf52dic(root[str(name)])
            break
#       f.close()
        gc.collect()

    for db_f in file_pointers:
        db_f.close()
    return a
    pass


def merge_experiments(file_list, mode = 'all_keys', import_nl=0, recompute_voisins = False,patient_cohort=False,hdf5=['/Volumes/hdd/BRIC/Processed_data/all_merged.hdf5']):
    """docstring for merge_experiments
        merges experiments into one workable dataset.

        mode :  all_keys imports all keys in found (generated by the normal normalization programs)
                basic_keys imports only basic keys [ 'data',  'platforms', 'rownames',  'cel_file', 'colnames', 'stemness']

        import_nl : 1 or 0 , appends the nl hierarchy at the end also.
    """
    import copy
    import gc
    from collections import OrderedDict
    entries = import_sample_data(file_list)
#   print 'Files to load:'
#   print '\n'.join(numpy.array(entries)[:,0])
    #put here list of all merged datasets.
    hdf5_files=[#'/Users/nrapin/BRIC/AML_NK_study/amls_NK_w_new_voisins_bis.hdf5',\
#               '/Users/nrapin/BRIC/MILE_study/amls_MILE_w_new_voisins_bis.hdf5',\
                '/Users/nrapin/BRIC/MILE_study/merged_Mile_w_p2.hdf5',\
                '/Users/nrapin/BRIC/BRIC_AML_Study/merged_bric.hdf5',\
                '/Users/nrapin/BRIC/BRIC_AML_Study/Kasumi/merged_Kasumi.hdf5',\
#               '/Users/nrapin/BRIC/OS_study/merged_os_plus2.hdf5',\
                '/Users/nrapin/BRIC/AML_NK_study/merged_AML_NK.hdf5',\
                '/Users/nrapin/BRIC/GSE6891/merged_GSE6891.hdf5'\
                ]
    if patient_cohort:
        hdf5_files=['/Volumes/hdd/BRIC/Processed_data/new_normalization/new_normalization_all_merged.hdf5']
    else:
        hdf5_files=['/Volumes/hdd/BRIC/Processed_data/all_merged.hdf5']

    import platform;
    if platform.uname()[1] == "SGI (Porse)":
            hdf5_files=['/home/nrapin/BRIC/MILE_study/merged_Mile_w_p2.hdf5',\
                        '/home/nrapin/BRIC/BRIC_AML_Study/merged_bric.hdf5',\
                        '/home/nrapin/BRIC/BRIC_AML_Study/Kasumi/merged_Kasumi.hdf5',\
        #               '/home/nrapin/BRIC/OS_study/merged_os_plus2.hdf5',\
                        '/home/nrapin/BRIC/AML_NK_study/merged_AML_NK.hdf5',\
                        '/home/nrapin/BRIC/GSE6891/merged_GSE6891.hdf5'\
                        ]


    #create a dic with all the merged_data for individual sample
    if hdf5 != ['/Volumes/hdd/BRIC/Processed_data/all_merged.hdf5']:
   	    hdf5_files = hdf5
   	    print 'Using custom exp. database'
    all_dic=OrderedDict()
    i=0
    file_pointers=[]
    for db_f in list(hdf5_files):
        print db_f
        file_pointers.append(h5py.File(db_f,'r'))
    for db_f in file_pointers:
        print db_f
#       db_f = '/Users/nrapin/BRIC/BRIC_AML_Study/merged_bric.hdf5'
        i+=1
        print 'Looking in' + str(db_f)
        root=db_f['root']
        for sample, l in entries:
            if sample in root:
                print sample + ' found in ' + db_f.filename.split('/')[-1]
                a= hdf52dic(root[str(sample)])
                print ' '.join(a.keys())
                if a=={}:
                    a= hdf52dic(root[str(sample)+'.CEL'])
                    if a == {}:
                        print sample + 'not working...'

                if recompute_voisins:
                    #a['data'] = filter_low_genes(a['data'],a['data'].mean())
                    probes_to_use  = select_probes_by_variance(a, var_thr=4)
                    data_to_use = range(len(a['colnames']))
                    [T,P,E] = run_PCA(a,probes=probes_to_use,interacive=False,fit_ellipse=False,Text_Labels=False)
                    [distance , data_to_use ] = PCA_neighbourg(copy.copy(a),T,E,a['cel_file'][-1],Variance_dims=6,mode='both',data_to_use=data_to_use)
                    probes_to_use  = select_probes_by_variance(a, var_thr=2, data_to_use=data_to_use)
                    [T,P,E] = run_PCA(a,data_to_use,probes_to_use,interacive=False,fit_ellipse=False)
                    [top_n,low_n, raw_fold_change,chimera_sample] = find_neighbour(copy.copy(a),T,E,a['cel_file'][-1],data_to_use=data_to_use,Variance_dims=6, n=500)
                    a['voisins']=[numpy.array(data_to_use),numpy.array(distance)]
                    a['top_500']=top_n
                    a['low_500']=low_n
                    a['chimera_sample']=chimera_sample
                    a['raw_fold_change']=raw_fold_change
                #print a
                all_dic[sample] = copy.copy( a )
#       f.close()
        gc.collect()

    for db_f in file_pointers:
        db_f.close()

#   now check that they have the basic keys:
    if mode != 'all_keys':
        all_keys=basic_keys=[ 'data',  'platforms', 'rownames',  'cel_file', 'colnames', 'stemness']
    else:
        all_keys=whole_set_of_keys=[ 'data', 'voisins', 'raw_fold_change', 'platforms', 'rownames', 'chimera_sample', 'cel_file', 'colnames',  'stemness']

    merged_data ={}
    for sample in all_dic:
        print sample
        #print all_dic[sample].keys()
        for key in all_keys:
            if key not in merged_data:
                merged_data[key]=[]
            if key in [ 'platforms',  'cel_file', 'colnames', 'stemness']:
                merged_data[key].append(all_dic[sample][key][-1])
#                   if key == 'colnames':
#                       print all_dic[sample][key][-1]
            elif key in ['data']:
                merged_data[key].append(all_dic[sample][key][:,-1])
            elif key in ['voisins','rownames']:
                merged_data[key].append(all_dic[sample][key])
            elif key in ['raw_fold_change' ,'chimera_sample' ]:
                merged_data[key].append(all_dic[sample][key])
            elif key in ['top_500','low_500']:
                for gene in all_dic[sample][key]:
                    if gene not in merged_data[key]:
#                           print gene
                        merged_data[key].append(gene)
    if merged_data =={}:
        print 'Something\'s wrong, no arrays found...'
        return
#check if there are by any chances some arrays that are from different platform
    platforms=OrderedDict() #number + length of rownames
    p_types=OrderedDict()  #number + rownames
    sample_platform=[]
#   print all_dic
    for index,sample in enumerate(merged_data['rownames']):
        ss = sample.shape
        if ss not in platforms.values():
            if platforms.keys() == []:
                platforms[1] = ss
                p_types[1] = merged_data['rownames'][index]
            else:
                platforms[len(platforms.keys())+1] = ss
                p_types[len(platforms.keys())+0] = merged_data['rownames'][index]

    sample_platform=[]
    for index,sample in enumerate(merged_data['rownames']):
        for p in platforms:
            if platforms[p][0] == sample.shape[0]:
                sample_platform.append(p)
#   print sample_platform
#       move to numpy arrays ..
    all_data=OrderedDict()
    if len(platforms) == 1:
        print 'Only one platform!'
        for key in all_keys:
            if key in [ 'top_n','low_n' ,'platforms', 'rownames',  'cel_file', 'colnames', 'stemness', 'raw_fold_change', 'chimera_sample']:
                all_data[key] = numpy.array(merged_data[key])
            elif key in ['data']:
                all_data[key] = numpy.array(merged_data[key]).T
            elif key in ['voisins']:
                all_data[key]=merged_data[key]
        all_data['rownames'] = all_dic[all_dic.keys()[0]]['rownames']
    else: #now we need to get to the bottom of things and merged based on names...that's dangerous...
    #find platform with smallest number of probe sets-
        #all_data=p_types
        shortest = -1
        tmpval=1e10
        for k,v in platforms.iteritems():
            print k,v[0]
            if v[0] < tmpval:
                tmpval=v[0]
                shortest=k
        listed=OrderedDict()
        for key in p_types:
            listed[key]=p_types[key].tolist()

        indexes=OrderedDict()
        all_platforms=listed.keys()
        all_platforms_minus_shortest = copy.copy(all_platforms)
        all_platforms_minus_shortest.pop(all_platforms.index(shortest))
        for larger_platform in all_platforms_minus_shortest:
            probes=[]
            for n in listed[shortest]:
                if n in listed[larger_platform]:
                    probes.append( listed[larger_platform].index(n) )
            indexes[larger_platform]= probes
        indexes[shortest]=range(len( listed[shortest] ))

        #now we're good to go.
        data=[]
        for ii in range(len(merged_data['colnames'])):
            cur_sample=numpy.array(merged_data['data'][ii])[indexes[sample_platform[ii]]]
            data.append(cur_sample)
        raw_fold_change=[]
        for ii in range(len(merged_data['colnames'])):
            cur_sample=numpy.array(merged_data['raw_fold_change'][ii])[indexes[sample_platform[ii]]]
            raw_fold_change.append(cur_sample)
        chimera_sample=[]
#       for ii in range(len(merged_data['colnames'])):
#           cur_sample=numpy.array(merged_data['chimera_sample'][ii])[indexes[sample_platform[ii]]]
#           chimera_sample.append(cur_sample)

        for key in all_keys:
            if key in [ 'top_n','low_n' ,'platforms',  'cel_file', 'colnames', 'stemness']:
                all_data[key] = numpy.array(merged_data[key])
            elif key in ['data']:
                all_data[key] = numpy.array(data).T
            elif key in ['raw_fold_change']:
                all_data[key] = numpy.array(raw_fold_change)
#           elif key in ['chimera_sample']:
#               all_data[key] = numpy.array(chimera_sample).T
            elif key in ['voisins']:
                all_data[key]=merged_data[key]
        all_data['rownames'] = p_types[shortest]
    print 'Merged %d experiments...'%len(all_data['colnames'])
    print '----------'
    missing = False
    for i in numpy.array(entries)[:,0]:
        name_o = i + '.CEL.gz'
        if not(name_o in  all_data['cel_file'] or i  in  all_data['cel_file']) :
            print 'Missing Sample : %s' %(i)
            missing = True
    if not missing:
        print 'All samples merged.'

    if import_nl == 1 and patient_cohort == False:
        print  'importing nl cells.'
        if platform.uname()[1] == "SGI (Porse)":
            nl = hdf_load('/home/nrapin/BRIC/normal_cells/nl_cells.hdf5')
        else:
            nl = hdf_load('/Users/nrapin/BRIC/normal_cells/nl_cells.hdf5')
        if recompute_voisins:
            nl['data'] = filter_low_genes(nl['data'],thr=nl['data'].mean())
        for key in ['platforms',  'cel_file', 'colnames',  'stemness', 'data']:
            all_data[key] = numpy.concatenate( [ all_data[key], nl[key] ] , axis=1)
    elif import_nl != 1 and import_nl != 0:
        print  'importing nl cells.'
        if platform.uname()[1] == "SGI (Porse)":
            nl = pickle_load('/home/nrapin/BRIC/normal_cells/nl_cells.pkl')
        else:
            nl = pickle_load('/Users/nrapin/BRIC/normal_cells/nl_cells.pkl')

        if recompute_voisins:
            nl['data'] = filter_low_genes(nl['data'],thr=nl['data'].mean())
        for key in ['platforms',  'cel_file', 'colnames',  'stemness', 'data']:
            all_data[key] = numpy.concatenate( [ all_data[key], nl[key] ] , axis=1)

    if import_nl == 1 and patient_cohort:
        print  'importing nl cells.'
        if platform.uname()[1] == "SGI (Porse)":
            print 'please update me'
            nl = hdf_load('/home/nrapin/BRIC/patient_cohort_nl/merged_data_batch_corrected_w_mono.pkl')
        else:
            nl = pickle_load('/Volumes/hdd/BRIC/patient_cohort_nl/merged_data_batch_corrected_w_mono.pkl')
            nl['stemness'] =  stemness_by_class(nl)
        if recompute_voisins:
            nl['data'] = filter_low_genes(nl['data'],thr=nl['data'].mean())
        all_data = merge_datasets_on_highest(nl,all_data)
        #for key in ['platforms',  'cel_file', 'colnames',  'stemness', 'data']:
            #all_data[key] = numpy.concatenate( [ all_data[key], nl[key] ] , axis=1)
    elif import_nl != 1 and import_nl != 0:
        print  'importing nl cells.'
        if platform.uname()[1] == "SGI (Porse)":
            nl = pickle_load('/home/nrapin/BRIC/patient_cohort_nl/merged_data_batch_corrected_w_mono.pkl')
        else:
            nl = pickle_load('/Users/nrapin/BRIC/patient_cohort_nl/merged_data_batch_corrected_w_mono.pkl')
            nl['stemness'] =  stemness_by_class(nl)
        if recompute_voisins:
            nl['data'] = filter_low_genes(nl['data'],thr=nl['data'].mean())
        for key in ['platforms',  'cel_file', 'colnames',  'stemness', 'data']:
            all_data[key] = numpy.concatenate( [ all_data[key], nl[key] ] , axis=1)


    return all_data
    pass

def run_ICA(data,samples=None,probes=None,out='ICA.pdf',interacive=False,fit_ellipse=True,Text_Labels=True,PCs=10,flow=None, no_plot=False):
    """docstring for run_ICA
        runs an ICA
        to reduce computationnal load, we run a PCA first.
    """
    import copy
    import mdp
    import pylab
#   from scikits.learn import fastica

    if probes == None:
        print 'addidng probes'
        probes = numpy.arange(data['data'].shape[0])
        print probes
    if samples == None:
        print 'addidng sample'
        samples = numpy.arange(data['data'].shape[1])
        print samples

#   ica = fastica.fastica

    markers= [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
    colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    X = numpy.matrix(data['data'][probes,:][:,samples])
    X = X - X.mean(0).repeat(X.shape[0], 0)
#   X = (X - numpy.mean(X, 0))/numpy.std(X, axis=0, ddof=0)
    X = numpy.array(X).T

    #center X
#   X = (X - numpy.mean(X, 0))

    #Rescale x to have zero mean and unit variance





    if flow==None:
#       flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.FastICANode(whitened = True)])
#       flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.FastICANode()

#       flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs)])
#       flow = mdp.nodes.PCANode(output_dim=PCs)
        flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.SFANode()])
        flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.SFANode()

#       flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.CuBICANode(whitened = True)])
#       flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.CuBICANode()

        flow.train(X)
        0
    #use first two components of ICA,

    if no_plot:
        return flow.execute(X)

    data_2d=numpy.array([flow.execute(X)[:,0], flow.execute(X)[:,1] ])
#   [K,W,S] = ica(X)
#   data_2d=numpy.array([K[0,:], K[1,:] ])
    maxi=data_2d.max()
    maxi=1
#   t1, t2 = K[:,0]/maxi, K[:,1]/maxi

    t1, t2 = flow.execute(X)[:,0]/maxi, flow.execute(X)[:,1]/maxi

    for i in range(len(data['colnames'][samples])):
#       print data['colnames'][samples][i] , data['stemness'][samples][i]

        marker=markers[ int(data['stemness'][samples][i]  ) ]
        color  = colors[ int(data['stemness'][samples][i]  ) ]

        if Text_Labels:
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=6)
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.5)
        else:
            if data.has_key('types')!=True:
                data['types']=data['colnames']
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['types'][samples][i], alpha=.5 )
#           pylab.legend([data['types'][samples][i]])
#           print data['colnames'][samples][i]
#           pylab.legend(data['colnames'][samples][i])
    if not Text_Labels: #then we need some sort of class detection an labelling.
        types = {} # determine types, based on names.
        if data.has_key('types'):
            for i,s in enumerate(data['types'][samples]):
                if s not in types.keys():
                    types[s]=int(data['stemness'][samples][i])
        else:
            types['data']=0
        #now we know which cathegories, and which number they have.
        used_types=[0]*len(types)
        types=types.keys() # temp hack.
        for i in range(len(data['colnames'][samples])):
            #find first of each type, and
            print data['types'][samples][i]
            print used_types
            print
            if used_types[ types.index(data['types'][samples][i]) ]==0 :#  and data['stemness'][samples][i] == '22' :
                pylab.text(t1[i],t2[i],data['types'][samples][i], fontsize=6)
                used_types[ types.index(data['types'][samples][i]) ]=1

    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors = ['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','m','b','c','y','k']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            try:
                if len(indexes) > 2:
                    eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
#                   eee.append(draw_pca_ellipse( numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]))
                    colors.pop(-1)
            except Exception, e:
                print e

    pylab.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return flow
    pass
#
def minimum_spanning_tree(data , probes_to_use=None, data_to_use=None, out = 'minimum_spanning_tree.pdf'):
    """docstring for minimum_spanning_tree
        makes a minimum spanning tree from correlation between samples-
    """
    import networkx as nx
    import pylab
    if probes_to_use == None:
        print 'addidng probes'
        probes_to_use = numpy.arange(data['data'].shape[0])
    if data_to_use == None:
        print 'addidng sample'
        data_to_use = numpy.arange(data['data'].shape[1])

    colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
    names = data['colnames'][data_to_use]

    cats={}
    new_names=[]
    for i in names:
        if i not in cats:
            cats[i]=1
            new_names.append(i)
        else:
            cats[i]+=1
            new_names.append(i+'_' +str(cats[i]))

    corr = fast_cor(data, probes_to_use=probes_to_use, data_to_use=data_to_use)
    G=nx.Graph()
    G.add_nodes_from(new_names)
    for index,i in enumerate(G.nodes()):
        G.node[i]['color'] = colors[numpy.array(data['stemness'], dtype = int)[ data_to_use [index]]]
#   return G
    for j,l in enumerate(corr):
        for k,i in enumerate(l):
            g1=new_names[j]
            g2=new_names[k]
#           G.add_edge(g1,g2)
            if i > .4:
                G.add_edge(g1,g2)

                G.edge[g1][g2]['weight'] = numpy.abs(i)*numpy.abs(i)
#           else:
#               G.edge[g1][g2]['weight'] = 0

    mst = nx.minimum_spanning_tree(G)

    pos =nx.graphviz_layout(mst)

    fig=pylab.figure(figsize=(20,12))

    nx.draw_networkx_edges(mst,pos,alpha=0.3, edge_color='m')
    node_color = [ mst.node[i]['color'] for i in mst.nodes()]
    node_size = [ 1 for i in mst.nodes()]
    nx.draw_networkx_nodes(mst,pos,node_size=100*numpy.array(node_size),node_color=node_color,alpha=0.4)
    edge_width = [ (mst.edge[i[0]][i[1]]['weight'])  for i in mst.edges()]
    nx.draw_networkx_edges(mst,pos,alpha=0.4,node_size=0,width=edge_width,edge_color='k')
    nx.draw_networkx_labels(mst,pos,fontsize=6)

    pylab.savefig(out)
    pass

def map_PCA(data,samples=None,probes=None, train=None, out='Mapped_PCA.pdf',interacive=False,fit_ellipse=False,Text_Labels=True,PCs=10,flow=None, no_plot=False, font_size=10):
    """docstring for run_ICA
        runs an PCA
        to reduce computationnal load, we run a PCA first.
    """
    import copy
    import mdp
    import pylab
#   from scikits.learn import fastica

    if probes == None:
        print 'addidng probes'
        probes = numpy.arange(data['data'].shape[0])
        print probes
    if samples == None:
        print 'addidng sample'
        samples = numpy.arange(data['data'].shape[1])
        print samples
    if train == None:
        print 'Training on all set'
        probes = numpy.arange(data['data'].shape[1])
        print samples


    markers= [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
    colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    X = numpy.matrix(data['data'][probes,:][:,train])
#   X = X - X.mean(0).repeat(X.shape[0], 0)

    #center X
#   X = (X - numpy.mean(X, 0))

    #Rescale x to have zero mean and unit variance
    X = (X - numpy.mean(X, 0))/numpy.std(X, axis=0, ddof=0)
    X = numpy.array(X).T



#   flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.SFANode()])
#   flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.SFANode()

    if flow==None:
#       flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.FastICANode(whitened = True)])
#       flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.FastICANode()

        flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs)])
        flow = mdp.nodes.PCANode(output_dim=PCs)

#       flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs), mdp.nodes.CuBICANode(whitened = True)])
#       flow = mdp.nodes.PCANode(output_dim=PCs) + mdp.nodes.CuBICANode()

        flow.train(X)
        0
    #use first two components of ICA,

    if no_plot:
        return flow.execute(X)

    Y = numpy.matrix(data['data'][probes,:][:,samples])
    Y = Y - Y.mean(0).repeat(Y.shape[0], 0)
    Y = numpy.array(Y).T

    #Rescale x to have zero mean and unit variance
    Y = (Y - numpy.mean(Y, 0))/numpy.std(Y, axis=0, ddof=0)

    data_2d=numpy.array([flow.execute(Y)[:,0], flow.execute(Y)[:,1] ])
#   [K,W,S] = ica(X)
#   data_2d=numpy.array([K[0,:], K[1,:] ])
    maxi=data_2d.max()
    maxi=1
#   t1, t2 = K[:,0]/maxi, K[:,1]/maxi

    t1, t2 = flow.execute(Y)[:,0]/maxi, flow.execute(Y)[:,1]/maxi

    for i in range(len(data['colnames'][samples])):
#       print data['colnames'][samples][i] , data['stemness'][samples][i]

        marker=markers[ int(data['stemness'][samples][i]  ) ]
        color  = colors[ int(data['stemness'][samples][i]  ) ]

        if Text_Labels:
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=font_size)
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.5)
        else:
            if data.has_key('types')!=True:
                data['types']=data['colnames']
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['types'][samples][i], alpha=.5 )
#           pylab.legend([data['types'][samples][i]])
#           print data['colnames'][samples][i]
#           pylab.legend(data['colnames'][samples][i])
    if not Text_Labels: #then we need some sort of class detection an labelling.
        types = {} # determine types, based on names.
        if data.has_key('types'):
            for i,s in enumerate(data['types'][samples]):
                if s not in types.keys():
                    types[s]=int(data['stemness'][samples][i])
        else:
            types['data']=0
        #now we know which cathegories, and which number they have.
        used_types=[0]*len(types)
        types=types.keys() # temp hack.
        for i in range(len(data['colnames'][samples])):
            #find first of each type, and
            print data['types'][samples][i]
            print used_types
            print
            if used_types[ types.index(data['types'][samples][i]) ]==0 :#  and data['stemness'][samples][i] == '22' :
                pylab.text(t1[i],t2[i],data['types'][samples][i], fontsize=font_size)
                used_types[ types.index(data['types'][samples][i]) ]=1

    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors = ['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','m','b','c','y','k']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            try:
                if len(indexes) > 2:
                    eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
#                   eee.append(draw_pca_ellipse( numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]))
                    colors.pop(-1)
            except Exception, e:
                print e

    pylab.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return flow
    pass
#

def find_AML_patient_nb(s):
	l=s.split('-')
	aml=''
	for i in l:
		if i.find('AML') != -1:
			aml=i
			break
	if aml == '':
		aml=s
	return aml
	pass

def make_correlation_mimimum_spanning_tree(data,probes=None,samples=None, min_corr=.9, max_node = 5, out='corr_tree.gml',fig_out='net.png'):
    """
    makes a minimmu, spanning tree from correlation winthin the data.

    """
    import networkx as nx
    import pylab
    if probes == None:
        print 'addidng probes'
        probes = numpy.arange(data['data'].shape[0])
    if samples == None:
        print 'addidng sample'
        samples = numpy.arange(data['data'].shape[1])

    a=fast_cor(data, probes_to_use=probes, data_to_use=samples, figure='corr_spanning_tree.pdf')
    en_G=nx.Graph()
    thr = min_corr
    indexes = numpy.arange(len(data['colnames'][samples]))
    for j,l in enumerate(a):
    #   print j
        b=numpy.sort(l)
        for i in b[-max_node:]:
            if i > thr and i <1:
                try:
                    k = indexes[l == i][0]
                except :
                    k = list(l).index(i)
                print data['colnames'][samples][j], '->' ,data['colnames'][samples][k], str(i)
                g1 = data['colnames'][samples][j]
                g2 = data['colnames'][samples][k]
                en_G.add_edge(g1,g2)
                en_G.edge[g1][g2]['weight'] = i
    en_tree = nx.minimum_spanning_tree(nx.Graph(en_G))
    nx.write_gml(en_tree,'corr_tree.gml')
    nx.write_gml(en_G,'corr_net.gml')
    #nx.draw_graphviz(en_G)
    nx.draw_networkx(en_G,pos=nx.spring_layout(en_G,iterations=500))
    pylab.savefig(fig_out)
    pylab.close()
    return a
    pass

def PCA_paper(data_,samples=None,probes=None, train=None, out='Mapped_PCA.pdf',interactive=False,Text_Labels=False,\
              PCs=10,flow=None, no_plot=False, font_size=3 , markers = 1, colors=0,loc ='best', show_neighbor=False , \
              several_ds=False, nvoisins = 2 , in_3d = False, color_code_gene=None ,recompute_voisins=False, show_nl=False,\
              marker_size=120):
    """docstring for PCA_paper
        runs a PCA
        also does a lot of other stuffs.
    """


    import numpy
    import copy
    import mdp
    import pylab
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse, Polygon
    data = copy.copy(data_)
#   from scikits.learn import fastica

    if probes == None:
        print 'addidng probes'
        probes = numpy.arange(data['data'].shape[0])
        print probes
    if samples == None:
        print 'addidng sample'
        samples = numpy.arange(data['data'].shape[1])
        print samples
    if train == None and samples != None:
        print 'Training on all set'
        train = numpy.arange(len(samples))
        print train
    if train == None:
        print 'Training on all set'
        train = numpy.arange(data['data'].shape[1])
        print samples

    print 'Using %d probes'%len(probes)
    print 'Using %d samples'%len(samples)
    print 'Using %d samples to train'%len(train)

    data['colnames'] = numpy.array([  i.strip('*').replace(' (P)','') for i in data['colnames'] ])

    if markers :
        markers = [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    else:
        markers = ['s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ]
    if colors :
        colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
    else:
        colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#faebd7','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' ,'#faebd7', '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
        colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    X = numpy.matrix(data['data'][probes,:][:,train])
#   X = X - X.mean(0).repeat(X.shape[0], 0)
#   X = (X - numpy.mean(X, 0))/numpy.std(X, axis=0, ddof=0)
    X = numpy.array(X).T

    if flow==None:
        flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs,svd=False)])
        flow = mdp.nodes.PCANode(output_dim=PCs,svd=False)
        flow.train(X)

    if no_plot:
        flow.stop_training()
        return flow #.execute(X)

    Y = numpy.matrix(data['data'][probes,:][:,samples])
#   Y = (Y - numpy.mean(Y, 0))
    Y=numpy.array(Y)
    #Rescale y to have zero mean and unit variance
#   Y = (Y - numpy.mean(Y, 0))/numpy.std(Y, axis=0, ddof=0)
    Y = numpy.array(Y).T

    print 'Y has shape', str(Y.shape)

    data_2d=numpy.array([flow.execute(Y)[:,0], flow.execute(Y)[:,1] ])
    maxi=data_2d.max()
#   maxi=1

    t1, t2 = flow.execute(Y)[:,0]/maxi, flow.execute(Y)[:,1]/maxi
    print t1.shape
    cats=[]
    cats_names=[]
    for j,i in enumerate(data['stemness'][samples]):
        if i not in cats:
            cats.append(i)
            cats_names.append(data['colnames'][samples][j])
    print cats
    print cats_names
    sd={}
    sd_gmp=0

#
    if show_neighbor:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        print 'to implement better..... but this fuction can\'t really be general anways ...'
        print 'l idee c est de mettre une ligne entre un aml et les n voisins les plus proches.'
        #first, we need to find the index of the normal sampmles.
        nl_index = numpy.array([j for j,i in enumerate(data['colnames']) if i.find('AML')==-1 ])
        for i,stem in enumerate(cats) : #(len(data['colnames'][samples])):
            color  = colors[ stem ]
            alpha = .5
            print cats_names[i]
            if cats_names[i].find('AML')!=-1 or cats_names[i].find('Non-leukemia')!=-1 or cats_names[i].find('NK')!=-1 or cats_names[i].find('MLL')!=-1 or cats_names[i] in ['U937', 'K562', 'K052', 'TF1', 'NB4', 'FKH-1', 'UCSD_AML1', 'THP1', 'NH', 'HL60', 'Kasumi-3', 'HNT-34', 'K051', 'HTB58', 'HEL', 'OIH1', 'U937_AML1_ETO', 'MOLM1']: # we have some aml.
                index = numpy.argwhere(data['colnames'][samples] == cats_names[i]).reshape(-1)
                print 'wWwwWwWwWwWwWwWwWwWwWwWwWw'
                print cats_names[i]
                print data['colnames'][samples][index]
                print 'wWwwWwWwWwWwWwWwWwWwWwWwWw'
                #now i can draw lines
                for aml in index:
                    print 'looking for '+ data['cel_file'][samples][aml]
                    sample = get_sample_from_hdf5(data['cel_file'][samples][aml])
                    if recompute_voisins:
                    #a['data'] = filter_low_genes(a['data'],a['data'].mean())
                        probes_to_use  = select_probes_by_variance(sample, var_thr=4)
                        data_to_use = range(len(sample['colnames']))
                        [T,P,E] = run_PCA(sample,probes=probes_to_use,interacive=False,fit_ellipse=False,Text_Labels=False)
                        [distance , data_to_use ] = PCA_neighbourg(copy.copy(sample),T,E,sample['cel_file'][-1],Variance_dims=6,mode='both',data_to_use=data_to_use)
                        probes_to_use  = select_probes_by_variance(sample, var_thr=2, data_to_use=data_to_use)
                        [T,P,E] = run_PCA(sample,data_to_use,probes_to_use,interacive=False,fit_ellipse=False)
                        [top_n,low_n, raw_fold_change,chimera_sample] = find_neighbour(copy.copy(sample),T,E,sample['cel_file'][-1],data_to_use=data_to_use,Variance_dims=6, n=500)
                        sample['voisins']=[numpy.array(data_to_use),numpy.array(distance)]

#########################  load nls from sample
                    nls=range(len(sample['colnames'])-1)
                    data['data'][:,-(len(sample['colnames'])-1):] = sample['data'][:,nls]


                    X = numpy.matrix(data['data'][probes,:][:,train])
                    X = numpy.array(X).T
                    if flow==None:
                        flow = mdp.Flow([mdp.nodes.PCANode(output_dim=PCs,svd=False)])
                        flow = mdp.nodes.PCANode(output_dim=PCs,svd=False)
                        flow.train(X)
                    if no_plot:
                        return flow.execute(X)
                    Y = numpy.matrix(data['data'][probes,:][:,samples])
                    Y=numpy.array(Y)
                    #Rescale y to have zero mean and unit variance
                    Y = numpy.array(Y).T
                    print 'Y has shape', str(Y.shape)
                    data_2d=numpy.array([flow.execute(Y)[:,0], flow.execute(Y)[:,1] ])
                    maxi=data_2d.max()
                    t1, t2 = flow.execute(Y)[:,0]/maxi, flow.execute(Y)[:,1]/maxi

#########################  load nls from sample




                    voisins = []
                    for ii in sample['cel_file'][sample['voisins'][0][1:1+nvoisins].tolist()]:
                        print '### v ->' + str(ii)
                        if ii+'.gz' in  data['cel_file'][samples]:
                            print 'cool'
                            voisins.append( data['cel_file'][samples].tolist().index(ii.split('.')[0]+'.CEL.gz')  )
                        elif ii in  data['cel_file'][samples]:
                            voisins.append( data['cel_file'][samples].tolist().index(ii.split('.')[0]+'.CEL.gz')  )
#                   [ data['cel_file'][samples].tolist().index(i.split('.')[0]+'.CEL.gz') for i in sample['cel_file'][sample['voisins'][0][1:4].tolist()]]
#here I make lines
                    if 1:
                        for v in voisins:
                            color =  colors[data['stemness'][samples][v] ]
                            pylab.plot( [t1[aml] , t1[v]] , [t2[aml] , t2[v]], '-',  color=color, alpha=.4 )

#or polygones.
                    else:
                        dots = [ [t1[v], t2[v]] for v in voisins ]
                        dots.append( [t1[aml],t2[aml]] )
                        ax1.add_patch(Polygon(dots, closed=True, fill=False,  color=color, alpha=1,hatch='+' ))#hatch='/'))



    if in_3d:
        import numpy as np
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        data_3d=numpy.array([flow.execute(Y)[:,0], flow.execute(Y)[:,1] ,flow.execute(Y)[:,2]])
        maxi=data_3d.max()
        maxi=1
        t1, t2 , t3 = flow.execute(Y)[:,0]/maxi, flow.execute(Y)[:,1]/maxi ,flow.execute(Y)[:,2]/maxi
        #t1, t2 , t3 = flow.execute(Y)[:,0]/flow.execute(Y)[:,0].max(), flow.execute(Y)[:,1]/flow.execute(Y)[:,1].max() ,flow.execute(Y)[:,2]/flow.execute(Y)[:,2].max()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.set_zlabel('PC 3')
        ax.w_xaxis.set_pane_color((.8, .8, .8, 1.0))
        ax.w_yaxis.set_pane_color((.7, .7, .7, 1.0))
        ax.w_zaxis.set_pane_color((.8, .8, .8, 1.0))
        MAX = 40
        for direction in (-1, 1):
            for point in np.diag(direction * MAX * np.array([1,1,1])):
                ax.plot([point[0]], [point[1]], [point[2]], 'w')


        for i,stem in enumerate(cats) : #(len(data['colnames'][samples])):
            marker = markers[ stem ]
            color  = colors[ stem ]
            print color
            alpha = .5
            size=2
            size=120
            index = numpy.argwhere(data['stemness'][samples] == stem).reshape(-1)
            x = t1[index]
            y = t2[index]
            z = t3[index]
            edgecolors='w'
            edgecolors=color
            ax.scatter(x, y, z, marker=marker, c=color, s=size/3., edgecolors = edgecolors , label=cats_names[i], alpha=alpha)
            #ax.set_xlabel('PC 1')
            #ax.set_ylabel('PC 2')
            #ax.set_zlabel('PC 3')
            #ax.text(x.mean(),y.mean()+y.std()*1.2, z.mean(),cats_names[i], color=color)
            #ax.w_xaxis.set_pane_color((.8, .8, .8, 1.0))
            #ax.w_yaxis.set_pane_color((.7, .7, .7, 1.0))
            #ax.w_zaxis.set_pane_color((.8, .8, .8, 1.0))
        pylab.show()
    else:
        for i,stem in enumerate(cats) : #(len(data['colnames'][samples])):

            marker = markers[ stem ]
            color  = colors[ stem ]
            alpha = .5
            size=2
            size=marker_size
            if cats_names[i].find('AML')!=-1 or cats_names[i].find('NK')!=-1:
                size=20
                alpha=.5

            if cats_names[i].find('AML')!=-1 and show_neighbor:
                size=70
                alpha=.8
            #size=120
#           print str(marker), str(color)
            index = numpy.argwhere(data['stemness'][samples] == stem).reshape(-1)
            x = t1[index]
            y = t2[index]

            if cats_names[i].find('GMP')!=-1:
                sd_gmp = numpy.array(flow.execute(Y)[index]).std(axis=0).mean()
            if cats_names[i].find('AML')!=-1:
                sd[ cats_names[i] ] = numpy.array(flow.execute(Y)[index]).std(axis=0).mean()
            if several_ds : # for merged verhaak and mile ds.
                numpy.argwhere(array([int(i.replace('GSM','').split('.')[0]) for i in data['cel_file']]) > 300000).reshape(-1)

            # show_computed nl
            if show_nl and numpy.all([ind in numpy.arange(data['chimera_sample'].shape[0]) for ind in samples[index]]):
                samples_without_chimera = numpy.arange(data['chimera_sample'].shape[0],data['data'].shape[1])
                print 'showing nl ------>  '+','.join(data['colnames'][samples[index]])
                Y_nl = numpy.matrix(data['chimera_sample'].T[probes,:][:,samples[index]])
                Y_nl=numpy.array(Y_nl)
               #Rescale y to have zero mean and unit variance
                Y_nl = numpy.array(Y_nl).T
                print 'chimera has shape ... check maxi nl value---', str(Y_nl.shape)
                #maxi_nl=data_2d_nl.max()
                t1_nl, t2_nl = flow.execute(Y_nl)[:,0]/maxi, flow.execute(Y_nl)[:,1]/maxi
                print index
                x_nl = t1_nl[numpy.arange((len(index)))]
                y_nl = t2[numpy.arange((len(index)))]
                edgecolors='k'
                pylab.scatter(x_nl,y_nl, marker=marker, color=color, s=size/2, edgecolors = edgecolors , label=cats_names[i]+'_nl', alpha=alpha)
                pylab.text(x_nl.mean(),y_nl.mean()+y_nl.std()*1.2,'_nl', color=color)#cats_names[i]+'_nl', color=color)


            if 0:#cats_names[i].find('NK')!=-1:
                from matplotlib import pyplot as P
                print 'et oui!'
                alphas = numpy.log(numpy.array(numpy.array(data['OS'],dtype=int),dtype=float)+1)/numpy.log(3575.)
                alphas[:7] = 0
                color = alphas
                print color
                cmap = P.matplotlib.cm.hot
                norm = P.matplotlib.colors.Normalize(vmin=0, vmax=1)
                pylab.scatter(x,y, marker=marker,c=color, s=size,  alpha=alpha, cmap = cmap, norm=norm )
                pylab.text(x.mean(),y.mean()+y.std()*1.2,cats_names[i], color='k')

            else:
#               size=2
#               alpha=.4
                edgecolors='k'
                if color_code_gene is not None:
                    cm = plt.get_cmap("winter")
                    index_gene = get_probe_index([color_code_gene],data['rownames'])
                    values_gene = data_['data'][index_gene,samples[index]]
                    col = [cm(float(iii)/(12)) for iii in values_gene]
                    pylab.scatter(x,y, marker=marker, c=col, s=size, edgecolors = edgecolors , label=cats_names[i], alpha=alpha)
                    print values_gene
                else:
                    pylab.scatter(x,y, marker=marker, color=color, s=size, edgecolors = edgecolors , label=cats_names[i], alpha=alpha)
                try:
                    pylab.text(x.mean(),y.mean()+y.std()*1.2,cats_names[i], color=color)
                except Exception, e:
                    print e
                    print data['colnames'][index]
    kwarg={'size':5 }
    pylab.legend(markerscale=.5,loc=loc, prop=kwarg)


    pylab.savefig(out,format='PDF')
    if interactive:
        pylab.show()
    else:
        pylab.close()


    for i in sd.keys():
        print 'std ', i , str(sd[i]/(sd_gmp))
    return flow
    pass

def tsne_paper(data_,samples=None,probes=None,  out='TSNE.pdf',interactive=False,Text_Labels=False,\
              PCs=3, font_size=3 , markers = 1, colors=0,loc ='best',  \
              in_3d = False, color_code_gene=None, marker_size=120, perplexity=30,verbose=1, plot=True, use_C_implementation=False):
    """docstring for tsne_paper
        runs a bh tsne
    """

    import bhtsne
    import numpy
    import copy
    import pylab
    import matplotlib.pyplot as plt
    data = copy.copy(data_)
    from sklearn.manifold import TSNE

#   from scikits.learn import fastica

    if probes == None:
        print 'addidng probes'
        probes = numpy.arange(data['data'].shape[0])
        #print probes
    if samples == None:
        print 'addidng sample'
        samples = numpy.arange(data['data'].shape[1])
        #print samples

    print 'Using %d probes'%len(probes)
    print 'Using %d samples'%len(samples)

    data['colnames'] = numpy.array([  i.strip('*').replace(' (P)','') for i in data['colnames'] ])

    if markers :
        markers = [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    else:
        markers = ['s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ]
        markers = ['o' , '^' , 'o' , '^' , 'o' , '^' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , '*','o' , '^' , 'o' , '^' , 'o' , '^' ,  'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ]
    if colors :
        colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
    else:
        colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#faebd7','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' ,'#faebd7', '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
        colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    X = data['data'][probes,:][:,samples].T

    if use_C_implementation:
        Y = numpy.array([i for i in bhtsne.bh_tsne(X, perplexity=perplexity,no_dims=PCs,verbose=verbose)])
    else:
        print 'Using sklearn TSNE'
        model = TSNE(n_components=3, random_state=999, verbose=1,perplexity=perplexity,init='pca', learning_rate = 1000,early_exaggeration = 40, n_iter=10000, n_iter_without_progress = 100)
        #model = TSNE(n_components=3, random_state=999, verbose=1, perplexity=perplexity, init='random', n_iter_without_progress = 1000)

        Y=model.fit_transform(X)


    print 'Y has shape', str(Y.shape)

    data_2d=numpy.array([Y[:,0], Y[:,1] ])
    maxi=data_2d.max()
#   maxi=1

    t1, t2 = Y[:,0]/maxi, Y[:,1]/maxi
    if plot == False:
        Y -= Y.min()
        Y /= Y.max()
        return Y
    cats=[]
    cats_names=[]
    for j,i in enumerate(data['stemness'][samples]):
        if i not in cats:
            cats.append(i)
            cats_names.append(data['colnames'][samples][j])
    print cats
    print cats_names
    sd={}
    sd_gmp=0

#

    if in_3d:
        import numpy as np
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        data_3d=numpy.array([Y[:,0], Y[:,1] ,Y[:,2]])
        maxi=data_3d.max()
        maxi=1
        t1, t2 , t3 = Y[:,0]/maxi, Y[:,1]/maxi ,Y[:,2]/maxi
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.set_zlabel('PC 3')
        ax.w_xaxis.set_pane_color((.8, .8, .8, 1.0))
        ax.w_yaxis.set_pane_color((.7, .7, .7, 1.0))
        ax.w_zaxis.set_pane_color((.8, .8, .8, 1.0))
        MAX = 40
        for direction in (-1, 1):
            for point in np.diag(direction * MAX * np.array([1,1,1])):
                ax.plot([point[0]], [point[1]], [point[2]], 'w')


        for i,stem in enumerate(cats) : #(len(data['colnames'][samples])):
            marker = markers[ stem ]
            color  = colors[ stem ]
            alpha = .5
            size=2
            size=120
            index = numpy.argwhere(data['stemness'][samples] == stem).reshape(-1)
            x = t1[index]
            y = t2[index]
            z = t3[index]
            edgecolors='k'
            ax.scatter(x, y, z, marker=marker, c=color, s=size/3., edgecolors = edgecolors , label=cats_names[i], alpha=alpha)
        pylab.show()
    else:
        for i,stem in enumerate(cats) : #(len(data['colnames'][samples])):
            marker = markers[ stem ]
            color  = colors[ stem ]
            alpha = .5
            size=2
            size=marker_size

            index = numpy.argwhere(data['stemness'][samples] == stem).reshape(-1)
            x = t1[index]
            y = t2[index]

            if cats_names[i].find('GMP')!=-1:
                sd_gmp = numpy.array(Y[index]).std(axis=0).mean()
            if cats_names[i].find('AML')!=-1:
                sd[ cats_names[i] ] = numpy.array(Y[index]).std(axis=0).mean()

            if 0:#cats_names[i].find('NK')!=-1:
                0
            else:
#               size=2
#               alpha=.4
                edgecolors='k'
                if color_code_gene is not None:
                    cm = plt.get_cmap("winter")
                    index_gene = get_probe_index([color_code_gene],data['rownames'])
                    values_gene = data_['data'][index_gene,samples[index]]
                    col = [cm(float(iii)/(12)) for iii in values_gene]
                    pylab.scatter(x,y, marker=marker, c=col, s=size, edgecolors = edgecolors , label=cats_names[i], alpha=alpha)
                    print values_gene
                else:
                    pylab.scatter(x,y, marker=marker, color=color, s=size, edgecolors = edgecolors , label=cats_names[i], alpha=alpha)
                try:
                    pylab.text(x.mean(),y.mean()+y.std()*1.2,cats_names[i], color=color)
                except Exception, e:
                    print e
                    print data['colnames'][index]
    kwarg={'size':5 }
    pylab.legend(markerscale=.5,loc=loc, prop=kwarg)


    pylab.savefig(out,format='PDF')
    if interactive:
        pylab.show()
    else:
        pylab.close()


    for i in sd.keys():
        print 'std ', i , str(sd[i]/(sd_gmp))

    if in_3d:
        return t1,t2,t3
    else:
        return t1,t2

    pass
def beeswarm(data, names,out='beeswarm.pdf',title='',plot_median=False, use_mean=False, order=None):
    """docstring for beeswarm
    plots a beeswram plot.
    data is an array or values
    names is the corresponding names for each data values.

    """
    from numpy import array
    import pylab
    import pandas as pa
    import seaborn as sns

    if (order !=None) and all((numpy.unique(order) == numpy.unique(names))):
        print 'Custom order of smaples.'
        orders =[]
        for i in order:
            orders.append(numpy.argwhere(names == i).ravel() )
        orders = numpy.concatenate(orders)
        data  =  data[orders]
        names = names[orders]

    dat=pa.DataFrame({'Expression':array(data) , 'cathegories':array(names)})
    sns.violinplot(x='Expression' , y='cathegories',data=dat, scale="width", linewidth = 1)
    pylab.savefig('violin_'+out)
    pylab.close()

    data = numpy.array(data)
    print "made violin plot!"
    if use_mean:
        meam_method = numpy.mean
    else:
        meam_method = numpy.median

    in_gene = title
    graphics=importr('graphics')
    beeswarm= importr('beeswarm')
    grdevices = importr('grDevices')
    all_classes_order=[]
    all_classes={}
    for index,classe in enumerate(names):
        if classe not in all_classes:
            all_classes_order.append(classe)
            all_classes[classe]=[]
            all_classes[classe].append(index)
        else:
            all_classes[classe].append(index)


    nn=robjects.StrVector(numpy.array(names))
    fcdf=robjects.FloatVector(numpy.array(data))
    cool=robjects.DataFrame({'Expression': fcdf , 'Class':nn})
    Classes= robjects.r.factor(nn, levels = robjects.r.unique(nn), ordered = 1) #robjects.r.factor(cool.rx2(2), levels = robjects.r.unique(cool.rx2(2)), ordered = 1)
    robjects.globalenv['Classes']=Classes

    a=beeswarm.beeswarm( robjects.r("Expression ~ Classes") , data = cool, method = "swarm",pch=1, bg=robjects.r("rainbow(6)") ,col = robjects.r("rainbow(8)") , las =2 )
    opt={'main':in_gene}
    graphics.title( **opt)
    grdevices.dev_off()

    x=numpy.array(a[0])
    y=numpy.array(a[1])
    c=numpy.array(a[3])

    fig=pylab.figure()
    for i,j in enumerate(x):
        pylab.plot(x[i],y[i],'o',color=c[i][0:-2],alpha=.7) # R adds FF at the end of the color string, which is bad.
    if plot_median:
        names =  numpy.array(names)
        cat=[]
        for i in names:
            if i not in cat:
                cat.append(i)
        for ii,i in enumerate(cat):
            z= numpy.argwhere(names == i).ravel()
            m=meam_method(data[z])
            pylab.plot([ii+1-.3, ii+1+.3],[m,m],'b',linewidth=1.5, alpha=.7)
#       pylab.plot([1,2],[0,0],'b',linewidth=1.5,alpha=.5)
    pylab.xticks(range(len(all_classes_order)+1), numpy.concatenate([[''],all_classes_order],axis=0), rotation=90,size=9)
    fig.autofmt_xdate(bottom=0.18)
#    pylab.axis((0,len(all_classes),4,14))
    pylab.ylabel('Median centered expression (log2)')
    pylab.title(in_gene)
    pylab.savefig(out)
    pylab.close()

    pass
def filter_low_genes(data,thr):
    """docstring for filter_low_genes
        takes a matrix, and put everything that's below the threshold at the thr value.
    """
    indexes=data < thr
    data[indexes] = thr

    return data
    pass

def stemness_by_class(data,patient_ds=0):
    """docstring for stemness_by_class
        put stemness values, based on labels.
    """
    colnames =  data['colnames']
    cat=[]
    for n in colnames:
        if n not in cat:
            print 'adding '+ n
            cat.append(n)
    cat=numpy.unique(colnames).tolist()
    stemness=[]
    for ii,n in enumerate(colnames):
        #print 'i' + str(n)
        stemness.append(cat.index(n))
    stemness =numpy.array(stemness)
    if patient_ds:
        for index,i in enumerate(colnames):
            print index
            stemness[index] = 10
            if i.find('HSC')!=-1 or i.find('MPP')!=-1:
                stemness[index] = 0
            if i.find('CMP')!=-1:
                stemness[index] = 1
            if i.find('GMP')!=-1:
                stemness[index] = 2
    return stemness
    pass

def get_date_from_cel_file(data):
    """ for i in `ls -1 | grep -v AML | grep -v KTh | grep -v APL` ; do echo $i "   "   `gunzip --stdout  $i | grep --text d.a.t.e` ; done
        or

        for i in `ls -1  | grep -v KTh ` ; do echo $i ; apt-cel-convert -f text -o . $i; done
        and
        for i in `ls -1 | grep -v AML | grep -v KTh | grep -v APL` ; do echo $i "    " `grep ^DatHeader  $i|awk '{print $8}'`; done

    """
    from subprocess import check_output
    dates=[]
    for i in data['cel_file']:
        string = "gunzip --stdout  %s  | grep --text d.a.t.e "%(i)
        print i
        dates.append(check_output(string, shell=True).replace('\x00','').replace('text/plain\x14affymetrix-scan-date(','')[:10])
    return numpy.array(dates)
    pass

def mayavi_PCA(T,data,data_to_use=None):
    """docstring for mayavi_PCA
        this takes the loadings from the pca (ie. the transposed coordinates.)
    """

    if data_to_use==None:
        data_to_use = range(len(data['colnames']))

    curdir=os.path.realpath('.')
    dic={}
    dic['points']=[]
    dic['name']=[]
    dic['stemness']=[]
    for j,i in enumerate(data_to_use):
        dic['points'].append(T[j][0:3])
        dic['name'].append(data['colnames'][i])
        dic['stemness'].append(int(data['stemness'][i]))
    dic['points']=numpy.array(dic['points'])
    os.chdir('/Users/nrapin/Dropbox/Python-stuff/')
    pickle_save(dic,'_pca3d.pkl')
    os.popen('/Library/Frameworks/Python.framework/Versions/Current/bin/python /Users/nrapin/Dropbox/Python-stuff/Mayavi_PCA.py')
    os.chdir(curdir)
#   return dic
    pass

def run_PCA(data,samples=None,probes=None,out='PCA.pdf',interacive=False,fit_ellipse=False,Text_Labels=True,PCs=10,colors = 0, fontsize=10, no_plot=False):
    """docstring for run_PCA
        runs a pca.
    """
    import copy
    import pca_module
    import pylab

    if probes == None:
        probes = range(data['data'].shape[0])
        print 'Adding probes...'
    if samples == None:
        samples = range(data['data'].shape[1])
        print 'Adding samples...'
    #    pylab.matplotlib.use('cairo.pdf')

#   markers= [ 'o', '*', '<', '>',  '+', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' , '+', '*', '<', '>',  'D', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' ]
#   colors=  ['r', 'b', 'g', 'r', 'k', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k' , 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k']
    markers= [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    if colors == 1:
        colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
    else:
        colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    X = copy.copy( numpy.matrix(data['data'][probes,:][:,samples]) )
    X = X - X.mean(0).repeat(X.shape[0], 0)
    X = numpy.array(X)
#   T, P, explained_var = pca_module.PCA_svd(merged_data['data'][probes_to_use,:][:,data_to_use].transpose(), standardize=True)
    T, P, E = pca_module.PCA_nipals2(X.transpose(), standardize=False, E_matrices=False,PCs=PCs)
    if no_plot:
        return T, P , E
    #use first two components
    data_2d=numpy.array([T[:,0], T[:,1]])
    maxi=data_2d.max()
#   maxi=1
    t1, t2 = T[:,0]/maxi, T[:,1]/maxi
    for i in range(len(data['colnames'][samples])):
#       print data[ samples[i] ][11]
        marker=markers[ int(data['stemness'][samples][i]  ) ]
        color  = colors[ int(data['stemness'][samples][i]  ) ]

        if Text_Labels:

            #pylab.text(t1[i],t2[i],data['cel_file'][samples][i], fontsize=6)
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=fontsize,withdash=False,)
            if 0:#data['cel_file'][samples][i].find('100531') == -1:
                pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.5)
            else:
                #print data['colnames'][samples][i]
                pylab.plot([t1[i]],[t2[i]], marker='o', color=color, markersize=5, label=unicode(data['colnames'][samples][i],errors='ignore'), alpha=.7)

        else:
            if data.has_key('types')!=True:
                data['types']=data['colnames']
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['types'][samples][i], alpha=.7 )
#           pylab.legend([data['types'][samples][i]])
#           print data['colnames'][samples][i]
#           pylab.legend(data['colnames'][samples][i])
    if not Text_Labels: #then we need some sort of class detection an labelling.
        types = {} # determine types, based on names.
        if data.has_key('types'):
            for i,s in enumerate(data['types'][samples]):
                if s not in types.keys():
                    types[s]=int(data['stemness'][samples][i])
        else:
            types['data']=0
        #now we know which cathegories, and which number they have.
        used_types=[0]*len(types)
        types=types.keys() # temp hack.
        for i in range(len(data['colnames'][samples])):
            #find first of each type, and
            #print data['types'][samples][i]
            #print used_types
            #print
            if used_types[ types.index(data['types'][samples][i]) ]==0:
                pylab.text(t1[i],t2[i],data['types'][samples][i], fontsize=fontsize)#6)
                used_types[ types.index(data['types'][samples][i]) ]=1



    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            try:
                if len(indexes) > 2:
                    eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
#                   eee.append(draw_pca_ellipse( numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]))
                    colors.pop(-1)
            except Exception, e:
                print e

    pylab.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return [T,P,E]
    pass


def run_PCA2(data,samples,probes,out='PCA.pdf',interacive=False,fit_ellipse=True,Text_Labels=True):
    """docstring for run_PCA
        runs a pca.
        here we use the incremental PCA method. instead of the nipals algorithm.
    """
    import pca_module
    import pylab
    #    pylab.matplotlib.use('cairo.pdf')
    markers= [ 'o', '*', '<', '>',  '+', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' , '+', '*', '<', '>',  'D', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' ]
    colors=  ['r', 'b', 'g', 'r', 'k', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k' , 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k']

#   here we use the incremental PCA method.
#   T, P, E = pca_module.PCA_nipals2(data['data'][probes,:][:,samples].transpose(), standardize=False, E_matrices=False)


#   pca = IPCA(data['data'].shape[1], 2)
#   for i in probes:
#       x = numpy.matrix(data['data'][i,:]).T
#       pca.update(x)
#   T = pca.components
    X=numpy.matrix(x[probes,:][:,samples])
    X = X - X.mean(0).repeat(X.shape[0], 0)
    [_, _, T] = numpy.linalg.svd(X)
    t2=numpy.array(T[1:2])[0]
    t1=numpy.array(T[0:1])[0]
#   return t1
    #use first two components
#   data_2d=numpy.array([T[:,0], T[:,1]])
#   maxi=data_2d.max()
#   t1, t2 = numpy.array(T[:,0]/maxi), numpy.array(T[:,1]/maxi)
    for i in range(len(data['colnames'][samples])):
        print i
#       print data[ samples[i] ][11]
        marker=markers[ int(data['stemness'][samples][i]  ) ]
        color  = colors[ int(data['stemness'][samples][i]  ) ]

        if Text_Labels:
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=6)
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.5)
        else:
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['types'][samples][i], alpha=.5)
#           pylab.legend([data['types'][samples][i]])
#           print data['colnames'][samples][i]
#           pylab.legend(data['colnames'][samples][i])
    if not Text_Labels: #then we need some sort of class detection an labelling.
        types = {} # determine types, based on names.
        if data.has_key('types'):
            for i,s in enumerate(data['types'][samples]):
                if s not in types.keys():
                    types[s]=int(data['stemness'][samples][i])
        else:
            types['data']=0
        #now we know which cathegories, and which number they have.
        used_types=[0]*len(types)
        types=types.keys() # temp hack.
        for i in range(len(data['colnames'][samples])):
            #find first of each type, and
            print data['types'][samples][i]
            print used_types
            print
            if used_types[ types.index(data['types'][samples][i]) ]==0:
                pylab.text(t1[i],t2[i],data['types'][samples][i], fontsize=6)
                used_types[ types.index(data['types'][samples][i]) ]=1

    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors = ['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','m','b','c','y','k']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            try:
                if len(indexes) > 2:
                    eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
                    colors.pop(-1)
            except Exception, e:
                print e

    pylab.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return [T]


    pass




def run_PCA3d(data,samples,probes,out='PCA.pdf',interacive=False,fit_ellipse=False,Text_Labels=True):
    """docstring for run_PCA
        runs a pca - plot in 3D.
    """

    import pca_module
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    markers= [ 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
    colors=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]

#   T, P, explained_var = pca_module.PCA_svd(merged_data['data'][probes_to_use,:][:,data_to_use].transpose(), standardize=True)
    T, P, E = pca_module.PCA_nipals2(data['data'][probes,:][:,samples].transpose(), standardize=False, E_matrices=False)
    #use first two components
    data_3d=numpy.array([T[:,0], T[:,1], T[:,2]] )
    maxi=data_3d.max()
    t1, t2 , t3 = T[:,0]/maxi, T[:,1]/maxi , T[:,2]/maxi

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(data['colnames'][samples])):
        marker=markers[ int(data['stemness'][samples][i]  ) ]
        color  = colors[ int(data['stemness'][samples][i]  ) ]

        if Text_Labels:
#           pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=6)
            print data['colnames'][samples][i].strip('_')
            ax.text(t1[i], t2[i], t3[i], data['colnames'][samples][i].strip('_'))
        else:
            ax.scatter([t1[i]], [t2[i]], [t3[i]], marker=marker, color=color, label=data['types'][samples][i].strip('_'), alpha=.5)
            ax.legend([data['types'][samples][i]])
#           print data['colnames'][samples][i]
#           pylab.legend(data['colnames'][samples][i])
    if not Text_Labels: #then we need some sort of class detection an labelling.
        types = {} # determine types, based on names.
        if data.has_key('types'):
            for i,s in enumerate(data['types'][samples]):
                if s not in types.keys():
                    types[s]=int(data['stemness'][samples][i])
        else:
            types['data']=0
        #now we know which cathegories, and which number they have.
        used_types=[0]*len(types)
        types=types.keys() # temp hack.
        for i in range(len(data['colnames'][samples])):
            #find first of each type, and
            print data['types'][samples][i]
            print used_types
            print
            if used_types[ types.index(data['types'][samples][i]) ]==0:
                ax.text(t1[i],t2[i],t3[i],data['types'][samples][i], fontsize=6)
                used_types[ types.index(data['types'][samples][i]) ]=1

    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors = ['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','m','b','c','y','k']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            try:
                if len(indexes) > 2:
                    eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
                    colors.pop(-1)
            except Exception, e:
                print e

    fig.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return [T,P,E]


    pass



def run_PCA_OS(data,samples,probes,OS,out='PCA.pdf',interacive=False,fit_ellipse=True,T=0,P=0,E=0):
    """docstring for run_PCA
        runs a pca.
    """
    a=[i[1] for i in OS]
    OS_d=[int(i) for i in a]
    maxOS=max(OS_d)
    color='black'
    import pca_module
    import pylab
    #    pylab.matplotlib.use('cairo.pdf')

    markers= [ 'o', '*', '<', '>',  '+', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' , '+', '*', '<', '>',  'D', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' ]
    colors=  ['r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k' , 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k']

#   T, P, explained_var = pca_module.PCA_svd(merged_data['data'][probes_to_use,:][:,data_to_use].transpose(), standardize=True)
    print type(T)
    if type(T) == type(0):
        print 'running PCA'
        T, P, E = pca_module.PCA_nipals2(data['data'][probes,:][:,samples].transpose(), standardize=False, E_matrices=False)
    #use first two components
    data_2d=numpy.array([T[:,0], T[:,1]])
    maxi=data_2d.max()
    t1, t2 = T[:,0]/maxi, T[:,1]/maxi
    for i in range(len(data['colnames'][samples])):
        if t1[i]<0:
#           print data['colnames'][samples][i]
            print i
#       print data[ samples[i] ][11]
        if int(data['stemness'][samples][i]) != 22:
            print 'pas cancer'
            marker=markers[ int(data['stemness'][samples][i]  ) ]
            color  = colors[ int(data['stemness'][samples][i]  ) ]
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.5)
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=6)
        elif (int(data['stemness'][samples][i]) == 22 ): # that's cancer data- with mapped OS!
            marker=markers[ int(data['stemness'][samples][i]  ) ]
            #find color.
            for index,patient in enumerate(OS):
                if str(patient[0])+'.CEL'== data['colnames'][samples][i]:
                    survie = int(patient[1])
                    print '->' + str(float(float(survie)/float(maxOS)))
                    print survie
                    print maxOS
                    color=(float(float(survie)/float(maxOS)),float(float(survie)/float(maxOS)),float(float(survie)/float(maxOS)))
                    if survie < 150:
                        color='r'
                    elif survie < 365:
                        color = 'y'
                    elif survie < 800:
                        color = 'g'
                    else:
                        color = 'b'
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.4)

    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors = ['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','m','b','c','y','k']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            if cat != '22':
                try:
                    if len(indexes) > 2:
                        eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
                        colors.pop(-1)
                except Exception, e:
                    print e

    pylab.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return [T,P,E]


    pass

def run_PCA_AML_NK(data,samples,probes,OS,out='PCA.pdf',interacive=False,fit_ellipse=True,T=0,P=0,E=0):
    """docstring for run_PCA
        runs a pca.
    """
    a=[i[1] for i in OS]
    OS_d=[int(i) for i in a]
    maxOS=max(OS_d)
    color='black'
    import pca_module
    import pylab
    #    pylab.matplotlib.use('cairo.pdf')

    markers= [ 'o', '*', '<', '>',  '+', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' , '+', '*', '<', '>',  'D', 'h', 'o', 'p', 's', 'v','+', '*', '<', '>' ]
    colors=  ['r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k' , 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b','b','k','k']

#   T, P, explained_var = pca_module.PCA_svd(merged_data['data'][probes_to_use,:][:,data_to_use].transpose(), standardize=True)
    print type(T)
    if type(T) == type(0):
        print 'running PCA'
        T, P, E = pca_module.PCA_nipals2(data['data'][probes,:][:,samples].transpose(), standardize=False, E_matrices=False)
    #use first two components
    data_2d=numpy.array([T[:,0], T[:,1]])
    maxi=data_2d.max()
    t1, t2 = T[:,0]/maxi, T[:,1]/maxi
    for i in range(len(data['colnames'][samples])):
        if t1[i]<0:
#           print data['colnames'][samples][i]
            print i
#       print data[ samples[i] ][11]
        if int(data['stemness'][samples][i]) != 22:
            print 'pas cancer'
            marker=markers[ int(data['stemness'][samples][i]  ) ]
            color  = colors[ int(data['stemness'][samples][i]  ) ]
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.5)
            pylab.text(t1[i],t2[i],data['colnames'][samples][i], fontsize=6)
        elif (int(data['stemness'][samples][i]) == 22 ): # that's cancer data- with mapped OS!
            marker=markers[ int(data['stemness'][samples][i]  ) ]
            #find color.
            for index,patient in enumerate(OS):
                if str(patient[0])+'.CEL'== data['colnames'][samples][i]:
                    survie = int(patient[1])
                    print '->' + str(float(float(survie)/float(maxOS)))
                    print survie
                    print maxOS
                    color=(float(float(survie)/float(maxOS)),float(float(survie)/float(maxOS)),float(float(survie)/float(maxOS)))
                    if survie < 150:
                        color='r'
                    elif survie < 365:
                        color = 'y'
                    elif survie < 800:
                        color = 'g'
                    else:
                        color = 'b'
            pylab.plot([t1[i]],[t2[i]], marker=marker, color=color, markersize=5, label=data['colnames'][samples][i], alpha=.4)

    if fit_ellipse:
        print 'Fitting ellipses.'
        cathegories=[]
        for s in data['stemness']:
            if s not in cathegories:
                cathegories.append(s)

#       cathegories=['11']
        colors = ['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','m','b','c','y','k']
        for cat in cathegories:
            indexes=[j for j,i in  enumerate(data['stemness']) if i == cat]
            print 'cat. %s as %d elements' %(cat,len(indexes))
            print indexes
            eee=[]
#           pylab.plot(numpy.array([t1[indexes]]), numpy.array([t2[indexes]]) )
            if cat != '22':
                try:
                    if len(indexes) > 2:
                        eee.append( fit_ellipse_2d(numpy.array([t1[indexes], t2[indexes]]).T,color=colors[-1]) )
                        colors.pop(-1)
                except Exception, e:
                    print e

    pylab.savefig(out,format='PDF')
    if interacive:
        pylab.show()
    else:
        pylab.close()

    return [T,P,E]


    pass

def _is_zero(x):
    """Return a boolean indicating whether the given vector is a zero vector up
    to a threshold.
    """
    return numpy.fabs(x).min() < 1.e-9

class IPCA(object):
    """Incremental PCA calculation object.

    General Parameters:
        m - Number of variables per observation
        n - Number of observations
        p - Dimension to which the data should be reduced
    """

    def __init__(self, m, p):
        """Creates an incremental PCA object for m-dimensional observations
        in order to reduce them to a p-dimensional subspace.

        @param m: Number of variables per observation.
        @param p: Number of principle components.

        @return: An IPCA object.
        """
        self._m = float(m)
        self._n = 0.0
        self._p = float(p)
        self._mean = numpy.matrix(numpy.zeros((m , 1), dtype=numpy.float64))
        self._covariance = numpy.matrix(numpy.zeros((m, m), dtype=numpy.float64))
        self._eigenvectors = numpy.matrix(numpy.zeros((m, p), dtype=numpy.float64))
        self._eigenvalues = numpy.matrix(numpy.zeros((1, p), dtype=numpy.float64))

    def update(self, x):
        """Updates with a new observation vector x.

        @param x: Next observation as a column vector (m x 1).
        """
        m = self._m
        n = self._n
        p = self._p
        mean = self._mean
        C = self._covariance
        U = self._eigenvectors
        E = self._eigenvalues

        if type(x) is not numpy.matrix or x.shape != (m, 1):
            raise TypeError('Input is not a matrix (%d, 1)' % int(m))

        # Update covariance matrix and mean vector and centralize input around
        # new mean
        oldmean = mean
        mean = (n*mean + x) / (n + 1.0)
        C = (n*C + x*x.T + n*oldmean*oldmean.T - (n+1)*mean*mean.T) / (n + 1.0)
        x -= mean

        # Project new input on current p-dimensional subspace and calculate
        # the normalized residual vector
        g = U.T*x
        r = x - (U*g)
        r = (r / numpy.linalg.norm(r)) if not _is_zero(r) else numpy.zeros_like(r)

        # Extend the transformation matrix with the residual vector and find
        # the rotation matrix by solving the eigenproblem DR=RE
        U = numpy.concatenate((U, r), 1)
        D = U.T*C*U
        (E, R) = numpy.linalg.eigh(D)

        # Sort eigenvalues and eigenvectors from largest to smallest to get the
        # rotation matrix R
        sorter = list(reversed(E.argsort(0)))
        E = E[sorter]
        R = R[:,sorter]

        # Apply the rotation matrix
        U = U*R

        # Select only p largest eigenvectors and values and update state
        self._n += 1.0
        self._mean = mean
        self._covariance = C
        self._eigenvectors = U[:, 0:p]
        self._eigenvalues = E[0:p]

    @property
    def components(self):
        """Returns a matrix with the current principal components as columns.
        """
        return self._eigenvectors

    @property
    def variances(self):
        """Returns a list with the appropriate variance along each principal
        component.
        """
        return self._eigenvalues


def fast_index(list,item):
    """docstring for fast_index
        reproduces the [...].index(x), fast.
    """
    for pos,elem in enumerate(list):
        if elem == item:
            return pos

    raise ValueError("list.index(x): x not in list")

    pass


def mahalanobis(g1,g2):
    """docstring for mahalanobis
    computes mahalanobis distance between two multivariate populations g1 and g2
    g1 and g2:

    rows: samples;
    cols: observations (genes, etc..);


    input:
    g1 = multivariate population 1.
    g2 = multivariate population 2.

    """
    n2=n1=1
    g1=numpy.matrix(g1)
    g2=numpy.matrix(g2)

    if g1.shape[1] != g2.shape[1]:
        print 'No way. matrix shapes don\'t agree (different number of cols)'
        return None

    if g2.shape[0] > 1:
        n2=g2.shape[0]
    if g1.shape[0] > 1:
        n1=g1.shape[0]

    variables=g1.shape[1]

    if g1.shape[0] > 1:
        m1=numpy.mean(g1,axis=0)
        cg1=numpy.array(numpy.subtract(g1,m1))
        cov1=numpy.cov(cg1,rowvar=0)
    else:
        cg1=numpy.array(g1)
        m1=numpy.array(g1)
        cov1=numpy.identity(variables)

    if g2.shape[0] > 1:
        m2=numpy.mean(g2,axis=0)
        cg2=numpy.array(numpy.subtract(g2,m2))
        cov2=numpy.cov(cg2,rowvar=0)
    else:
        cg2=numpy.array(g2)
        m2=numpy.array(g2)
        cov2=numpy.identity(variables)

    COV= numpy.matrix( (cov1*n1+cov2*n2)/(n1+n2) )
    COVI=  COV.I

    mean_diff=numpy.matrix(m1 - m2)

    return float(mean_diff*COVI*mean_diff.T)

    pass


def make_phylip_dist_file(data,names,probes=None,dist='c',number_names=False):
    """docstring for make_phylip_dist_file
    computes distmatrix, and out put the file as infile for neighbourgs.
    """

    if number_names:
        a=[o.replace(' ','0') for o in     ['%4d'%(i+1) for i in range(len(names))]]
        names=a
    if probes==None:
        probes=range( numpy.array(data).shape[0] )


    nb_entries=len(names)
    matrix= numpy.zeros([nb_entries,nb_entries])
    f=open('infile','w')
    f.write('%s\n'%str(nb_entries))


    print 'Computing distance:'
    prog = ProgressBar(0, nb_entries , 50, mode='fixed')
#   Bio.Cluster.distancematrix(merged_data['data'].transpose(),dist='s')
    for i,index_i in enumerate(range(nb_entries)):
        for j,index_j in enumerate(range(nb_entries)):
            if matrix[j,i]==0 and matrix[i,j]==0 :
#               matrix[i,j]= ( hcluster.correlation(data[probes,index_i],data[probes,index_j]) )
                matrix[i,j]= (Bio.Cluster.distancematrix(( ( data[probes,index_i] ) ,  (data[probes,index_j])  ) , dist=dist)[1][0])
                matrix[j,i]=matrix[i,j]
            else:
                continue
                continue
        #progress bar code..
        oldprog = str(prog)
        prog.update_amount(i)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)
#make list of names (uniques!):
    for i in range(nb_entries):
        names[i]= ''.join([ c for c in names[i] if c not in ( '?' )])

    unique_names=[]
    for i in names:
        if i not in unique_names:
            unique_names.append(i)
    dup_key_index=[0]*len(unique_names)

    for i,name in enumerate(names):
        if name in unique_names:
            index=unique_names.index(name)
            dup_key_index[index]+=1

    mono_names=[]
    for i in names:
        index=unique_names.index(i)
        left_index=dup_key_index[index]
        if left_index - 1:
            dup_key_index[index]-=1
            mono_names.append(str(left_index)+'_'+str(i))
        else:
            mono_names.append(i)

    names=mono_names
    for i in range(nb_entries):
        names[i]= ''.join([ c for c in names[i] if c not in ( '(' , ')' , ':' , ';' , ',' , '[' , ']')])
        l=len(names[i])
        x = 10 - l
        if x < 0:
            x=0
        f.write('%s '%str(  '_'.join([o for o in names[i].split()]).strip('():;,[]')[0:10]   +' '*x ) )
        for j in range(nb_entries):
            f.write('%5f '%(matrix[i,j]) )
        f.write('\n')
    f.close()
    pass



def PCA_neighbourg(merged_data,T,E,sample,Variance_dims=3,data_to_use=None,n=15,mode='index',Sort=True):
    """docstring for PCA_neighbourg
        gives back list of neighbour
        mode is : index, ditance or both.

    """
    import copy
    c_merged_data=copy.deepcopy(merged_data)
    if data_to_use==None: #if parameter not given as argument, then use everything
        data_to_use=range(len(merged_data['stemness']))
        print 'hehe'
    else:
        for key in [ 'platforms', 'cel_file','colnames', 'data', 'stemness']:
            print key
            try:
                c_merged_data[key]=copy.copy(merged_data[key][data_to_use])
            except Exception,e:
                c_merged_data[key]=copy.copy(merged_data[key][:,data_to_use])
    #makes whiegted geometric distance
    #find where the sample is-


    sample_index=None
    for i,pl in enumerate( c_merged_data['stemness']):
        print pl
        if  c_merged_data['cel_file'][i]==sample :
            sample_index=i
            print 'found! -> '+ c_merged_data['cel_file'][i]
            break
#       if  pl==str(19) or pl==str(13) or merged_data['cel_file'][i]==sample :
#           sample_index=i
#           print 'found! -> '+ merged_data['cel_file'][i]


    if sample_index==None:
        print 'Sample not found...'
        return -1

    #compute distance
    dist=[]
    for i,array in enumerate(c_merged_data['colnames']):
        print array,i
        sum_v=0
        for j in range(Variance_dims):
            sum_v+= numpy.square(T[sample_index][j]-T[i][j])# * E[j]
        dist.append(numpy.sqrt(sum_v))

    #sort data according to that.
    x=stats.rankdata(dist)-1
    ranks=x.tolist()
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()

    if Sort:
        order=[x.index(i) for i in sorted( x , key=lambda i:( int(i[0]),i[1] ) )]
    else:
        order = range(len(x))

    dist=numpy.array(dist)

    print 'Clotest samples to %s:' % sample
    print '\n'.join([ '%d\t%s\t%f'%(i, c_merged_data['colnames'][order][i],  dist[order][i]) for i in range(len(c_merged_data['colnames'])) if (i < n+1) and (i != 0 )])

    d=[ dist[order[i]] for i in range(len(c_merged_data['colnames'])) if (i < n+1) ]
    ind=[ order[i] for i in range(len(c_merged_data['colnames'])) if (i < n+1) ]
    if mode=='index':
        return ind
    elif mode=='distance':
        return d
    elif mode =='both':
        return d,ind
    pass
#
def mapped_neighbourgs(data,pca_node,sample,PCA_dims=3,probes_to_use=None,data_to_use=None,n=15,mode='index',Sort=True):
    """docstring for PCA_neighbourg
        gives back list of neighbour
        mode is : index, ditance or both.
    """

    if data_to_use==None: #if parameter not given as argument, then use everything
        data_to_use=range(len(data['colnames']))
    if probes_to_use==None: #if parameter not given as argument, then use everything
        probes_to_use=range(len(data['rownames']))
    if pca_node.input_dim != len(probes_to_use) :
        print 'PCA Node trained on different number of probes.'
        print ' provide correct probes_to_use.'
        return -1
    if n == None:
        n = len(data_to_use)

    data_to_use.sort()
    #makes whiegted geometric distance
    #find where the sample is-


    sample_index=None
    for cat in ['cel_file', 'stemness','colnames']:
        try:
#           print 'looking in '+ cat
            if sample in data[cat]:
#               print 'it\'s in there'
                sample_index=data[cat][data_to_use].tolist().index(sample)
                print sample + ' has index ' + sample_index
        except Exception, e:
            pass
#           print 'missing ' + cat
    if sample_index==None:
        print 'Sample not found...'
        return -1

    #compute distance
    T = pca_node.execute(data['data'][probes_to_use,:][:,data_to_use].T)
    dist=[]
    for i,array in enumerate(data_to_use):
        sum_v=0
        for j in range(PCA_dims):
#           print j,i
            sum_v+= numpy.square(T[sample_index][j]-T[i][j])
        dist.append(numpy.sqrt(sum_v))

    #sort data according to that.
    x=stats.rankdata(dist)-1
    ranks=x.tolist()
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()

    if Sort:
        order=[x.index(i) for i in sorted( x , key=lambda i:( int(i[0]),i[1] ) )]
    else:
        order = range(len(x))

    dist=numpy.array(dist)

    print 'Clotest samples to %s:' % sample
    print '\n'.join([ '%d\t%s\t%f'%(i, data['colnames'][data_to_use][order][i],  dist[order][i]) for i in range(len(data['colnames'][data_to_use])) if (i < n+1) and (i != 0 )])

    d = [ dist[order[i]] for i,j in enumerate(data_to_use) if (i < n+1) ]
    ind = numpy.array(data_to_use)[order[:n]]

    if mode=='index':
        return ind
    elif mode=='distance':
        return d
    elif mode =='both':
        return d,ind
    pass

def find_neighbour(merged_data,T,E,sample,data_to_use=None,Variance_dims=3,disp=0,n=100):
    """docstring for find_neighbour
    returns a list of n, up and down regulated genes. from a merged sample.
    the function takes care fo finding the neighbour samples.
    """

    if data_to_use==None: #if parameter not given as argument, then use everything
        data_to_use=range(len(merged_data['stemness']))
        print 'hehe'
    else:
        for key in [  'cel_file','colnames', 'stemness']:
            print key
            print data_to_use
            merged_data[key]=merged_data[key][data_to_use]
        merged_data['data']=merged_data['data'][:,data_to_use]

#makes weihgted geometric distance
    #find where the sample is-
    print 'looking for %s:'%sample
    sample_index=None
    for i,pl in enumerate( merged_data['stemness']):
        print pl
        if  merged_data['cel_file'][i]==sample :
            sample_index=i
            print 'found! -> '+ merged_data['cel_file'][i]
            break
        if  pl==str(19) or pl==str(13) or merged_data['cel_file'][i]==sample :
            sample_index=i
            print 'found! -> '+ merged_data['cel_file'][i]


    if sample_index==None:
        print 'Sample not found...'
        return [None,None,None,None]

    #compute distance
    dist=[]
    for i,array in enumerate(merged_data['colnames']):
        sum_v=0
        for j in range(Variance_dims):
            sum_v+= numpy.square(T[sample_index][j]-T[i][j]) * E[j]
        dist.append(numpy.sqrt(sum_v))

    #sort data according to that.
    x=stats.rankdata(dist)-1
    ranks=x.tolist()
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
    order=[x.index(i) for i in sorted( x , key=lambda i:( int(i[0]),i[1] ) )]
    dist=numpy.array(dist)
    print dist[order]

    print 'Clostest samples to %s:' % sample
    print '\n'.join([ '%d\t%s\t%f'%(i, merged_data['colnames'][order][i],  dist[order][i]) for i in range(len(merged_data['colnames'])) if (i < 5) and (i != 0 )])

    liste_array=[]
    for i in [1,2,3,4]:
        liste_array.append([ int(order[i]) ])
    w=[]
    for i in [1,2,3,4]:
        print  numpy.exp( -1.0 * dist[order][i])
        w.append([ numpy.exp( -1.0* dist[order][i] /5) ])
    print w
#   constructs an expression profile based on neighbours and weitghts:
    constructed_profile=numpy.zeros([len(merged_data['data'][:,0]), len(w) ])
    for i,s in enumerate(liste_array):
        constructed_profile[:,i]=merged_data['data'][:,s].reshape(-1)*w[i]
    chimera_sample=numpy.sum(constructed_profile,axis=1) / numpy.sum(w)

    # computes fold change, raw:
    raw_fold_change=numpy.subtract(merged_data['data'][:,sample_index], chimera_sample)
    #filter low expressed genes..
    fold_change=numpy.zeros(raw_fold_change.shape)
    for i in range(len(fold_change)):
        if (chimera_sample[i] < 8.0) and (merged_data['data'][:,sample_index][i] < 8.0):
            fold_change[i]=0
        else:
            fold_change[i]=raw_fold_change[i]

    ranked = stats.rankdata(fold_change).tolist()

    x=numpy.array([ranked,range(len(ranked))]).transpose().tolist()
    v=sorted( x , key=lambda i:( int(i[0]),i[1] ) )
    y={}
    for i,item in enumerate(x):
        y[str(item[0])]=i
    order_probes=[y[str(i[0])] for i in v]

    #find top ranked and bottom ranked:
    top_100=order_probes[0:n]
    low_100= order_probes[-1*(n+1):-1]


    if disp:
        import pylab
        pylab.close()
        pylab.plot(numpy.sort(fold_change[top_100]))
        pylab.plot(numpy.sort(fold_change[low_100]))

    return merged_data['rownames'][top_100],merged_data['rownames'][low_100], raw_fold_change , chimera_sample
    pass






def plot_survival(all_data,data_to_use, cluster_to_use,key_to_use='stemness',no_plot=False):
    """docstring for plot_survival
    take a dic, with OS and event OS plus other attribute.
    plot according to the key.
    """
# Some R packages...
    survival=importr('survival')
    utils=importr('utils')
    if (no_plot == False):
        grdevices = importr('grDevices')

    cat = []
    for i in all_data[key_to_use][data_to_use]:
        if i not in cat:
            cat.append(i)
    nclusters=len(cat)

    print 'found %d clusters.'%nclusters
    print 'using %s'%(key_to_use)
#    print cat
    OS_dic={}
    for key in ['OS','EventOS',key_to_use]:
#   for key in ['OS',key_to_use]:
#        print key
        OS_dic[key]=numpy.array([[]])
        for index in data_to_use:
            if  all_data[key_to_use][index] in cluster_to_use:
                val=all_data[key][index]
#                print key, index
                try:
                    OS_dic[key]=numpy.concatenate((OS_dic[key],numpy.array([[int(val)]])), axis=1)
                except Exception, e:
                    if val=='na':
                        OS_dic[key]=numpy.concatenate((OS_dic[key],numpy.array([[  numpy.nan  ]])), axis=1)
                    else:
                        OS_dic[key]=numpy.concatenate((OS_dic[key],numpy.array([[val]])), axis=1)
        if type(OS_dic[key][0]) != type(' ') :
#            print 'chiffres'
 #           print robjects.Vector(OS_dic[key].T)
            OS_dic[key]=robjects.Vector(OS_dic[key].T)
        else:
 #           print 'lettres'
 #           print robjects.StrVector(OS_dic[key].T)
            OS_dic[key]=robjects.StrVector(OS_dic[key].T)
 #           print OS_dic[key]
    OS_dic['colnames'] = OS_dic[key_to_use]
    R_OS=robjects.DataFrame(OS_dic)
 #   print R_OS
#    print 'Using %d clusters.' %len(cluster_to_use)

    fit=survival.survfit(robjects.r("Surv(OS, EventOS) ~ %s"%key_to_use) ,data=R_OS)
#   fit=survival.survfit(robjects.r("Surv(OS) ~ %s"%key_to_use) ,data=R_OS)
    if no_plot == False:
        figure='OS_fc_color_many_probes.pdf'
        grdevices.pdf(figure)
    colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
    #colors=robjects.StrVector([c.upper() +'FF'  for i,c in enumerate(colors) if i <= nclusters])
    the_colors={}
    for i in range(len(all_data[key_to_use][data_to_use])):
        if all_data['stemness'][data_to_use][i] not in the_colors.keys():
#           print data['stemness'][data_to_use][order1][i], colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
#           print the_colors.keys()
            the_colors[all_data['stemness'][data_to_use][i]] = colors[ int(all_data['stemness'][data_to_use][i]) % len(colors)]
            print str(all_data['stemness'][data_to_use][i]) , str(colors[ int(all_data['stemness'][data_to_use][i]) % len(colors)])
   #print the_colors

    R_colors = colors =robjects.StrVector([c.upper() +'FF'  for c in  the_colors.values() ])
    #print R_colors
    #colors=robjects.r("rainbow(%d)"%nclusters)
#   robjects.r.plot(fit, col=colors )

#   return(fit,R_OS)
    #print no_plot
    pval = ggkm(fit,R_colors,no_plot=no_plot)
    pval = numpy.array(pval)[0]
    if no_plot == False:
        grdevices.dev_off()

    return pval


#   cl_names=robjects.StrVector(numpy.array(cat))
#   robjects.r.legend("topright",cl_names, col=colors , pch=robjects.Vector(numpy.array([3]*nclusters ) ) )

#####   pass

#


def bigcorPar(x, nblocks = 10, verbose = True, ncore=30):
    """
    runs a correlation analysis ... multicore.
    """

    robjects.r('''
    bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore=30, ...){
  library(ff, quietly = TRUE)
  require(doMC)
    registerDoMC(cores = ncore)

    NCOL <- ncol(x)

    ## test if ncol(x) %% nblocks gives remainder 0
    if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    #corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
    corMAT <- matrix(data=0,nrow=NCOL,ncol=NCOL)


    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

    ## create all unique combinations of blocks
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)

    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    results <- foreach(i = 1:nrow(COMBS)) %dopar% {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        COR <- cor(MAT[, G1], MAT[, G2], ...)
        corMAT[G1, G2] <- COR
        corMAT[G2, G1] <- t(COR)
        COR <- NULL
    }

    gc()
    return(corMAT)
}
    ''')
    r_bigcorpar =robjects.globalenv['bigcorPar']

    corr_mat = r_bigcorpar(x, nblocks = nblocks, verbose = verbose, ncore=ncore)
    return numpy.array(corr_mat)
    pass
#

    pass

def ggkm(sfit,R_colors, no_plot=False):
    """docstring for ggkm
    """

    robjects.r('''
    # Create a Kaplan-Meier plot using ggplot2
    #http://mcfromnz.wordpress.com/2012/05/05/kaplan-meier-survival-plot-with-at-risk-table-by-sub-groups/
    ggkm <- function(sfit,
                     R_colors,
                     no_plot,
                     table = TRUE,
                     returns = TRUE,
                     xlabs = "Time",
                     ylabs = "Survival Probability",
                     xlims = c(0,max(sfit$time)),
                     ylims = c(0,1),
                     ystratalabs = NULL,
                     ystrataname = NULL,
                     timeby = 100,
                     main = "Kaplan-Meier Plot",
                     pval = TRUE,
                     subs = NULL,
                     ...) {

        #############
        # libraries #
        #############

        require(ggplot2)
        require(survival)
        require(gridExtra)
        require(reshape)
        require(plyr)
        #################################
        # sorting the use of subsetting #
        #################################

        if (no_plot == FALSE){

        #print("hdhd")
        times <- seq(0, max(sfit$time), by = timeby)

        if(is.null(subs)){
            subs1 <- 1:length(levels(summary(sfit)$strata))
            subs2 <- 1:length(summary(sfit,censored=T)$strata)
            subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$strata)
        } else{
            for(i in 1:length(subs)){
                if(i==1){
                    ssvar <- paste("(?=.*\\b=",subs[i],sep="")
                }
                if(i==length(subs)){
                    ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
                }
                if(!i %in% c(1, length(subs))){
                    ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
                }
                if(i==1 & i==length(subs)){
                    ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
                }
            }
            subs1 <- which(regexpr(ssvar,levels(summary(sfit)$strata), perl=T)!=-1)
            subs2 <- which(regexpr(ssvar,summary(sfit,censored=T)$strata, perl=T)!=-1)
            subs3 <- which(regexpr(ssvar,summary(sfit,times = times,extend = TRUE)$strata, perl=T)!=-1)
        }

        ##################################
        # data manipulation pre-plotting #
        ##################################

        if(is.null(ystratalabs)) ystratalabs <- as.character(levels(summary(sfit)$strata)[subs1])
        if(is.null(ystrataname)) ystrataname <- "strata"
        m <- max(nchar(ystratalabs))
        times <- seq(0, max(sfit$time), by = timeby)

        .df <- data.frame(                      # data to be used in the survival plot
            time = sfit$time[subs2],
            n.risk = sfit$n.risk[subs2],
            n.event = sfit$n.event[subs2],
            surv = sfit$surv[subs2],
            strata = factor(summary(sfit, censored = T)$strata[subs2]),
            upper = sfit$upper[subs2],
            lower = sfit$lower[subs2]
    )
		#print("Cnahge me!!!! ystratalabs")
		#print("Cnahge me!!!! ystratalabs")
		#ystratalabs <- c("above median","below median")
		#print("Cnahge me!!!! ystratalabs")
		#print("Cnahge me!!!! ystratalabs")
        levels(.df$strata) <- ystratalabs       # final changes to data for survival plot
        zeros <- data.frame(time = 0, surv = 1,
                            strata = factor(ystratalabs, levels=levels(.df$strata)),
                            upper = 1, lower = 1)
        .df <- rbind.fill(zeros, .df)
        d <- length(levels(.df$strata))
        print(ystratalabs)
        ###################################
        # specifying plot parameteres etc #
        ###################################
        eee <- as.array(R_colors)
        p <- ggplot( .df,  aes(time, surv)) +
            geom_step(aes(linetype = strata, colour  = strata )  , size = 0.7  ) +
            theme_bw() +
            theme(axis.title.x = element_text(vjust = 0.5)) +
            scale_x_continuous(xlabs, breaks = times, limits = xlims) +
            scale_y_continuous(ylabs, limits = ylims) +
            theme(panel.grid.minor = element_blank()) +
            theme(legend.position = c(ifelse(m < 10, .28, .35),ifelse(d < 4, .25, .35))) +    # MOVE LEGEND HERE [first is x dim, second is y dim]
            theme(legend.key = element_rect(colour = NA)) +
            labs(linetype = ystrataname) +
            theme(plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 2.5)),"lines"))  +
            labs(title = main)


        ## Create a blank plot for place-holding
        ## .df <- data.frame()
        blank.pic <- ggplot(.df, aes(time, surv)) +
            geom_blank() + theme_bw() +
            theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
                 axis.title.x = element_blank(),axis.title.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid.major = element_blank(),panel.border = element_blank())

        #####################
        # p-value placement #
        #####################a

        if(pval) {
            sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
            pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
            pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
            p <- p + annotate("text",x = 0.6 * max(sfit$time),y = 0.1,label = pvaltxt)
            print(pvaltxt)
        }

        ###################################################
        # Create table graphic to include at-risk numbers #
        ###################################################

        if(table) {
            risk.data <- data.frame(
                strata = factor(summary(sfit,times = times,extend = TRUE)$strata[subs3]),
                time = summary(sfit,times = times,extend = TRUE)$time[subs3],
                n.risk = summary(sfit,times = times,extend = TRUE)$n.risk[subs3]
                )
            risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))

            data.table <- ggplot(risk.data,aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
                #, color = strata)) +
                geom_text(size = 3.5) + theme_bw() +
                scale_y_discrete(breaks = as.character(levels(risk.data$strata)),
                                 labels = rev(ystratalabs)) +
                                     # scale_y_discrete(#format1ter = abbreviate,
                                     # breaks = 1:3,
                                     # labels = ystratalabs) +
                                     scale_x_continuous("Numbers at risk", limits = xlims) +
                                     theme(axis.title.x = element_text(size = 10, vjust = 1),
                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.border = element_blank(),axis.text.x = element_blank(),
                                          axis.ticks = element_blank(),axis.text.y = element_text(face = "bold",hjust = 1))

                                     #opts(axis.title.x = theme_text(size = 10, vjust = 1),
                                     #    panel.grid.major = theme_blank(), panel.grid.minor = theme_blank(),
                                     #    panel.border = theme_blank(),axis.text.x = theme_blank(),
                                     #    axis.ticks = element_blank(),axis.text.y = theme_text(face = "bold",hjust = 1))

            data.table <- data.table +
                theme(legend.position = "none") + xlab(NULL) + ylab(NULL)

            data.table <- data.table +
                theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.28 * m), "lines")) # ADJUST POSITION OF TABLE FOR AT RISK

            #######################
            # Plotting the graphs #
            #######################

            ## p <- ggplotGrob(p)
            ## p <- addGrob(p, textGrob(x = unit(.8, "npc"), y = unit(.25, "npc"), label = pvaltxt,
            ## gp = gpar(fontsize = 12)))
            grid.arrange(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                         ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))

            if(returns) {
                a <- arrangeGrob(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                                 ncol = 1, heights = unit(c(2, .1, .25), c("null", "null", "null")))
                print(pval)
                return(pval)
            }
        } else {
            ## p <- ggplotGrob(p)
            ## p <- addGrob(p, textGrob(x = unit(0.5, "npc"), y = unit(0.23, "npc"),
            ## label = pvaltxt, gp = gpar(fontsize = 12)))

        }
    }else{ # no plot
            sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
            pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
            pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
        #    p <- p + annotate("text",x = 0.6 * max(sfit$time),y = 0.1,label = pvaltxt)
            print(pvaltxt)
            return(pval)

    }


    }   ''')
    r_ggkm =robjects.globalenv['ggkm']

    pval = r_ggkm(sfit,R_colors,no_plot)
#   print pval
#   print type(sss)
    return pval
    pass
#

def heatmap_map(data, train ,probes_to_use=None, data_to_use=None, mode = None, order1=None,order2=None,center_g=False,center_g2 = False, Corder = False, rescale_data = False, Text_labels= False, prefix  = '', ncolors=1000):
    """docstring for heatmap
        makes a heatmap of input data. cluster genes and samples.
        data is the whole data
        train is the indexes of the mapped data.
        train is the data used to train the data. (so genes are ordered with the train an)
    """
    import pylab
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    import matplotlib.ticker as ticker

    colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
    pylab.close()
    nullfmt   = plt.NullFormatter()         # no labels
    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
        print 'added data_to_use'
    if probes_to_use == None:
        probes_to_use=range(data['data'].shape[0])
        print 'added probes_to_use'

    cat=[]
    for i in data['stemness']:
        if i not in cat:
            cat.append(i)
    names=['-'*(cat.index(i)+1) for i in data['stemness']]
    values=data['data'][probes_to_use,:][:,data_to_use].T

    #train_samples have to be sorted..
    cat_train =[]
    for i in data['colnames'][train]:
        if i not in cat_train:
            cat_train.append(i)

    train_sorted = numpy.concatenate([numpy.argwhere(data['colnames'][train] == str(i)).reshape(-1) for i in cat_train])
    train = train[train_sorted]
    print "train"
    print train

    if center_g or center_g2:
        data_pca = (values - numpy.mean(values, 0))/numpy.std(values, axis=0, ddof=1)
    #and put betweem -1 and 1 genewise
        data_pca = (data_pca-numpy.array([data_pca.min(axis=1)]*data_pca.shape[1]).T)
        data_pca = numpy.array([data_pca[:,i] / data_pca[:,i].max() for i in range(data_pca.shape[1])]).T
        data_pca = data_pca*2-1
        if center_g2:
            #and put betweem -1 and 1 sample wise
            data_pca = data_pca.T
            data_pca = (data_pca-numpy.array([data_pca.min(axis=1)]*data_pca.shape[1]).T)
            data_pca = numpy.array([data_pca[:,i] / data_pca[:,i].max() for i in range(data_pca.shape[1])]).T
            data_pca = data_pca*2-1
            data_pca = data_pca.T
    else:
        data_pca = values

    data_pcaT = values.T

    if 1: #uses orange optimal leaf ordering clustering
        import orange, orngClustering
        print 'clustering samples...'
        if order1==None:
            root = orngClustering.hierarchicalClustering(orange.ExampleTable(data_pca),\
                    distanceConstructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    linkage=orange.HierarchicalClustering.Ward,\
                    order=Corder)
            order1 = numpy.array(root.mapping)
            orngClustering.dendrogram_draw(prefix+"samples-dendrogram.eps", root )
            print 'number of samples', order1.shape
        print 'clustering genes...'#[:,order1]),\
        if order2==None:
            root = orngClustering.hierarchicalClustering(orange.ExampleTable(data_pcaT),\
                    distanceConstructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    linkage=orange.HierarchicalClustering.Ward,\
                    order=Corder)
            order2 = numpy.array(root.mapping)
            print 'number of genes', order2.shape
            print 'done'
            orngClustering.dendrogram_draw(prefix+"genes-dendrogram.eps", root )
#   else:
#       a=hcluster.pdist(data_pca,'correlation')
#       b=hcluster.linkage(a,'complete')
##      b=hcluster.linkage(data_pca,'ward')
#       c=hcluster.dendrogram(b, orientation='right' ,get_leaves = True, color_threshold=0., labels = names , no_plot=False)
#       order1=c['leaves']
#       pylab.savefig('sample_tree.pdf')
#       pylab.close()
#       a=hcluster.pdist(data_pcaT[:,order1],'correlation')
##      b=hcluster.linkage(a,'weighted')
#       b=hcluster.linkage(data_pcaT[:,order1],'ward')
#       c=hcluster.dendrogram(b,orientation='bottom',no_labels=True, no_plot=True)
#
#       order2= c['leaves']
#       order1.reverse()


    fig = plt.figure(facecolor='w')
    ax3 = SubplotZero(fig, 121)
    fig.add_subplot(ax3)

    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.0, 1.0),
                      (0.5, 1.0, 1.0),
                      (1.0, 0.0, 0.0))}
    my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,ncolors)
    print order1
    print order2
    X=data_pca[:,order2][order1,:]
    mapped_X = data['data'][probes_to_use,:][:,train].T
    mapped_X = mapped_X[:,order2]
    print X.shape
    print mapped_X.shape
    X = numpy.concatenate([X,mapped_X],axis=0)

    if 0:
        print data['colnames'][data_to_use][order1]
        print data['stemness'][data_to_use][order1]

    if rescale_data:
        X = (X - numpy.mean(X, 0))/numpy.std(X, axis=0, ddof=1)
    #and put betweem -1 and 1 genewise
        X = (X-numpy.array([X.min(axis=1)]*X.shape[1]).T)
        X = numpy.array([X[:,i] / X[:,i].max() for i in range(X.shape[1])]).T
        X = X*2-1

    ccc = ax3.imshow( X, interpolation='nearest', cmap= my_cmap, aspect='auto',  )
#   ccc = ax3.imshow(X , interpolation='nearest', cmap= 'hot', aspect='auto',  )
#   ccc = ax3.imshow( X, interpolation='gaussian', cmap= my_cmap, aspect='auto',  )
    f = open(prefix+'heatmap.txt','w')
    f.write('\t')
    for i in order1:
        f.write('%s\t'% (data['colnames'][data_to_use][i]))
    f.write('\n')
    for i in order2:
        f.write('%s\t'% (data['rownames'][probes_to_use][i]) )
        for j in order1:
            f.write('%f\t'% values[j,i])
        f.write('\n')
    f.close()


    kwargs= {'cax':ax3,'ax':ax3}
#
    ax4 = SubplotZero(fig, 122, sharey=ax3)
    fig.add_subplot(ax4)

#   print len (order1)
    used_labels =[]
    for i,order_i in enumerate(order1):
        print i
        label = numpy.array(data['colnames'][data_to_use])[order_i]
        if label not in used_labels:
            used_labels.append(label)
            ax4.plot([-1,1],[i,i], '-',linewidth=1.,color=colors[ int(data['stemness'][data_to_use][order_i]) % len(colors)],label=label )
        else:
            ax4.plot([-1,1],[i,i], '-',linewidth=1.,color=colors[ int(data['stemness'][data_to_use][order_i]) % len(colors)] )
#
    for j,order_j in enumerate(train):
        print j
        label = numpy.array(data['colnames'][train])[j]
        if label not in used_labels:
            used_labels.append(label)
            ax4.plot([-1,1],[i+j,i+j], '-',linewidth=1.,color=colors[ int(data['stemness'][train][j]) % len(colors)],label=label )
        else:
            ax4.plot([-1,1],[i+j,i+j], '-',linewidth=1.,color=colors[ int(data['stemness'][train][j]) % len(colors)] )


    if Text_labels:
#       print numpy.array(data['colnames'][data_to_use])[order1]
        n,p = merge_labels( numpy.array(data['colnames'][data_to_use])[order1] , range(len([order1])) )
        for i,j in enumerate(n):
            ax4.text(1.5, p[i],  j , size=5)
        maxx = max(p)
        n,p = merge_labels( numpy.array(data['colnames'][train]) , range(len([train])) )
        for i,j in enumerate(n):
            ax4.text(1.5, p[i]+maxx,  j , size=5)


    kwarg={'size':5 }
    pylab.legend(markerscale=.5, prop=kwarg)

#marker='s', linewidths=0,
#       print 'There is something wrong here!!'
    if mode != None: #this is used to report the AML NK genes mutations.
        for index,key in enumerate(['cluster','cebpa' , 'npm1', 'flt3']):
            ax4.text(-1+index*2,-10,key,rotation ='vertical', size=6)
        print order1
        for i in range(len(data['stemness'][data_to_use])):
            for index,key in enumerate(['cebpa' , 'npm1', 'flt3']):
                #print 'There is something wrong here!!'
#               print i
                if data[key][data_to_use][order1][i] == 1:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='k')
                elif data[key][data_to_use][order1][i] == -1:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='y')
                else:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='w')
    else:
        0
#       ax3.set_yticks(range(len(order1)))
#       kwargs = {'size':'xx-small'}
#       ax3.set_yticklabels(range(len(data['colnames'][data_to_use])),numpy.array(data['colnames'][data_to_use])[order1] , **kwargs)



    ax4.set_xlim( -1, 75 ) #hack!!
    for ax in [ax3,ax4]:
        for direction in ["left", "right", "bottom", "top"]:
            ax.axis[direction].set_visible(False)



    pylab.colorbar(ccc, ax=ax4 ,fraction= .05, orientation ='vertical', )

#   ax3.set_xticks(ccc,range(len(order2)))
#   ax3.set_xticklabels(ccc,range(len(probes_to_use)),numpy.array(data['rownames'][probes_to_use])[order2] )


    fig.savefig(prefix+'heatmap_genes.pdf',dpi=600)
    fig.show()
    pylab.close()

    the_colors={}
    cat=[]
    for i in range(len(data['stemness'][data_to_use])):
        if data['stemness'][data_to_use][order1][i] not in the_colors.keys():
#           print data['stemness'][data_to_use][order1][i], colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
#           print the_colors.keys()
            the_colors[data['stemness'][data_to_use][order1][i]] = colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
    print the_colors
    return order1, order2

    pass

def normal_vs_cancer(data, data_to_use = None,probes_to_use=None,p = .005, fc = 1, limit=None , dir = None ):
    qval=importr('qvalue')

    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
        print 'added data_to_use'
    if probes_to_use == None:
        probes_to_use=range(data['data'].shape[0])
        print 'added probes_to_use'

    if not data.has_key('chimera_sample'):
        print 'chimera_sample'
    dat = {}
    dat['data'] = numpy.concatenate([data['data'][probes_to_use,:][:,data_to_use] ,data['chimera_sample'].T[probes_to_use,:][:,data_to_use] ],axis=1)
    dat['colnames'] = numpy.concatenate([ ['cancer']*len(data_to_use) , ['nl']*len(data_to_use)] )
    dat['rownames'] = data['rownames'][probes_to_use]
    dat['stemness']  = stemness_by_class(dat)
    e=numpy.array(limma_multiclass(dat, p = p, fc=fc, return_raw_r =1))
    id_indexed = numpy.array(e[0,:])

    if dir is None:
        fc_indexes = numpy.argwhere(numpy.abs(numpy.array(e[1,:],dtype = 'f')) > fc).reshape(-1)
    elif numpy.sign(dir) == 1:
        fc_indexes = numpy.argwhere(numpy.array(e[1,:],dtype = 'f') > fc).reshape(-1)
    else:
        fc_indexes = numpy.argwhere(numpy.array(e[1,:],dtype = 'f') < fc).reshape(-1)

    p_indexes = numpy.argwhere(numpy.abs(numpy.array(e[4,:],dtype = 'f')) < p).reshape(-1)

    qvalues=qval.qvalue(numpy.array(e[4,:],dtype = 'f'))
    qvals= numpy.array(qvalues.rx2(3))

    print 'found %d genes that have fc above %f'%(len(fc_indexes),fc)
    print 'found %d genes that have p bellow %f'%(len(p_indexes),p)
    overlap  = set( fc_indexes ).intersection(set(p_indexes))
    if limit is None:
        return  get_probe_index(id_indexed[list(overlap)],data['rownames'])
    else:
    #   print fc_indexes[:limit]
        overlap  = fc_indexes[:limit]
        return  get_probe_index(id_indexed[overlap],data['rownames'])

    pass

def transpose_data_dic(data):
    data['data'] = data['data'].T
    data['colnames'], data['rownames'] = data['rownames'], data['colnames']
    return data
    pass
def trim_data_dic(data, data_to_use=None,probes_to_use=None):
    import copy
    _data = copy.copy(data)
    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
        print 'added data_to_use'
    if probes_to_use == None:
        probes_to_use=range(data['data'].shape[0])
        print 'added probes_to_use'

    _data['data']  = _data['data'][probes_to_use,:][:,data_to_use]
    _data['colnames'] =  _data['colnames'][data_to_use]
    _data['rownames'] = _data['rownames'][probes_to_use]
    try:
        _data['stemness'] =  _data['stemness'][data_to_use]
    except LookupError,e:
        print e
    return _data
    pass


def idoify_heatmap(dat,probes=None):
    '''
    make heatmap the way ido amit does..
    '''
    import scipy
    from copy import deepcopy
    from numpy import *
    import matplotlib.colors as mcolors

    if probes == None:
        probes=range(data['data'].shape[0])
        print 'added probes_to_use'

    colconv = mcolors.ColorConverter().to_rgb

    p2u =  probes

    genes = dat['rownames'][p2u]


    my_cmap =  make_colormap( [colconv('white'), colconv('white'), .2, \
                                  colconv('white'), colconv('#edecf2'),.5,\
                                  colconv('#edecf2'), colconv('#9792b3'),.9,\
                                  colconv('#9792b3'), colconv('#5F1A4D'),.99,\
                                  colconv('#5F1A4D')],ncolors=60)

    dat['data']  = array(dat['data'] /dat['data'].sum(axis=0))
    #Then, each gene was scaled by dividing expression by the 98 percentile for that gene.
    percentile98 =  array([nanpercentile(dat['data'][i],98,interpolation='lower') for i in range(len(dat['data'][:,0]))])  # check here     the nan percentile
    zereos =[]
    for index,v in enumerate(percentile98):
        dat['data'][index]  = dat['data'][index] /v
    dat['data'] = nan_to_num(dat['data'])
    #for index,v in enumerate(percentile98):
    #   if v >0:
    #       dat['data'][index]  = dat['data'][index] /v
    #   else:
    #       print index
    #       zereos.append(index)
    #       #a = dat['data'][index,:]
    #       #v= min(a[a!=0])
    #       dat['data'][index]  = dat['data'][index] / nanmax(dat['data'][index,:])# percentile(dat['data'][index,:],80.9)
    for i in genes:
        if i in dat['rownames'][zereos]:
            print i + ' has only zeros'
    #The resulting score was converted to the log scale,
    a = dat['data'].ravel()
    print a.shape
    minimum = numpy.min(a[a!=0])
    maximum = numpy.max(a[a!=0])
    dat['data'][dat['data']==0] = minimum / 100
    dat['data'] = log(dat['data'])
    #a rolling mean with window size of 5 cells was calculated, to account for the visualization of such a large number of cells.
    #dat['data'] = array([smoothTriangle(dat['data'][i],5,mode=False) for i in  range(len(dat['data'][:,0]))  ])

    copydat = deepcopy(dat)
    copydat['rownames'] = dat['rownames'][p2u]
    copydat['data'] = dat['data'][p2u,:]
    copydat['data'] = array([smoothTriangle(copydat['data'][i],7,mode=False) for i in  range(len(copydat['data'][:,0]))  ])
    a = copydat['data'].ravel()
    minimum = numpy.min(a[a!= minimum])
    maximum = sort(nan_to_num(unique(a[a!=0])))[-2] # i m not kidding.

    a,b = heatmap(copydat,prefix='Ido_heatmap_',gene_labels= 1, order1= range(len(dat['colnames']))[::-1],order2=range(len(p2u)), colorscale=my_cmap,Corder=1,Text_labels=False,center_g2=1)
    pass


def heatmap(data, probes_to_use=None, data_to_use=None, mode = None, order1=None,order2=None,\
            center_g=False,center_g2 = False, Corder = False, rescale_data = False, \
            Text_labels= False, prefix  = '', ncolors=100,rescale_gene_wise = False, rescale_sample_wise =False , \
            colorscale=None,gene_labels=False,cmax=None,cmin=None, no_plot=False, \
            which_text_key=None, pdf_output=True,interpolation='nearest', gene_text_color = None):
    """docstring for heatmap
        makes a heatmap of input data. cluster genes and samples.
        NB: gene_text_color should be a list of colors with the same length of the gene lenght.
    """
    print 'Try biclust R pacakge.http://cran.r-project.org/web/packages/biclust/index.html'
    import pylab
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    import matplotlib.ticker as ticker

    colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
    pylab.close()
    if gene_labels == False:
        nullfmt   = plt.NullFormatter()         # no labels

    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
        print 'added data_to_use'
    if probes_to_use == None:
        probes_to_use=range(data['data'].shape[0])
        print 'added probes_to_use'

    cat=[]
    for i in data['stemness']:
        if i not in cat:
            cat.append(i)
    names=['-'*(cat.index(i)+1) for i in data['stemness']]
    values=data['data'][probes_to_use,:][:,data_to_use].T

    if center_g or center_g2 :
        data_pca = (values - numpy.mean(values, 0))/numpy.std(values, axis=0, ddof=1)
    #and put betweem -1 and 1 genewise
        data_pca = (data_pca-numpy.array([data_pca.min(axis=1)]*data_pca.shape[1]).T)
        data_pca = numpy.array([data_pca[:,i] / data_pca[:,i].max() for i in range(data_pca.shape[1])]).T
        data_pca = data_pca*2-1
        if center_g2:
            #and put betweem -1 and 1 sample wise
            data_pca = numpy.array(values-numpy.array([values.min(axis=0)]*values.shape[0]))
            data_pca = numpy.array([data_pca[:,i] / data_pca[:,i].max() for i in range(data_pca.shape[1])])
            data_pca = data_pca*2-1
            data_pca = data_pca.T
    else:
        data_pca = values
        if rescale_gene_wise:
            print 'Rescaling gene wise'
            data_pca = values.T
            arr = []
            for index,i in enumerate(data_pca):
                data_min =  i.min()
                x= (i-data_min)
                x /=  max(x)
                x = x*2 -1
                arr.append(x)
            data_pca = numpy.array(arr).T
        if rescale_sample_wise:
            print 'Rescaling sample wise'
            data_pca = values
            arr = []
            for index,i in enumerate(data_pca.T):
                data_min =  i.min()
                x= (i-data_min)
                x /=  max(x)
                x = x*2 -1
                arr.append(x)
            data_pca = numpy.array(arr).T


    data_pcaT = values.T

    if 1: #uses orange optimal leaf ordering clustering
        import orange, orngClustering
        print 'clustering samples...'
        if order1==None:
            root = orngClustering.hierarchicalClustering(orange.ExampleTable(data_pca),\
                    distance_constructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    #linkage=orange.HierarchicalClustering.Complete,\
                    linkage=orange.HierarchicalClustering.Ward,\
                    order=Corder)
            order1 = numpy.array(root.mapping)
            orngClustering.dendrogram_draw(prefix+"samples-dendrogram.eps", root )
            print 'number of samples', order1.shape
        print 'clustering genes...'#[:,order1]),\
        if order2==None:
            root = orngClustering.hierarchicalClustering(orange.ExampleTable(data_pcaT),\
                    distance_constructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    linkage=orange.HierarchicalClustering.Ward,\
                    order=Corder)
            order2 = numpy.array(root.mapping)
            print 'number of genes', order2.shape
            print 'done'
            orngClustering.dendrogram_draw(prefix+"genes-dendrogram.eps", root )
#   else:
#       a=hcluster.pdist(data_pca,'correlation')
#       b=hcluster.linkage(a,'complete')
##      b=hcluster.linkage(data_pca,'ward')
#       c=hcluster.dendrogram(b, orientation='right' ,get_leaves = True, color_threshold=0., labels = names , no_plot=False)
#       order1=c['leaves']
#       pylab.savefig('sample_tree.pdf')
#       pylab.close()
#       a=hcluster.pdist(data_pcaT[:,order1],'correlation')
##      b=hcluster.linkage(a,'weighted')
#       b=hcluster.linkage(data_pcaT[:,order1],'ward')
#       c=hcluster.dendrogram(b,orientation='bottom',no_labels=True, no_plot=True)
#
#       order2= c['leaves']
#       order1.reverse()

    if no_plot:
        return order1, order2

    fig = plt.figure(facecolor='w')
    ax3 = plt.Subplot(fig, 121)
    fig.add_subplot(ax3)
    print('-')
    cdict = {'red':   ((0.0, 0.4, 0.0),
                    (0.2, 1.0, 1.0),
                   (0.8, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 1.0, 1.0)
                  ,(0.2, 1.0, 1.0),
                   (0.8, 0.0, 0.0),
                   (1.0, 0.4, 0.0)),

         'blue':  ((0.0, 0.7, 0.0),
                   (0.2, 1.0, 1.0),
                   (0.8, 0.0, 0.0),
                   (1.0, 0.7, 0.0))
        }


### cdict = {'red':   ((0.0, 1.0, 1.0),
###                (0.8, 0.0, 0.1),
###                (1.0, 1.0, 1.0)),
###
###      'green': ((0.0, 1.0, 1.0),
###                (0.8, 0.0, 0.0),
###                (1.0, 0.4, 0.0)),
###
###      'blue':  ((0.0, 1.0, 1.0),
###                (0.8, 0.0, 0.0),
###                (1.0, 0.7, 0.0))
###     }
### cdict = {'red':   ((0.0, 1.0, 1.0),
###                   (0.0, 1.0, 1.0)),
###
###         'green': ((0.0, 0.0, 1.0),
###                   (0.0, 0.0, 1.0)),
###
###         'blue':  ((0.0, 0.0, 1.0),
###                   (1.0, 1.0, 1.0))
###        }

    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.0, 1.0),
                      (0.5, 1.0, 1.0),
                      (1.0, 0.0, 0.0))}

#   cdict = {'red': ((0.0, 0.0, 0.0),
#                    (0.35, 1.0, 1.0),
#                    (0.65, 1.0, 1.0),
#                    (1.0, 1.0, 1.0)),
#            'green': ((0.0, 0.0, 0.0),
#                      (0.35, 1.0, 1.0),
#                      (0.65, 1.0, 1.0),
#                      (1.0, 0.0, 0.0)),
#            'blue': ((0.0, 0.0, 1.0),
#                     (0.35, 1.0, 1.0),
#                     (0.65, 1.0, 1.0),
#                     (1.0, 0.0, 0.0))}

    my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,ncolors)
#   print order1
#   print order2
    X=data_pca[:,order2][order1,:]
    if 0:
        print data['colnames'][data_to_use][order1]
        print data['stemness'][data_to_use][order1]

    if rescale_data:
        X = (X - numpy.mean(X, 0))/numpy.std(X, axis=0, ddof=1)
    #and put betweem -1 and 1 genewise
        X = (X-numpy.array([X.min(axis=1)]*X.shape[1]).T)
        X = numpy.array([X[:,i] / X[:,i].max() for i in range(X.shape[1])]).T
        X = X*2-1

    if gene_labels == False:
        plt.axis('off')
        print('--')
        ax3.set_axis_off()
        ax = plt.gca()

#    my_cmap =  pylab.get_cmap('pink')
    if cmin==None:
        if colorscale == None:
            ccc = ax3.imshow( X, interpolation=interpolation, cmap=  my_cmap, aspect='auto',  ) #my_cmap, aspect='auto',  )
        else:
            ccc = ax3.imshow(X , interpolation=interpolation, cmap= colorscale, aspect='auto',  )
    else:
        if colorscale == None:
            ccc = ax3.imshow( X, interpolation=interpolation, cmap=  my_cmap, aspect='auto',vmin=cmin, vmax=cmax ) #my_cmap, aspect='auto',  )
        else:
            ccc = ax3.imshow(X , interpolation=interpolation, cmap= colorscale, aspect='auto', vmin=cmin, vmax=cmax )

#   ccc = ax3.imshow( X, interpolation='sinc', cmap= my_cmap, aspect='auto',  )
    if gene_labels:
        print 'gene labels!'
        ddd= plt.xticks(numpy.arange(len(data['rownames'][probes_to_use][order2])),data['rownames'][probes_to_use][order2])
        ddd= fig.autofmt_xdate()
        for i,j in enumerate(data['rownames'][probes_to_use][order2]):
            print data['rownames'][probes_to_use][order2][i]
            if gene_text_color == None:
                label_color ='k'
            else:
                try:
                    gene_text_color = numpy.array(gene_text_color)
                    #find the color in stemness:
                    stemness_i  = numpy.argwhere(data['colnames'] == gene_text_color[order2][i]).ravel()[0]
                    label_color = colors[data['stemness'][data_to_use][int(stemness_i)]]
                except Exception,e:
                    unique_colors=[]
                    gene_text_color = numpy.array(gene_text_color)
                    for ii in gene_text_color:
                        if ii not in unique_colors:
                            unique_colors.append(ii)
                    unique_colors = numpy.array(unique_colors)
                    stemness_i  = numpy.argwhere(unique_colors == gene_text_color[order2][i]).ravel()[0]
                    label_color = colors[data['stemness'][data_to_use][int(stemness_i)]]
            ax3.text( i+.75, 50+len(order1), j , size=5, rotation='vertical',ha = 'right',color=label_color)
        #plt.axes.set_xticklabels(data['rownames'][probes_to_use][order2], fontdict=None, minor=False)


    f = open(prefix+'heatmap.txt','w')
    f.write('\t')
    for i in order1:
        f.write('%s\t'% (data['colnames'][data_to_use][i]))
    f.write('\n')
    for i in order2:
        f.write('%s\t'% (data['rownames'][probes_to_use][i]) )
        for j in order1:
            f.write('%f\t'% values[j,i])
        f.write('\n')
    f.close()
    #remove Axis
    if gene_labels == False:
        plt.gca().axison = False
        pylab.gca().axes.get_yaxis().set_visible(False)
        pylab.gca().axes.get_xaxis().set_visible(False)
        for direction in ["left", "right", "bottom", "top"]:
            ax3.spines[direction].set_visible(False)
            plt.gca().axison = False
            pylab.gca().axison = False

    kwargs= {'cax':ax3,'ax':ax3}
#
    ax4 = pylab.Subplot(fig, 122, sharey=ax3)
    fig.add_subplot(ax4)

#   print len (order1)
    used_labels =[]
    for i,order_i in enumerate(order1):
#       print i
        label = numpy.array(data['colnames'][data_to_use])[order_i]
        #print 'ZZZZZZZZZZZis is why the labels are strange.'
        #label = numpy.array(data['cel_file'][data_to_use])[order_i]
        if label not in used_labels:
            used_labels.append(label)
            ax4.plot([-1,1],[i,i], '-',linewidth=1.,color=colors[ int(data['stemness'][data_to_use][order_i]) % len(colors)],label=label )
        else:
            ax4.plot([-1,1],[i,i], '-',linewidth=1.,color=colors[ int(data['stemness'][data_to_use][order_i]) % len(colors)] )
    if Text_labels:
#       print numpy.array(data['colnames'][data_to_use])[order1]
        try:
            n,p = merge_labels( numpy.array(data['colnames'][data_to_use])[order1] , range(len([order1])) )
            #n,p = merge_labels( numpy.array(data['cel_file'][data_to_use])[order1] , range(len([order1])) )

            for i,j in enumerate(n):
                ax4.text(1.5, p[i],  j , size=5)
            pass
        except Exception, e:
            print e
            for i,j in enumerate(numpy.array(data['colnames'][data_to_use])[order1]):
                ax4.text(1.5, i ,  j , size=5)

        else:
            if which_text_key == None:
                for i,j in enumerate(numpy.array(data['colnames'][data_to_use])[order1]):
                       ax4.text(1.5, i ,  j , size=5)
            else:
                for i,j in enumerate(numpy.array(data[which_text_key][data_to_use])[order1]):
                       ax4.text(1.5, i ,  j , size=5)


    kwarg={'size':5 }
    pylab.legend(markerscale=.5, prop=kwarg)

#marker='s', linewidths=0,
#       print 'There is something wrong here!!'
    if mode != None: #this is used to report the AML NK genes mutations.
        for index,key in enumerate(['cluster','cebpa' , 'npm1', 'flt3']):
            ax4.text(-1+index*2,-10,key,rotation ='vertical', size=6)
#       print order1
        for i in range(len(data['stemness'][data_to_use])):
            for index,key in enumerate(['cebpa' , 'npm1', 'flt3']):
                #print 'There is something wrong here!!'
#               print i
                if data[key][data_to_use][order1][i] == 1:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='k')
                elif data[key][data_to_use][order1][i] == -1:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='y')
                else:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='w')
    else:
        0
#       ax3.set_yticks(range(len(order1)))
#       kwargs = {'size':'xx-small'}
#       ax3.set_yticklabels(range(len(data['colnames'][data_to_use])),numpy.array(data['colnames'][data_to_use])[order1] , **kwargs)



    ax4.set_xlim( -1, 75 ) #hack!!
    pylab.colorbar(ccc, ax=ax4 ,fraction= .05, orientation ='vertical', )


#   return [ax3,ax4]
    for ax in [ax3,ax4]:
        pylab.gca().axes.get_yaxis().set_visible(False)
        pylab.gca().axes.get_xaxis().set_visible(False)
        for direction in ["left", "right", "bottom", "top"]:
            ax.spines[direction].set_visible(False)
            plt.gca().axison = False
            pylab.gca().axison = False


    ax3.set_xticks([])
    ax3.set_yticks([])
#   ax3.set_xticks(ccc,range(len(order2)))
    ax3.set_xticklabels(['',''] )
    ax3.set_yticklabels(['',''] )
    ax3.get_xaxis().set_visible(False)

    print('---')

    if pdf_output:
        fig.savefig(prefix+'heatmap_genes.pdf',dpi=600)
    else:
        fig.savefig(prefix+'heatmap_genes.eps',dpi=600)
    #fig.show()
    #pylab.close()

    if 0:
        the_colors={}
        cat=[]
        for i in range(len(data['stemness'][data_to_use])):
            if data['stemness'][data_to_use][order1][i] not in the_colors.keys():
    #           print data['stemness'][data_to_use][order1][i], colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
    #           print the_colors.keys()
                the_colors[data['stemness'][data_to_use][order1][i]] = colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
#       print the_colors
    if rescale_data:
        return order1, order2 ,X
    else:
        return order1, order2
    pass

def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    #cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))

    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]

def show_colormap(cmap):
    im = np.outer(np.ones(10), np.arange(100))
    fig, ax = plt.subplots(2, figsize=(6, 1.5), subplot_kw=dict(xticks=[], yticks=[]))
    fig.subplots_adjust(hspace=0.1)
    ax[0].imshow(im, cmap=cmap)
    #ax[1].imshow(im, cmap=grayify_cmap(cmap))
    pass
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

def plot_cluster(data, probes_to_use=None, data_to_use=None, order1=None,order2=None,\
            center_g=False,center_g2 = False, Corder = False, \
            cluster_samples = True,\
            out = 'clustering.eps',\
            labels = None):
    """docstring for plot_cluster
        makes a clustering of the data. cluster genes or samples.
    """
    import pylab
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    import matplotlib.ticker as ticker
    import orange, orngClustering

    colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    pylab.close()

    if out.find('.eps') == -1:
        out = out + '.eps'


    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
        print 'added data_to_use'
    if probes_to_use == None:
        probes_to_use=range(data['data'].shape[0])
        print 'added probes_to_use'

    if labels == None:
        if cluster_samples:
            labels = data['colnames'][data_to_use]
        else:
            labels = data['colnames'][probes_to_use]

    cat=[]
    for i in data['stemness'][data_to_use]:
        if i not in cat:
            cat.append(i)
#    names=['-'*(cat.index(i)+1) for i in data['stemness']]
    values=data['data'][probes_to_use,:][:,data_to_use].T

    if center_g or center_g2 :
        data_pca = (values - numpy.mean(values, 0))/numpy.std(values, axis=0, ddof=1)
    #and put betweem -1 and 1 genewise
        data_pca = (data_pca-numpy.array([data_pca.min(axis=1)]*data_pca.shape[1]).T)
        data_pca = numpy.array([data_pca[:,i] / data_pca[:,i].max() for i in range(data_pca.shape[1])]).T
        data_pca = data_pca*2-1
        if center_g2:
            #and put betweem -1 and 1 sample wise
            data_pca = data_pca.T
            data_pca = (data_pca-numpy.array([data_pca.min(axis=1)]*data_pca.shape[1]).T)
            data_pca = numpy.array([data_pca[:,i] / data_pca[:,i].max() for i in range(data_pca.shape[1])]).T
            data_pca = data_pca*2-1
            data_pca = data_pca.T
    else:
        data_pca = values

    data_pcaT = values.T
    if cluster_samples :
        print 'clustering samples...'

        color_palette =  numpy.array([hex_to_rgb(i) for i in colors])[numpy.unique(cat)]
        print color_palette
        if order1==None:
            stem_data =  orange.ExampleTable( numpy.array([data['stemness'][data_to_use]]).T )
            root = orngClustering.hierarchicalClustering(orange.ExampleTable(data_pca),\
                    distance_constructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    #linkage=orange.HierarchicalClustering.Complete,\
                    linkage=orange.HierarchicalClustering.Ward,\
                    order=Corder)
            order1 = numpy.array(root.mapping)
            orngClustering.dendrogram_draw(out, root ,data = stem_data,  labels=list(labels), color_palette=color_palette  )#(,color_palette=color_palette ) # data=orange.ExampleTable(data['stemness'][data_to_use]),
            print 'number of samples', order1.shape
    else:
        print 'clustering genes...'#[:,order1]),\
        if order2==None:
            root = orngClustering.hierarchicalClustering(orange.ExampleTable(data_pcaT),\
                    distance_constructor=orange.ExamplesDistanceConstructor_Euclidean, \
                    linkage=orange.HierarchicalClustering.Ward,\
                    order=Corder)
            order2 = numpy.array(root.mapping)
            print 'number of genes', order2.shape
            print 'done'
            orngClustering.dendrogram_draw(out, root )

    pass

def get_sorted_index(x):
    """
    x is an array of numbers, returns the index order so that x[order] is sorted.

    input x
    output order (x[order] is sorted)
    """
    if type(x) == type([]):
        ranks=x
    elif type(x) == type(numpy.array([])):
        ranks=x.tolist()
    else:
        print 'wrong type of array as input.'
        print type(x)
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
    order=numpy.array([x.index(i) for i in sorted( x , key=lambda i:( float(i[0]),i[1] ) )])
    return order
    pass

def merge_labels(labels,positions):
    """docstring for merge_labels
        try to smartly merge lables on the side of the heatmap-

    """
    mini=0
    maxi=0
    names=[]
    new_positions=[]

    names=[labels[0]]
#   new_positions=[positions[0]]
    for index,name in enumerate(labels):
        if index!=0:
            if (labels[index -1] != labels[index]) :
                if (maxi - mini < 1 ):
                    names.pop( -1 )
                    names.append(labels[index+1])
                    mini=index
                    maxi=index
                else:
                    names.append( labels[index] )
                    new_positions.append( (mini+maxi+1)/2 )
                    mini=index
                    maxi=index
            else:
                maxi+=1
    new_positions.append( (mini+ len(labels))/2 )
    return names, new_positions
    pass


def quad_venn(dic,out='venn.pdf'):
    """docstring for quad_venn

        uses R VennDiagram lib to make the 4 way venn

        library(VennDiagram)
        venn.plot <- draw.quad.venn(
        area1 = 720,
        area2 = 86,
        area3 = 50,
        area4 = 0,
        n12 = 44,
        n13 = 27,
        n14 = 32,
        n23 = 38,
        n24 = 32,
        n34 = 20,
        n123 = 18,
        n124 = 17,
        n134 = 11,
        n234 = 13,
        n1234 = 6,
        category = c("First", "Second", "Third", "Fourth"),
        fill = c("orange", "red", "green", "blue"),
        lty = "dashed",
        cex = 2,
        cat.cex = 2,
        cat.col = c("orange", "red", "green", "blue")
        );
        tiff(filename = "Four-set Venn diagram.tiff", compression = "lzw");
        grid.draw(venn.plot);
        dev.off();
    """
    from collections import OrderedDict
    if type(dic) != type(dict()):
        print "please provide the correct data as input (a dic, with genes as values and classes as key)"
        #return
    if not (len(dic.keys()) == 4 or len(dic.keys()) == 5):
        print "four classes ... "
        return
    venn = importr('VennDiagram')
    base = importr('base')
    grdevices = importr('grDevices')
    grdevices.pdf(out)

    dic = OrderedDict(dic)

    sets=OrderedDict()
    sets['len'] = []
    sets['set'] = []
    sets['names'] = []
    for i in dic:
        sets['len'].append( len(dic[i]))
        sets['set'].append( set(dic[i]))
        sets['names'].append(i)

    print sets

    print robjects.StrVector( dic.keys() )
    area1 = sets['len'][0]
    area2 = sets['len'][1]
    area3 = sets['len'][2]
    area4 = sets['len'][3]
    n12 = len ( sets['set'][0].intersection(sets['set'][1]) )
    n13 = len ( sets['set'][0].intersection(sets['set'][2]) )
    n14 = len ( sets['set'][0].intersection(sets['set'][3]) )
    n23 = len ( sets['set'][1].intersection(sets['set'][2]) )
    n24 = len ( sets['set'][1].intersection(sets['set'][3]) )
    n34 = len ( sets['set'][2].intersection(sets['set'][3]) )
    n123 = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]) )
    n124 = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][3]) )
    n134 = len ( sets['set'][0].intersection(sets['set'][2]).intersection(sets['set'][3]) )
    n234 = len ( sets['set'][1].intersection(sets['set'][2]).intersection(sets['set'][3]) )
    n1234 =len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]).intersection(sets['set'][3]) )
    colors  = robjects.StrVector([ "yellow", "green","red", "blue" , ])

    base.cat.cex = 2
    base.cat.col = colors

    if len(dic.keys()) == 4 :
        venn.draw_quad_venn( area1 = area1, area2 = area2, area3 = area3, area4 = area4,\
        n12 = n12, n13 = n13,   n14 = n14, n23 = n23,   n24 = n24, n34 = n34, \
        n123 = n123, n124 = n124, n134 = n134, n234 = n234, n1234 = n1234,\
        category =  robjects.StrVector( dic.keys() ) ,\
        fill = colors,\
#       lty = "dashed",\
        cex = 2)
    else:
        colors  = robjects.StrVector([ "yellow", "green","red", "blue" , "white"])
        area5  = sets['len'][4]
        n12    = len ( sets['set'][0].intersection(sets['set'][1]) )
        n13    = len ( sets['set'][0].intersection(sets['set'][2]) )
        n14    = len ( sets['set'][0].intersection(sets['set'][3]) )
        n15    = len ( sets['set'][0].intersection(sets['set'][4]) )
        n25    = len ( sets['set'][1].intersection(sets['set'][4]) )
        n23    = len ( sets['set'][1].intersection(sets['set'][2]) )
        n24    = len ( sets['set'][1].intersection(sets['set'][3]) )
        n34    = len ( sets['set'][2].intersection(sets['set'][3]) )
        n35    = len ( sets['set'][2].intersection(sets['set'][4]) )
        n45    = len ( sets['set'][3].intersection(sets['set'][4]) )
        n123   = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]) )
        n124   = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][3]) )
        n125   = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][4]) )
        n134   = len ( sets['set'][0].intersection(sets['set'][2]).intersection(sets['set'][3]) )
        n135   = len ( sets['set'][0].intersection(sets['set'][2]).intersection(sets['set'][4]) )
        n145   = len ( sets['set'][0].intersection(sets['set'][3]).intersection(sets['set'][4]) )
        n235   = len ( sets['set'][1].intersection(sets['set'][2]).intersection(sets['set'][4]) )
        n234   = len ( sets['set'][1].intersection(sets['set'][2]).intersection(sets['set'][3]) )
        n245   = len ( sets['set'][1].intersection(sets['set'][3]).intersection(sets['set'][4]) )
        n345   = len ( sets['set'][2].intersection(sets['set'][3]).intersection(sets['set'][4]) )
        n1234  = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]).intersection(sets['set'][3]) )
        n1235  = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]).intersection(sets['set'][4]) )
        n2345  = len ( sets['set'][1].intersection(sets['set'][2]).intersection(sets['set'][3]).intersection(sets['set'][4]) )
        n1245  = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][3]).intersection(sets['set'][4]) )
        n1345  = len ( sets['set'][0].intersection(sets['set'][2]).intersection(sets['set'][3]).intersection(sets['set'][4]) )
        n12345 = len ( sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]).intersection(sets['set'][3]).intersection(sets['set'][4]) )
        print area1
        print area2
        print area3
        print area4
        print area5
        print n12
        print n13
        print n14
        print n15
        print n25
        print n23
        print n24
        print n34
        print n35
        print n45
        print n123
        print n124
        print n125
        print n134
        print n135
        print n145
        print n234
        print n235
        print n245
        print n345
        print n1234
        print n1235
        print n2345
        print n1245
        print n1345
        print n12345
        print sets['set'][0].intersection(sets['set'][1]).intersection(sets['set'][2]).intersection(sets['set'][3]).intersection(sets['set'][4])
#        return area1, area2 , area3 , area4,area5,\
#        n12 , n13 , n14 , n15 , n23 , n24 , n25 , n34 , n35 , n45 ,\
#        n123 , n124, n125 , n134 , n135 , n145, n234 , n235 , n245, n345 ,\
#        n1234 ,  n1235, n1245, n1345, n2345, n1234
        colors  = robjects.StrVector(["dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"])
        base.cat.cex = 2
        base.cat.col = colors
        venn.draw_quintuple_venn( area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5,\
        n12 = n12, n13 = n13, n14 = n14, n15 = n15, n23 = n23, n24 = n24, n25 = n25, n34 = n34, n35 = n35, n45 = n45,\
        n123 =n123, n124 =n124, n125 =n125, n134 =n134, n135 =n135, n145 =n145, n234 =n234, n235 =n235, n245 =n245, n345 =n345,\
        n1234 = n1234, n1235 = n1235, n1245 = n1245, n1345 = n1345, n2345 = n2345, n12345 = n1234,\
        category =  robjects.StrVector( dic.keys() ) ,\
        fill = colors,\
        #lty = "dashed",\
        cex = 2)

    grdevices.dev_off()

pass


def Venn_diag(label1, label2, label1_txt='A',label2_txt='B',overlap_txt='AnB',out='venn.pdf'):
    """docstring for Venn_diag
    """
    numpy.seterr(all='ignore')

    import pylab
    pylab.close()

    label1=set(label1)
    label2=set(label2)
    n1= len(label1)
    n2= len(label2)

    A= len(label1.intersection(label2))*1./ (n1+n2-len(label1.intersection(label2))*1.)*(n1+n2)#./ numpy.min([n1,n2])
#   A=A*A*3.14

    print n1, n2 , A

    r1=numpy.sqrt((n1)/3.14)
    R2= numpy.sqrt((n2)/3.14)
    r=numpy.sqrt(r1**2/(r1**2+R2**2))
    R=numpy.sqrt(R2**2/(r1**2+R2**2))
#   r=r1/(r1+R2)
#   R=R2/(r1+R2)
    A = numpy.sqrt(A / (r1**2+R2**2) )
#   r=r1
#   R=R2
    print r,R,A
    dd=[]
    for d in numpy.arange(0,R+r,.0001):
        diff  = abs(A - (r*r*numpy.arccos((d*d+r*r-R*R)/(2*d*r)) + R*R*numpy.arccos((d*d-r*r+R*R)/(2*d*R)) \
                        - 1/2*numpy.sqrt( (-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R) ) ) )
        dd.append(diff)
        if diff < .001:
            print 'sol =' + str(d)
#           break

    dd=numpy.array(dd)
    x = numpy.arange(0,R+r,.0001)
    d = x[numpy.argwhere(dd > 0).reshape(-1)][ numpy.argmin(dd[numpy.argwhere(dd > 0).reshape(-1)])]

    print d

    pylab.axes()
    cir = pylab.Circle((0,0), radius=r, alpha =.5, fc='#ed2921' ) # 'y')
    pylab.gca().add_patch(cir)
    if d - R < 0:
        pylab.text((-r-d+R),0,label1_txt)
        pylab.text( ((-r-d+R) + d)/2.  ,0,overlap_txt)
    else:
        pylab.text(0,0,label1_txt)
        pylab.text(d/2.,0,overlap_txt)

    cir = pylab.Circle((d,0), radius=R, alpha =.5, fc='#5f14fa' ) #'b')
    pylab.gca().add_patch(cir)
    pylab.text(d,0,label2_txt)
    pylab.axis('scaled')
    pylab.axis('off')
    pylab.savefig(out)

    pass

def heatmap_OS(data, probes_to_use=None, data_to_use=None, mode = None,labelx=None,labely=None):
    """docstring for heatmap
        makes a heatmap of input data. cluster genes and samples.
    """
    import pylab
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    import matplotlib.ticker as ticker

    colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']


    pylab.close()

#   matplotlib.use('Agg')


    nullfmt   = plt.NullFormatter()         # no labels

    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
        print 'added data_to_use'
    if probes_to_use == None:
        probes_to_use=range(data['data'].shape[0])
        print 'added probes_to_use'

    x=numpy.array([data['OS'][data_to_use],range(len(data_to_use))]).transpose().tolist()
    order_OS=[x.index(i) for i in sorted( x , key=lambda i:( int(i[0]),i[1] ) )]


    fig = plt.figure(facecolor='w')

#   ax1 = SubplotZero(fig, 221)
#   fig.add_subplot(ax1)

    cat=[]
    for i in data['stemness']:
        if i not in cat:
            cat.append(i)
    names=['-'*(cat.index(i)+1) for i in data['stemness']]
#   print names

    data_pca = data['data'][probes_to_use,:][:,data_to_use][:,order_OS].T
    data_pcaT = data['data'][probes_to_use,:][:,data_to_use][:,order_OS]
#   print data_pcaT

    order1=order_OS

#   order2=range(len(probes_to_use))

    a=hcluster.pdist(data_pcaT[:,order1],'correlation')
    b=hcluster.linkage(a,'complete')
    b=hcluster.linkage(data_pcaT[:,order1],'ward')
    c=hcluster.dendrogram(b,orientation='bottom',no_labels=True, no_plot=True)
    order2=c['leaves']

    order1.reverse()


    ax3 = SubplotZero(fig, 121)
    fig.add_subplot(ax3)
#   ax3= fig.add_subplot(121)
    #now run Gage on gene clusters!yeah!
#   mapping=[(c['color_list'][i], c['leaves'][i]) for i in range(len(c['color_list']))]
#   clusters=[]
#   for color in c['color_list']:
#       if color not in clusters:
#           clusters.append(color)

    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.0, 1.0),
                      (0.5, 1.0, 1.0),
                      (1.0, 0.0, 0.0))}
    my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,64)

    ccc = ax3.imshow( data_pca[order1,:][:,order2], interpolation='nearest', cmap=my_cmap, aspect='auto',  )
    ddd = ax3.set_yticks(range(len(order2)))
    eee = ax3.set_yticklabels(range(len(data['colnames'][data_to_use])),numpy.array(data['colnames'][data_to_use])[order1] )

    kwargs= {'cax':ax3,'ax':ax3}
#
    ax4 = SubplotZero(fig, 122, sharey=ax3)
    fig.add_subplot(ax4)

#   print len (order1)
    for i in range(len(data['stemness'][data_to_use])):
#       if int(data['stemness'][data_to_use][order1][i]) == 5:
        ax4.plot([-1,1],[i,i], '-',linewidth=1.,color=colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)])
#marker='s', linewidths=0,
    if mode != None: #this is used to report the AML NK genes mutations.
        for index,key in enumerate(['cluster','cebpa' , 'npm1', 'flt3']):
            ax4.text(-1+index*2,-10,key,rotation ='vertical', size=6)
        for i in range(len(data['stemness'][data_to_use])):
            for index,key in enumerate(['cebpa' , 'npm1', 'flt3']):
                if data[key][data_to_use][i] == 1:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='k')
                elif data[key][data_to_use][i] == -1:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='y')
                else:
                    ax4.plot([1+index*2,3+index*2],[i,i], '-',linewidth=1.,color='w')


    for i in range(len(x)):
        ax4.scatter(data['OS'][data_to_use][order_OS][i],i)



#   ax4.set_xlim( -1, int(data['OS'][data_to_use][order_OS][-1]) ) #hack!!
    for ax in [ax3,ax4]:
        for direction in ["left", "right", "bottom", "top"]:
            ax.axis[direction].set_visible(True)



    pylab.colorbar(ccc, ax=ax4 ,fraction= .05, orientation ='vertical', )

    fig.savefig('heatmap_genes.pdf',dpi=600)
#   fig.show()
#   pylab.close()

    the_colors={}
    cat=[]
    for i in range(len(data['stemness'][data_to_use])):
        if data['stemness'][data_to_use][order1][i] not in the_colors.keys():
#           print data['stemness'][data_to_use][order1][i], colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
#           print the_colors.keys()
            the_colors[data['stemness'][data_to_use][order1][i]] = colors[ int(data['stemness'][data_to_use][order1][i]) % len(colors)]
    print the_colors
    pass

def cox_analysis(data,probes_to_use=None, data_to_use=None, cluster_to_use=None ,mode = None, seed = None, granular_data=None,GS=None, GS_to_use=None):
    """Cox analysis.
    runs cox analysis. for ALL variables.
    Note that the variables are transformed into binary data.

    """
    import copy
    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
    if probes_to_use == None:
        data_to_use=range(data['data'].shape[0])
    if seed==None:
        seed=999
    # row names have to be changed, as affy numbering sucks big time in any language to have variable names.

    names=[]
    convertion_names={}
    for i,n in enumerate(probes_to_use):
        names.append('x'+str(i))
        convertion_names['x'+str(i)] = data['rownames'][n]

    OS_dic=dict()
    for i,cancer in enumerate(data_to_use):
#       print '__'+ str(data['colnames'][cancer])
        if data['stemness'][cancer] in cluster_to_use:
            for j,p_index in enumerate(probes_to_use):
#           print p_index
#               OS_dic[data['rownames'][p_index]]=[]
                if data['OS'][cancer] !='0':
                    if granular_data is None:
                        OS_dic[names[j]]  = numpy.array([data['data'][p_index,data_to_use]])
                    else:
                        OS_dic[names[j]] = numpy.round(numpy.array([data['data'][p_index,data_to_use]]),decimals=granular_data)
#                       print OS_dic[names[j]]
    for j,p_index in enumerate(probes_to_use):
        OS_dic[names[j]]=robjects.Vector(numpy.array(OS_dic[names[j]]).T)
#now add the remaining keys;
    keys = [ 'OS', 'cebpa',  'EventOS', 'gender', 'age', 'npm1', 'blast_cells', 'flt3']#,'hsc_score']
    keys = [ 'OS',  'EventOS', 'age' , 'p_FISH_pos_HSC']
    for key in keys:
        print 'Processing : ' + key
        if key not in  [ 'OS',  'EventOS','cebpa',   'gender', 'npm1',  'flt3', 'age' , 'FISH_Status','p_FISH_pos_HSC' ]:
            c_array = copy.copy(numpy.array(data[key][data_to_use]))
            print c_array
            u_i = numpy.argwhere(c_array< numpy.mean(c_array)).reshape(-1)
            d_i = numpy.argwhere(c_array >= numpy.mean(c_array)).reshape(-1)
            c_array[u_i] = 0
            c_array[d_i] = 1
            OS_dic[key]=robjects.IntVector(c_array.T)
        else:
            OS_dic[key]=robjects.Vector(numpy.array(data[key][data_to_use]).T)
    #Gage results are put here in the dic: (from Gage_analysis.py)
    if GS:
        print 'adding %s to the genes for rsf'%GS
        gs_names=data['Gage_gs_order']
        matrix=data[GS]
        for i,name in enumerate(numpy.array(gs_names)[GS_to_use]):
            print name
            if numpy.sum(numpy.isnan(matrix[data_to_use,GS_to_use[i]])) < 1:
                OS_dic[name] = numpy.round( matrix[data_to_use,GS_to_use[i]] ,decimals=granular_data)

    #put clusters in binary form and as variables. so we see if they have a meaning.
    for cluster in cluster_to_use:
        print cluster
        print cluster_to_use
        c=numpy.zeros(len(data_to_use))
        c[ numpy.argwhere(data['stemness'] == str(cluster)) [0].T] = 1.
        for ii,s in enumerate(data['stemness'][data_to_use]):
            if int(s) == int(cluster):
                c[ii] = 1
            else:
                c[ii] = 0
#       print numpy.array(c)
        OS_dic['cluster_'+str(cluster)]=robjects.IntVector(numpy.array(c).T)
#       print OS_dic['cluster_'+str(cluster)]
    R_OS=robjects.DataFrame(OS_dic)
    print R_OS
    robjects.r('''
    cox_rsf <- function(pbc){
    library(survcomp)

    if (library("survival", logical.return = TRUE)
        & library("Hmisc", logical.return = TRUE))
    {
    cox.weights <- function(rsf.f, rsf.data) {
        event.names <- all.vars(rsf.f)[1:2]
        print(event.names)
        p <- ncol(rsf.data) - 2
        event.pt <- match(event.names, names(rsf.data))
        predictor.pt <- setdiff(1:ncol(rsf.data), event.pt)
        sapply(1:p, function(j) {
          cox.out <- coxph(rsf.f, rsf.data[, c(event.pt, predictor.pt[j])])
          pvalue <- summary(cox.out)$coef[5]
          print(summary(cox.out))
          print("-----------------------------------")
          if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
        })
    }

    }
    coxph.control(iter.max = 200000 )
    #cox.out <- coxph( Surv(OS,EventOS) ~  cebpa +blast_cells+  age + gender + npm1 + flt3+cluster_2, data=pbc ,  method="efron" )
    cox.out <- coxph( Surv(OS,EventOS) ~   age +  cluster_1 +age:cluster_1 , data=pbc ,  method="efron" )
    print(summary(cox.out))
    cox.zph(cox.out)
    cox.out <- coxph( Surv(OS,EventOS) ~   age + gender + npm1 + flt3 + cluster_0 +  flt3:cluster_2, data=pbc ,  method="efron" )
    print(summary(cox.out))

    print("coin coin")
    cox.wts <- cox.weights(Surv(OS,EventOS) ~ ., pbc)
#   print(pbc)
#   cox.out <- coxph( Surv(OS,EventOS) ~ cebpa+hsc_sr_5+cluster_core+gender+age+cluster_1+cluster_2+cluster_3+cluster_0+cluster_4+flt3+blast_cells+npm1 , data=pbc ,method="efron" )
#   print(summary(cox.out))
    cox.out <- coxph( Surv(OS,EventOS) ~  cebpa + blast_cells + flt3 + age + gender + npm1 + cluster_1+cluster_2+cluster_3+cluster_0+cluster_4+cluster_5 , data=pbc ,method="exact" )
    cox.out <- coxph( Surv(OS,EventOS) ~  cluster_1+cluster_2+cluster_3+cluster_0+cluster_4 , data=pbc ,method="exact" )
    print(summary(cox.out))

    return(cox.wts)
    }
    ''')
    cox = robjects.globalenv['cox_rsf']


    cox(R_OS)

    return
    pass



def rsf(data,probes_to_use=None, data_to_use=None, cluster_to_use=None ,mode = None, seed = None, granular_data=None,GS=None):
    """docstring for rsf
    seed is

    """
    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
    if probes_to_use == None:
        data_to_use=range(data['data'].shape[0])
    if seed==None:
        seed=999
    # row names have to be changed, as affy numbering sucks big time in any language to have variable names.

    names=[]
    convertion_names={}
    for i,n in enumerate(probes_to_use):
        names.append('x'+str(i))
        convertion_names['x'+str(i)] = data['rownames'][n]

    OS_dic=dict()
    for i,cancer in enumerate(data_to_use):
        print '__'+ str(data['colnames'][cancer])
        # uses specific clusters...
        if cluster_to_use:
            if data['stemness'][cancer] in cluster_to_use:
                for j,p_index in enumerate(probes_to_use):
#               print p_index
#                   OS_dic[data['rownames'][p_index]]=[]
                    if data['OS'][cancer] !='0':
                        if granular_data is None:
                            OS_dic[names[j]]  = numpy.array([data['data'][p_index,data_to_use]])
                        else:
                            OS_dic[names[j]] = numpy.round(numpy.array([data['data'][p_index,data_to_use]]),decimals=granular_data)
#                           print OS_dic[names[j]]
        #uses everything..
        else:
            for j,p_index in enumerate(probes_to_use):
#               print p_index
#                   OS_dic[data['rownames'][p_index]]=[]
                if data['OS'][cancer] !='0':
                    if granular_data is None:
                        OS_dic[names[j]]  = numpy.array([data['data'][p_index,data_to_use]])
                    else:
                        OS_dic[names[j]] = numpy.round(numpy.array([data['data'][p_index,data_to_use]]),decimals=granular_data)

    for j,p_index in enumerate(probes_to_use):
        OS_dic[names[j]]=robjects.Vector(numpy.array(OS_dic[names[j]]).T)
#now add the remaining keys;
    keys = [ 'OS', 'cebpa', 'stemness', 'EventOS', 'gender', 'age', 'npm1', 'blast_cells', 'flt3']
    for key in keys:
        OS_dic[key]=robjects.Vector(numpy.array(data[key][data_to_use]).T)

    #Gage results are put here in the dic: (from Gage_analysis.py)
    if GS:
        print 'adding %s to the genes for rsf'%GS
        gs_names=data['Gage_gs_order']
        matrix=data[GS]
        for i,name in enumerate(gs_names):
            print name
            if numpy.sum(numpy.isnan(matrix[data_to_use,i])) < 1:
                OS_dic[name] = numpy.round( matrix[data_to_use,i] ,decimals=granular_data)


    #put clusters in binary form and as variables. so we see if they have a meaning.
    if cluster_to_use:
        for cluster in cluster_to_use:
            print cluster
            c=[]
            for s in data['stemness']:
                if s == cluster:
                    c.append(1)
                else:
                    c.append(0)
            OS_dic['cluster_'+str(cluster)]=robjects.Vector(numpy.array(c)[data_to_use].T)
#       print OS_dic['cluster_'+str(cluster)]
    R_OS=robjects.DataFrame(OS_dic)

#   print R_OS


#   return R_OS
    robjects.r('''
            function_rsf <- function(pbc,ntree, seed){
            library(randomSurvivalForest)
            set.seed(seed)
            mtry <- (dim(pbc)[2])^(3/4)
            print(mtry)
            nrep=300
#           print(pbc)
            #Weights calculated with random forest-
#           vs3 <- varSel(Surv(OS,EventOS) ~ ., pbc, method= "vh", nodesize=2, mtry=mtry,big.data=TRUE, nsplit=3 , nrep=nrep, ntree=ntree)

            #Weights calculted with a numivariate cox model
            if (library("survival", logical.return = TRUE)
                & library("Hmisc", logical.return = TRUE))
            {
              cox.weights <- function(rsf.f, rsf.data) {
                event.names <- all.vars(rsf.f)[1:2]
                p <- ncol(rsf.data) - 2
                event.pt <- match(event.names, names(rsf.data))
                predictor.pt <- setdiff(1:ncol(rsf.data), event.pt)
                sapply(1:p, function(j) {
                  cox.out <- coxph(rsf.f, rsf.data[, c(event.pt, predictor.pt[j])])
                  pvalue <- summary(cox.out)$coef[5]
                  if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
                })
              }

              rsf.f <- as.formula(Surv(Time, Censoring) ~ .)
              cox.wts <- cox.weights(Surv(OS,EventOS) ~ ., pbc)
              vs3 <- varSel(Surv(OS,EventOS) ~ ., pbc, method= "vh", nodesize=2, mtry=mtry,big.data=TRUE, nsplit=3 , nrep=nrep, ntree=ntree, predictorWt = cox.wts)

            }
#
            return(vs3$varselect)
            }
            ''')

    r_rsf_out =robjects.globalenv['function_rsf']
#   print robjects.globalenv['function_rsf']

    aaa=r_rsf_out(R_OS,3000,seed)
#   importance=robjects.r.plot(aaa)
    imp_variables=[]
    for i in numpy.array(aaa.rownames):
        if i in convertion_names.keys():
            imp_variables.append(convertion_names[i])
        else:
            imp_variables.append(i)

    return(aaa, imp_variables)
    pass


def Go_analysis(probelist):
    """docstring for Go_analysis
    this is an R function
    """

    robjects.r('''
    hyperGOinlist = function(probelist){

        library(GO.db); library(Category);library(hgu133plus2.db);library(GOstats)
        foundGOterms = c("GO:0000082","GO:0000216","GO:0000226","GO:0000278","GO:0000819","GO:0000910","GO:0001709","GO:0002367","GO:0002437","GO:0005976","GO:0006081","GO:0006082","GO:0006091","GO:0006342","GO:0006403","GO:0006446","GO:0006518","GO:0006520","GO:0006730","GO:0006790","GO:0006944","GO:0006984","GO:0007006","GO:0007050","GO:0007091","GO:0007618","GO:0008366","GO:0008610","GO:0009057","GO:0009306","GO:0009308","GO:0009896","GO:0010564","GO:0015911","GO:0015931","GO:0015980","GO:0016052","GO:0016441","GO:0019884","GO:0022403","GO:0022900","GO:0031047","GO:0031345","GO:0032507","GO:0032623","GO:0033554","GO:0034440","GO:0035383","GO:0035821","GO:0035966","GO:0036075","GO:0042180","GO:0042246","GO:0042769","GO:0043242","GO:0043244","GO:0043414","GO:0043487","GO:0043647","GO:0043933","GO:0044248","GO:0044262","GO:0044282","GO:0045786","GO:0045913","GO:0046782","GO:0046907","GO:0048002","GO:0048278","GO:0048524","GO:0050434","GO:0050658","GO:0050709","GO:0050777","GO:0051186","GO:0051236","GO:0051439","GO:0051651","GO:0051817","GO:0071843","GO:0072593","GO:0072594","GO:0090068","GO:0090342","GO:0097190","GO:2001022","GO:0000082","GO:0000216","GO:0000910","GO:0001881","GO:0002263","GO:0002366","GO:0002437","GO:0002562","GO:0006446","GO:0006518","GO:0006520","GO:0006944","GO:0006984","GO:0007006","GO:0007091","GO:0009308","GO:0009791","GO:0009895","GO:0016445","GO:0019884","GO:0021695","GO:0021766","GO:0021987","GO:0022900","GO:0030307","GO:0031018","GO:0031345","GO:0032506","GO:0033013","GO:0034440","GO:0042168","GO:0042246","GO:0042558","GO:0043299","GO:0045454","GO:0045767","GO:0045912","GO:0046148","GO:0046782","GO:0046902","GO:0048002","GO:0048278","GO:0048524","GO:0050434","GO:0050792","GO:0051186","GO:0051439","GO:0051705","GO:0071843","GO:0072594","GO:0090068","GO:2000243","GO:2001022","GO:0000216","GO:0006446","GO:0006730","GO:0006984","GO:0007006","GO:0007062","GO:0009994","GO:0019047","GO:0019048","GO:0019059","GO:0019884","GO:0030069","GO:0031069","GO:0032273","GO:0032507","GO:0034440","GO:0036075","GO:0043487","GO:0045454","GO:0045767","GO:0045768","GO:0046796","GO:0048002","GO:0048806","GO:0051439","GO:0051651","GO:0051701","GO:0051983","GO:0071843","GO:0090068","GO:2001022","GO:0000082","GO:0000083","GO:0000086","GO:0000216","GO:0000226","GO:0000278","GO:0000710","GO:0000819","GO:0001824","GO:0001894","GO:0002200","GO:0002250","GO:0002260","GO:0002263","GO:0002366","GO:0002367","GO:0002377","GO:0002443","GO:0002562","GO:0002566","GO:0002686","GO:0002700","GO:0002702","GO:0002905","GO:0005975","GO:0006066","GO:0006069","GO:0006082","GO:0006091","GO:0006139","GO:0006342","GO:0006403","GO:0006446","GO:0006518","GO:0006520","GO:0006629","GO:0006725","GO:0006730","GO:0006766","GO:0006790","GO:0006793","GO:0006805","GO:0006818","GO:0006970","GO:0006979","GO:0006984","GO:0007006","GO:0007050","GO:0007052","GO:0007062","GO:0007080","GO:0007091","GO:0007131","GO:0007569","GO:0007618","GO:0007622","GO:0008104","GO:0008610","GO:0009057","GO:0009059","GO:0009308","GO:0009310","GO:0009314","GO:0009408","GO:0009410","GO:0009890","GO:0009892","GO:0009893","GO:0009896","GO:0010033","GO:0010243","GO:0010259","GO:0010467","GO:0010564","GO:0010604","GO:0010605","GO:0010638","GO:0010941","GO:0010948","GO:0012501","GO:0015031","GO:0015931","GO:0015980","GO:0016042","GO:0016051","GO:0016052","GO:0016358","GO:0016445","GO:0017144","GO:0018904","GO:0019047","GO:0019048","GO:0019058","GO:0019059","GO:0019079","GO:0019080","GO:0019083","GO:0019222","GO:0019538","GO:0019835","GO:0019884","GO:0021544","GO:0021756","GO:0021766","GO:0022402","GO:0022403","GO:0022411","GO:0022415","GO:0022602","GO:0022607","GO:0022900","GO:0030069","GO:0030307","GO:0030522","GO:0030901","GO:0031047","GO:0031323","GO:0031324","GO:0031325","GO:0031345","GO:0031348","GO:0031647","GO:0032102","GO:0032387","GO:0032507","GO:0032602","GO:0032606","GO:0032612","GO:0032886","GO:0032922","GO:0032943","GO:0033013","GO:0033043","GO:0033327","GO:0033554","GO:0033627","GO:0034440","GO:0034641","GO:0034763","GO:0035176","GO:0035383","GO:0040008","GO:0040014","GO:0042023","GO:0042092","GO:0042168","GO:0042180","GO:0042558","GO:0042592","GO:0042752","GO:0042769","GO:0043412","GO:0043414","GO:0043449","GO:0043487","GO:0043647","GO:0043933","GO:0044092","GO:0044093","GO:0044248","GO:0044249","GO:0044255","GO:0044260","GO:0044262","GO:0044282","GO:0044283","GO:0045132","GO:0045143","GO:0045184","GO:0045454","GO:0045767","GO:0045786","GO:0045823","GO:0045833","GO:0045933","GO:0046148","GO:0046209","GO:0046434","GO:0046483","GO:0046649","GO:0046661","GO:0046677","GO:0046782","GO:0046796","GO:0046902","GO:0046907","GO:0046950","GO:0048002","GO:0048145","GO:0048146","GO:0048147","GO:0048512","GO:0048519","GO:0048522","GO:0048523","GO:0048524","GO:0048610","GO:0050434","GO:0050658","GO:0050777","GO:0050790","GO:0050792","GO:0051090","GO:0051171","GO:0051172","GO:0051186","GO:0051224","GO:0051236","GO:0051303","GO:0051304","GO:0051310","GO:0051321","GO:0051439","GO:0051640","GO:0051649","GO:0051651","GO:0051656","GO:0051701","GO:0051705","GO:0051726","GO:0055086","GO:0060255","GO:0060759","GO:0060761","GO:0060996","GO:0070192","GO:0070271","GO:0070663","GO:0070665","GO:0070727","GO:0070997","GO:0071216","GO:0071842","GO:0071843","GO:0072593","GO:0072594","GO:0072657","GO:0080090","GO:0080135","GO:0090068","GO:0090342","GO:0090398","GO:0097190","GO:2000177","GO:2000178","GO:2000241","GO:2000242","GO:2000243","GO:2000772","GO:2001021","GO:2001022","GO:2001235","GO:0000082","GO:0000083","GO:0000086","GO:0000216","GO:0000226","GO:0000278","GO:0000819","GO:0000910","GO:0001541","GO:0002260","GO:0005975","GO:0005976","GO:0006081","GO:0006082","GO:0006091","GO:0006139","GO:0006403","GO:0006446","GO:0006518","GO:0006520","GO:0006629","GO:0006644","GO:0006725","GO:0006730","GO:0006766","GO:0006790","GO:0006818","GO:0006979","GO:0007006","GO:0007050","GO:0007052","GO:0007062","GO:0007076","GO:0007091","GO:0007098","GO:0007099","GO:0007131","GO:0007440","GO:0007566","GO:0007569","GO:0007617","GO:0007618","GO:0007625","GO:0008104","GO:0008366","GO:0008610","GO:0009057","GO:0009266","GO:0009308","GO:0009408","GO:0009791","GO:0009892","GO:0010035","GO:0010171","GO:0010243","GO:0010467","GO:0010564","GO:0010604","GO:0010605","GO:0010638","GO:0010639","GO:0010770","GO:0010948","GO:0015031","GO:0015931","GO:0015980","GO:0016042","GO:0016051","GO:0016441","GO:0017144","GO:0019048","GO:0019098","GO:0019538","GO:0019835","GO:0019884","GO:0022402","GO:0022403","GO:0022415","GO:0022602","GO:0022607","GO:0022900","GO:0030307","GO:0031047","GO:0031324","GO:0031345","GO:0032231","GO:0032506","GO:0032507","GO:0032612","GO:0032846","GO:0032886","GO:0033013","GO:0033043","GO:0033554","GO:0034440","GO:0034641","GO:0040014","GO:0042023","GO:0042168","GO:0042180","GO:0042558","GO:0042769","GO:0043242","GO:0043412","GO:0043414","GO:0043449","GO:0043487","GO:0043933","GO:0044092","GO:0044093","GO:0044248","GO:0044249","GO:0044255","GO:0044260","GO:0044262","GO:0044282","GO:0044283","GO:0044403","GO:0045104","GO:0045184","GO:0045454","GO:0045767","GO:0045785","GO:0045786","GO:0045927","GO:0045932","GO:0046148","GO:0046209","GO:0046483","GO:0046796","GO:0046907","GO:0048002","GO:0048546","GO:0048565","GO:0050658","GO:0050803","GO:0051186","GO:0051236","GO:0051439","GO:0051651","GO:0051701","GO:0051726","GO:0055086","GO:0060021","GO:0060135","GO:0060322","GO:0060323","GO:0060324","GO:0060325","GO:0070271","GO:0070727","GO:0070997","GO:0071842","GO:0071843","GO:0072593","GO:0072594","GO:0080135","GO:0090068","GO:0090398","GO:0097190","GO:2001022","GO:0000216","GO:0006446","GO:0019884","GO:0034440","GO:0048002","GO:0048806","GO:0051439","GO:0051701","GO:2001022","GO:0001569","GO:0001881","GO:0001919","GO:0001921","GO:0002053","GO:0002757","GO:0006900","GO:0008360","GO:0009595","GO:0010464","GO:0010829","GO:0016441","GO:0022604","GO:0030100","GO:0031047","GO:0032273","GO:0032612","GO:0033205","GO:0040014","GO:0040018","GO:0043414","GO:0045807","GO:0046677","GO:0050854","GO:0051650","GO:0060491","GO:0071216","GO:0071604","GO:1900076","GO:1900077","GO:0002088","GO:0002089","GO:0006342","GO:0007062","GO:0010171","GO:0014010","GO:0031333","GO:0032272","GO:0035821","GO:0043242","GO:0051702","GO:0051817","GO:0060251","GO:0060253","GO:0060322","GO:0060323","GO:0001759","GO:0002053","GO:0002089","GO:0003151","GO:0003279","GO:0006342","GO:0007405","GO:0007569","GO:0009994","GO:0010171","GO:0010464","GO:0010522","GO:0017145","GO:0021772","GO:0021846","GO:0021988","GO:0030901","GO:0033627","GO:0034381","GO:0045168","GO:0045766","GO:0048645","GO:0051209","GO:0051282","GO:0051283","GO:0060322","GO:0060323","GO:0060324","GO:0060325","GO:0060349","GO:0060411","GO:0070252","GO:0072676","GO:0090342","GO:0090343","GO:0090398","GO:0097006","GO:2000772","GO:0000819","GO:0006342","GO:0007062","GO:0010522","GO:0014010","GO:0051209","GO:0051282","GO:0051283","GO:0060711","GO:0000086","GO:0000819","GO:0001541","GO:0001542","GO:0001824","GO:0001825","GO:0002053","GO:0002088","GO:0002089","GO:0002347","GO:0002418","GO:0003205","GO:0003206","GO:0003279","GO:0006342","GO:0006446","GO:0006482","GO:0007062","GO:0007091","GO:0007098","GO:0007405","GO:0007569","GO:0007589","GO:0007595","GO:0009895","GO:0009994","GO:0010464","GO:0010522","GO:0010883","GO:0010884","GO:0010886","GO:0010948","GO:0014047","GO:0014855","GO:0019048","GO:0019059","GO:0019835","GO:0021846","GO:0021915","GO:0022612","GO:0030168","GO:0030522","GO:0030728","GO:0031018","GO:0031344","GO:0031345","GO:0031346","GO:0031647","GO:0031668","GO:0032271","GO:0032272","GO:0032273","GO:0032350","GO:0032845","GO:0033627","GO:0034086","GO:0034329","GO:0034330","GO:0034405","GO:0035148","GO:0035821","GO:0042446","GO:0042594","GO:0043242","GO:0043414","GO:0043449","GO:0043576","GO:0043647","GO:0044403","GO:0045454","GO:0045596","GO:0045637","GO:0045646","GO:0045787","GO:0045807","GO:0045833","GO:0046677","GO:0046718","GO:0046782","GO:0046888","GO:0048146","GO:0048286","GO:0048469","GO:0048524","GO:0050434","GO:0050792","GO:0051172","GO:0051271","GO:0051602","GO:0051653","GO:0051701","GO:0051702","GO:0051817","GO:0051828","GO:0052126","GO:0060008","GO:0060038","GO:0060041","GO:0060349","GO:0060411","GO:0060419","GO:0060688","GO:0060760","GO:0071216","GO:0071496","GO:0090287","GO:0090398","GO:0097190","GO:2000241","GO:2000243","GO:2001021","GO:2001235","GO:0001525","GO:0001759","GO:0001824","GO:0001832","GO:0002053","GO:0002088","GO:0002089","GO:0002200","GO:0002377","GO:0002562","GO:0003170","GO:0003179","GO:0003205","GO:0003206","GO:0003279","GO:0006482","GO:0007162","GO:0007405","GO:0007632","GO:0008542","GO:0009895","GO:0009994","GO:0010171","GO:0010464","GO:0010469","GO:0010522","GO:0010524","GO:0014855","GO:0016445","GO:0021510","GO:0021517","GO:0021772","GO:0021846","GO:0021915","GO:0021988","GO:0031018","GO:0031345","GO:0032411","GO:0033627","GO:0034405","GO:0035019","GO:0035148","GO:0042246","GO:0043270","GO:0043299","GO:0043414","GO:0045136","GO:0045168","GO:0045596","GO:0045684","GO:0045766","GO:0045787","GO:0045833","GO:0045912","GO:0048145","GO:0048521","GO:0048645","GO:0050432","GO:0050920","GO:0050922","GO:0055021","GO:0060038","GO:0060043","GO:0060322","GO:0060323","GO:0060349","GO:0060411","GO:0060419","GO:0060688","GO:0060744","GO:0061383","GO:0070849","GO:0072657","GO:0090150","GO:0090342","GO:0097190","GO:0000082","GO:0000083","GO:0000086","GO:0000216","GO:0000226","GO:0000278","GO:0000819","GO:0000910","GO:0001556","GO:0002562","GO:0007050","GO:0007052","GO:0007076","GO:0007080","GO:0007091","GO:0007098","GO:0007569","GO:0009994","GO:0010564","GO:0010639","GO:0010948","GO:0014855","GO:0022402","GO:0022403","GO:0032465","GO:0032507","GO:0032886","GO:0033205","GO:0033598","GO:0035821","GO:0040001","GO:0045185","GO:0045786","GO:0045787","GO:0051293","GO:0051303","GO:0051310","GO:0051313","GO:0051439","GO:0051640","GO:0051651","GO:0051653","GO:0051656","GO:0051726","GO:0051817","GO:0051983","GO:0051984","GO:0055021","GO:0060038","GO:0060043","GO:0060045","GO:0060419")

        allg = get("hgu133plus2ENTREZID")
        allg = as.data.frame(unlist(as.list(allg)))

        entrez.ids <- unique(allg[probelist,])

        params <- new("GOHyperGParams", geneIds=entrez.ids, annotation=c("hgu133plus2"), ontology="BP", pvalueCutoff=0.05, conditional=FALSE, testDirection="over")

        hgOver = hyperGTest(params)
        t = summary(hgOver)

        return = t[t[,"GOBPID"] %in% foundGOterms,][c("Term", "Pvalue")]

    }''')

    r_go_analysis =robjects.globalenv['hyperGOinlist']

    results = r_go_analysis(probelist)

    return results

    pass


def incremental_rsf(data,probes_to_use=None, data_to_use=None, cluster_to_use=None ,mode = None, seed = None, \
                    granular_data=None,GS=None, GS_to_use=None,binarize =False):
    """docstring for incremental_rsf

    this is run after rsf, when the variable selection has been done.
    this would tippically run fast, with a few important variables as input.

    """
    import copy
    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
    if probes_to_use == None:
        data_to_use=range(data['data'].shape[0])
    if seed==None:
        seed=999
    # row names have to be changed, as affy numbering sucks big time in any language to have variable names.

    names=[]
    convertion_names={}
    for i,n in enumerate(probes_to_use):
        names.append('x'+str(i))
        convertion_names['x'+str(i)] = data['rownames'][n]

    OS_dic=dict()
    for i,cancer in enumerate(data_to_use):
#       print '__'+ str(data['colnames'][cancer])
        if data['stemness'][cancer] in cluster_to_use:
            for j,p_index in enumerate(probes_to_use):
#           print p_index
#               OS_dic[data['rownames'][p_index]]=[]
                if data['OS'][cancer] !='0':
                    if granular_data is None:
                        OS_dic[names[j]]  = numpy.array([data['data'][p_index,data_to_use]])
                    else:
                        OS_dic[names[j]] = numpy.round(numpy.array([data['data'][p_index,data_to_use]]),decimals=granular_data)
#                       print OS_dic[names[j]]
    for j,p_index in enumerate(probes_to_use):
        OS_dic[names[j]]=robjects.Vector(numpy.array(OS_dic[names[j]]).T)
#now add the remaining keys;
    keys = [ 'OS', 'cebpa',  'EventOS', 'gender', 'age', 'npm1', 'blast_cells', 'flt3']#,'hsc_score']
#   keys = [ 'OS',  'EventOS']
    if binarize :
        for key in keys:
            if key not in  [ 'OS',  'EventOS']:
                c_array = copy.copy(numpy.array(data[key][data_to_use]))
                u_i = numpy.argwhere(c_array<c_array.mean()).reshape(-1)
                d_i = numpy.argwhere(c_array >= c_array.mean()).reshape(-1)
                c_array[u_i] = 0
                c_array[d_i] = 1
                OS_dic[key]=robjects.Vector(c_array.T)
            else:
                OS_dic[key]=robjects.Vector(numpy.array(data[key][data_to_use]).T)
    else:
        for key in keys:
            print key
            OS_dic[key]=robjects.Vector(numpy.array(data[key][data_to_use]).T)

    #Gage results are put here in the dic: (from Gage_analysis.py)
    if GS:
        print 'adding %s to the genes for rsf'%GS
        gs_names=data['Gage_gs_order']
        matrix=data[GS]
        for i,name in enumerate(numpy.array(gs_names)[GS_to_use]):
            print name
            if numpy.sum(numpy.isnan(matrix[data_to_use,GS_to_use[i]])) < 1:
                OS_dic[name] = numpy.round( matrix[data_to_use,GS_to_use[i]] ,decimals=granular_data)


    #put clusters in binary form and as variables. so we see if they have a meaning.
    for cluster in cluster_to_use:
        print cluster
        print cluster_to_use
        c=numpy.zeros(len(data_to_use))
        c[ numpy.argwhere(data['stemness'][data_to_use] == cluster ).reshape(-1)] = 1.
        print numpy.array(c)
        OS_dic['cluster_'+str(cluster)]=robjects.Vector(numpy.array(c).T)
#       print OS_dic['cluster_'+str(cluster)]
    R_OS=robjects.DataFrame(OS_dic)

#   print R_OS


#   return R_OS
    robjects.r('''
            function_rsf <- function(pbc,ntree, seed){
            library(randomSurvivalForest)
            library(grDevices)
            set.seed(seed)
#           print(pbc$EventOS)
#           print(dim(pbc))
            mtry <- (dim(pbc)[2])^2#(3/4)
            pbc.out <- rsf(as.formula("Survrsf(OS,EventOS)~.") , pbc, mtry=mtry, ntree=ntree, splitrule = "logrank", forest=T, big.data=FALSE)
            print(pbc.out)
            print(colnames(pbc))
            imp <- pbc.out$importance
            pnames <- pbc.out$predictorNames
#           print(pnames)
            pnames.order <- pnames[rev(order(imp))]
            pdf('variable_importance.pdf')
#           plot.variable(pbc.out,3,type="rel.freq",partial=T,n.pred=length(pnames))
            plot.variable(pbc.out,3,type="time",partial=T,n.pred=length(pnames))
            dev.off()
            pdf('error_rsf.pdf')
            plot.error(pbc.out)

            dev.off()
            n.pred <- length(pnames)
            pbc.err <- rep(0, n.pred)
            max_var <- 10
            for (k in 1:n.pred){
#           for (k in 1:max_var){
#               formula <- as.formula(paste("Survrsf(OS,EventOS)~",paste(pnames.order[1:k],collapse="+")))
                formula <- as.formula(paste("Survrsf(OS)~",paste(pnames.order[1:k],collapse="+")))
                print(paste(k,"/",n.pred))
                pbc.err[k] <- rsf( formula , pbc, ntree=ntree, mtry=mtry, splitrule="logrank", big.data=TRUE)$err.rate[ntree]
                }

            pbc.imp.out <- as.data.frame(cbind(round(rev(sort(imp)),4),round(pbc.err,4),round(-diff(c(0.5,pbc.err)),4)),row.names=pnames.order)

            colnames(pbc.imp.out) <-c("Imp","Err","Drop Err")

            pdf('error_drop.pdf')
            plot(pbc.err)
            dev.off()

            print(pbc.imp.out)
            print(pbc.err)
            return(pbc.imp.out)
            }
            ''')

    r_rsf_out =robjects.globalenv['function_rsf']

    robjects.r('''
            function_final_rsf <- function(pbc,ntree, seed, n.pred){
            library(randomSurvivalForest)
            library(grDevices)
            set.seed(seed)
            mtry <- (dim(pbc)[2])^(3/4)
#           pbc.out <- rsf(as.formula("Survrsf(OS,EventOS)~.") , pbc, mtry=mtry, ntree=ntree, splitrule = "logrank", forest=T, big.data=FALSE)
            pbc.out <- rsf(as.formula("Survrsf(OS)~.") , pbc, mtry=mtry, ntree=ntree, splitrule = "logrank", forest=T, big.data=FALSE)
            imp <- pbc.out$importance
            pnames <- pbc.out$predictorNames
            pnames.order <- pnames[rev(order(imp))]

            #pbc.err <- rep(0, n.pred) # this one we already know...

            formula <- as.formula(paste("Survrsf(OS,EventOS)~",paste(pnames.order[1:n.pred],collapse="+")))
            pbc.final.out <- rsf( formula , pbc, ntree=ntree, mtry=mtry, splitrule="logrank", big.data=TRUE)

            return(pbc.final.out)
            }
            ''')

    grow_final_forest = robjects.globalenv['function_final_rsf']


    robjects.r('''
    cox_rsf <- function(pbc){
    library(randomSurvivalForest)

    if (library("survival", logical.return = TRUE)
        & library("Hmisc", logical.return = TRUE))
    {
    cox.weights <- function(rsf.f, rsf.data) {
        event.names <- all.vars(rsf.f)[1:2]
        print(event.names)
        p <- ncol(rsf.data) - 2
        event.pt <- match(event.names, names(rsf.data))
        predictor.pt <- setdiff(1:ncol(rsf.data), event.pt)
        sapply(1:p, function(j) {
          cox.out <- coxph(rsf.f, rsf.data[, c(event.pt, predictor.pt[j])])
          pvalue <- summary(cox.out)$coef[5]
          print(summary(cox.out))
          print("-----------------------------------")
          if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
        })
    }

    }
    print("coin coin")
    cox.wts <- cox.weights(Surv(OS,EventOS) ~ ., pbc)
#   print(pbc)
#   cox.out <- coxph( Surv(OS,EventOS) ~ cebpa+hsc_sr_5+cluster_core+gender+age+cluste1+cluster_2+cluster_3+cluster_0+cluster_4+flt3+blast_cells+npm1 , data=pbc ,method="breslow" )
#   print(summary(cox.out))

    print("-----------------------------------")
#   cox.out <- coxph( Surv(OS,EventOS) ~  hsc_score + cebpa + blast_cells + flt3 + age + gender + npm1  , data=pbc ,method="breslow" )
    cox.out <- coxph( Surv(OS) ~  hsc_score + cebpa + blast_cells + flt3 + age + gender + npm1  , data=pbc ,method="breslow" )
    print(summary(cox.out))
    print("-----------------------------------")
#   cox.out <- coxph( Surv(OS,EventOS) ~   cluster_0 + cluster_1+cluster_2+cluster_3+cluster_4+cluster_5 +cebpa + blast_cells+  flt3 + age + gender + npm1   , data=pbc ,method="breslow" )
    cox.out <- coxph( Surv(OS) ~   cluster_0 + cluster_1+cluster_2+cluster_3+cluster_4+cluster_5 +cebpa + blast_cells+  flt3 + age + gender + npm1   , data=pbc ,method="breslow" )
    print(summary(cox.out))


    return(cox.wts)
    }
    ''')
    cox = robjects.globalenv['cox_rsf']



#   print robjects.globalenv['function_rsf']

    aaa=r_rsf_out(R_OS,500,seed)

    n_pred=numpy.argmin(numpy.array(aaa)[1])

    Final_Forest  =  grow_final_forest(R_OS,5000,seed,n_pred)
    cox(R_OS)
    importance=robjects.r.plot(aaa)
    imp_variables=[]
    for i in numpy.array(aaa.rownames):
        if i in convertion_names.keys():
            imp_variables.append(convertion_names[i])
        else:
            imp_variables.append(i)

    return(aaa, imp_variables, Final_Forest)
    pass

def multicorr(data):

    robjects.r('''
    bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){
    library(ff, quietly = TRUE)
      require(doMC)
    if(ncore=="all"){
        ncore = multicore:::detectCores()
        registerDoMC(cores = ncore)
    } else{
        registerDoMC(cores = ncore)
    }

    NCOL <- ncol(x)
    print('this is a test0')
    ## test if ncol(x) %% nblocks gives remainder 0
    while (NCOL %% nblocks != 0){
    print('this is a test')
    ##stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
    nblocks <- nblocks + 1
    print('this is a test 1')
    }

    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    corMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))


    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

    ## create all unique combinations of blocks
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)

    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    results <- foreach(i = 1:nrow(COMBS)) %dopar% {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        #COR <- cor(MAT[, G1], MAT[, G2], ...)
        COR <- cor(x[, G1], x[, G2], ...)
        corMAT[G1, G2] <- COR
        corMAT[G2, G1] <- t(COR)
        COR <- NULL
    }

    gc()
    print(corMAT)
    corMAT2 <- matrix(corMAT)
    print(corMAT2)
    return(as.matrix(corMAT2))
    }''')
    bigcorPar_r = robjects.globalenv['bigcorPar']
    cor=numpy.array(bigcorPar_r(data))
    return numpy.array(cor)

#
def test_rsf(data,forest,probes_to_use=None, data_to_use=None, cluster_to_use=None ,mode = None, seed = None, granular_data=None,GS=None, GS_to_use=None ):
    """docstring for incremental_rsf

    this is run after rsf, when the variable selection has been done.
    This test the error rate of the grown forest.

    """
    if data_to_use == None:
        data_to_use=range(data['data'].shape[1])
    if probes_to_use == None:
        data_to_use=range(data['data'].shape[0])
    if seed==None:
        seed=999
    # row names have to be changed, as affy numbering sucks big time in any language to have variable names.

    names=[]
    convertion_names={}
    for i,n in enumerate(probes_to_use):
        names.append('x'+str(i))
        convertion_names['x'+str(i)] = data['rownames'][n]

    OS_dic=dict()
    print data_to_use
    for i,cancer in enumerate(data_to_use):
#       print '__'+ str(data['colnames'][cancer])
        print 'TEST DATA DOESN T HAVE  A CLUSTER NUMBER! TRY TO DO SOMETHIG!!'
        if data['stemness'][cancer] in cluster_to_use:
            for j,p_index in enumerate(probes_to_use):
                print p_index
#               OS_dic[data['rownames'][p_index]]=[]
                if data['OS'][cancer] !='0':
                    if granular_data is None:
                        OS_dic[names[j]]  = numpy.array([data['data'][p_index,data_to_use]])
                    else:
                        OS_dic[names[j]] = numpy.round(numpy.array([data['data'][p_index,data_to_use]]),decimals=granular_data)
#                       print OS_dic[names[j]]
    for j,p_index in enumerate(probes_to_use):
        OS_dic[names[j]]=robjects.Vector(numpy.array(OS_dic[names[j]]).T)
#now add the remaining keys;
    keys = [ 'OS', 'cebpa', 'stemness', 'EventOS', 'gender', 'age', 'npm1', 'blast_cells', 'flt3']
    for key in keys:
        OS_dic[key]=robjects.Vector(numpy.array(data[key][data_to_use]).T)

    #Gage results are put here in the dic: (from Gage_analysis.py)
    if GS:
        print 'adding %s to the genes for rsf'%GS
        gs_names=data['Gage_gs_order']
        matrix=data[GS]
        for i,name in enumerate(numpy.array(gs_names)[GS_to_use]):
            print name
            if numpy.sum(numpy.isnan(matrix[data_to_use,GS_to_use[i]])) < 1:
                OS_dic[name] = numpy.round( matrix[data_to_use,GS_to_use[i]] ,decimals=granular_data)


    #put clusters in binary form and as variables. so we see if they have a meaning.
    for cluster in cluster_to_use:
        print cluster
        c=[]
        for s in data['stemness']:
            if s == cluster:
                c.append(1)
            else:
                c.append(0)
        OS_dic['cluster_'+str(cluster)]=robjects.Vector(numpy.array(c)[data_to_use].T)
#       print OS_dic['cluster_'+str(cluster)]
    R_OS=robjects.DataFrame(OS_dic)

#   print R_OS

    rsf=importr('randomSurvivalForest')
    pred = rsf.predict_rsf(forest,R_OS, split_depth = 1, outcome ="test")

    return pred


    pass

def select_similar_samples(data , data_to_use = None,  p2u=None, max_pop=3):
    """
    selects samples based on correlation.
    """
    import numpy
    if data_to_use == None:
        data_to_use = range(data['data'].shape[1])
    if p2u == None:
        p2u = range(data['data'].shape[0])
#   print p2u
#   print data_to_use
    classes  = {}
    for index,i in enumerate(data['colnames'][data_to_use]):
        if i not in classes:
            classes[i] = 1
        else:
            classes[i] += 1

    print 'Found %d different cell types'%len(classes)
    pop_indexes={}
    class_corrs = {}
    for celltype in classes:
        print celltype
        samples  = numpy.argwhere(data['colnames'] == celltype).reshape(-1)
        pop_indexes[celltype] = samples
#       print samples
        pop_corr = fast_cor(data, data_to_use = samples ,probes_to_use=p2u)
        x = pop_corr.mean(axis=0)
        ranks=x.tolist()
        x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
        order_pop=numpy.array([x.index(i) for i in sorted( x , key=lambda i:( float(i[0]),i[1] ) )])[::-1]
        class_corrs[celltype] = order_pop
    similar_samples = numpy.corrcoef(data['data'][:,numpy.concatenate([ pop_indexes[pop][class_corrs[pop]][:max_pop] for pop in class_corrs])] [ p2u ].T)


    return numpy.concatenate([ pop_indexes[pop][class_corrs[pop]][:max_pop] for pop in class_corrs])


def Kim_paper_heatmap(data ,  test , train = None ,select  ='P.Value', \
                      nprobes = 1000 , out = 'heat_kim_',disp_correlation =False,\
                      cap= 500  , sort_amls =  False ,which_text_key=None)  :
    """
    This makes the heatmap that will make us famous

    """

    import numpy
    hsc =[]
    gmp = []
    try:
        hsc  = numpy.argwhere(data['colnames'] == 'HSC').reshape(-1)
        gmp  = numpy.argwhere(data['colnames'] == 'GMP').reshape(-1)
        lpm  = numpy.argwhere(data['colnames'] == 'late_PM').reshape(-1)
        mpp =  numpy.argwhere(data['colnames'] == 'MPP').reshape(-1)
        pass
    except Exception, e:
        print e


    if train == None:
        train = numpy.concatenate([hsc,gmp])
        #train = numpy.concatenate([hsc,lpm,mpp,gmp])
#    p2u = select_probes_by_variance(data ,  data_to_use  = numpy.concatenate([test , train]), var_thr=4 )
    p2u = limma_multiclass(data, p= 1e-3, pval=select , limit = nprobes, data_to_use =train, one_against_all=0 )['indexes']
#   p2u = MINiML.limma_multiclass(data, p= 5e-2, pval='adj.P.Val',limit = 1000, data_to_use =concatenate([hsc,gmp]), one_against_all=0 )['indexes']

    # easy
#   for i in [t1517, t821, inv16, del7q,inv3]:
#       MINiML.heatmap(copy_aml, probes_to_use= array(p2u)[order_rankdiff],center_g=0,Corder=01, Text_labels =0, prefix =  '%f_%f_gmp_%d_non_centered_filtered_%s'%(p,fc,amount,i), \
#                           data_to_use=  concatenate([gmp, i, hsc,]) , rescale_data=0, ncolors=30,center_g2=0 , \
#                           order1= range(len(concatenate([gmp, i, hsc,])) )    , order2=range(len( order_rankdiff) ),rescale_gene_wise=1)
    # selection of 3 gmp and 3 hsc

    corrrr = numpy.corrcoef(data['data'][:,numpy.concatenate([gmp,  hsc,])][ p2u ].T)
    hsc_corr=corrrr[len(gmp):len(gmp)+len(hsc),len(gmp):len(gmp)+len(hsc)]
    x = hsc_corr.mean(axis=0)
    ranks=x.tolist()
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
    order_hsc=numpy.array([x.index(i) for i in sorted( x , key=lambda i:( float(i[0]),i[1] ) )])[::-1]

    gmp_corr=corrrr[0:len(gmp),0:len(gmp)]
    x = gmp_corr.mean(axis=0)
    ranks=x.tolist()
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
    order_gmp=numpy.array([x.index(i) for i in sorted( x , key=lambda i:( float(i[0]),i[1] ) )])[::-1]
    if disp_correlation:
        heatmap(data,probes_to_use=p2u,data_to_use=numpy.concatenate([gmp[order_gmp[:3]] , \
            hsc[order_hsc[:3]]]),prefix = 'hsc_gmp_centered_',center_g = 1,Corder=1)
        heatmap(data,probes_to_use=p2u,data_to_use=numpy.concatenate([gmp[order_gmp[:3]] , \
            hsc[order_hsc[:3]]]),prefix = 'hsc_gmp_',Corder=1)

    similar_samples = numpy.corrcoef(data['data'][:,numpy.concatenate([gmp[order_gmp[:3]] ,  hsc[order_hsc[:3]]])][ p2u ].T)

#   data['data'][:,concatenate([gmp[order_gmp[:3]] ,  hsc[order_hsc[:3]]])][ p2u ].T
#   close()
#   imshow(similar_samples ,aspect = 'auto',interpolation='nearest')
#   colorbar()

    # now make heatmap of test.
    x_rankdiff = [(data['data'][i,gmp].mean()-data['data'][i,hsc].mean() ) for  i in p2u  ]
    ranks=x_rankdiff
    x_rankdiff=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
    order_rankdiff=[x_rankdiff.index(i) for i in sorted( x_rankdiff , key=lambda i:( float(i[0]),i[1] ) )]

    #and with the same amout of genes:
    x1 = order_rankdiff[:cap]
    x2 = order_rankdiff[-cap:]
    final_probes = numpy.array(p2u)[numpy.concatenate([x1[::-1],x2[::-1]])]

    print ' Using : '+ str(len(final_probes))
    if sort_amls:
        print 'Sorting AMLs'
        final_probes1 = numpy.array(p2u)[x1]
        amls_rankdiff1 = data['data'][final_probes1,:][:,test].mean(axis = 1 ).tolist()
        print len(amls_rankdiff1)
        ranks=amls_rankdiff1
        amls_rankdiff1=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
        order_rankdiff1=[amls_rankdiff1.index(i) for i in sorted( amls_rankdiff1 , key=lambda i:( float(i[0]),i[1] ) )]
        order_rankdiff1

        final_probes2 = numpy.array(p2u)[x2]
        amls_rankdiff2 = data['data'][final_probes2,:][:,test].mean(axis = 1 ).tolist()
        print len(amls_rankdiff2)
        ranks=amls_rankdiff2
        amls_rankdiff2=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
        order_rankdiff2=[amls_rankdiff2.index(i) for i in sorted( amls_rankdiff2 , key=lambda i:( float(i[0]),i[1] ) )]
        final_probes = numpy.array(p2u)[numpy.concatenate([numpy.array(x1)[order_rankdiff1][::-1],numpy.array(x2)[order_rankdiff2]])]
        order_rankdiff2


    print 'Rescaling genewise'
    data_pca = data['data'][final_probes,:][:,numpy.concatenate([gmp[order_gmp[:3]], test, hsc[order_hsc[:3]],])]
    data_pca = data_pca - data_pca.min(axis=0)
    arr = []
    for i in data_pca.T:
        if max(i) > 0:
            arr.append( i/max(i) )
        else:
            arr.append(i)
    data_pca = numpy.array(arr)
    data_pca = data_pca*2-1
    data_pca = data_pca.T


    score_n=data['colnames'][numpy.concatenate([gmp[order_gmp[:3]], test, hsc[order_hsc[:3]] ])].tolist()
    score_n.append('up_dn')
    score_n.append('up_dn')
    score_gmp=[]
    score_hsc=[]
    EE=0
    for i in data_pca.T:
        print score_n[EE],
        EE+=1
        score_gmp.append(i[:cap].mean())
        score_hsc.append(i[-cap:].mean())
        print 'score 1 -->  \t'+ str( sum((i[:cap] > .2))) +'\t',
        print  str( sum((i[-cap:] > .20))) + '\t  <-- score 2'

#       print 'score 1 -->  \t'+ str(i[:cap].mean()) +'\t',
#       print  str(i[-cap:].mean()) + '\t  <-- score 2'
    score_gmp.append(-1)
    score_hsc.append(-1)
    score_gmp.append(1)
    score_hsc.append(1)

    beeswarm(score_gmp, score_n, out=out+'gmp_beeswarm.pdf', plot_median=False, use_mean=False)
    beeswarm(score_hsc, score_n, out=out+'hsc_beeswarm.pdf', title='', plot_median=False, use_mean=False)

    print len(final_probes)
    print numpy.concatenate([gmp[order_gmp[:3]], test, hsc[order_hsc[:4]],])
    print range(len(numpy.concatenate([gmp[order_gmp[:3]], test, hsc[order_hsc[:3]],])) )

    test=select_similar_samples(data,test,max_pop=len(test))

    heatmap(data,
          probes_to_use= final_probes, \
##            probes_to_use= p2u, \
            center_g=0, \
            Corder=1, \
            Text_labels =1,  \
            prefix =  out, \
##            data_to_use= numpy.concatenate([train,test]),\
          data_to_use=  numpy.concatenate([gmp[order_gmp[:3]], test, hsc[order_hsc[:3]],]) , \
            rescale_data=0,  \
#           ncolors=30, \
            center_g2=0 , \
          order1= range(len(numpy.concatenate([gmp[order_gmp[:3]], test, hsc[order_hsc[:3]],])) ),  \
#           order2=[]\
          order2= range(len(numpy.array(p2u)[numpy.concatenate([x1[::-1],x2[::-1]])]  )), \
            rescale_gene_wise=1,\
           which_text_key=which_text_key
            )


    for i in final_probes:
        print data['rownames'][i] , data['data'][i,gmp[order_gmp[:3]]].mean(),   data['data'][i,hsc[order_hsc[:3]]].mean()

    return p2u
    pass


    pass
def limma_fc(data,data_to_use, p = 0.05, pval='adj_p_val', fc = 0):
    """docstring for limma_fc
    groups are determined by stemness!
    """
    classes=[]
    classes_names=[]

    if 'stemness' not in data.keys():
        print 'No stemness key in dataset!'
        return

    for index,i in enumerate(data['stemness'][data_to_use]):
        if i not in classes:
            classes.append(i)
            classes_names.append(data['colnames'][data_to_use][index])

    if len(classes)!=2:
        print 'number of classes != 2 (%d)'%len(classes)
        return
    comparission_name  = '%s_vs_%s' % (classes_names[0],classes_names[1])
    print 'comparing %s'%comparission_name

    stemness = numpy.array(data['stemness'], dtype=int)
    classes_indexes={}
    for i in classes:
        classes_indexes[i]=numpy.array(range(len(stemness)))[stemness == i]
#   print data_to_use
#   print data['stemness'][data_to_use]
#   print 'classes ' + str(classes)
#   for i,j in enumerate(data['stemness'][data_to_use]):
#       print classes_indexes
#       print classes.index(int(j))
#       classes_indexes[str(classes.index(j))].append(i)
#   print classes_indexes

    limma=importr('limma')
    base=importr('base')

    #produce design matrix:
    #make 0s and 1s instead of stemness values
    stemness = numpy.array( data['stemness'][data_to_use], dtype=int)
    stemness = numpy.array(  numpy.floor( stemness / stemness.max()),dtype=int)
#   design = numpy.array([numpy.ones(len(data_to_use)),stemness ],dtype=int)
#   design = design.T
#   design = robjects.Matrix(design)
#   colnames = numpy.array(['vide', comparission_name])
#   base.colnames(design) = robjects.Vector(colnames)

    design_dic={}
    design_dic['vide'] = robjects.Vector ( numpy.ones(len(data_to_use)) )
    design_dic[comparission_name] = robjects.Vector ( stemness )
    r_d_d = robjects.DataFrame( design_dic )
    design = base.as_matrix(r_d_d)
#   print design

    e_matrix = robjects.Matrix(data['data'][:,data_to_use])
#   print comparission_name.replace(' ','.')
    fit = limma.lmFit(e_matrix, design)
    corr_fit =  limma.eBayes(fit)
#   print corr_fit
#   print data['rownames'].shape[0]
#   results = limma.topTable(corr_fit,  adjust="BH",number=data['rownames'].shape[0])
    results = limma.topTable(corr_fit, coef="%s"%comparission_name.replace(' ','.'), adjust="BH",number=data['rownames'].shape[0])
#   print results
    r=numpy.array(results)
    out={}
    out['index'] = numpy.array([ int(i)-1 for i in numpy.array(results.rownames)])
    out['logFC'] = r[0]
    out['av_expr'] = r[1]
    out['p_val'] = r[3]
    out['adj_p_val'] = r[4]
#   out['B'] = r[5]

    #this is the threshold for significantly +- genes
    indexes = range(sum(out[pval] < p))

    sig_genes={}
    sig_genes['indexes'] = out['index'][indexes]
    sig_genes['p_val'] = out[pval][indexes]
    sig_genes['fc'] = []
    for i in sig_genes['indexes']:
        c1=data['data'][i,:][:,classes_indexes[classes[0]]].mean()
        c2=data['data'][i,:][:,classes_indexes[classes[1]]].mean()
        sig_genes['fc'].append(c1-c2)
    sig_genes['fc'] = numpy.array(sig_genes['fc'])

    sig_genes['up'] = sig_genes['indexes'][sig_genes['fc'] > fc]
    sig_genes['dn'] = sig_genes['indexes'][sig_genes['fc'] < -1.*fc]

#   print len(sig_genes['fc'])


    return sig_genes
    pass
#

def smoothTriangle(data,degree=.5,mode=False):
        """performs moving triangle smoothing with a variable degree."""
        """note that if dropVals is False, output length will be identical
        to input length, but with copies of data at the flanking regions"""
        if mode:
            triangle=numpy.array(range(degree)+[degree]+range(degree)[::-1])+1
            smoothed=[]
            for i in range(degree,len(data)-degree*2):
                point=data[i:i+len(triangle)]*triangle
                smoothed.append(sum(point)/sum(triangle))
            smoothed=[smoothed[0]]*(degree+degree/2)+smoothed
            while len(smoothed)<len(data):smoothed.append(smoothed[-1])
            return numpy.array(smoothed)
        else:
            return pa.Series(data).ewm(com=degree).mean()
#
def plot_exprs_histo(data,samples=None, probes = None,smoothdegree = 200, scale= .01):
    """docstring for plot_exprs_histo"""
    import pylab
    if samples == None:
        print 'addidng sample'
        samples = numpy.arange(data['data'].shape[1])
    if probes == None:
        print 'addidng probes'
        probes = numpy.arange(data['data'].shape[0])



    print 'Using %d samples'%len(samples)

    data['colnames'] = numpy.array([  i.strip('*').replace(' (P)','') for i in data['colnames'] ])
    colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#faebd7','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' ,'#faebd7', '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
    bins= numpy.arange(0,18,scale)
    if  data['data'][probes,:].min() < 0 :
        bins= numpy.arange(data['data'][probes,:].min() ,data['data'][probes,:].max(), scale)

    for s in samples:
        a=numpy.histogram(data['data'][probes,s],bins=bins)[0]
        #print smoothTriangle(a,smoothdegree).shape
        #print bins.shape
        pylab.plot(bins[:-1],smoothTriangle(a,smoothdegree),c= colors[data['stemness'][s]],alpha=.3)

    pass


def permute_genes_in_model(xx):
	"""

	index=xx[0]
	model = xx[1]
	sig_len = xx[2]
	iterations = xx[3]

	"""
	from numpy.random import randint
	from numpy.random import RandomState
	index=xx[0]
	model = xx[1]
	sig_len = xx[2]
	iterations = xx[3]
	prng = RandomState(index)

	if sig == None:
		to_run = model.keys()
	else:
		to_run = [sig_len]




	pass

def _pperm_test(xx):
    from numpy.random import randint
    from numpy.random import RandomState
    from numpy import *
    from MINiML import *
    import gc
    index = xx[0]
    nl= xx[1]
    converted_s= xx[2]

    prng = RandomState(index)

    all_pca_models = {}
    len_u =len(unique([len(converted_s[j]) for j in converted_s]))
    for index_,i in enumerate(unique([len(converted_s[j]) for j in converted_s])):
        print '###### perm: ' + str(index) +  ' sig No: ' +str(float(index_)/float(len_u)) + ' ######'
        p2u = prng.randint(len(nl['rownames']), size=i)
        all_pca_models[i] = PCA_paper(nl,probes=p2u,no_plot=1)
    pickle_save(all_pca_models,'perm_%d_%d.pkl'%(index,prng.randint(len(nl['rownames'])*2, size=1)))
    del all_pca_models
    gc.collect()
    return index
    pass

def _pgeom_test(xx):
    #zip(data_to_use,range(len(data_to_use)),[signatures]*len(data_to_use), [fc]*len(data_to_use),[data]*len(data_to_use),[data_to_use]*len(data_to_use))
    from numpy import argwhere,ravel
    from scipy.stats  import hypergeom
    from numpy import concatenate, array, argwhere, ravel
    ii= xx[1]
    index = xx[0]
    signatures = xx[2]
    fc = xx[3]
    rownames=xx[4]
    data_to_use = xx[5]
    top_n = xx[6]
    probe_to_gene=xx[7]

    universe = set(concatenate(array(signatures.values())).ravel())
    print str(float(ii) / float(len(data_to_use))*100.) + ' % done.'
    up_mat = numpy.zeros([1,len(signatures.keys())])
    dn_mat = numpy.zeros([1,len(signatures.keys())])
    mini_fc = 3
    p2u = argwhere(abs(fc) > mini_fc).ravel()
    while len(p2u) < top_n*2:
        mini_fc -= .2
        #print str(index)+'.',
        p2u = argwhere(abs(fc) > mini_fc).ravel()
    sorted_i = get_sorted_index(fc[p2u])
    dn_index = (p2u[sorted_i])[0:top_n]
    up_index = (p2u[sorted_i])[-top_n:]
    dn_genes = concatenate([probe_to_gene.get(i) for i in rownames[dn_index]  if probe_to_gene.get(i)!=None])
    dn_genes = set(dn_genes[dn_genes != 'NA'])
    up_genes = concatenate([probe_to_gene.get(i) for i in rownames[up_index]  if probe_to_gene.get(i)!=None])
    #print str(index)+' :\t'+'\t'.join(up_genes[0:5])
    up_genes = set(up_genes[up_genes != 'NA'])
    r=[[],[],[]]
    for k in [1,2]:
        if k ==1:
            for jj, sig in enumerate(signatures):
                x= len(up_genes.intersection(set(signatures[sig]))) -1
                M= len(universe)+len(signatures[sig])
                n= len(signatures[sig])
                N= 500-x
                #print (x,M,n,N)
                up_mat[0,jj] =hypergeom.sf(x,M,n,N)
                #print hypergeom.sf(x,M,n,N)
        else:
            for jj, sig in enumerate(signatures):
                x=len(dn_genes.intersection(set(signatures[sig]))) -1
                M=len(universe)+len(signatures[sig])
                n=len(signatures[sig])
                N=500-x
                dn_mat[0,jj] = hypergeom.sf(x,M,n,N)
    r[0]=up_mat[0]
    r[1]=dn_mat[0]
    r[2]= ii,index
    return r
    pass

def p_sig_median(xx):

    dat_c= xx[0]
    s=xx[1]
    signature = xx[2]
    s_i=get_probe_index(signature,dat_c['rownames'],verbose=False)
    if len(s_i)>0:
        s_val=numpy.median(dat_c['data'][s_i],axis=0)
    else:
        s_val=numpy.array([0]*len(dat_c['colnames']))
        print s_val
    return  s_val

def p_sig_mean(xx):

    dat_c= xx[0]
    s=xx[1]
    signature = xx[2]
    s_i=get_probe_index(signature,dat_c['rownames'],verbose=False)
    if len(s_i)>0:
        s_val=numpy.mean(dat_c['data'][s_i],axis=0)
        print s_val
    else:
        s_val=numpy.array([0]*len(dat_c['colnames']))
        print s_val
    return  s_val

def compute_sig_median(dat,sig, array_type= None,vs_nl=False,probes_to_use = None,compute = 'median'):
    '''computes median signature for a list of signature'''
    from multiprocessing import Pool
    from itertools import izip, repeat
    import copy
    pool = Pool(processes=7)


    print 'Array type : '+ str(array_type)

    if array_type==None:
        dat_c  = copy.copy(convert_affy_2_genes_in_data(dat,compute_fc=0))
    else:
        if array_type == 'as_is':
            dat_c = copy.copy(dat)
            print "not converting probenames"
        else:
            dat_c = copy.copy(convert_affy_2_genes_in_data(dat, array_type=array_type, compute_fc=False))
            dat_c['rownames'] = numpy.array([i.upper() for i in dat_c['rownames']])

    if vs_nl:
        dat_c['data'] = dat['raw_fold_change'].T
        dat_c['colnames'] = dat['colnames'][probes_to_use]
    sig_data=[]
    sig_data_names=[]



    l = len(sig.keys()) # finish me here...
    data_for_paralell = izip(repeat(dat_c,l),sig.keys(),sig.values())
    if compute == 'median':
        results =  pool.map(p_sig_median, data_for_paralell)
    else:
        results =  pool.map(p_sig_mean, data_for_paralell)
    return results



    for s in signatures:
        s_i=get_probe_index(signatures[s],dat_c['rownames'])
        if len(s_i)>0:
            s_val=median(dat['data'][s_i],axis=0)
            sig_data.append(s_val)
            sig_data_names.append(s)
    dat['data'] = concatenate([dat['data'],sig_data],axis=0)
    dat['rownames'] = concatenate([dat['rownames'],sig_data_names],axis=0)







def score_patients_vs_signatures(data, signature, data_to_use, top_n = 500,run_in_parallel=True):
    AnnotationDbi=importr('AnnotationDbi')
    hg133=importr('hgu133plus2.db')
    import collections
    from collections import OrderedDict
    from scipy.stats  import hypergeom
    from numpy import concatenate, array, argwhere, ravel
    from multiprocessing import Pool
    from itertools import izip, repeat
    pool = Pool(processes=7)
    signatures = OrderedDict(read_signatures(signature))
    universe = set(concatenate(array(signatures.values())).ravel())
    p2go=AnnotationDbi.as_list(hg133.hgu133plus2SYMBOL)
    probe_to_gene = collections.defaultdict(list)
    gene_to_probe = collections.defaultdict(list)
    p2go_names=numpy.array(p2go.names)
    for i,names in enumerate(p2go_names):
        try:
            for go in  p2go[i]:
                go_term=str(go)
                probe_to_gene[names].append(go_term)
                gene_to_probe[go_term].append(names)
        except Exception,e:
            pass
    #concatenate([gene_to_probe.get(i) for i in a[a.keys()[0]]  if gene_to_probe.get(i)!=None])

    up_mat = numpy.zeros([len(data_to_use),len(signatures.keys())])
    dn_mat = numpy.zeros([len(data_to_use),len(signatures.keys())])
    fc = data['raw_fold_change'].T

    if run_in_parallel == False: #run in serie...
        for ii,index in enumerate(data_to_use):
            print str(float(ii) / float(len(data_to_use))*100.) + ' % done.'
            mini_fc = 3
            p2u = argwhere(abs(fc[:,index]) > mini_fc).ravel()
            while len(p2u) < top_n*2:
                mini_fc -= .2
                #print '.',
                p2u = argwhere(abs(fc[:,index]) > mini_fc).ravel()
            sorted_i = get_sorted_index(fc[p2u,index])
            dn_index = (p2u[sorted_i])[0:top_n]
            up_index = (p2u[sorted_i])[-top_n:]
            dn_genes = concatenate([probe_to_gene.get(i) for i in data['rownames'][dn_index]  if probe_to_gene.get(i)!=None])
            dn_genes = set(dn_genes[dn_genes != 'NA'])
            up_genes = concatenate([probe_to_gene.get(i) for i in data['rownames'][up_index]  if probe_to_gene.get(i)!=None])
            #print str(index)+' :\t'+'\t'.join(up_genes[0:5])
            up_genes = set(up_genes[up_genes != 'NA'])

            for k in [1,2]:
                if k ==1:
                    for jj, sig in enumerate(signatures):
                        x= len(up_genes.intersection(set(signatures[sig]))) -1
                        M= len(universe)+len(signatures[sig])
                        n= len(signatures[sig])
                        N= 500-x
                        #print (x,M,n,N)
                        up_mat[ii,jj] = hypergeom.sf(x,M,n,N)
                        #print hypergeom.sf(x,M,n,N)
                else:
                    for jj, sig in enumerate(signatures):
                        x=len(dn_genes.intersection(set(signatures[sig]))) -1
                        M=len(universe)+len(signatures[sig])
                        n=len(signatures[sig])
                        N=500-x
                        dn_mat[ii,jj] = hypergeom.sf(x,M,n,N)
        return dn_mat,up_mat
    else:
        print 'running in //'
        l = len(data_to_use)
        data_for_paralell = izip(data_to_use,range(l),\
                                    repeat(signatures,l), \
                                    fc[:,data_to_use].T,\
                                    repeat(data['rownames'],l),\
                                    repeat(data_to_use,l),\
                                    repeat(top_n,l)  ,\
                                    repeat(probe_to_gene,l))
        results =  pool.map(_pgeom_test, data_for_paralell)
        return array([ array(results[i])[1] for i in range(len(data_to_use )) ]), array([ array(results[i])[0] for i in range(len(data_to_use )) ])



        #hypergeom.sf(x,M,n,N),
        # x : number of selected items having the property of interest (number of genes in signature)
        # M : M is total number of objects, (universe length)
        # n : n is total number of Type I objects (length of signature)
        # N : total number of selected items

    pass

def _pgeom_test2(xx):
    #zip(data_to_use,range(len(data_to_use)),[signatures]*len(data_to_use), [fc]*len(data_to_use),[data]*len(data_to_use),[data_to_use]*len(data_to_use))
    from numpy import argwhere,ravel ,sort, unique
    from scipy.stats  import hypergeom, rankdata
    from numpy import concatenate, array, argwhere, ravel
    ii= xx[1]
    index = xx[0]
    signatures = xx[2]
    fc = xx[3]
    rownames=xx[4]
    data_to_use = xx[5]
    top_n = xx[6]
    probe_to_gene=xx[7]

    universe = set(concatenate(array(signatures.values())).ravel())
    print str(float(ii) / float(len(data_to_use))*100.) + ' % done.'
    up_mat = numpy.zeros([1,len(signatures.keys())])
    ranks = rankdata(fc)
    sorted_ranks = sort(ranks)[-100:]
    p2u = []
    for r in unique(sorted_ranks):
        p2u.append(argwhere(ranks == r).ravel())
    p2u = concatenate(p2u)
    up_genes = rownames[p2u]

    #print str(index)+' :\t'+'\t'.join(up_genes[0:5])
    up_genes = set(up_genes[up_genes != 'NA'])
    r=[[],[]]
    for jj, sig in enumerate(signatures):
        x= len(up_genes.intersection(set(signatures[sig]))) -1
        M= len(universe)+len(signatures[sig])
        n= len(signatures[sig])
        N= 500-x
        #print (x,M,n,N)
        up_mat[0,jj] =hypergeom.sf(x,M,n,N)
        #print hypergeom.sf(x,M,n,N)
    r[0]=up_mat[0]
    r[1]= ii,index
    return r
    pass


def score_data_vs_signatures(data, signature, data_to_use, top_n = 500, run_in_parallel=True):
    '''
    score signature in the dataset (hypergeometric test). rownames have to match the genes in the signatures...
    '''
    import collections
    from collections import OrderedDict
    from scipy.stats  import hypergeom
    from numpy import concatenate, array, argwhere, ravel
    from multiprocessing import Pool
    from itertools import izip, repeat
    pool = Pool(processes=7)
    signatures = OrderedDict(read_signatures(signature))
    universe = set(concatenate(array(signatures.values())).ravel())

    up_mat = numpy.zeros([len(data_to_use),len(signatures.keys())])
    fc = data['data']

    print 'running in //'
    l = len(data_to_use)
    data_for_paralell = izip(data_to_use,range(l),\
                                repeat(signatures,l), \
                                fc[:,data_to_use].T,\
                                repeat(data['rownames'],l),\
                                repeat(data_to_use,l),\
                                repeat(top_n,l)  ,\
                                repeat(None,l))
    results =  pool.map(_pgeom_test2, data_for_paralell)
    return array([ array(results[i])[1] for i in range(len(data_to_use )) ]), array([ array(results[i])[0] for i in range(len(data_to_use )) ])



        #hypergeom.sf(x,M,n,N),
        # x : number of selected items having the property of interest (number of genes in signature)
        # M : M is total number of objects, (universe length)
        # n : n is total number of Type I objects (length of signature)
        # N : total number of selected items

    pass


def reduce_cluster_number(data_, probes = None, samples = None, cluster_numbers=10, min_corr = .8):
    ''' take the data, look into the stemness varable, and try to reduce the cluster numbers,
        by looking at those clusters whose correlation is high.
        cluster number is the number of clusters to achive.
        min corr is the minimum correlation between all the clusters. is any pair of clusters has a correlation aboce that value,
        they will be merged.
    '''
    from numpy import unique, argwhere
    import copy
    data = copy.copy(data_)
    if probes == None:
        probes= arange(len(data['rownames']))
    if samples != None:
        data['data'] = data['data'][:,samples]
        data['colnames'] = data['colnames'][samples]
        data['stemness'] = data['stemness'][samples]
    ave = average_cells_by_stemness(data)
    c=fast_cor(ave,probes_to_use=probes)
    c[c==1.]
    c[c==1.] = -1

    while ( (len(unique(data['stemness'])) > cluster_numbers) or (min_corr > c.max() ) ) :
        #find max
        argwhere(c == c.max())
        # find names of the two most close clusters...
        c1,c2 = ave['colnames'][argwhere(c == c.max())[1]]
        data['stemness'][argwhere(data['colnames'] == c2).ravel()] = c1
        data['colnames'][argwhere(data['colnames'] == c2).ravel()] = c1
        #re-run  ad libitum;)
        ave = average_cells_by_stemness(data)
        c=fast_cor(ave,probes_to_use=probes)
        c[c==1.]
        c[c==1.] = -1
        print str(len(unique(data['stemness']))) + ' cluster left!'
    data['stemness'] = stemness_by_class(data)
    return data
    pass


def limma_multiclass(data, data_to_use=None, p =None, pval='adj.P.Val', fc = None, limit = None, one_against_all = False , return_raw_r = False):
    """docstring for limma_fc
    groups are determined by stemness!
    one_against_all option does what is says. ;)

    """
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import pandas as pa
    import scipy.stats  as stats
    import copy
    import gc
    classes=[]
    classes_names=[]

    if 'stemness' not in data.keys():
        print 'No stemness key in dataset!'
        return
    if data_to_use == None:
        data_to_use = range(len(data['colnames']))

    data['colnames'] = numpy.array([  i.strip('*').replace(' (P)','') for i in data['colnames'] ])
    colnames_orig = numpy.array([i for i in data['colnames']],dtype='str')
    stemness_orig = numpy.array([i for i in data['stemness']],dtype='i')
    for index,i in enumerate(data['stemness'][data_to_use]):
        if i not in classes:
            classes.append(i)
            classes_names.append(data['colnames'][data_to_use][index])

    print 'number of classes is %d'%len(classes)

    comparission_name = '"' + '", " '.join(classes_names) + '"'
    comparission_name = comparission_name.replace(' ','')

    print 'comparing %s'%comparission_name


    limma=importr('limma')
    base=importr('base')
    f = base.factor(robjects.StrVector(data['colnames'][data_to_use] ), levels=robjects.StrVector(classes_names))
    design =  robjects.Matrix(numpy.array([data['colnames'][data_to_use] == i for i in classes_names], dtype=int).T)
    design.colnames = robjects.StrVector(classes_names)

#   print design

    e_matrix = robjects.Matrix(data['data'][:,data_to_use])
    e_matrix.colnames = robjects.StrVector(data['colnames'][data_to_use])
    e_matrix.rownames = robjects.StrVector(data['rownames'])
#   return e_matrix
    fit = limma.lmFit(e_matrix, design)
    comparisson = []

    #this is some initialisation of the variable.
    data_copy = copy.deepcopy(data)

#do one vs all comp.
    if one_against_all:
        classes = numpy.unique(data_copy['colnames'][data_to_use])

        #print type(stemness_orig)
        probes_to_use = []
        for c in classes:
            data_copy['colnames'] = colnames_orig
            data_copy['stemness'] = stemness_orig
            c1 = numpy.argwhere(data['colnames'] == c).ravel()
            c2 = numpy.argwhere(data['colnames'] != c).ravel()
            data_copy['colnames'][c1] = '1'
            data_copy['colnames'][c2] = '9'
            data_copy['stemness'][c1] = '1'
            data_copy['stemness'][c2] = '9'
            probes_to_use.append( limma_multiclass(data_copy, data_to_use=data_to_use, p = p, pval=pval, fc = fc, limit = limit, one_against_all = False , return_raw_r = False) )
            gc.collect()
        data_copy['colnames'] = colnames_orig
        data_copy['stemness'] = stemness_orig


        if len(probes_to_use) ==0:
            return nuumpy.array([])

        else:
            print probes_to_use
            return numpy.unique(numpy.concatenate(  [list(i) for i in probes_to_use if len(i) >0]  ))
#           limma:
    else:
# do pairwise cmp.
        nb_comp =(len(classes_names) * (len(classes_names)-1) ) / 2
        contrast_matrix = numpy.zeros([nb_comp , len(classes_names)]).T
        iterator=0
        for i in range(len(classes_names)):
            for j in range(i,len(classes_names)):
#           for j in range(len(classes_names)): #USED FOR ALL PAIRWISE COMPARISSON, nb_comp = (len(classes_names) * (len(classes_names)-1) ) here
                if i != j:
                    comparisson.append(classes_names[i]+'-'+classes_names[j])
#                   print i,j
                    contrast_matrix[i,iterator] =1
                    contrast_matrix[j,iterator] = -1
                    iterator+=1
#   print iterator , nb_comp
    contrast_matrix = robjects.Matrix(contrast_matrix)
    contrast_matrix.colnames = robjects.StrVector(comparisson)
    contrast_matrix.rownames = robjects.StrVector(classes_names)

    #print contrast_matrix
    fit2 = limma.contrasts_fit(fit, contrast_matrix)
    corr_fit =  limma.eBayes(fit2,0.01)



    #print corr_fit
#   print data['rownames'].shape[0]
#   results = limma.topTable(corr_fit,  adjust="BH",number=data['rownames'].shape[0])

    probes_to_use = set([])
    for i in range(1,nb_comp+1):
        nb = data['rownames'].shape[0]
        results = limma.topTable(corr_fit, coef=i, adjust="fdr",number=nb ,sort="none")#adjust="BH",number=nb)
        #print limma.topTable(corr_fit, coef=0, adjust="fdr",number=10)#adjust="BH",number=nb)
        if return_raw_r:
            return pandas2ri.ri2py(results)
#       print results
        out =  pandas2ri.ri2py(results)
        out.index = numpy.arange(out.shape[0])
        #out['ID'] = numpy.array(robjects.r.rownames(results))
        #print out['ID']
        if limit is None:
                limit = out.shape[0]
        if fc == None:
            fc =  0
        if p == None:
            p = .05
        out['logFC'] = numpy.abs(out['logFC'])
        #print fc
        #print p
        new_probes = numpy.array(out[ numpy.logical_and ( out['logFC'] > fc , out[pval] < p )].sort_values(by=pval).index.values,dtype='i')[:limit]
        out = out[ numpy.logical_and ( out['logFC'] > fc , out[pval] < p )].sort_values(by=pval)
        print 'Found', out.shape[0], 'genes with p <', p , 'and FC >', fc
        print out
###
###
###        rn=base.rownames
###
####       return results
###        indexes_p=[]
###        indexes_fc=[]
###        ###a = [numpy.asarray(r.rx2(ii)) for ii in r.names]
###        a = [numpy.array(r)[ii] for ii in range(len(r.names))]
###        #a.append( numpy.array(get_probe_index( rn(r), data['rownames'])) )
###        a.append(  rn(r) )
###        #return rn(r)
###        cnames = numpy.array(r.names).tolist()
###        cnames.append('index')
###        #return a
###        try:
###            out =  pandas2ri.ri2pandas(r)
###        except Exception,e:
###            out =  pandas2ri.ri2py(r)
###        out['ID'] = numpy.array(robjects.r.rownames(r))
###        #out= numpy.recarray(len (a[0]) , dtype= [ (cnames[ii],j.dtype) for ii,j in enumerate(a)])
###        #for ii in range(len(cnames)):
###        #    out[cnames[ii]] = a[ii]
###        #del a
####       return out
###        #this is the threshold for significantly +- genes
###
###        if fc != None:
###            #find the most differentiallyly expressed genes (+ and -), sort them
###            #find those diff above threshold
###            indexes_fc_up = out[out['logFC'] > fc].sort('logFC',ascending=False)['ID']
###            indexes_fc_dn = out[out['logFC'] < -fc].sort('logFC',ascending=True)['ID'][:nb]
###            indexes_fc = numpy.concatenate([indexes_fc_up[:nb] ,indexes_fc_dn[:nb] ])
###            print 'Found %d genes with FC > %f in %s'%(len(indexes_fc), fc,comparisson[i-1])
###        else:
###            indexes_fc =  out['ID']
###
###        if p != None:
###            #p values come ordered... but we still recompute, in case.
###            # 'P.Value', 'adj.P.Val'
###            indexes_p =  out[out[pval] < p].sort(pval)['ID'][:nb]
###            print 'Found %d genes with p < %f in %s'%(len(indexes_p), p,comparisson[i-1])
###            if one_against_all:
###                print indexes_p
###
###        else:
###            indexes_p = out['ID']
###
###
###        indexes_p     = set(indexes_p)
###        indexes_fc    = set(indexes_fc)
###
###        #probes_to_use = probes_to_use.union( indexes_p.intersection(indexes_fc)     )
###        new_probes = indexes_p.intersection(indexes_fc)
        #print new_probes
        if 0:#one_against_all:
            probes_to_use = probes_to_use.union( new_probes )
        else:
            probes_to_use = probes_to_use.union(new_probes).difference(probes_to_use.intersection(new_probes))

    #sig_genes={}
    #sig_genes['indexes'] = get_probe_index(list(probes_to_use),data['rownames'])
#       sig_genes['p_val'] = out[pval][probes_to_use]

#   print len(sig_genes['fc'])

    if one_against_all:
        print corr_fit
    #probes_to_use = get_probe_index(list(probes_to_use),data['rownames'])
    probes_to_use = numpy.array(list(probes_to_use),dtype='i')
    return  probes_to_use #,sig_genes
    pass

#

def go_mouse_analysis(genes,html_out='example.html',use_go_slim=False):
    import mygene
    import goenrich
    import subprocess

    genes = [i.upper() for i in numpy.array(genes).tolist()]

    mg = mygene.MyGeneInfo()

    # build the ontology
    if use_go_slim:
        #O = goenrich.obo.ontology('/Users/nrapin/Dropbox/BRIC/People/Felicia/Single_cells_poolled/db/goslim_generic.obo')
        O = goenrich.obo.ontology('/Users/nrapin/Dropbox/BRIC/People/Felicia/Single_cells_poolled/db/goslim_metagenomics.obo')
    else:
        O = goenrich.obo.ontology('/Users/nrapin/Dropbox/BRIC/People/Felicia/Single_cells_poolled/db/go-basic.obo')
    # use all entrez geneid associations form gene2go as background
    # use goenrich.read.goa('db/gene_association.goa_ref_human.gz') for uniprot
    gene2go = goenrich.read.gene2go('/Users/nrapin/Dropbox/BRIC/People/Felicia/Single_cells_poolled/db/gene2go.gz')
    values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
    # propagate the background through the ontology
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(O, values, background_attribute)
    # extract some list of entries as example query
    query = [i['entrezgene'] for i in mg.querymany(genes ,  scopes='symbol', species='human')  if i.has_key('entrezgene')]
    print numpy.array(query)
     #gene2go['GeneID'].unique()[:20]
    # for additional export to graphviz just specify the gvfile argument
    # the show argument keeps the graph reasonably small
    df = goenrich.enrich.analyze(O, query, background_attribute,gvfile='test.dot')
    # generate html
    df[df.q < .05 ].dropna().dropna().sort_values(by='p').head(100).to_html(html_out)
    #make graph
    subprocess.check_call(['dot', '-Tpng', 'test.dot', '-o', html_out.replace('html','png')])
    subprocess.check_call(['rm', 'test.dot' ])
    return df
    pass



def hist2d(x,y,bins=50, interpolation='nearest'):
    """docstring for hist2d
        makes a heatmap of x,y scatter
    """
    import pylab
    heatmap, yedges, xedges = numpy.histogram2d(y, x, bins=bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    pylab.clf()
    pylab.imshow(numpy.flipud(numpy.log(heatmap)), extent=extent, aspect='auto',interpolation=interpolation)
    pylab.colorbar()
#   plt.show()
    pass

def select_probes_by_ANOVA(data, p=0.01, data_to_use=None,mode='p'):
    """docstring for select_probes_by_ANOVA
        looks into the stemness variable to generate the vectors and perform one way Ftest to check wether or not to include probe
    """
    from scipy import stats
    StatsR = importr('stats')

    if data_to_use==None:
        data_to_use=range(len(data['colnames']))
        data_to_use=numpy.array(data_to_use)

    probes_to_use=[]
    matrix= data['data'][:,data_to_use]
    print matrix.shape
    unique_types=[]
    for n in data['stemness']:
        if n not in unique_types:
            unique_types.append(n)
    #loop through probes
    pvals=numpy.zeros(matrix.shape[0])
    for i,probe in enumerate(range(matrix.shape[0])):
        #build array based on type
        cat={}
        for j,s in enumerate(data['stemness'][numpy.array(data_to_use)]):
            if cat.has_key(s):
                cat[s].append(matrix[i,j])
            else:
                cat[s]=[]
                cat[s].append(matrix[i,j])
        pvals[i]=stats.f_oneway(*cat.values())[1]

#   print pvals
    pvals_corrected=numpy.array(StatsR.p_adjust(pvals))


    for i,probe in enumerate(range(matrix.shape[0])):
        if (pvals_corrected[i] < p): #and (numpy.sort(matrix[i])[-2] > 7)  :
            probes_to_use.append(i)

#   print qvalues
    if mode=='p':
        return probes_to_use
    else:
        qval=importr('qvalue')
        qvalues=qval.qvalue(pvals_corrected)
        qvals= numpy.array(qvalues.rx2(3))
        p2u=[i for i,j in enumerate(qvals) if j < p ]
        return p2u, qvals[p2u]
    pass

def select_probes_by_variance_single_cell(all_data, var_thr=1.0 ,expr_thr=0.0, at_least=None , data_to_use =  None, disp = None):
    """docstring for select_probes_by_variance

    return probes_to_use.
    """



    if data_to_use==None:
        data_to_use=range(len(all_data['colnames']))



    probes_to_use=[]
    matrix= all_data['data'][:,data_to_use]
    nbsamples=matrix.shape[1]
    sd=numpy.std(matrix,axis=1)
    median_sd=numpy.median(sd)

    if at_least == None:
        for i,probe in enumerate(range(matrix.shape[0])):
            if sd[i] > var_thr*median_sd and numpy.abs( numpy.mean(matrix[i]) ) > expr_thr :
                probes_to_use.append(i)
                if disp:
                    print all_data['rownames'][i], sd[i]/median_sd
    else:
#       numpy.arange(matrix.shape[0] ) [numpy.std(matrix,axis=1) > numpy.std(data['data']) * 1]

        for i,probe in enumerate(range(matrix.shape[0])):
            sd[i] = numpy.std(matrix[probe][matrix[probe] != 0 ])
        median_sd=numpy.median(sd)
        for i,probe in enumerate(range(matrix.shape[0])):
            if sd[i] > var_thr*median_sd and numpy.abs( numpy.mean(matrix[i]) ) > expr_thr :
                probes_to_use.append(i)

    print 'Found '+str(len(probes_to_use))
    return probes_to_use
    pass


def select_probes_by_variance(all_data, var_thr=1.0 ,expr_thr=0.0, at_least=None , data_to_use =  None, disp = None):
    """docstring for select_probes_by_variance

    return probes_to_use.
    """



    if data_to_use==None:
        data_to_use=range(len(all_data['colnames']))



    probes_to_use=[]
    matrix= all_data['data'][:,data_to_use]
    nbsamples=matrix.shape[1]
    sd=numpy.std(matrix,axis=1)
    median_sd=numpy.median(sd)

    if at_least == None:
        for i,probe in enumerate(range(matrix.shape[0])):
            if sd[i] > var_thr*median_sd and numpy.abs( numpy.mean(matrix[i]) ) > expr_thr :
                probes_to_use.append(i)
                if disp:
                    print all_data['rownames'][i], sd[i]/median_sd
    else:
#       numpy.arange(matrix.shape[0] ) [numpy.std(matrix,axis=1) > numpy.std(data['data']) * 1]

        for i,probe in enumerate(range(matrix.shape[0])):
            if sd[i] > var_thr*median_sd and numpy.sort(matrix[i])[at_least] > expr_thr :
                probes_to_use.append(i)

    print 'Found '+str(len(probes_to_use))
    return probes_to_use
    pass


def find_neighbour_by_class(merged_data,T,E,sample,data_to_use=None,Variance_dims=3,disp=0,n=100):
    """docstring for find_neighbour
    returns a list of n, up and down regulated genes. from a merged sample.

    The query sample is projected onto a PCA space of normal cells.
    Normal Cells are first averaged

    the function takes care fo finding the neighbour samples.
    The input of a normal PCA is optionnal, as we rerun the PCA if the E matrix is not given.

    """

    if data_to_use==None: #if parameter not given as argument, then use everything
        data_to_use=range(len(merged_data['stemness']))
        print 'hehe'
    else:
        for key in [ 'platforms', 'cel_file','colnames', 'data', 'stemness']:
            merged_data[key]=merged_data[key][:,data_to_use]


#makes weihgted geometric distance
    #find where the sample is-
    sample_index=None
    for i,pl in enumerate( merged_data['stemness']):
        print pl
        if  merged_data['cel_file'][i]==sample :
            sample_index=i
            print 'found! -> '+ merged_data['cel_file'][i]
            break
        if  pl==str(19) or pl==str(13) or merged_data['cel_file'][i]==sample :
            sample_index=i
            print 'found! -> '+ merged_data['cel_file'][i]


    if sample_index==None:
        print 'Sample not found...'
        return -1

    #compute distance
    dist=[]
    for i,array in enumerate(merged_data['colnames']):
        sum_v=0
        for j in range(Variance_dims):
            sum_v+= numpy.square(T[sample_index][j]-T[i][j]) * E[j]
        dist.append(numpy.sqrt(sum_v))

    #sort data according to that.
    x=stats.rankdata(dist)-1
    ranks=x.tolist()
    x=numpy.array([ranks,range(len(ranks))]).transpose().tolist()
    order=[x.index(i) for i in sorted( x , key=lambda i:( int(i[0]),i[1] ) )]
    dist=numpy.array(dist)

    print 'Clostest samples to %s:' % sample
    print '\n'.join([ '%d\t%s\t%f'%(i, merged_data['colnames'][order][i],  dist[order][i]) for i in range(len(merged_data['colnames'])) if (i < 5) and (i != 0 )])

    liste_array=[]
    for i in [1,2,3,4]:
        liste_array.append([ int(order[i]) ])
    w=[]
    for i in [1,2,3,4]:
        w.append([ numpy.exp( -1.0* dist[order][i] /5) ])

#   constructs an expression profile based on neighbours and weitghts:
    constructed_profile=numpy.zeros([len(merged_data['data'][:,0]), len(w) ])
    for i,s in enumerate(liste_array):
        constructed_profile[:,i]=merged_data['data'][:,s].reshape(-1)*w[i]
    chimera_sample=numpy.sum(constructed_profile,axis=1) / numpy.sum(w)

    # computes fold change, raw:
    raw_fold_change=numpy.subtract(merged_data['data'][:,sample_index], chimera_sample)
    #filter low expressed genes..
    fold_change=numpy.zeros(raw_fold_change.shape)
    for i in range(len(fold_change)):
        if (chimera_sample[i] < 8.0) and (merged_data['data'][:,sample_index][i] < 8.0):
            fold_change[i]=0
        else:
            fold_change[i]=raw_fold_change[i]

    ranked = stats.rankdata(fold_change).tolist()

    x=numpy.array([ranked,range(len(ranked))]).transpose().tolist()
    v=sorted( x , key=lambda i:( int(i[0]),i[1] ) )
    y={}
    for i,item in enumerate(x):
        y[str(item[0])]=i
    order_probes=[y[str(i[0])] for i in v]

    #find top ranked and bottom ranked:
    top_100=order_probes[0:n]
    low_100= order_probes[-1*(n+1):-1]


    if disp:
        import pylab
        pylab.close()
        pylab.plot(numpy.sort(fold_change[top_100]))
        pylab.plot(numpy.sort(fold_change[low_100]))

    return merged_data['rownames'][top_100],merged_data['rownames'][low_100], raw_fold_change , chimera_sample
    pass



def KEGG_analysis(top_n,low_n,use_SOAP=False):
    """docstring for KEGG_analysis
        input:  gene_symbols: up_genes, down_genes

        for each gene in top/low n positions,
        find it in KEGG,
        then find pathway(s),
        get all genes in the pathway.

      *****************************************************
      * report as much as possible in terms of affyprobes *
      *****************************************************

    Two modes of oeration are available:
    SOAP    : gets all the data directly from KEGG in japan. (very slow.)
    pyKEGG  : get the data from locally stored files generated thanks to the KEGG.py module. (very fast)

        save picture of pathway(s) (optionnal an commented here.)
        find drugs associated with these pathway(s).  (optionnal an commented here.)
    """

    [gene2affy,affy2gene] = load_affy2gene()



    if use_SOAP:
        from ZSI.ServiceProxy import ServiceProxy
        converted_rownames_high=[]
        converted_rownames_low=[]

        for i in low_n:
            try:
                converted_rownames_low.append(affy2gene[i])
            except Exception,e:
                converted_rownames_low.append(i)
        for i in top_n:
            try:
                converted_rownames_high.append(affy2gene[i])
            except Exception,e:
                converted_rownames_high.append(i)

# some display that' not displayed...
    #   a=[[converted_rownames_low[i],converted_rownames_high[i]] for i in range(len(converted_rownames_high)) ]
    #   for i in a:
    #       print i


        #load SOAP KEGG WSDL
        service= ServiceProxy(wsdl='http://soap.genome.jp/KEGG.wsdl')#, tracefile=sys.stdout)
        # find gene symbol in kegg-
        all_kegg_genes_up=get_KEGG_genes(converted_rownames_high)
        all_kegg_genes_down=get_KEGG_genes(converted_rownames_low)
            # get all pathways:
        all_kegg_pathways_up=get_KEGG_pathways_from_genes(all_kegg_genes_up)
        all_kegg_pathways_down=get_KEGG_pathways_from_genes(all_kegg_genes_down)
        probe_list_by_pathway_up=built_probe_list_from_patways(all_kegg_pathways_up,gene2affy)
        probe_list_by_pathway_down=built_probe_list_from_patways(all_kegg_pathways_down,gene2affy)

        all_kegg_drugs=[]
#       for pathways in all_kegg_pathways:
#           try:
#               pathway=service.bget(string=pathways)
#               print pathway['return']
#               url=service.mark_pathway_by_objects(pathway_id=pathways, object_id_list=all_kegg_genes) #all_kegg_genes used to be [entry] before
#               download(str(url['return']))
#               drug=service.get_drugs_by_pathway(pathway_id=pathways)
#               if drug != []:
#                   print drug['return']
#                   if len(drug['return'])>10:
#                       print 'too much compound to be true...'
#                       drug['return']=[]
#                   for compound in drug['return']:
#                       if compound not in all_kegg_drugs:
#                           all_kegg_drugs.append(compound)
#           except Exception,e :
#               print e
#

    else: #use local stuffs from KEGG.py
        [kegg2affy,affy2kegg] = KEGG.generate_kegg2affy()
        [p2g,g2p] = KEGG.generate_kegg_gene2pathway()


        all_kegg_genes_up=[]
        for i in top_n:
            if affy2kegg.has_key(i):
                all_kegg_genes_up.append(affy2kegg[i])
            else:
                print '%s is not reported in KEGG..'%(i)

        all_kegg_genes_down=[]
        for i in low_n:
            if affy2kegg.has_key(i):
                all_kegg_genes_down.append(affy2kegg[i])
            else:
                print '%s is not reported in KEGG..'%(i)


        # get all pathways:
        all_kegg_pathways_up=[]
        for gene in all_kegg_genes_up:
            if g2p.has_key(gene):
                p_list=g2p[gene]
                for p in p_list:
                    all_kegg_pathways_up.append(p)
            else:
                print '%s is not part of a pathway..' %(gene)
        all_kegg_pathways_down=[]
        for gene in all_kegg_genes_down:
            if g2p.has_key(gene):
                p_list=g2p[gene]
                for p in p_list:
                    all_kegg_pathways_down.append(p)
            else:
                print '%s is not part of a pathway..' %(gene)


        probe_list_by_pathway_up={}
        for p in all_kegg_pathways_up:
            glist=[]
            for g in p2g[p]:
                if kegg2affy.has_key(g):
                    glist.append(kegg2affy[g])
                else:
                    print '%s not in KEGG / affy list' %(str(g))
            probe_list_by_pathway_up[p]=glist

        probe_list_by_pathway_down={}
        for p in all_kegg_pathways_down:
            glist=[]
            for g in p2g[p]:
                if kegg2affy.has_key(g):
                    glist.append(kegg2affy[g])
                else:
                    print '%s not in KEGG / affy list' %(str(g))
            probe_list_by_pathway_down[p]=glist

        all_kegg_drugs=[]


    print '--------------------> Done.'




    results={}
    results['genes']={}
    results['genes']['up']=all_kegg_genes_up
    results['genes']['down']=all_kegg_genes_down
    results['pathways']={}
#   results['pathways']['list']=all_kegg_pathways
    results['pathways']['up']=probe_list_by_pathway_up
    results['pathways']['down']=probe_list_by_pathway_down
    results['drugs']=all_kegg_drugs

    for d in results['drugs']:
        aaa=service.bget(string=d)
        for i in aaa['return'].splitlines():
            if i.find('NAME')!=-1:
                print i

#   a=g[0].split()[0]
#   c=service.get_pathways_by_genes(genes_id_list=[a])
#   pathway=c['return'][0]
#   service.bget(string=patway)
#   service.get_drugs_by_pathway(pathway_id='path:eco00020')
#   service.mark_pathway_by_objects(pathway_id='path:eco00970', object_id_list=obj_list)
#   download(...)

    return results
    pass


def load_affy2gene():
    """docstring for load_affy2gene"""

    AnnotationDbi=importr('AnnotationDbi')
    hg133=importr('hgu133plus2.db')
    import collections
    p2go=AnnotationDbi.as_list(hg133.hgu133plus2SYMBOL)
    probe_to_gene = collections.defaultdict(list)
    gene_to_probe = collections.defaultdict(list)
    p2go_names=numpy.array(p2go.names)
    for i,names in enumerate(p2go_names):
        try:
            for go in  p2go[i]:
                go_term=str(go)
                probe_to_gene[names].append(go_term)
                gene_to_probe[go_term].append(names)
        except Exception,e:
            pass
    return probe_to_gene,gene_to_probe,

    pass




def built_probe_list_from_patways(all_kegg_pathways,gene2affy):
    """docstring for built_gene_list_from_patways
    input : patways list from get_kegg_pathways
    out: {'pathways1: [affy_probe_1,affy_probe_2,...,affy_probe_n,'}
    """

    result={}
    service= ServiceProxy(wsdl='http://soap.genome.jp/KEGG.wsdl')#, tracefile=sys.stdout)
    for pathways in all_kegg_pathways:
        p_name= str(pathways.split(':')[-1])
        print 'Looking for affy probes in ' + p_name
        result[p_name]=[]
        try:
#           pathway=service.bget(string=pathways) #ONLY INFOS here.
            all_genes=service.get_genes_by_pathway(pathway_id=pathways)
            for gene in all_genes['return']:
                name=convert_kegg_gene(gene)
                if len(name) > 0:
                    for entry in name:
                        try:
                            affyname=gene2affy[entry]
                            result[p_name].append(affyname)
                            print affyname+' -> '+entry
                        except Exception, e:
                            print str(entry) + ' not on chip..'
#           get_pathway_picture(pathways,upgenes=all_kegg_genes_up,downgenes=all_kegg_genes_down)
        except Exception,e :
            print e
    return result
    pass


def convert_kegg_gene(gene):
    """docstring for convert_kegg_gene
        to index index in array.
    """
    refseq2gene, gene2refseq=load_ref_seq_to_gene_names()
    name=[]
    service= ServiceProxy(wsdl='http://soap.genome.jp/KEGG.wsdl')#, tracefile=sys.stdout)
    z=service.bconv(string=gene)
    text=str(z['return'])
    lines=text.splitlines()
    for line in lines:
        split=line.split('\t')
        kegg=split[0]
        other=split[1]
        equivalent=split[2]
        if other.split(':')[0]=='rs':
            try:
                name.append(refseq2gene[other.split(':')[-1]])
                print name
            except Exception, e:
                xxx=0
#               print str(other.split(':')[-1]) + ' not in db...'
    return name
    pass

def get_pathway_picture(path,upgenes=[],downgenes=[],color_up='red',color_down='blue'):
    """docstring for get_pathway_pucture
        get a KEGG path way picture.
        it is possible to provide a list of down /up regulated genes too.
    """

    service= ServiceProxy(wsdl='http://soap.genome.jp/KEGG.wsdl')#, tracefile=sys.stdout)
    if len(upgenes+downgenes) == 0:
        ulr_r=service.mark_pathway_by_objects(pathway_id=path, object_id_list=[])
    else:
        bg_list=[]
        fg_list= len(upgenes+downgenes)*['black'] #black text..
        for i in upgenes:
            bg_list.append(color_up)
        for i in downgenes:
            bg_list.append(color_down)
        obj_list=upgenes+downgenes
        ulr_r=service.color_pathway_by_objects(pathway_id=path,object_id_list=obj_list, fg_color_list= fg_list, bg_color_list= bg_list)


    download(ulr_r['return'])
    pass

def get_KEGG_genes(converted_rownames):
    """docstring for get_KEGG_genes
    from gene names.
    """
    service= ServiceProxy(wsdl='http://soap.genome.jp/KEGG.wsdl')#, tracefile=sys.stdout)
    #now the strategy is easy:
    #   for each gene in top/low 100,
    #   find it in KEGG,
    #   then find pathway(s),
    #   save picture of pathway(s)
    #   find drugs associated with these pathway(s).
    # find gene symbol in kegg-
    all_kegg_genes=[]
    for gene in converted_rownames:
        request=unicode('genes %s'%gene)
        print 'Looking for symbol for %s' % gene
        #finds all genes with a given name.
        r=service.bfind(string=request)
        #split text of resutls
        t=r['return'].split(';')
        #get all genes from human (the 'hsa' thing, means human in KEGG)
        g=[str(i) for i  in t if str(i).find('hsa')!=-1]
        #there might be a lot of them... so
#       print g[-1] #this is a human readable name
            #we need to filter out junk here..
        if len(g)> 10:
            print 'hum fishy query, trying with exact name...'
            g=[unicode('hsa:'+gene)] #trying with exact name.
        for entry in g:
            print entry.split()[-1] #this is a human readable name
            for item in entry.split():
                if item.find('hsa')!=-1:
                    if item  not in all_kegg_genes:
                        all_kegg_genes.append(item)
    return all_kegg_genes
    pass

def get_KEGG_pathways_from_genes(all_kegg_genes):
    """docstring for get_KEGG_pathways_from_genes"""
    service= ServiceProxy(wsdl='http://soap.genome.jp/KEGG.wsdl')
    # get all pathways:
    all_kegg_pathways=[]
    for entry in all_kegg_genes:
        pathway_Response=service.get_pathways_by_genes(genes_id_list=[entry])
        if pathway_Response is not []:
            print 'patways found for ' +entry
            for pathways in pathway_Response['return']:
                if pathways not in all_kegg_pathways:
                    all_kegg_pathways.append(pathways)
    return all_kegg_pathways
    pass


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def make_gene_plots_UMIs( genes_list, dat,prefix = 'genes_',order1=None, col_wrap=6,sharey=False):
    if order1 == None:
        order1 = numpy.unique(dat['colnames'])
    #split so that we see  50 genes / page..
    from numpy import *
    from pylab import close, savefig
    import seaborn as sns
    nsplits = numpy.array(genes_list).shape[0] /30
    if nsplits ==0:
        nsplits =1
    for iiii, iii in enumerate(array_split(genes_list,nsplits)):
        selected_genes =pa.DataFrame()
        #for gene in  concatenate(differenciation_genes.values()):
        for gene in  iii:
            print gene
            index = get_probe_index([gene],dat['rownames']).ravel()
            if len(index) != 0:
                index = index[0]
                aa = pa.DataFrame([dat['colnames'],dat['data'][index,:]],index=['cluster',gene]).T
                #aa = pa.DataFrame([cnames,dat['data'][index,:]],index=['cluster',gene]).T
                #values = [aa[aa.cluster == i][gene].astype('f').values  for i in ['l','m','e','12'] ]
                values = [aa[aa.cluster == i][gene].astype('f').values  for i in array(arange(12)+1,dtype='str') ]
                aa.columns=['cluster','val']
                aa['gene'] = [gene]*aa.shape[0]
                if selected_genes.shape == (0,0):
                    selected_genes = aa
                else:
                    selected_genes = selected_genes.append(aa)
        selected_genes['val'] =  selected_genes['val'].astype('f')
        #add the percentage of expressed genes:
        mean_df = []
        selected_genes['height'] = 0
        for gene in selected_genes.gene.unique():
            for c in selected_genes.cluster.unique():
                df = selected_genes[numpy.logical_and(selected_genes.gene == gene,selected_genes.cluster == c )].val
                not_zeros =  sum(df != 0)
                zeros =  sum(df == 0)
                mean_val = df.mean() #df[df != 0].mean()
                #print gene , zeros, not_zeros, not_zeros*1./(not_zeros+zeros),mean_val
                mean_df.append( (gene , c,zeros, not_zeros, not_zeros*1./(not_zeros+zeros) , mean_val ,not_zeros*1./(not_zeros+zeros) *     mean_val))
                selected_genes.iloc[argwhere(numpy.logical_and(selected_genes.gene == gene,selected_genes.cluster == c )).ravel(),3] =  not_zeros*1./ (not_zeros+zeros) * mean_val
                print not_zeros*1./ (not_zeros+zeros) * mean_val
        mean_df  = pa.DataFrame(mean_df,columns=['gene' , 'cluster','zeros', 'not_zeros', 'pc_not_null' , 'mean', 'height'])

        #plot all in one go...
        close()
        #g = sns.FacetGrid(selected_genes[selected_genes.val != 0], col="gene", col_wrap=6, size=2, hue="cluster",sharey=False)
        sns.set(font_scale=.55)
        g = sns.FacetGrid(selected_genes, col="gene", col_wrap=col_wrap, size=2, hue="cluster",sharey=sharey)
        ####g.map(sns.pointplot, "solutions", "score", color=".3", ci=None);
        g.map( sns.barplot ,  'cluster','val', order=order1, ci = None ,saturation=0.3)
        ###g.map( sns.boxplot ,  'cluster','val', order=order1 ,saturation=0.3)
    ##g.map( sns.swarmplot , 'cluster','val',order=order1, edgecolor='gray',color="black",size = 7,marker='o')
        ####g.map(sns.pointplot, "solutions", "score", color=".3", ci=None);
        g.map( sns.barplot ,  'cluster','height', order=order1, ci = None, hatch='////' , saturation=0.999)
        g.fig.subplots_adjust( wspace=.2, hspace=.2);
        g.set_xticklabels(rotation=45)
        savefig(prefix+ str(iiii) +'.pdf')
        close()
        return selected_genes
        ## biggger_clusters = {'e':[1,2,3,4],
        ##         'l':[5,6,7,8],
        ##         'm':[9,10,11]
        ##         }
        ##
        ## biggger_clusters = {
        ##     'EM':[2,9],
        ##     'EL':[1,3,5],
        ##      'E':[4],
        ##      'M':[11],
        ##      'LM':[6,8,10],
        ##      'L':[7]
        ##         }
        ## biggger_clusters = {
        ##       'HSC':[18],
        ##       'E':[8,7,6,5,4,3],
        ##       'Em':[12,2],
        ##       'EL':[1,20,9],
        ##       'L':[11],
        ##       'eLM':[10,19,16],
        ##       'Me':[15,13,14],
        ##       'M':[17],
        ##          }
##
        ## for i in biggger_clusters:
        ##     for cl in biggger_clusters[i]:
        ##         indexes = argwhere(selected_genes.cluster == str(cl)).ravel()
        ##         selected_genes.cluster.iloc[indexes] = i
        ##
        ## mean_df = []
        ## selected_genes['height'] = 0
        ## for gene in selected_genes.gene.unique():
        ##     for c in selected_genes.cluster.unique():
        ##         df = selected_genes[numpy.logical_and(selected_genes.gene == gene,selected_genes.cluster == c )].val
        ##         not_zeros =  sum(df != 0)
        ##         zeros =  sum(df == 0)
        ##         mean_val = df.mean()
        ##         #print gene , zeros, not_zeros, not_zeros*1./(not_zeros+zeros),mean_val
        ##         mean_df.append( (gene , c,zeros, not_zeros, not_zeros*1./(not_zeros+zeros) , mean_val ,not_zeros*1./(not_zeros+zeros) *     mean_val))
        ##         selected_genes.iloc[argwhere(numpy.logical_and(selected_genes.gene == gene,selected_genes.cluster == c )).ravel(),3] =  not_zeros*1./ (## not_zeros+zeros) * mean_val
        ##         print not_zeros*1./ (not_zeros+zeros) * mean_val
        ## mean_df  = pa.DataFrame(mean_df,columns=['gene' , 'cluster','zeros', 'not_zeros', 'pc_not_null' , 'mean', 'height'])
        ##
        ## cluster_order = ['e','l','m','12']
        ## cluster_order = ['E', 'M', 'L', 'EL', 'EM', 'LM' ,'12']
        ## close()
        ## g = sns.FacetGrid(selected_genes, col="gene", col_wrap=6, size=2, hue="cluster",sharey=False)
        ## g.map( sns.barplot ,  'cluster','val', order=cluster_order, ci = None ,saturation=0.3)
        ## #g.map( sns.swarmplot , 'cluster','val',order=['l','m','e','12'], edgecolor='gray',color="black",size = 1.5,marker='o')
        ## g.map( sns.barplot ,  'cluster','height', order=cluster_order, ci = None, hatch='////' , saturation=0.999)
        ## g.fig.subplots_adjust( wspace=.2, hspace=.2);
        ## #savefig('sigs_genes_big_clusters.pdf')
        ## savefig(prefix+'big_clusters_'+ str(iiii) +'.pdf')
        ## close()

    pass


def plot_clustered_stacked(dfall, labels=None, title="multiple stacked bar plot",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot.
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    import pandas as pd
    import matplotlib.cm as cm
    import numpy as np
    import matplotlib.pyplot as plt

    n_df = len(dfall)
    n_col = len(dfall[0].columns)
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)
    colors=['#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#faebd7','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' ,'#faebd7', '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#faebd7','#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']

    for index,df in enumerate(dfall) : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      color = colors[index],
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * (i / n_col))
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]
    for index,i in enumerate(range(n_df)):
        n.append(axe.bar(0, 0, color=colors[index], hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1])
    axe.add_artist(l1)
    return axe


def pickle_save_gz(var,out_file='file.pkl'):
    """docstring for pickle_save
        Pickle a variable, saves as var_name.pkl
    """
    import cPickle, gzip
    print 'Saving '+out_file
    with gzip.open(out_file, 'wb') as f:
        cPickle.dump(var, f,-1)
    pass

def pickle_load_gz(in_file, inplace=False):
    """docstring for pickle_load"""
    import cPickle, marshal, gzip
    x=None
    if inplace:
        if os.path.exists(in_file):
            with gzip.open(in_file, 'rb') as f:
                loaded_object = cPickle.load(f)
            return loaded_object
        else:
            print 'File %s could not be found ! ' % (in_file)
    else:
        if '/' in in_file:
            filename = in_file.split('/')[-1]
            path = '/'.join(in_file.split('/')[:-1])
        else:
            filename = in_file
            path = ''
        run_sh('cp %s %s/temp_%s'%(in_file, path,filename) ,text='')
        run_sh('gunzip %s/temp_%s'%(path,filename) ,text='')
        x=pickle_load('%s/temp_%s'%(path,filename.replace('.gz','') ) )
        run_sh('rm -f %s/temp_%s'%(path,filename.replace('.gz','')) ,text='')
        run_sh('rm -f %s/temp_%s'%(path,filename) ,text='')
    return x
    pass

def run_sh(command,text='Im running a command'):
    """docstring for run_sh
        runs a command with ps.Popen, wait until sthis is done.
    """
    import subprocess as sp
    print command.split(' ')
    p = sp.Popen(command.split(' '),stdout = sp.PIPE,  stdin=sp.PIPE,  stderr=sp.PIPE)

    output, err = p.communicate() #now wait
    print output
    return output,err
    pass

def pickle_save(var,out_file='file.pkl'):
    """docstring for pickle_save
        Pickle a variable, saves as var_name.pkl
    """
    import cPickle
    output = open(out_file, 'wb')
    print 'Saving '+out_file
    cPickle.dump(var, output,-1)
    output.close()

    pass

def pickle_load(in_file):
    """docstring for pickle_load"""
    import cPickle
    x=None
    if os.path.exists(in_file):
        pkl_file = open(in_file, 'rb')
        x=cPickle.load(pkl_file)
        pkl_file.close()
    else:
        print 'File %s could not be found ! ' % (in_file)

    return x
    pass

def multisplit(string,separators, mode='str'):
    """docstring for multisplit
        string.split, but with more than one separator.
        separators as a string (mode =str) ',\t;'
                      a list of strings [',' , 'www' ]
    """
    if mode == 'str':
        seps=[]
        for i in range(len(separators)):
            seps.append(separators[i])
    else:
        seps=separators

    splited_s=[string]
    for sep in seps:
        for index,s in enumerate(splited_s):
#           print s, index
            splited=s.split(sep)
            splited_s.pop(index)
            for s_i in splited:
                splited_s.insert(index,s_i)
    splited_s.reverse()
    return splited_s
    pass

def get_probe_index(names,rownames,verbose=True):
    """docstring for get_probe_index
        get smth like [209302_at, 32625_at],
        and produces their index in rownames.
    """
    pdic={}
    for i,probe in enumerate(rownames):
        pdic[probe]=i

    indexes=[]
    for i in names:
        if pdic.has_key(i):
            indexes.append(pdic[i])
        else:
            if verbose:
                print 'Missing probe %s in list...' %(i)
    return numpy.array(indexes)
    pass

def download(url):
    """Copy the contents of a file from a given URL
    to a local file.
    """
    import urllib
    webFile = urllib.urlopen(url)
    localFile = open(url.split('/')[-1], 'w')
    localFile.write(webFile.read())
    webFile.close()
    localFile.close()



def make_heatmap(data,colnames=None,rownames=None,figure='correlation.pdf'):
    """docstring for make_heatmap
    gets a matrix/nd array,
    and export a pdf.
    """
    dims=data.shape
    if colnames==None:
        colnames=range(dims[0])
    if rownames==None:
        rownames=range(dims[1])
    TRUE=1
    FALSE=0
    # R Libraries to load:
    packages=['affy','hgu133plus2.db','hgu133a.db','hgu133b.db',  'gplots','genefilter' , 'Biobase',  'squash',  'amap']
    affy, hgu133plus2, hgu133a, hgu133b, gplots, genefilter, Biobase , squash, amap = import_rlib(packages)


    sm=rebuild_R_object(data=data, colnames=colnames ,rownames=rownames)

    grdevices = importr('grDevices')
    grdevices.pdf(figure)
    nul = robjects.r('NULL')
    scaling=robjects.r('\"none\"') # , scale=scaling
    dendrogramme = robjects.r('\"row\"') #c("both","row","column","none"),
#   fclust = as_dendrogramme(amap.hcluster(sm, method = "pea",link = "ward") )
#   fdist = robjects.r.cor()
    gplots.heatmap_2(sm, scale=scaling, dendrogram=dendrogramme , trace="none", col= gplots.bluered(20),  symm=TRUE, margins = robjects.r('c(10,10)'), cexCol = 0.5,cexRow = 0.5, main = 'Correlation of sample expression data')
#   gplots.heatmap_2(sm,  col= gplots.bluered(20),   trace="none" ,margins = robjects.r('c(10,10)'), cexCol = 0.5,cexRow = 0.5, main = 'Correlation of sample expression data')
    grdevices.dev_off()
    pass

#
def fast_cor(data,probes_to_use=None, data_to_use  = None ,figure='correlation.pdf',dist="s",plot=False):
    """docstring for compute_cor
    computes correlation,
    and export a pdf.
    """
    TRUE=1
    FALSE=0
    if probes_to_use == None:
        probes_to_use = range(data['data'].shape[0])
    if data_to_use == None:
        data_to_use = range(data['data'].shape[1])


    # R Libraries to load:
#   packages=['affy','hgu133plus2.db','hgu133a.db','hgu133b.db',  'gplots','genefilter' , 'Biobase',  'squash',  'amap']
#   affy, hgu133plus2, hgu133a, hgu133b, gplots, genefilter, Biobase ,  squash, amap = import_rlib(packages)

    if plot:
    	gplots = importr('gplots')
        base =importr('base')
    nb_s = len(data_to_use)
    similarity_matrix=numpy.zeros( [nb_s, nb_s] )

    print 'Computing Correlation:'
    distance_matrix = Bio.Cluster.distancematrix(data['data'][probes_to_use,:][:,data_to_use].T , dist=dist)

    #reconstruct the similarity matrix out of this triangular matrix for R to be happy.

    for i in range(nb_s):
        similarity_matrix[i,i] = 1.
        for j in range(i):
            similarity_matrix[i,j]=similarity_matrix[j,i] = 1. - distance_matrix[i][j]


    if plot:
        #return similarity_matrix
    	sm= robjects.Matrix(similarity_matrix) #rebuild_R_object(data=similarity_matrix, colnames= data['colnames'][data_to_use],rownames=data['colnames'][data_to_use])
        sm.colnames =  robjects.StrVector(data['colnames'][data_to_use])
        sm.rownames =  robjects.StrVector(data['colnames'][data_to_use])
        grdevices = importr('grDevices')
    	grdevices.pdf(figure)
    	nul = robjects.r('NULL')
    	scaling=robjects.r('\"none\"') # , scale=scaling
    	dendrogramme = robjects.r('\"row\"') #c("both","row","column","none"),
        cat=numpy.unique(data['colnames'][data_to_use]).tolist()
        c=['#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF','#800000','#FF0000','#800080','#FF00FF','#008000','#00FF00','#808000','#FFFF00','#000080','#0000FF','#008080','#00FFFF' ]
        colors= robjects.StrVector([c[cat.index(i)] for i in  data['colnames'][data_to_use]    ] )
#   fclust = as_dendrogramme(amap.hcluster(sm, method = "pea",link = "ward") )
#   fdist = robjects.r.cor()
    	gplots.heatmap_2(sm, scale=scaling, dendrogram=dendrogramme , trace="none", col= gplots.bluered(20),  symm=TRUE, margins = robjects.r('c(10,10)'), cexCol = 0.5,cexRow = 0.5, main = 'Correlation of sample expression data', ColSideColors = colors )
#   gplots.heatmap_2(sm,  col= gplots.bluered(20),   trace="none" ,margins = robjects.r('c(10,10)'), cexCol = 0.5,cexRow = 0.5, main = 'Correlation of sample expression data')
    	grdevices.dev_off()

    return similarity_matrix


    pass






def compute_cor(data,data_to_use,probes_to_use,figure='correlation.pdf'):
    """docstring for compute_cor
    computes correlation,
    and export a pdf.
    """
    TRUE=1
    FALSE=0
    # R Libraries to load:
    packages=['affy','hgu133plus2.db','hgu133a.db','hgu133b.db',  'gplots','genefilter' , 'Biobase', 'squash',  'amap']
    affy, hgu133plus2, hgu133a, hgu133b, gplots, genefilter, Biobase , bias, squash, amap = import_rlib(packages)

    nb_s = len(data['colnames'][data_to_use])
    similarity_matrix=numpy.zeros( [nb_s, nb_s] )

    print 'Computing Correlation:'
    prog = ProgressBar(0, len(data_to_use) , 50, mode='fixed')
#   Bio.Cluster.distancematrix(merged_data['data'].transpose(),dist='s')
    for i,index_i in enumerate(data_to_use):
        for j,index_j in enumerate(data_to_use):
            if similarity_matrix[j,i]==0 and similarity_matrix[i,j]==0 :
#               similarity_matrix[i,j]= (1 - hcluster.correlation(data['data'][probes_to_use,index_i],data['data'][probes_to_use,index_j]) )
                similarity_matrix[i,j]= (1 - Bio.Cluster.distancematrix(( ( data['data'][probes_to_use,index_i] ) ,  (data['data'][probes_to_use,index_j])  ) , dist="c")[1][0])
#               similarity_matrix[i,j]=  cosine_dist(data['data'][probes_to_use,index_i],data['data'][probes_to_use,index_j] )
                similarity_matrix[j,i]=similarity_matrix[i,j]
#               print data['data'][probes_to_use,index_i].shape
            else:
                continue
                continue

    #progress bar code..
        oldprog = str(prog)
        prog.update_amount(i)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)
    print

    sm=rebuild_R_object(data=similarity_matrix, colnames= data['colnames'][data_to_use],rownames=data['colnames'][data_to_use])
    print sm
    grdevices = importr('grDevices')
    grdevices.pdf(figure)
    nul = robjects.r('NULL')
    scaling=robjects.r('\"none\"') # , scale=scaling
    dendrogramme = robjects.r('\"row\"') #c("both","row","column","none"),
#   fclust = as_dendrogramme(amap.hcluster(sm, method = "pea",link = "ward") )
#   fdist = robjects.r.cor()
    gplots.heatmap_2(sm, scale=scaling, dendrogram=dendrogramme , trace="none", col= gplots.bluered(20),  symm=TRUE, margins = robjects.r('c(10,10)'), cexCol = 0.5,cexRow = 0.5, main = 'Correlation of sample expression data')
#   gplots.heatmap_2(sm,  col= gplots.bluered(20),   trace="none" ,margins = robjects.r('c(10,10)'), cexCol = 0.5,cexRow = 0.5, main = 'Correlation of sample expression data')
    grdevices.dev_off()
    pass

def run_cv(data, data_to_use = None, probes_to_use = None, nchunks=None,variance_filter=True):
    """docstring for run_cv
        the data
    """
    from mvpa2.suite import Dataset, LinearCSVMC, MappedClassifier, CrossValidation, SVDMapper, NFoldPartitioner

    if probes_to_use== None:
        print 'Using all features ...'
        probes_to_use= range(len(data['rownames']))
        if variance_filter:
            print 'Variance filtering..'
            probes_to_use = select_probes_by_variance(data, var_thr=4)
            print 'Using ', str(len(probes_to_use)), ' features.'
    if data_to_use== None:
        print 'Using all Examples ...'
        data_to_use= range(len(data['colnames']))

    ds = Dataset(data['data'][probes_to_use,:][:,data_to_use].T)
    ds.sa['colnames'] = data['colnames'][data_to_use]
    ds.sa['targets'] = ds.sa.colnames

    #divide into chuncks for cross valitaiton. #easy and stupid..

    if nchunks == None:
        nchunks = len(ds.sa['colnames']) # LOO cross validation

    chunks = []

    for i,j in enumerate(ds.sa['colnames']):
        chunks.append(i%nchunks)

    ds.sa['chunks'] = chunks

    print
    # a bit of display.
    print 'Found %d classes:\n'%(ds.sa['colnames'].unique.shape[0]) + '\n'.join(ds.sa['colnames'].unique)

    baseclf = LinearCSVMC()
#   metaclf = MappedClassifier(baseclf, SVDMapper())
    metaclf = baseclf
    cvte = CrossValidation(metaclf, NFoldPartitioner( cvtype = 3),  enable_ca=['stats'])
    cv_results = cvte(ds)
    print
    print '---------------------'
    print 'Error is : %.3f ' %(numpy.mean(cv_results)*100.), '% |'
    print '---------------------'

    print
    print '---------------------'
    print numpy.round(cvte.ca.stats.stats['ACC%'], 1)
    print '---------------------'

    return cvte,metaclf
    pass


#
def run_classifier_sweep(data, data_to_use = None, probes_to_use = None, nchunks=None,variance_filter=True):
    """docstring for run_cv
        the data
    """
    from mvpa2.suite import *

    if probes_to_use== None:
        print 'Using all features ...'
        probes_to_use= range(len(data['rownames']))
        if variance_filter:
            print 'Variance filtering..'
            probes_to_use = select_probes_by_variance(data, var_thr=4)
            print 'Using ', str(len(probes_to_use)), ' features.'
    if data_to_use== None:
        print 'Using all Examples ...'
        data_to_use= range(len(data['colnames']))

    ds = Dataset(data['data'][probes_to_use,:][:,data_to_use].T)
    ds.sa['colnames'] = data['colnames'][data_to_use]
    ds.sa['targets'] = ds.sa.colnames

    #divide into chuncks for cross valitaiton. #easy and stupid..

    if nchunks == None:
        nchunks = len(ds.sa['colnames']) # LOO cross validation

    chunks = []

    for i,j in enumerate(ds.sa['colnames']):
        chunks.append(i%nchunks)

    ds.sa['chunks'] = chunks

    print
    # a bit of display.
    print 'Found %d classes:\n'%(ds.sa['colnames'].unique.shape[0]) + '\n'.join(ds.sa['colnames'].unique)


    #disable warnings:
    old_warnings = numpy.geterr()
    numpy.seterr(all='ignore')
    for (dataset, datasetdescr), clfs_ in \
        [
        ((ds,
          "AML Dataset"),
          clfswh[:] ),#clfswh['multiclass']),
        ]:
        # XXX put back whenever there is summary() again
        #print "%s\n %s" % (datasetdescr, dataset.summary(idhash=False))
        print " Classifier on %s\n" \
                "                                          :   %%corr   " \
                "#features\t train  predict full" % datasetdescr
        for clf in clfs_:
            print "  %-40s: "  % clf.descr,
            # Let's prevent failures of the entire script if some particular
            # classifier is not appropriate for the data
            try:
                # Change to False if you want to use CrossValidation
                # helper, instead of going through splits manually to
                # track training/prediction time of the classifiers
                do_explicit_splitting = True
                if not do_explicit_splitting:
                    cv = CrossValidation(
                        clf, NFoldPartitioner(), enable_ca=['stats', 'calling_time'])
                    error = cv(dataset)
                    # print cv.ca.stats
                    print "%5.1f%%      -    \t   -       -    %.2fs" \
                          % (cv.ca.stats.percent_correct, cv.ca.calling_time)
                    continue

                # To report transfer error (and possibly some other metrics)
                confusion = ConfusionMatrix()
                times = []
                nf = []
                t0 = time.time()
                #TODO clf.ca.enable('nfeatures')
                partitioner = NFoldPartitioner()
                for nfold, ds in enumerate(partitioner.generate(dataset)):
                    (training_ds, validation_ds) = tuple(
                        Splitter(attr=partitioner.space).generate(ds))
                    clf.train(training_ds)
                    #TODO nf.append(clf.ca.nfeatures)
                    predictions = clf.predict(validation_ds.samples)
                    confusion.add(validation_ds.targets, predictions)
                    times.append([clf.ca.training_time, clf.ca.predicting_time])

                tfull = time.time() - t0
                times = np.mean(times, axis=0)
                #TODO nf = np.mean(nf)
                # print confusion
                #TODO print "%5.1f%%   %-4d\t %.2fs  %.2fs   %.2fs" % \
                print "%5.1f%%       -   \t %.2fs  %.2fs   %.2fs" % \
                      (confusion.percent_correct, times[0], times[1], tfull)
                #TODO      (confusion.percent_correct, nf, times[0], times[1], tfull)
            except LearnerError, e:
                print " skipped due to '%s'" % e



    pass

def select_best_clf(dataset_, clfs):

    """Select best model according to CVTE

    Helper function which we will use twice -- once for proper nested
    cross-validation, and once to see how big an optimistic bias due
    to model selection could be if we simply provide an entire dataset.

    Parameters
    ----------
    dataset_ : Dataset
    clfs : list of Classifiers
      Which classifiers to explore

    Returns
    -------
    best_clf, best_error
    """
    best_error = None
    for clf in clfs:
        cv = CrossValidation(clf, NFoldPartitioner())
        # unfortunately we don't have ability to reassign clf atm
        # cv.transerror.clf = clf
        try:
            error = np.mean(cv(dataset_))
        except LearnerError, e:
            # skip the classifier if data was not appropriate and it
            # failed to learn/predict at all
            continue
        if best_error is None or error < best_error:
            best_clf = clf
            best_error = error
        verbose(4, "Classifier %s cv error=%.2f" % (clf.descr, error))
    verbose(3, "Selected the best out of %i classifiers %s with error %.2f"
            % (len(clfs), best_clf.descr, best_error))
    return best_clf, best_error

def average_cells_by_stemness(merged_data,exp2=False, disp=True):
    """docstring for average_cells_by_stemness"""
    #create averaged data from all categories of cell-
    averaged_data = {}
    celltypes=[]
    for i in merged_data['stemness']:
#       i=str(i)
        if i not in celltypes:
            celltypes.append(i)
    if disp:
        print celltypes
    for key in [ 'colnames', 'stemness']:
        averaged_data[key]=[]
    averaged_data['data']=numpy.zeros([merged_data['data'].shape[0],len(celltypes)])
    averaged_data['std']=numpy.zeros([merged_data['data'].shape[0],len(celltypes)])
    for index,celltype in enumerate(celltypes):
        try:
            samples = numpy.argwhere(merged_data['stemness'] == celltype).reshape(-1)
        except:
            samples = [i for i,j in enumerate(merged_data['stemness']) if j==celltype ]
            print 'fishy fishy'
            if disp:
                 print str(merged_data['colnames'][samples[0]]),len(samples) , '\t',samples
        for key in [ 'colnames', 'stemness']:
            averaged_data[key].append(merged_data[key][samples][0])
        if exp2:
            averaged_data['data'][:,index]=numpy.median(numpy.exp2(merged_data['data'][:,samples]), axis=1)
            averaged_data['std'][:,index]=numpy.std(numpy.exp2(merged_data['data'][:,samples]), axis=1)
        else:
            averaged_data['data'][:,index]=numpy.median(merged_data['data'][:,samples], axis=1)
            averaged_data['std'][:,index]=numpy.std(merged_data['data'][:,samples], axis=1)
#       print averaged_data['data'][:,index].shape
    for key in [ 'colnames', 'stemness']:
        averaged_data[key]=numpy.array(averaged_data[key])
    #done
    #print '\n'.join([ '(%f , %s)'%(numpy.median(averaged_data['data'][:,i]), averaged_data['colnames'][i]) for i in range(len(averaged_data['colnames'])) ])

    return averaged_data
    pass

def gcRMA(_file,Saved_data=False,Bias_correct=False, celfilepath='/Volumes/B52/Converted/', platforms=['HG-U133A','HG-U133B','HG-U133_Plus_2']):
    """docstring for RMA
    run RMA on samples.
    """
    print 'Je suis bugee!!! check how data is loaded!!'
    #load RMA from R as txt file...



    TRUE=1
    FALSE=0
    packages=['affy','hgu133plus2.db','hgu133a.db','hgu133b.db',  'gplots','genefilter' , 'Biobase', 'bias' , 'squash',  'amap', 'gcrma']
    affy, hgu133plus2, hgu133a, hgu133b, gplots, genefilter, Biobase , bias, squash, amap ,gcrma = import_rlib(packages)

    f=open(_file)
    t=f.read()
    f.close()
    data=[]
    i=0
    for sample in t.splitlines():
        #skip first line
        if i==0:
            i+=1
        else:
            data.append( sample.split('\t') )
            for j,field in enumerate(data[-1]):
                data[-1][j] = field.strip() #clears white spaces before and after..
            i+=1

    #separate the data according to the platform and normalize them separatelly:
    #possible platform are  : HG-U133A, HG-U133B, HG-U133_Plus_2 (all probes in A+B = probes in Plus_2)
    index_by_platform={}
    norm_data={}
#   platforms = ['HG-U133A','HG-U133B','HG-U133_Plus_2','HG-U133_Plus_2_new' ]
#   platforms = ['HG-U133A','HG-U133B','HG-U133_Plus_2' ]
    for platform in platforms:
        print 'Processing samples from %s' % platform
        index_by_platform[platform] = [i for i,j in enumerate(data) if data[i][7]==platform ]

        #Run RMA Normalization from R.
        norm_data[platform]={}
        #prepare object for justRMA.
        CEL_list=[]
        CEL_list_name=[]
        for sample_index in index_by_platform[platform]:
            CEL_list.append( celfilepath +data[sample_index][1])
            CEL_list_name.append( data[sample_index][1])
        cellist=robjects.StrVector(CEL_list)
        cellistname=robjects.StrVector(CEL_list_name)
        norm_data[platform]['list']=CEL_list
        samplenames=robjects.StrVector(CEL_list_name) #should be changed.. -- SEEMS OK NOW.

        print cellist
        print cellistname

        dname='HSB_Norm_data_%s.dat' % platform
        rname='HSB_Norm_rownames_%s.dat' % platform
        cname='HSB_Norm_colnames_%s.dat' % platform

        if (not Saved_data) and (not os.path.isfile(dname)):
            print 'Running gcRMA analysis...'
            if Bias_correct:
                affybatch = affy.read_affybatch(filenames = cellist)
                affybatch_RMA = affy.rma(affybatch)
                affybatch_bias_matrix = bias.getBiasMetrics(affybatch, affybatch_RMA)
                affybatch_RMA_Bias_corrected = bias.biasCorrection(affybatch_RMA, affybatch_bias_matrix)
                HSB_Norm_exp=Biobase.exprs(affybatch_RMA_Bias_corrected)
            else:

                affybatch = affy.read_affybatch(filenames = cellist)
                HSB_Norm  = gcrma.gcrma(affybatch)
                HSB_Norm_exp=Biobase.exprs(HSB_Norm)

            HSB_Norm_colnames = numpy.array(samplenames) #remove the .CEL extention from the name..
            HSB_Norm_rownames = numpy.array(HSB_Norm_exp.rownames)
            HSB_Norm_data = numpy.array(HSB_Norm_exp)

            HSB_Norm_data.dump(dname)
            HSB_Norm_rownames.dump(rname)
            HSB_Norm_colnames.dump(cname)
        else:
            print 'Loading from numpy...'
            HSB_Norm_colnames = numpy.load(cname)
            HSB_Norm_rownames = numpy.load(rname)
            HSB_Norm_data = numpy.load(dname)

            #data loaded, now recreate R object.
#           HSB_Norm = rebuild_R_object(data=HSB_Norm_data,rownames=HSB_Norm_rownames,colnames=HSB_Norm_colnames)
            print 'Done.'
        #save in  dictionarry.
        norm_data[platform]['data'] = HSB_Norm_data
        norm_data[platform]['colnames'] = HSB_Norm_colnames
        norm_data[platform]['rownames'] = HSB_Norm_rownames

    #now remove data from HG-133_plus2 that is not present in HG133A or B, merge raw normalized data.

    #merge new HG-U133_Plus_2 together with new experiments.
    try:
        norm_data['HG-U133_Plus_2']['data']= numpy.concatenate( (norm_data['HG-U133_Plus_2']['data'], norm_data['HG-U133_Plus_2_new']['data']), axis=1)
        norm_data['HG-U133_Plus_2']['colnames']= numpy.concatenate( (norm_data['HG-U133_Plus_2']['colnames'], norm_data['HG-U133_Plus_2_new']['colnames']), axis=1)
        norm_data.pop('HG-U133_Plus_2_new')
    except:
        print 'Bric experiments already in db..'
    return [data, norm_data]

    pass


robjects.r('''
    qn <- function(celfiles) {
    library(affy)
    print(celfiles)
    rawdata<-ReadAffy(filenames=celfiles)
    PM<-probes(rawdata,which="pm")
    AffyInfo<-dimnames(PM)[[1]]
    cutpos<-regexpr("\\\d+$",AffyInfo,perl=T)
    AffyID<-substr(AffyInfo,1,cutpos-1)
    probe<-as.numeric(substr(AffyInfo,cutpos,nchar(AffyInfo)))
    data.bgc<-bg.correct(rawdata,method="rma")
    data.bgc.q<-normalize.AffyBatch.quantiles(data.bgc,type="pmonly")
    pm.bgc.q<-probes(data.bgc.q,which="pm")
    normalized<-cbind(AffyID,probe,pm.bgc.q)
    return (normalized)
    }
    ''')

r_quantile_normalize = robjects.globalenv['qn']
def quantile_norm(self):
    """docstring for quantile_norm
        run basic quantile norm. almost like RMA.
    """
    TRUE=1
    FALSE=0
    packages=['affy','hgu133plus2.db','hgu133a.db','hgu133b.db',  'gplots','genefilter' , 'Biobase', 'bias' , 'squash',  'amap']
    affy, hgu133plus2, hgu133a, hgu133b, gplots, genefilter, Biobase , bias, squash, amap = import_rlib(packages)

    f=open(_file)
    t=f.read()
    f.close()
    data=[]
    i=0
    for sample in t.splitlines():
        #skip first line
        if i==0:
            i+=1
        else:
            data.append( sample.split('\t') )
            for j,field in enumerate(data[-1]):
                data[-1][j] = field.strip() #clears white spaces before and after..
            i+=1

    #separate the data according to the platform and normalize them separatelly:
    #possible platform are  : HG-U133A, HG-U133B, HG-U133_Plus_2 (all probes in A+B = probes in Plus_2)
    index_by_platform={}
    norm_data={}
#   platforms = ['HG-U133A','HG-U133B','HG-U133_Plus_2','HG-U133_Plus_2_new' ]
#   platforms = ['HG-U133A','HG-U133B','HG-U133_Plus_2' ]
    for platform in platforms:
        print 'Processing samples from %s' % platform
        index_by_platform[platform] = [i for i,j in enumerate(data) if data[i][7]==platform ]

        #Run RMA Normalization from R.
        norm_data[platform]={}
        #prepare object for justRMA.
        CEL_list=[]
        CEL_list_name=[]
        for sample_index in index_by_platform[platform]:
            CEL_list.append( celfilepath +data[sample_index][1])
            CEL_list_name.append( data[sample_index][1])
        cellist=robjects.StrVector(CEL_list)
        cellistname=robjects.StrVector(CEL_list_name)
        norm_data[platform]['list']=CEL_list
        samplenames=robjects.StrVector(CEL_list_name) #should be changed.. -- SEEMS OK NOW.

        print cellist
        print cellistname

        dname='HSB_Norm_data_%s.dat' % platform
        rname='HSB_Norm_rownames_%s.dat' % platform
        cname='HSB_Norm_colnames_%s.dat' % platform

        if (not Saved_data) and (not os.path.isfile(dname)):
            print 'Running quantile analysis...'
            if Bias_correct:
                affybatch = affy.read_affybatch(filenames = cellist)
                affybatch_RMA = affy.rma(affybatch)
                affybatch_bias_matrix = bias.getBiasMetrics(affybatch, affybatch_RMA)
                affybatch_RMA_Bias_corrected = bias.biasCorrection(affybatch_RMA, affybatch_bias_matrix)
                HSB_Norm_exp=Biobase.exprs(affybatch_RMA_Bias_corrected)
            else:

                HSB_Norm =  r_quantile_normalize(cellistname)
                # now check everything...
                print 'BOAGG'
                HSB_Norm_exp=Biobase.exprs(HSB_Norm)

            HSB_Norm_colnames = numpy.array(samplenames) #remove the .CEL extention from the name..
            HSB_Norm_rownames = numpy.array(HSB_Norm_exp.rownames)
            HSB_Norm_data = numpy.array(HSB_Norm_exp)

            HSB_Norm_data.dump(dname)
            HSB_Norm_rownames.dump(rname)
            HSB_Norm_colnames.dump(cname)
        else:
            print 'Loading from numpy...'
            HSB_Norm_colnames = numpy.load(cname)
            HSB_Norm_rownames = numpy.load(rname)
            HSB_Norm_data = numpy.load(dname)

            #data loaded, now recreate R object.
#           HSB_Norm = rebuild_R_object(data=HSB_Norm_data,rownames=HSB_Norm_rownames,colnames=HSB_Norm_colnames)
            print 'Done.'
        #save in  dictionarry.
        norm_data[platform]['data'] = HSB_Norm_data
        norm_data[platform]['colnames'] = HSB_Norm_colnames
        norm_data[platform]['rownames'] = HSB_Norm_rownames

    #now remove data from HG-133_plus2 that is not present in HG133A or B, merge raw normalized data.

    #merge new HG-U133_Plus_2 together with new experiments.
    try:
        norm_data['HG-U133_Plus_2']['data']= numpy.concatenate( (norm_data['HG-U133_Plus_2']['data'], norm_data['HG-U133_Plus_2_new']['data']), axis=1)
        norm_data['HG-U133_Plus_2']['colnames']= numpy.concatenate( (norm_data['HG-U133_Plus_2']['colnames'], norm_data['HG-U133_Plus_2_new']['colnames']), axis=1)
        norm_data.pop('HG-U133_Plus_2_new')
    except:
        print 'Bric experiments already in db..'
    return [data, norm_data]



    pass

def import_sample_data(_file,sep='\t',verbose=False, line_to_skip=1):
    """docstring for import_sample_data
    import data from tab separated file containing sample description.
    """
    if _file.split('.')[-1] == 'gz':
        f=gzip.open(_file,'r')
    else:
        f=open(_file,'r')
    t=f.readlines()
    f.close()
    data=[]
    for sample in t[line_to_skip:]:
        data.append( [i.strip() for i in sample.strip().split(sep)] )
        if verbose:
            print '#%s#'% data[-1]
    return data
    pass

def load_ref_seq_to_gene_names():
    """docstring for load_ref_seq_to_gene_names
        load ref correspondance between refseq and official gene names used in pubmed.
    """
    refseq2gene={}
    gene2refseq={}

    f=open('/Users/nrapin/Dropbox/Python-stuff/gene_name_to_refseq.txt')
    t=f.read()
    f.close()
    data=[]
    i=0
    for sample in t.splitlines():
        #skip first line
        if (i==0) or (sample.split()[0].find('withdrawn')!=-1):
            i+=1
        else:
            data.append( sample.split() )
            for j,field in enumerate(data[-1]):
                data[-1][j] = field.strip() #clears white spaces before and after..
            i+=1
    for i in data:
#       print i
        try:
            refseq2gene[i[1]]=i[0]
            gene2refseq[i[0]]=i[1]
        except Exception,e:
            s=0
#           print str(i) + 'not referenced in refseq.'
    return refseq2gene, gene2refseq
    pass


def radviz(frame, class_column, ax=None, **kwds):
    """RadViz - a multivariate data visualization algorithm

    Parameters:
    -----------
    frame: DataFrame object
    class_column: Column name that contains information about class membership
    ax: Matplotlib axis object, optional
    kwds: Matplotlib scatter method keyword arguments, optional

    Returns:
    --------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.text as text
    import random
    import pandas.core.common as com

    def random_color(column):
        random.seed(column)
        return [random.random() for _ in range(3)]

    def normalize(series):
        a = min(series)
        b = max(series)
        return (series - a) / (b - a)

    column_names = [column_name for column_name in frame.columns
                    if column_name != class_column]

    df = frame[column_names].apply(normalize)

    if ax is None:
        ax = plt.gca(xlim=[-1, 1], ylim=[-1, 1])

    classes = set(frame[class_column])
    to_plot = {}

    for class_ in classes:
        to_plot[class_] = [[], []]

    n = len(frame.columns) - 1
    s = numpy.array([(numpy.cos(t), numpy.sin(t))
                  for t in [2.0 * numpy.pi * (i / float(n))
                            for i in range(n)]])

    for i in range(len(frame)):
        row = df.irow(i).values
        row_ = numpy.repeat(numpy.expand_dims(row, axis=1), 2, axis=1)
        y = (s * row_).sum(axis=0) / row.sum()
        class_name = frame[class_column].iget(i)
        to_plot[class_name][0].append(y[0])
        to_plot[class_name][1].append(y[1])

    for class_ in classes:
        line = ax.scatter(to_plot[class_][0],
                          to_plot[class_][1],
                          color=random_color(class_),
                          label=com.pprint_thing(class_), **kwds)
    ax.legend()

    ax.add_patch(patches.Circle((0.0, 0.0), radius=1.0, facecolor='none'))

    for xy, name in zip(s, column_names):

        ax.add_patch(patches.Circle(xy, radius=0.025, facecolor='gray'))

        if xy[0] < 0.0 and xy[1] < 0.0:
            ax.text(xy[0] - 0.025, xy[1] - 0.025, name,
                    ha='right', va='top', size='small')
        elif xy[0] < 0.0 and xy[1] >= 0.0:
            ax.text(xy[0] - 0.025, xy[1] + 0.025, name,
                    ha='right', va='bottom', size='small')
        elif xy[0] >= 0.0 and xy[1] < 0.0:
            ax.text(xy[0] + 0.025, xy[1] - 0.025, name,
                    ha='left', va='top', size='small')
        elif xy[0] >= 0.0 and xy[1] >= 0.0:
            ax.text(xy[0] + 0.025, xy[1] + 0.025, name,
                    ha='left', va='bottom', size='small')

    ax.axis('equal')
    return ax


def RMA(_file,sample,Sample_list='/Users/nrapin/BRIC/MILE_study/Mile_samples.txt',Saved_data=False,Bias_correct=False, \
	celfilepath='/Volumes/B52/Converted/', platforms=['HG-U133A','HG-U133B','HG-U133_Plus_2'], mergeplus2=False):
    """docstring for RMA
    run RMA on samples.
    """
    #load RMA from R as txt file...

    TRUE=1
    FALSE=0

    #packages=['affy','hgu133plus2.db','hgu133a.db','hgu133b.db',  'gplots','genefilter' , 'Biobase',  'squash',  'amap']
    #affy, hgu133plus2, hgu133a, hgu133b, gplots, genefilter, Biobase ,  squash, amap = import_rlib(packages)

    packages = ['affy','Biobase']
    affy , Biobase = import_rlib(packages)


    data=import_sample_data(_file)
    all_samples=import_sample_data(Sample_list)
    print all_samples
    #find sample and append it to the list. (now the sample might be from a U133A platform, and maybe there might be a B platform with it)
#   print data
#   print _file
    if sample != '':
        print 'Im here'
        for j in all_samples:
            if j[0]=='Sample':
                print 'OOOOOOOOOOOOOOOOOOOOOOOO'
                sys.exit(2)
            for i in j:
                if i==sample:
                    data.append(j)
                    print 'and there too!'
                    break
        if data[-1][7]=='HG-U133A':
            print 'Trying to locate B platform associated with sample'
            try:
                for j in all_samples:
                    if (j[2]==data[-1][2]) and (j[7]!=data[-1][7]): #same smple name, but different platform.
                        data.append(j)
                        print 'B sample is ' + j[0]
                        break
            except Exception, e:
                print '%s does not have a HG-U133B counterpart. '%(data[-1][0])
        if data[-1][7]=='HG-U133B':
            print 'Trying to locate A platform associated with sample'
            found=0
            try:
                for j in all_samples:
                    if (j[2]==data[-1][2]) and (j[7]!=data[-1][7]): #same smple name, but different platform.
                        data.append(j)
                        found=1
                        print 'A sample is ' + j[0]
                        print 'Now checking if the file has been run already...'
                        if os.path.isfile('%s.CEL.gz_merged_data.pkl'%(j[0])):
                            print 'Already run..'
                            raise SystemExit
                        break
                if found == 0:
                    print 'hum...'
                    raise SystemExit

            except Exception, e:
                print '%s does not have a HG-U133A counterpart. '%(data[-1][0])


        else:
            print '---->'+ data[-1][7]

    #separate the data according to the platform and normalize them separatelly:
    #possible platform are  : HG-U133A, HG-U133B, HG-U133_Plus_2 (all probes in A+B = probes in Plus_2)
    index_by_platform={}
    norm_data={}
#   platforms = ['HG-U133A','HG-U133B','HG-U133_Plus_2','HG-U133_Plus_2_new' ]
#   platforms = ['HG-U133A','HG-U133B','HG-U133_Plus_2' ]
    for platform in platforms:

        print 'Processing samples from %s' % platform
        index_by_platform[platform] = [i for i,j in enumerate(data) if data[i][7]==platform ]
        print numpy.array(data)[index_by_platform[platform]]
        #Run RMA Normalization from R.
        norm_data[platform]={}
        #prepare object for justRMA.
        CEL_list=[]
        CEL_list_name=[]
        for sample_index in index_by_platform[platform]:
            #if celfilepath +data[sample_index][1]+'' not in CEL_list:
            CEL_list.append( celfilepath +data[sample_index][1]+'')
            #if data[sample_index][1]+'' not in CEL_list_name:
            CEL_list_name.append( data[sample_index][1]+'')
        cellist=robjects.StrVector(CEL_list)
        cellistname=robjects.StrVector(CEL_list_name)
        norm_data[platform]['list']=CEL_list
        samplenames=robjects.StrVector(CEL_list_name) #should be changed.. -- SEEMS OK NOW.

        print cellist
        print cellistname

        dname='HSB_Norm_data_%s.dat' % platform
        rname='HSB_Norm_rownames_%s.dat' % platform
        cname='HSB_Norm_colnames_%s.dat' % platform

        if (not Saved_data) or (not os.path.isfile(dname)):
            print 'Running RMA analysis...'
            if Bias_correct:
                affybatch = affy.read_affybatch(filenames = cellist)
                affybatch_RMA = affy.rma(affybatch)
                affybatch_bias_matrix = bias.getBiasMetrics(affybatch, affybatch_RMA)
                affybatch_RMA_Bias_corrected = bias.biasCorrection(affybatch_RMA, affybatch_bias_matrix)
                HSB_Norm_exp=Biobase.exprs(affybatch_RMA_Bias_corrected)
            else:
                print 'len celllist name ' + str( len(cellistname))
#               print cellistname
                HSB_Norm =  affy.justRMA(filenames = cellistname, celfile_path=celfilepath , background=TRUE, normalize=TRUE,  bgversion=2, sampleNames=cellistname)
                HSB_Norm_exp=Biobase.exprs(HSB_Norm)
#               print HSB_Norm_exp

            cnames=robjects.r('colnames')
            rnames=robjects.r('rownames')
            HSB_Norm_colnames = cnames(HSB_Norm) #remove the .CEL extention from the name..
            HSB_Norm_rownames = rnames(HSB_Norm)#_exp.rownames)
            HSB_Norm_data = numpy.array(HSB_Norm_exp)
            print 'Not saving data on disk.'
            #HSB_Norm_data.dump(dname)
            #HSB_Norm_rownames.dump(rname)
            #HSB_Norm_colnames.dump(cname)
        else:
            print 'Loading from numpy...'
            HSB_Norm_colnames = numpy.load(cname)
            HSB_Norm_rownames = numpy.load(rname)
            HSB_Norm_data = numpy.load(dname)

            #data loaded, now recreate R object.
#           HSB_Norm = rebuild_R_object(data=HSB_Norm_data,rownames=HSB_Norm_rownames,colnames=HSB_Norm_colnames)
            print 'Done.'
        #save in  dictionarry.
        norm_data[platform]['data'] = HSB_Norm_data
        norm_data[platform]['colnames'] = HSB_Norm_colnames
        norm_data[platform]['rownames'] = HSB_Norm_rownames
        print HSB_Norm_data
    if mergeplus2:
    #merge new HG-U133_Plus_2 together with new experiments.
        try:
            norm_data['HG-U133_Plus_2']['data']= numpy.concatenate( (norm_data['HG-U133_Plus_2']['data'], norm_data['HG-U133_Plus_2_new']['data']), axis=1)
            norm_data['HG-U133_Plus_2']['colnames']= numpy.concatenate( (norm_data['HG-U133_Plus_2']['colnames'], norm_data['HG-U133_Plus_2_new']['colnames']), axis=1)
            norm_data.pop('HG-U133_Plus_2_new')
        except:
            print 'Bric experiments already in db..'

    return [data, norm_data]

    pass
#
def score_msig_db(data, sig=None):
    """docstring for score_msig_db
        uses the expression data to score all genesin array. return the score of each signature for each sample, and a rownames (the order fo the signatures.)
    """
    if sig == None:
        print '!! Using the signature file generated by GSEA on ver. 3 of MSIG and plus2 platform. C2 and C5'
        sig = import_sample_data('/Users/nrapin/BRIC/databases/Cmap2/msig_C2_C5_converted_to_plus_2.gmt',line_to_skip=0)
    else:
        sig = import_sample_data(sig,line_to_skip=0)

    # first, transfor sgi into a dic wih indexes:
    print 'Processing signature db:'
    prog = ProgressBar(0, len(sig) , 50, mode='fixed')
    sig_dic= {}
    for index, i in enumerate(sig):
        sig_dic[i[0]] = get_probe_index(i[2:], data['rownames'], verbose=False)
        oldprog = str(prog)
        prog.update_amount(index)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)

#   del sig
#   return sig_dic
    print 'Computing scores for samples.'

    prog = ProgressBar(0, len(sig) , 50, mode='fixed')
    sig_names = numpy.array(sig_dic.keys())
    sig_mat = numpy.zeros([len(sig), len(data['colnames'])])
    for index_i, i in enumerate(sig_dic):
        for index_j, j in enumerate(data['data'].T):
            try:
                sig_mat[index_i,index_j] =  numpy.median(numpy.array(j[sig_dic[i]]))
            except Exception, e:
                sig_mat[index_i,index_j] = 0
        oldprog = str(prog)
        prog.update_amount(index_i)
        if oldprog != str(prog):
            print prog, "\r",
            sys.stdout.flush()
            oldprog=str(prog)

    return sig_names , sig_mat
    pass

def center_array_values(matrix, mode='median'):
    """docstring for center_array_values
        works by removing the median expression accross genes to the actual expression value.
        return centered matrix

        options:
            - median
            - median/sd
            - median/var
    """
    centered=numpy.zeros(matrix.shape)
    if mode == 'median':
        nbsamples=matrix.shape[1]
        median=numpy.median(matrix,axis=1)
        median_matrix= numpy.array([median]*nbsamples).transpose()
        centered = numpy.subtract(matrix, median_matrix)

    if mode == 'median/sd':
        nbsamples=matrix.shape[1]
        median=numpy.median(matrix,axis=1)
        sd=numpy.std(matrix,axis=1)
        median_matrix= numpy.array([median]*nbsamples).transpose()
        sd_matrix= numpy.array([sd]*nbsamples).transpose()
        centered = numpy.divide(numpy.subtract(matrix, median_matrix),sd_matrix)

    if mode == 'median/var':
        nbsamples=matrix.shape[1]
        median=numpy.median(matrix,axis=1)
        var=numpy.var(matrix,axis=1)
        median_matrix= numpy.array([median]*nbsamples).transpose()
        var_matrix= numpy.array([var]*nbsamples).transpose()
        centered = numpy.divide(numpy.subtract(matrix, median_matrix),var_matrix)

    return centered
    pass

def cosine_dist(x,y):
    """docstring for cosine_dist
    calculates cosine distance of two vector.
    """
    magnitude_x=numpy.sqrt( numpy.sum( numpy.square(x) ) )
    magnitude_y=numpy.sqrt( numpy.sum( numpy.square(y) ) )
    dot_product=numpy.dot(x,y)

    return dot_product/(magnitude_x*magnitude_y)

    pass

def merge_array_data(norm_data,data,center_data=0,use_AB=True,mergeplu2=False):
    """merge_array_data:

    merge data from A, B and PLUS_2 affimetrix platforms.

    input is a dictionary with keys as platform, each platform has 3 keys: rownames, colnames and the actual data.

    returns a dictionarry with 3 keys { merged_rownames, colnames and the actual merged_data}
    """
    lp2=norm_data['HG-U133_Plus_2']['colnames'].shape[0]
    if use_AB:
        la=norm_data['HG-U133A']['colnames'].shape[0]
        lb=norm_data['HG-U133B']['colnames'].shape[0]
        if (la != lb):
            print 'samples in A are not the same as in B---\n prepare for Crash!'
    if use_AB:
        b=norm_data['HG-U133B']['rownames'].tolist()
        a=norm_data['HG-U133A']['rownames'].tolist()
        b_C=norm_data['HG-U133B']['colnames'].tolist()
        a_C=norm_data['HG-U133A']['colnames'].tolist()
    p2=norm_data['HG-U133_Plus_2']['rownames'].tolist()

    #find probes uniques to p2 compared to a or b:
    if use_AB:
        p2_dic={}
        for i,probe in enumerate(p2):
            p2_dic[probe]=i
        ab_dic={}
        aplusb=a+b
        for i,probe in enumerate(aplusb):
            ab_dic[probe]=i
        uniques= [p2_dic[probe] for probe in p2 if not ( ab_dic.has_key(probe) )  ]
#old slow version....of the 8 lines aboves
#       uniques= [p2.index(probe) for probe in p2 if not ( (probe in b ) or  (probe in a) )  ]

#   print 'done'
    #merge expression data

        ran=range(len(p2))
        index_present_in_P2 = [i for i in ran if i not in uniques]
        rows=len(index_present_in_P2)
        cols=sum([norm_data[i]['colnames'].shape[0] for i in norm_data.keys()])
        cols-=la #la and lb should be same, as they represent the same samples-jut on != platforms. therefore code line above takes into account la too many namples.
        m=numpy.zeros([rows,cols])
        #map data from A nd B array based on the Annotation table in data
        annotation= [[i[0].replace('.CEL','').replace('.gz',''),i[2], i[7]] for i in data if  i[7] == 'HG-U133B' or i[7] == 'HG-U133A']
        pairs=[]
        while annotation != []: # While loop is necessary , as the break statement exits all for loops for some reason...
            for sample_name_1 in annotation:
                name_A=sample_name_1[1]
                p_1=sample_name_1[2]
                for sample_name_2 in annotation:
                    name_B=sample_name_2[1]
                    p_2=sample_name_2[2]
                    if (name_A==name_B) and (p_1!= p_2) :
                        print '(%s(%s), %s(%s))' %(sample_name_1[1],p_1,sample_name_2[1],p_2 )
                        # add to pairs (in ordered platform, A then B)
                        if p_2=='HG-U133B':
                            pairs.append( (sample_name_1[0],sample_name_2[0] ) )
                        else:
                             pairs.append( (sample_name_2[0],sample_name_1[0] ) )
                        # remove from annotation.
                        annotation.remove(sample_name_2)
                        annotation.remove(sample_name_1)
                        break
        #check that the code is not doing zouaveries.
        if len(pairs)!=la:
            print 'Fishy number of smaple pairs...'
            raise


        A_B_rownames= a+b

    # create matrix with rows = samples , cols = expression data in platform A + B
        A_B_matrix=numpy.zeros([len(A_B_rownames),la])
        P2_matrix=numpy.asarray(norm_data['HG-U133_Plus_2']['data'])
        print pairs
        for i,pair in enumerate(pairs):
            a_index=a_C.index(str(pair[0].replace('.CEL','').replace('.gz',''))  + '.CEL.gz')
            b_index=b_C.index(str(pair[1].replace('.CEL','').replace('.gz','')) + '.CEL.gz')
            a_data=norm_data['HG-U133A']['data'][:,a_index]
            b_data=norm_data['HG-U133B']['data'][:,b_index]
            A_B_matrix[:,i]=numpy.concatenate((a_data,b_data),axis=1)

        labr = range(int(lp2), int(lp2+la) )
########################## end of part specific to A+B
    else:
        index_present_in_P2 = range(len(p2))
        P2_matrix=numpy.asarray(norm_data['HG-U133_Plus_2']['data'])

    lp2r = range(int(lp2))

    if center_data:
        print 'Using centered data.'
        if use_AB:
            A_B_centered = center_array_values(A_B_matrix,mode='median')
        P_2_centered = center_array_values(P2_matrix,mode='median')
        #build matrix.
        if use_AB:
            for (i,probe_i) in enumerate(index_present_in_P2):
                r_index= int(A_B_rownames.index(p2[probe_i]))
                m[i][lp2r]=P_2_centered[probe_i,:]
                m[i][labr]=A_B_centered[r_index,:]
        else:
            m=numpy.zeros(numpy.asarray(norm_data['HG-U133_Plus_2']['data']).shape)
            m=P_2_centered

    else:
        print 'Using non-centered data'

    #build matrix.
    if use_AB:
        for (i,probe_i) in enumerate(index_present_in_P2):
            r_index= int(A_B_rownames.index(p2[probe_i]))
            m[i][lp2r]=P2_matrix[probe_i,:]
            m[i][labr]=A_B_matrix[r_index,:]
    else:
        m=numpy.zeros(numpy.asarray(norm_data['HG-U133_Plus_2']['data']).shape)
        m=P2_matrix


    # Wrap up everyhing!
    result={}

    final_rownames=[p2[probe_i] for probe_i in index_present_in_P2 ]
    final_platforms =[]
    final_colnames=[]
    final_stemness=[]
    final_celfile=[]
    final_X1=[]
    final_X2=[]
    for celfile in norm_data['HG-U133_Plus_2']['colnames']:
        for sample_description in data:
            if celfile == sample_description[1]:
                final_colnames.append(sample_description[5] )# + ' '  + sample_description[1]) #+sample_description[7][-1])+sample_description[7][-1])
                final_stemness.append(sample_description[11])
                final_celfile.append(sample_description[1])
                final_X1.append(sample_description[3])
                final_X2.append(sample_description[4])
                if sample_description[7]=='HG-U133_Plus_2':
                    final_platforms.append(sample_description[7])
                else:
                    if mergeplu2:
                        final_platforms.append('HG-U133_Plus_2')
                    else:
                        final_platforms.append('HG-U133_Plus_2_new')
#               print sample_description[2], celfile
    if use_AB:
        for i,pair in enumerate(pairs):
            for sample_description in data:
                if str(pair[0])+'.CEL.gz' == sample_description[1]:
                    print 1
                    final_stemness.append(sample_description[11])
                    final_colnames.append(sample_description[5] )#  + ' ' + sample_description[1]) #+sample_description[7][-1])
                    final_platforms.append(sample_description[7])
                    final_celfile.append(sample_description[1])
                    final_X1.append(sample_description[3])
                    final_X2.append(sample_description[4])
    #               print sample_description[2]


    #Sort order of cel files according to stemness: (this is already done in the excel file)
#   [(data.index(i),i[5],i[11],i[1]) for i in sorted(data, key=lambda i: int(i[11] )  )]
# this is why people should program in C!
    x=numpy.array([final_colnames,final_stemness,final_celfile]).transpose().tolist()
    print x
#   [(x.index(i),i[0],i[1]) for i in sorted( x , key=lambda i:( int(i[1]),i[0],i[2] ) )]
    order=[x.index(i) for i in sorted( x , key=lambda i:( int(i[1]),i[0],i[2] ) )]
    print len(order)
    print len(final_celfile)
    result['cel_file']=numpy.array(final_celfile)[order]
    result['stemness']=numpy.array(final_stemness)[order]
    result['platforms']=numpy.array(final_platforms)[order]
    result['rownames']=numpy.array(final_rownames)
    result['colnames']=numpy.array(final_colnames)[order]
    result['data']=numpy.copy(m[:,order])
    result['X1'] =  numpy.array(final_X1)[order]
    result['X2'] =  numpy.array(final_X2)[order]
    print

    return result
    pass




def plot_corr(x,y):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import pylab
    sns.regplot(numpy.array(x,dtype='double'),numpy.array(y,dtype='double'))
    main, x_marg, y_marg = plt.gcf().axes
    sns.despine(ax=main)
    sns.despine(ax=x_marg, left=True)
    sns.despine(ax=y_marg, bottom=True)
    pylab.savefig('correlaiton.pdf')

    pass


def merge_array_data_only_A(norm_data,data,center_data=0,use_AB=True,mergeplu2=False):
    """merge_array_data:

    merge data from A, and PLUS_2 affimetrix platforms.

    input is a dictionary with keys as platform, each platform has 3 keys: rownames, colnames and the actual data.

    returns a dictionarry with 3 keys { merged_rownames, colnames and the actual merged_data}
    """
    lp2=norm_data['HG-U133_Plus_2']['colnames'].shape[0]
    if use_AB:
        la=norm_data['HG-U133A']['colnames'].shape[0]
    if use_AB:
        a=norm_data['HG-U133A']['rownames'].tolist()
        a_C=norm_data['HG-U133A']['colnames'].tolist()
    p2=norm_data['HG-U133_Plus_2']['rownames'].tolist()

    #find probes uniques to p2 compared to a or b:
    if use_AB:
        p2_dic={}
        for i,probe in enumerate(p2):
            p2_dic[probe]=i
        ab_dic={}
        aplusb=a
        for i,probe in enumerate(aplusb):
            ab_dic[probe]=i
        uniques= [p2_dic[probe] for probe in p2 if not ( ab_dic.has_key(probe) )  ]
#old slow version....of the 8 lines aboves
#       uniques= [p2.index(probe) for probe in p2 if not ( (probe in b ) or  (probe in a) )  ]

#   print 'done'
    #merge expression data

        ran=range(len(p2))
        index_present_in_P2 = [i for i in ran if i not in uniques]
        rows=len(index_present_in_P2)
        cols=sum([norm_data[i]['colnames'].shape[0] for i in norm_data.keys()])

        m=numpy.zeros([rows,cols])
        #map data from A nd B array based on the Annotation table in data
        annotation= [[i[0].strip('.CEL'),i[2], i[7]] for i in data if i[7] == 'HG-U133A']

        pairs=annotation

        A_B_rownames= a

    # create matrix with rows = samples , cols = expression data in platform A + B
        A_B_matrix=numpy.zeros([len(A_B_rownames),la])
        P2_matrix=numpy.asarray(norm_data['HG-U133_Plus_2']['data'])
        for i,pair in enumerate(pairs):
            a_index=a_C.index(str(pair[0])  + '.CEL')
            a_data=norm_data['HG-U133A']['data'][:,a_index]
            A_B_matrix[:,i]=a_data

        labr = range(int(lp2), int(lp2+la) )
########################## end of part specific to A+B
    else:
        index_present_in_P2 = range(len(p2))
        P2_matrix=numpy.asarray(norm_data['HG-U133_Plus_2']['data'])

    lp2r = range(int(lp2))

    if center_data:
        print 'Using centered data.'
        if use_AB:
            A_B_centered = center_array_values(A_B_matrix,mode='median')
        P_2_centered = center_array_values(P2_matrix,mode='median')
        #build matrix.
        if use_AB:
            for (i,probe_i) in enumerate(index_present_in_P2):
                r_index= int(A_B_rownames.index(p2[probe_i]))
                m[i][lp2r]=P_2_centered[probe_i,:]
                m[i][labr]=A_B_centered[r_index,:]
        else:
            m=numpy.zeros(numpy.asarray(norm_data['HG-U133_Plus_2']['data']).shape)
            m=P_2_centered

    else:
        print 'Using non-centered data'

    #build matrix.
    if use_AB:
        for (i,probe_i) in enumerate(index_present_in_P2):
            r_index= int(A_B_rownames.index(p2[probe_i]))
            m[i][lp2r]=P2_matrix[probe_i,:]
            m[i][labr]=A_B_matrix[r_index,:]
    else:
        m=numpy.zeros(numpy.asarray(norm_data['HG-U133_Plus_2']['data']).shape)
        m=P2_matrix


    # Wrap up everyhing!
    result={}

    final_rownames=[p2[probe_i] for probe_i in index_present_in_P2 ]
    final_platforms =[]
    final_colnames=[]
    final_stemness=[]
    final_celfile=[]
    for celfile in norm_data['HG-U133_Plus_2']['colnames']:
        for sample_description in data:
            if celfile == sample_description[1]:
                final_colnames.append(sample_description[5] )# + ' '  + sample_description[1]) #+sample_description[7][-1])+sample_description[7][-1])
                final_stemness.append(sample_description[11])
                final_celfile.append(sample_description[1])
                if sample_description[7]=='HG-U133_Plus_2':
                    final_platforms.append(sample_description[7])
                else:
                    if mergeplu2:
                        final_platforms.append('HG-U133_Plus_2')
                    else:
                        final_platforms.append('HG-U133_Plus_2_new')
#               print sample_description[2], celfile
    if use_AB:
        for i,pair in enumerate(pairs):
            for sample_description in data:
                if str(pair[0])+'.CEL' == sample_description[1]:
                    final_stemness.append(sample_description[11])
                    final_colnames.append(sample_description[5] )#  + ' ' + sample_description[1]) #+sample_description[7][-1])
                    final_platforms.append(sample_description[7])
                    final_celfile.append(sample_description[1])
    #               print sample_description[2]


    #Sort order of cel files according to stemness: (this is already done in the excel file)
#   [(data.index(i),i[5],i[11],i[1]) for i in sorted(data, key=lambda i: int(i[11] )  )]
# this is why people should program in C!
    x=numpy.array([final_colnames,final_stemness,final_celfile]).transpose().tolist()
#   [(x.index(i),i[0],i[1]) for i in sorted( x , key=lambda i:( int(i[1]),i[0],i[2] ) )]
    order=[x.index(i) for i in sorted( x , key=lambda i:( int(i[1]),i[0],i[2] ) )]

    result['cel_file']=numpy.array(final_celfile)[order]
    result['stemness']=numpy.array(final_stemness)[order]
    result['platforms']=numpy.array(final_platforms)[order]
    result['rownames']=numpy.array(final_rownames)
    result['colnames']=numpy.array(final_colnames)[order]
    result['data']=numpy.copy(m[:,order])
    print

    return result
    pass







def corr(x,y):
    """docstring for corr
    pearson correlation...
    """
    xv=numpy.asarray(x)
    yv=numpy.asarray(y)
    return 1- Bio.Cluster.distancematrix( (xv ,yv), dist="c")[1][0]
    pass

def clean_signatures(sig,gene2affy):
    """docstring for clean_signatures
    clean gene names that are not found in affy name probes
    """
    for gene_set in sig.iterkeys():
        for regulation in sig[gene_set].iterkeys():
            toremove=[]
            for gene in sig[gene_set][regulation]:
                if not (gene in gene2affy):
                    print 'missing gene\t%s in\t%s %s'%(gene,gene_set, regulation)
                    toremove.append(gene)
#           print toremove
            for gene in toremove:
                sig[gene_set][regulation].remove(gene)
    return sig
    pass

def activation_scores(signatures,exp_data,row_names,col_names,samples2use,gene2affy):
    """docstring for activation scores"""
    rows=len(signatures)
    print 'Signatures:',rows
    cols=len(samples2use) # for all : GSE13204_tr_exp_data.shape[1]
    print 'Samples:',cols
    probe_list=row_names.tolist()
    m=numpy.zeros([rows,cols])
    #find indexes of samples to be used in the whole dataset:
    allsamples=col_names.tolist()
#   print 'all samples\n',allsamples
    samples_index=[allsamples.index(i) for i in samples2use]
#   print samples_index

    #loop over signatures
    j=0
    for s in signatures.iterkeys():
        #look up down,up or both regulated genes
        regulation=['up','dn','both']
        for reg_type in regulation:
            if reg_type == 'dn':
                reg_sign=-1
            else:
                reg_sign=1

            if reg_type in signatures[s]:
                # find index of different regulated genes in the affy exp table
                regulated_indexes = [probe_list.index(gene2affy[i]) for i in signatures[s][reg_type]]
                #not regulated is the rest. we just make a list from 0..n and remove the up_indexes.
                not_reg  = range(exp_data.shape[0])
                for i in regulated_indexes:
                    try:
                        not_reg.remove(i)
                    except Exception, e:
                        print e
                        pass
                # now we have the list of regulated genes. compute z-score-thinguy across all samples to be used.
                for i in range(cols):
                    n1=len(regulated_indexes)
                    n2=len(not_reg)
                    sd= numpy.sqrt(n1*n2*(n1+n2+1)/12.0)
                    zscore=reg_sign*mannwhitneyu_1(exp_data[:,samples_index[i]][regulated_indexes,:], exp_data[:,samples_index[i]][not_reg,:],sd)[1]
                    print zscore
                    if reg_type == 'dn':
                        m[j,i]-=zscore
                    else:
                        m[j,i]+= zscore
                    if m[j,i] < 0:
                        print 'signature %s (%s) in sample %s (%s)is negative %f' %(s, reg_type, allsamples[samples_index[i]] , samples2use[i] ,m[j,i])
                    else:
                        print 'signature %s (%s) in sample %s (%s)is positive %f' %(s, reg_type, allsamples[samples_index[i]] , samples2use[i] ,m[j,i])
#                   print 'Bogumil : ', m[j,i]
#                   print 'Nico    : ',reg_sign*mannwhitneyu(exp_data[:,samples_index[i]][regulated_indexes,:], exp_data[:,samples_index[i]][not_reg,:],use_continuity=False)[1]
#                   print
                #a faire en C!
        j+=1
    return m
    pass


def hdf_load(f):
    """docstring for hdf_load
        load a data in hdf5 format to a dictionary- this is be faster than pickle_load, but less compatible.
        this can also handle larger filesself.
        """
    hdf_file=h5py.File(f,'r')

    rootnode=hdf_file[list(hdf_file)[0]]

    dic = hdf52dic(rootnode)

    hdf_file.close()
    return dic

    pass

def hdf52dic(g,dic={}):
    """docstring for hdf52doc

    gets an hdf5 group and makes a dictionarry out of it.
    """

    if dic==None:
        dic={}

    for key in g.iterkeys():
#       print key
        if isinstance(g[key], h5py.highlevel.Group):
#           print '->'
            dic[key]={}
            dic[key]=hdf52dic(g[key],dic[key])

        elif isinstance(g[key], h5py.highlevel.Dataset):
            dic[key]=g[key][...]
#           print '-**'

    return dic


    pass

def hdf_save(d,f,compression_type='lzf', force_write=False):
    """docstring for hdf_save
    try to save the data d into the file f in hdf5 format.
    """
    exit=0
    if os.path.isfile(f) and force_write == False:
        print "file already exist..."
        var = raw_input("Overwrite? [y/n] ")
        if var == 'y':
            exit = 0
            os.remove(f)
            hdf_file=h5py.File(f)
            hdf_file.create_group('root')
        else:
            exit=1
    else:
        try:
            os.remove(f)
        except:
            pass
        hdf_file=h5py.File(f)
        hdf_file.create_group('root')


    if exit:
        return

    root_node=hdf_file['root']
    dict2hdf5(d,root_node, compression_type='gzip')

    hdf_file.close()
    pass


def dict2hdf5( node , sg, compression_type='lzf'):
    """docstring for dict2hdf5
        get a dictionnary, go through children, make appropriate groups/datasets in hdf5 file
        inputs:
            node  : Dictionnary.
            sg    : HDF5 file.
    """
#   compression_type='gzip'
    #put attributes into description

    dup_key_index={}

    print 'begining'
        ######################################
        #subnode with children become groups, they are of type dict, isinstance( node ,dict) == True
        #subnodes without any become datasets.  they are of type list, isinstance( node ,list)  == True
        #######################################

    if isinstance( node , dict): #end of the path, this is therefore a group
        for subnode in node:
#           print subnode
            key = subnode # subnode has children, so let's go recursive!
#           print key
            key = key.replace('-','_').replace(' ','_')
            if isinstance(node[subnode],list) or isinstance(node[subnode],numpy.ndarray):
#               print 'This is a list'
                try:
                    val = node[subnode]
                except Exception,e:
                    print e
                    print val
                    pass
                lval=len(val)
                if (lval == 0) or (val == None):
                    print 'Empty Field'
                else:
                    if not(key in sg):
#                       print val
                        try:
                            sg.create_dataset(key, data=numpy.array(val), compression=compression_type)
                            dup_key_index[key]=0
                        except Exception,e:
                            print e
                            print val
                            pass

                    else:
                        try:
                            dup_key_index[key]+=1
                            sg.create_dataset(key+'_'+str(dup_key_index[key]), data=numpy.array(val), compression=compression_type)
                        except Exception,e:
                            print e
                            print val
                            pass
            elif isinstance(node[subnode],dict):
#               print key,'>'
                if not (key in sg): #and this key is not in the group yet...
                    subsg=sg.create_group(key)
                    dup_key_index[key]=0
                else: #we put a number..
                    if dup_key_index.has_key(key):
                        dup_key_index[key]+=1
                    else:
                        dup_key_index[key]=0
                    subsg=sg.create_group(key+'_'+str(dup_key_index[key]))

                dict2hdf5(node[subnode],subsg)
                    #check atrributes
    else: #if isinstance( subnode ,list)
        print 'OULALA!'

    return sg
    pass




def histo(data, bins=30):
    """docstring for histo"""
    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    data=np.array(data)
    mu = np.median(data)
    sigma = np.std(data)

    # the histogram of the data
    n, bins, patches = plt.hist(data, bins ,normed=1, facecolor='white', alpha=0.75)

    print np.sum(n*np.diff(bins))

    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)

    plt.xlabel('data')
    plt.ylabel('nb (normed)')
#   plt.grid(True)

    plt.show()

    pass


def rebuild_R_object(data=numpy.zeros((1,1)), rownames='row names',colnames='col names',nrows=None,ncols=None):
    """docstring for rebuild_R_object
        rebuild R table from numpy arrays.
    """
    from pandas import DataFrame
    import pandas.rpy.common

    print colnames
    print { colnames[0]:data[0]}
    df = DataFrame({ j+'_'+str(i):data[i] for i,j in enumerate(colnames) }, index=[j+'_'+str(i) for i,j in enumerate(colnames)])
    r_dataframe =  pandas.rpy.common.convert_to_r_dataframe(df)
    #return df
    return r_dataframe
    ##### old code....

    if nrows == None:
        nrows=numpy.matrix(data).shape[1]
    if ncols == None:
        ncols=numpy.matrix(data).shape[0]
    if ncols > 1:
        m = robjects.r.matrix(data,nrow=nrows, ncol=ncols,dimnames= robjects.r.list(robjects.StrVector(rownames), robjects.StrVector(colnames)  ) )
        m = robjects.r.matrix(data,nrow=nrows, ncol=ncols,dimnames= robjects.r.list(robjects.StrVector(rownames), robjects.StrVector(colnames)  ) )
    else:
        m = robjects.r.matrix(data,nrow=nrows, ncol=ncols, dimnames= robjects.r.list(robjects.StrVector(rownames)  ) )

    return m
    pass

def import_rlib(packages):
    """
    Import rlibs so that they can be used.
    Install them if not present ...

    Works as
    """
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr

    libs=[]
    for lib in packages:
        try:
            libs.append(importr(lib))
        except Exception,e:
            print e
            print 'Installing missing package...'
            base = importr('base')
            base.source("http://www.bioconductor.org/biocLite.R")
            bioclite = robjects.globalenv['biocLite']
            bioclite(lib)
            importr(lib)
    return libs
    pass

def find_samples(text,descr):
    samples=[]
    for i in descr.splitlines():
            if i.find(text)!=-1:
#               print i.split(',')[0].strip()
                samples.append(i.split(',')[0].strip())
    return samples

def get_data_from_name(name='',file_name='',verbose=0):
    """docstring for get_data_from_name
        you give a name (and a file as path or h5py file_object), it gives the data.
        the verbose option prints out the description of the columns .
    """
    if name=='':
        print 'Give a sample name...'
        return
    if file_name=='':
        print 'Give a file name...'
        return
    try:
        if str(file_name.__class__) == '<class \'h5py.highlevel.File\'>':
            h5f=file_name
        else:
            try:
                h5f=h5py.File(file_name,'r')
            except Exception,e:
                print e
                return
    except Exception,e:
        print e
        return

    #GSM495861
    for i in h5f.keys():
        if name in h5f[  i  ].keys():
#           if verbose:
#               h5f[ '/'.join([ i , name , 'Data-Table']
#               print  'eee'
            return  h5f[ '/'.join([ i , name , 'Data-Table',name+'-tbl-1'])  ][...]
            break

    pass


def reconnect_to_ncbi_server(f,curdir,email,currrootdir='/pub/geo/DATA/MiniMl'):
    """
    tries to reconnect in case NCBI ftp is down or something..
    """
    try:
        f.close()
    except Exception, e:
        print e
        pass


    #ftp address for GEO datasets:
    GEO_home_address="ftp.ncbi.nih.gov"#/pub/geo/DATA/"
    print 'Reconnecting to ftp [%s] as (%s) ...' % (GEO_home_address, email)
    #connect to the FTP-
    try:
        f=ftplib.FTP(GEO_home_address,"anonymous",email)
        print 'Connected!'
    except Exception, e:
        print e
        return -1

    #Now we want to go to the GEO dir. we don't want to import the whole NCBI
    #... or do we?
    print curdir

    ftp_dir=curdir
    print 'going to '+ ftp_dir
    for dir in currootdir.split('/').pop()+curdir.split('/'):
        if (dir != ''):
            print 'going to %s' % dir
            f.cwd(dir)

    return f

def append_xml_branch(subroot,f, dir, email,depth=0,currrootdir='/pub/geo/DATA/MiniMl'):
    """
    Crawl the ftp server from the 'dir' dicrectory and goes everywhere.
    """
    import time
    import ftplib
    from lxml import etree


    interval_between_requests=.1 #in seconds..
    verbosity=1

    verbose = verbosity


    f.set_debuglevel(verbosity)

    char = ['']
    for i in range(int(depth)):
        tmpnode=subroot.getparent()
        char.append(tmpnode.tag)
#       print char
    char.reverse()
    tmp_dir_char=''
    for i in char:
        tmp_dir_char+=i+'/'
#   print 'finished char' + str(char)

    try:
        f.cwd(dir)
        #PWD on locla server...
        rpwd=f.pwd()
    except Exception, e:
        print e
        print 'Waiting 2 minutes...'
        time.sleep(20)
        no_connection=1
        while no_connection:
            print 'Trying to reconnect'
            try:
                time.sleep(200)
                print 2
                f=reconnect_to_ncbi_server(f,tmp_dir_char,email,currrootdir)
                if f is not -1:
                    no_connection=0
            except Exception, e:
                print e
                no_connection=1
                pass


#   spacechar=''
#   for i in range(depth-1):
#       spacechar+='--#'
#   print spacechar + dir
    subdirs = []
    listing = []
    sizes   = []
    if verbose > 0 : print 'Listing remote directory %r...' % (dir,)

    try:
        f.retrlines('LIST', listing.append)
    except Exception, e:
        print e
        print 'Waiting 2 minutes...'
        time.sleep(2)
        no_connection=1
        while no_connection:
            print 'Trying to reconnect'
            try:
                time.sleep(2)
                f=reconnect_to_ncbi_server(f,tmp_dir_char,email,currootdir)
                no_connection=0
            except Exception, e:
                print e
                no_connection=1
                pass


    filesfound = []
    for line in listing:
        if verbose > 1: print '-->', repr(line)

        words = line.split(None, 8)
        if len(words) < 6:
            if verbose > 1: print 'Skipping short line'
            continue
        filename = words[-1].lstrip()
        infostuff = words[-5:-1]
        mode = words[0]
#       if verbose: print 'Mode ',mode, 'info', infostuff, 'filename' , filename
        if mode[0] == 'd': # filename  is a directory
            if verbose > 1:
                print 'Remembering subdirectory', repr(filename)
            subdirs.append(filename)
            continue
        if filename[0] != '.': # filename  is not a hidden file/directory, we try to ignore thoses..
            filesfound.append(str(filename))
            sizes.append(pretty_filesize(int(infostuff[0])))

############## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#####################
            try: #file names can't have numbers in the start of their file names to be listed in the xml..!!!
                etree.SubElement(subroot, str(filename))
                attrib=subroot[-1].attrib
                attrib['Size']=sizes[-1]
            except Exception,e:
                print e
            continue
############## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#####################

    if verbose > 0 and filesfound:
        print 'FILES:'
        for files in filesfound:
            print files, sizes[filesfound.index(str(files))]
    if verbose > 0 and subdirs:
        print 'DIRECTORY:'+ str(subdirs)

    for subdir in subdirs:
        try:
#           nobody want to make a ddos attack...
            time.sleep(interval_between_requests)

            subsubdir = etree.SubElement(subroot, str(subdir))
#           print(etree.tostring(subsubdir, pretty_print=True))
            append_xml_branch(subsubdir,f,subdir,email,depth+1)
        except Exception,e:
            print e
#           print(etree.tostring(subroot, pretty_print=True))
            print 'Quelque chose de poissoneux vient de se passer.....'
            print 'waiting 2 minutes...'
            time.sleep(2)
            no_connection=1
            while no_connection:
                print 'trying to reconnect'
                try:
                    f=reconnect_to_ncbi_server(f,subdir,email)
                    no_connection=0
                except Exception, e:
                        print e
                        no_connection=1
                        pass
            continue
    #print 'going back from .. ' + dir
    f.cwd('..')

def make_GEOftp_xml_tree(dir,email,ftp_dir='/pub/geo/DATA/MINiML',outfile='out.xml'):
    """
    List all directories and files from ftp.ncbi.nih.gov/pub/geo/DATA/ and co
    and generate an xml file with that information.
    this will be used to access quickly to any data on the site.
    """

    #this is the xml ogject that contains infos about the ftp site.
    root = etree.Element(ftp_dir.split('/')[-1])

    # Define the local directory name to put data in
    curdir  = os.curdir
    ddir    = curdir+dir

    # If directory doesn't exist make it
    if not os.path.isdir(ddir):
        os.mkdir(ddir)

    # go to the local directory where you want to put the data
    os.chdir(ddir)
    if 1: #ugly hack!comme il est vilain!
        #ftp address for GEO datasets:
        GEO_home_address="ftp.ncbi.nih.gov"#/pub/geo/DATA/"
        print 'Connecting to ftp [%s] as (%s) ...' % (GEO_home_address, email)
        #connect to the FTP-
        f=ftplib.FTP(GEO_home_address,"anonymous",email)
        print 'Connected!'

        #Now we want to go to the GEO dir. we don't want to import the whole NCBI
        #... or do we?
#       ftp_dir='/pub/geo/DATA/'
        print 'going to '+ ftp_dir
        for dir in ftp_dir.split("/"):
            if (dir != ''):
                print 'going to %s' % dir
                f.cwd(dir)
    dir = ''
    print 'listing all in '+ dir
    append_xml_branch(root,f,dir,email,currrootdir=ftp_dir)

#   open a file
    out = open(outfile,'w')
#   dump the xml into it..
    out.write("%s\t" % (etree.tostring(root, encoding='utf-8',pretty_print=1)))
#   done
    f.close()


def decode_xml(html_string):
    converted = UnicodeDammit(html_string,  smartQuotesTo='xml',isHTML=0)
    if not converted.unicode:
        raise unicodeDecodeError(
            "Failed to detect encoding, tried [%s]",
             ', '.join(converted.triedEncodings))
# print converted.originalEncoding
    return converted

def print_all(key,root):
    """Print various infos from the XML file"""
    if type(root).__name__ != '_Element':
        print "Expecting an etree Element .."
        print "was given  %s" % type(root).__name__
        raise

    if key=='Description' or key=='description':
        for node in root.findall('*//'+MINiML_header+str(key)):
            try:
                tag = str(node.getparent().tag).split(MINiML_header)[-1]
                txt = node.getparent().text
                print tag , txt , node.text
            except Exception, e:
                print e
            pass
    else:
        if len(root) == 0:
            print 'YEAH'
            print('\n'.join([str(e.getparent().tag).split(MINiML_header)[-1] + ':\t'+ str(e.text) for e in root.findall('*//'+MINiML_header+str(key))]))
        else:
            #print '%s' % str(root.getparent().tag).split(MINiML_header)[-1]
            for node in root.getchildren():
                print '    ' + node.tag.split(MINiML_header)[-1]
                print_all(node.tag.split(MINiML_header)[-1] , node)


def filesize(size,unit):
    """docstring for filesize:

        produce a pretty print
        from a size in byte.
    """
    if unit == 'bytes':
        return int(size)
    elif unit == 'KB':
        return int(size)*1024
    elif unit == 'MB':
        return int(size)*1024*1024
    elif unit == 'GB':
        return int(size)*1024*1024*1024

def pretty_filesize(bytes):
    if bytes >= 1073741824:
        return str(bytes / 1024 / 1024 / 1024) + ' GB'
    elif bytes >= 1048576:
        return str(bytes / 1024 / 1024) + ' MB'
    elif bytes >= 1024:
        return str(bytes / 1024) + ' KB'
    elif bytes < 1024:
        return str(bytes) + ' bytes'

def cool():
    """docstring for cool"""
    print 'cool'
    pass

def get_platform_database_size(infile):
    """get database size. from xml generated by maketree() """
    parser = etree.XMLParser(recover=1,encoding='utf-8')
    tree = etree.parse(infile,parser)
    root = tree.getroot()
    size=0
    for e in root:
        for n in e:
            a=n.values()
            a=a[0].split()
            print a
            size+=filesize(a[0],a[1])
    return pretty_filesize(size)



def guess_data_structure(filename,  col_nb=0, records=10, text_option=0,separator=''):
    """ Read a file with an known number of columns.
        the number of column comes from the xml.
        then try to determine for each column the type of data (number or string)
        If 0 is given as col_nb, will try to get the number of col itself...
    """
    import numpy as N

    if separator is '':
        separator_list=['\t',' ',',',';','    ','/','MUHAHA']
    else:
        separator_list=separator

    data=[]
    datatype=[]
    compute=1

    if col_nb != 0: # we know the col number..
#       print 'getting file records..'
        for separator in separator_list:
            data=[]
            nb_read_record=[]
            NR=0 #awk style line numbering ref.
            f=open(filename, 'r')
            for line in f:
                fields = map(str,line.strip().split(separator))
#               print fields
#               print 'line %s:'%NR, len(fields), fields
                nb_read_record.append(len(fields))
                #add the data, making sure biologists did not put white space before their text entries
                if len(fields) != records:
                    for i in range(len(fields) , col_nb):
                        fields.append('')

                data.append(fields)
#               print len(fields), '## ', fields

                NR+=1
                if NR > records: #stops at the 10th line, that's a very small sample of the file
                    break
            f.close()
            if max(nb_read_record) == col_nb:
#               print 'done #%s#' %separator
                break
            if separator==separator_list[-1]:
                print 'Can\'t guess data structure from all separator types...'
                compute = 0
                datatype=None
                break

    else: # we don't know the col number..
        abort=0
        for separator in separator_list:
            for col_nb in range(900): #thats already a lot of column...
            #print 'getting file records..'
                data=[]
                nb_read_record=[]
                NR=0 #awk style line numbering ref.
                f=open(filename, 'r')
                for line in f:
                    fields = line.strip().split(separator)
#                   print 'line %s:'%NR, len(fields), col_nb
#                   print fields, len(fields)
                    nb_read_record.append(len(fields))
                    data.append(fields)
                    NR+=1
                    if NR > records: #stops at the 10th line, that's a very small sample of the file
                        break
                f.close()
                if (N.mean(nb_read_record) == col_nb) and (len(fields[0]) < 10):
#                   print 'done #%s# %d' %(separator,N.mean(nb_read_record))
#                   print fields
                    abort = 1
                    break
                if separator==separator_list[-1]:
                    print 'Can\'t guess data structure from all separator types...'
                    compute = 0
                    datatype=None
                    break
            if abort:
                break

    if compute:
        #now get the type of data in there...
        datatbl=N.array(data[0:records-1][0:col_nb-1]).transpose() #to work colum by column
        #clear data to release a bit of memory
        data=[]
        data=None
        for i, line in enumerate(datatbl):
#           print line
            max_str_len=0
            if text_option:
                for j, entry in enumerate(line):
                    tmp=len(entry)
                    if tmp > max_str_len:
                        max_str_len=tmp
                strl='a'+str((max_str_len)*5+1) #makes a lot of room for the data in case..
                datatype.append(strl)
            else:
                u={'a15':0, '<f4':0}
                for j, entry in enumerate(line):
#                   print i , j, '!', entry, '!', line
#                   print i, j
                    try:
                        a=float(entry)
                        u['<f4']+=1
                    except:
                        u['a15']+=1
                        tmp=len(entry)
                        if tmp > max_str_len:
                            max_str_len=tmp
#                       print 'str length %d, max %d, %s' %(tmp,max_str_len,entry)
                        pass
                #ugly ... (here we assumethat the highest count is the datatype)
                if u['<f4']<u['a15']:
                    strl='a'+str(max_str_len+1)
#                   print strl
                    datatype.append(strl)
                else:
                    datatype.append('<f8')
#               print datatype
        return N.dtype(','.join(datatype))
    else:
        return None

def read_array(filename, dtype, separator='\t',text_option=0,rows=999999):
    """ Read a file with an arbitrary number of columns.
The type of data in each column is arbitrary
It will be cast to the given dtype at runtime
    """
    import numpy as N
    cast = N.cast
    data = [[] for dummy in xrange(len(dtype))]
    f=open(filename, 'r')
    j=0
    for line in f.readlines():
        j+=1
        if j>rows:
            print "MAX  ROW!"
            break
#       print '#'+line+'#'
        fields = line.strip().split(separator)
        if len(fields) < len(dtype): #for missing data
            for k in range(len(fields) , len(dtype)):
                fields.append('')
        for i, number in enumerate(fields):
            if text_option:
#               number=str(number)
                if number=='':
#                   print "33333"
                    number='Empty'
#               print number
#               print fields
            else:
                try:
#                   print "!!!!"

                    float(number)

                except:
                    if number=='':
#                       print "########"
                        number=N.nan
#                   number=str(number)
                    pass
            data[i].append(number)
    f.close()
    for i in xrange(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
#   for i in range(len(data[0])):
#       print str(data[0][i])+' \t '+str(data[1][i])+' \t '+str(data[2][i])+' \t '+str(data[3][i])+' \t '+str(data[4][i])
    return N.rec.array(data, dtype=dtype)



# h5py.File * read_xml_contributor(etree.Element * node){}
def read_xml_node(node):
    """node is a branch of a MINiML xml file : node argument is an 'Element tree' from lxml
        the idea is to create a hdf5 file in memory, return it, and then append it to a
        bigger file. the hdf5 file mirrors the xml.
    """

    print 'Parsing',node.tag.split(MINiML_header)[-1], 'information.'
    contrib_ghost_file = h5py.File(tmp_file_path+'/CT.hdf5', mode='w') #the 'core' driver means 'file in RAM'
    #find iid
#   Entry_id=node.get('iid')
    Entry_id=node.tag.split(MINiML_header)[-1]

    print Entry_id
    #create a dir in the ghost file with that iid name
    if not Entry_id in contrib_ghost_file:
        sg=contrib_ghost_file.create_group(Entry_id)
    else:
        sg=contrib_ghost_file[Entry_id]
#############################
    mirror_xml(node,sg,compression_type=6)     # The magic happens here.
#############################

    return contrib_ghost_file


def mirror_xml( node , sg, compression_type='lzf'):
    """docstring for mirror_xml
        get a node from lxml, go through children, make appropriate groups/datasets in hdf5 file
        inputs:
            node    : lxml etree Element.
            sg      : pointer to h5py file group.
    """
#   compression_type='gzip'
    #put attributes into description
    for i in node.attrib:
#       print node.attrib
        sg.attrs[i]=node.attrib[i]

    dup_key_index={}
    for subnode in node.iterchildren():

        ######################################
        #subnode with children become groups,
        #subnodes without any become datasets.
        #######################################

        if (subnode.getchildren() == []): #end of the Xpath, this is therefore a group
            #check atrributes
            for i in subnode.attrib:
                print subnode.attrib, i
                sg.attrs[i]=subnode.attrib[i]

            key=subnode.tag.split(MINiML_header)[-1]
#           print '     #', subnode.sourceline
            try:
                val=str(subnode.text).strip().encode('ascii', 'ignore')
            except Exception,e:
                print e
                val=subnode.text.encode('ascii', 'ignore')
                print '--> converted successfuly!'
#               print val
                pass
            lval=len(val)
            if (lval == 0) or (val == 'None'):
#               print 'Empty Field'
                if key.lower()=='sample-ref':
                    val=subnode.attrib['ref']
#                   print val
                else:
                    val='Empty Field'
                lval=len(val)
            print key,val,'*'

############ Special case: data-table, where there is  atext file to load-
            if key.lower()=='external-data': # data file!!
                print 'NOT NOT NOT Importing data-table:'
                break
                nb_col=len(node.findall(MINiML_header+'Column'))
                datatbl=subnode.text.strip().encode('ascii', 'ignore')
#               platform_id=node.getparent().getparent().find(MINiML_header+'Platform').attrib['iid']
                platform_id=node.getparent().getparent().find(MINiML_header+'Series').attrib['iid']
                print 'the data %s as %d cols' % (datatbl,nb_col)
                print GEO_data_dl_dir+platform_id+'_family/'  + datatbl
#               print 'gessing file structure from XML...'
                datas=None
                datatype=None
                try:
                    if datatbl.find('GPL') != -1:

                        print 'Platform description:', nb_col
                        datatype = guess_data_structure(GEO_data_dl_dir+platform_id+'_family/'  + datatbl,nb_col,5000,text_option=1)
                        print 'Succesfully guessed data structure...'
                    else:
                        datatype = guess_data_structure(GEO_data_dl_dir+platform_id+'_family/'  + datatbl,nb_col)
                except Exception, e:
                    print e,
                    print 'Impossible to guess data structure...'
                    pass
#               print 'geting file'
                if datatype is not None:
                    try:
                        if datatbl.find('GPL') != -1:
                            if subnode.attrib['rows']:
                                lignes=subnode.attrib['rows']
                            else:
                                lignes=999999
                            datas = read_array(GEO_data_dl_dir+ platform_id+'_family/' + datatbl,datatype,text_option=1,rows=lignes)
                        else:
                            datas = read_array(GEO_data_dl_dir+ platform_id+'_family/' + datatbl,datatype)
                            print datas
#                       print 'importing into hdf5'
                        try:
                            dset = sg.create_dataset(datatbl.split('.')[0], data=datas, compression=compression_type)
                        except Exception,e:
                            print e
                            pass
                    except Exception,e:
                        print e
                        pass
                else:
                    print 'trying to guess the data structure now...'
                    try:
                        datatype = guess_data_structure(GEO_data_dl_dir+platform_id+'_family/' + datatbl)
                    except Exception,e:
                        print e
                        print 'putting raw file in DB'
                        print datatbl
                        print GEO_data_dl_dir+platform_id+'_family/'  + val.split('/')[-1]
                        datas=open(GEO_data_dl_dir+platform_id+'_family/'+datatbl , 'r').read()
                        ######
                        pass
####### End of datatable loading. ####################################
            else:
                if not(key in sg):
                    sg.create_dataset(key, data=numpy.array([val],dtype='a'+str(lval)), compression=compression_type)
                    dup_key_index[key]=0
                else:
                    dup_key_index[key]+=1
                    sg.create_dataset(key+'_'+str(dup_key_index[key]), data=numpy.array([val],dtype='a'+str(lval)), compression=compression_type)
        else:
            key=subnode.tag.split(MINiML_header)[-1] # subnode has children, so let's go recursive!
            for i in subnode.attrib: #replace key with iid in xml instead of kind..
                print subnode.attrib, i
                if 'iid' in subnode.attrib.keys():
                    key=subnode.attrib['iid']
#           print key,'>'
            if not (key in sg): #and this key is not in the group yet...
                subsg=sg.create_group(key)
                dup_key_index[key]=0
            else: #we put a number..
                dup_key_index[key]+=1
                subsg=sg.create_group(key+'_'+str(dup_key_index[key]))
            mirror_xml(subnode,subsg)
    pass




def read_xml_empty(node):
    """node is unknown"""
    print 'Unknown Node'

def converter(x):
    if x == 'N/A':
        return numpy.nan
    else:
        return float(x)
    pass

def parse_node(node, node_type=''):
    print 'Parsing a %s' % node_type
#   import_func = getattr(MINiML, 'read_xml_%s' % node_type, read_xml_empty)
    try:
        return read_xml_node(node)
    except Exception,e:
        print e
        read_xml_empty(node)

    #   print ["%s" % s.tag.split('{http://www.ncbi.nlm.nih.gov/projects/geo/info/MINiML}')[-1] for s in node.getchildren() ]
#   return import_func(node) the one with getattr


def parse_GEO_xml(infile):
    #open dbfile
#   db=h5py.File(tmp_file_path,'a')

    try:
#       parser = etree.XMLParser(recover=1,encoding='utf-8')
#       tree = etree.parse(infile,parser)
        root = tree.getroot()
    except (IOError, OSError,Exception):
        print 'Now that\'s one ugly XML you got there..'
        print 'saving reformated xml...'
        xmlfile=open(infile,'r')
        print 'file opened %s' % infile
        texte=xmlfile.read()
        print 'text read...'
        converted=decode_xml(texte)
        print 'text converted'
        parser = etree.XMLParser(huge_tree=1,recover=1,encoding=converted.originalEncoding)
        root = etree.XML(texte,parser)
        print 'xml regenerated'
        texte=None
        converted=None
        # move funny xml somewhere else.
        shutil.move(infile, infile+'_orig')
        f=open(infile,'w+')
        f.write(unicode(etree.tostring(root,pretty_print=1)))
        f.close()
        pass

    print 'Loaded %s entries' % len(root)

    try:
        read_xml_node(root)
#       shutil.rmtree('../CT.hdf5')
#       shutil.move('CT.hdf5','..')
    except Exception,e :
        print e

#   for i , node in enumerate(root):
#       entry_type=(node.tag.split(MINiML_header))[-1]
#
#       try:
#           tmpf=parse_node(node,entry_type.lower())
#           if (tmpf is not None) and (list(tmpf)[0] not in db):
#               tmpf.copy(list(tmpf)[0], db)
#               tmpf.close()
#       except KeyError:
#           if list(tmpf)[0] in db:
#               print 'Record already present'
#               tmpf.close()
#           else:
#               print 'Something else'
#               tmpf.close()

def untar(filename, destination):
    import tarfile
    f = tarfile.open(filename, mode='r:gz')
    f.extractall(destination)
    f.close()



def check_fs(dir):
    """docstring for check_fs:
    this function creates all directories needed to have GEO in your computer.
    """
    #check if GEO dirs are in the home dir
    home = os.environ['HOME']
    #home = '/Volumes/Donnees/'
    os.chdir(home)
    if not os.path.isdir(home+'/GEO'):
        os.mkdir(home+'/GEO')
    os.chdir(home+'/GEO')
    #then make the other
    dir= dir.split(home)[-1]
    for subdir in dir.split('/'):
        if subdir is not '':
            if not os.path.isdir(subdir):
                os.mkdir(subdir)
    pass


def get_GEO_file(params,mail):
    """
        get file in params from ftp, decompress it and add it in database.
    """

    #dataset to be retrieved:
    dset = params

    # Define the local directory name to put data in
    curdir= os.curdir
    ddir = GEO_data_dl_dir

    # If directory doesn't exist make it
    if not os.path.isdir(ddir):
        check_fs(ddir)
    # Change the local directory to where you want to put the data
    os.chdir(ddir)


    #ftp address for GOE datasets:
    GEO_home_address="ftp.ncbi.nih.gov"#"/pub/geo/DATA/"
    print 'Connecting to ftp [%s] as (%s) ...' % (GEO_home_address, mail)
    #connect to the FTP-
    f=ftplib.FTP(GEO_home_address,"anonymous",mail)
    print 'Connected!'

    #Return the welcome message sent by the server in reply to the initial connection.
    #(This message sometimes contains disclaimers or help information that may be relevant to the user.)
#   print f.getwelcome()
    #go to the specific directory:
    #we do this by going in each directory sequentialy.

    if dset[1] == 'S':
        ftp_dir='/pub/geo/DATA/MINiML/by_series/'+dset+'/'
    if dset[1] == 'P':
        ftp_dir='/pub/geo/DATA/MINiML/by_platform/'+dset+'/'

    try:
        for dir in ftp_dir.split("/"):
            if (dir != ''):
                print 'going to %s' % dir
                not_passed=1
                while not_passed:
                    try:
                        f.cwd(dir)
                        not_passed=0
                    except Exception,e:
                        print e
                        print 'ftp dir', f.pwd()
                        print 'going to -> ', dir
                        time.sleep(10)
    except Exception, e:
        print e
    pass


    # change the remote directory
#   print 'Changing directory to %s' % ftp_dir
#   f.cwd('ftp_dir')
    print 'List of files in the FTP'
    print f.nlst()
    # define filename
    file_list=f.nlst()

    for file_names in  file_list:
        if not os.path.isfile(file_names):
            try:
                print 'getting file %s' % file_names
#               print pretty_filesize( f.size(file_names))
                # get the remote file to the local directory
                f.retrbinary("RETR %s" % file_names, open(file_names,"wb").write)
            except Exception, e:
                print e

    for file_names in  file_list:
        entryname=file_names.split('.')[0]
        print entryname
        if  (entryname != "GPL570_family"):
            print 'Uncompressing file...', file_names
            print os.getcwd()
            untar(file_names,entryname)

    print 'Adding %s to local database...' % entryname
    os.chdir(file_names.split('.')[0])
    print 'please wait...'
    parse_GEO_xml(entryname+'.xml')
    print 'cleaning'
#   print 'To be done.'
        #file=open(file_names,"wb")
        #f.retrbinary("RETR %s" % file_names,   handleDownload)

    # Close FTP connection
    f.close()
    os.chdir('..')
