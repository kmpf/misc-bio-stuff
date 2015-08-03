#!/usr/bin/python
#collect_metrics_seqXTo2d.py 

import sys
import re
import argparse
from argparse import RawTextHelpFormatter
import logging
import pprint
import yaml
import json
from collections import OrderedDict 
import os
import numpy
#from Bio import SeqIO
import ast
pp = pprint.PrettyPrinter(indent=4, width=10)










def read_arguments():
    parser = argparse.ArgumentParser(description="parses log files from various SeqX programs and turns them into a 2d table")
    parser.add_argument('infiles', nargs='+', type=argparse.FileType('r'),
                    default=sys.stdin, help ="Infiles: default=stdin")

    parser.add_argument('--secondary-files', nargs='+', type=argparse.FileType('r'),
                    default=None, help ="for paired end with two output files")

    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stdout, help ="Outfile: default=stdout")

    parser.add_argument('--program', nargs='?', type=str, choices=["post_sawdust", "cutadapt", "htseq-count" ],
                    default=None, help ="which program generated the logfile")

    parser.add_argument('--program-version', nargs=1, type=str, 
                    default=any, help ="which program generated the logfile")

    parser.add_argument('--htseq-strip-dot-gencode', nargs='?', type=bool, 
                        default=True, help ="strip .dot form gencode identifier ")


    parser.add_argument('--step-id', nargs='?', type=str,  
                         default="None", help ="include a column with name ")

    parser.add_argument('--split', nargs='?', type=str, 
                    default=None, help ="split basename by SPLIT")

                        
    parser.add_argument('--logfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stderr, help ="logfile  default=stderr ")
       

    parser.add_argument('-d', '--debug',  help="Print lots of debugging statements",
                        action="store_const", dest="loglevel", const=logging.DEBUG,
                        default=logging.WARNING)


    parser.add_argument('-v', '--verbose',  help="Be verbose",
        action="store_const", dest="loglevel", const=logging.INFO    )


    parser.add_argument( '--version', action='version', version='%(prog)s ' + __version__)

    return parser.parse_args()


def eval_args(args):

    logging.basicConfig(stream=args.logfile,  level=args.loglevel)
    logger = logging.getLogger()

#    logger.debug('This message should go to the log file')
#    logger.info('So should this')
#    logger.warning('And this, too')

    args.fnames = []
    args.fnames_sec = []

    for f in args.infiles:
        basename = os.path.basename(f.name)
        if args.split:
            basename = basename.split(args.split)[0]
        args.fnames.append (basename)

    if args.secondary_files:
        for f in args.secondary_files:
            basename = os.path.basename(f.name)
            if args.split:
                basename = basename.split(args.split)[0]
                args.fnames_sec.append (basename)


    if args.program == "cutadapt":
        fname = "-".join([args.step_id, 'clipping.txt'])
        args.outfile2 = open(fname, "w")

    if logger.isEnabledFor(logging.DEBUG):
        for k, v in args.__dict__.iteritems():
            pp.pprint([k,v])

    return (args)
 


## called from main                     
def init_metrics():
    """Returns a dict like thingy for counting"""
    #new test add 
    metrics = []

    indict =     OrderedDict()



    indict['sub_multiple']['r2']  = 0
   
    
    indict['unique']  = 0
    indict['sub_unique']  = {}
    indict['sub_unique']['r1']  = 0
    indict['sub_unique']['r2']  = 0
    indict['sub_unique']['splits']  = 0
    indict['sub_unique']['split_types'] = {}

    for types in split_types:
        indict['sub_unique']['split_types'][types]  = 0


    outdict =     OrderedDict()
    outdict['templates']  = 0
    outdict['unmapped']  = 0
    outdict['multiple']  = 0
    outdict['sub_multiple']  = {}
    outdict['sub_multiple']['r1']  = 0
    outdict['sub_multiple']['r2']  = 0
    outdict['unique']  = 0
    outdict['sub_unique']  = {}
    outdict['sub_unique']['r1']  = 0
    outdict['sub_unique']['r2']  = 0
    outdict['sub_unique']['splits']  = 0

    metrics.append(indict)
    metrics.append(outdict)
    return metrics


class OrderedDictYAMLLoader(yaml.Loader):
    """
    A YAML loader that loads mappings into ordered dictionaries.
    https://gist.github.com/enaeseth/844388
    """
 
    def __init__(self, *args, **kwargs):
        yaml.Loader.__init__(self, *args, **kwargs)
 
        self.add_constructor(u'tag:yaml.org,2002:map', type(self).construct_yaml_map)
        self.add_constructor(u'tag:yaml.org,2002:omap', type(self).construct_yaml_map)
 
    def construct_yaml_map(self, node):
        data = OrderedDict()
        yield data
        value = self.construct_mapping(node)
        data.update(value)
 
    def construct_mapping(self, node, deep=False):
        if isinstance(node, yaml.MappingNode):
            self.flatten_mapping(node)
        else:
            raise yaml.constructor.ConstructorError(None, None,
                'expected a mapping node, but found %s' % node.id, node.start_mark)
 
        mapping = OrderedDict()
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            try:
                hash(key)
            except TypeError, exc:
                raise yaml.constructor.ConstructorError('while constructing a mapping',
                    node.start_mark, 'found unacceptable key (%s)' % exc, key_node.start_mark)
            value = self.construct_object(value_node, deep=deep)
            mapping[key] = value
        return mapping


def iterme(d, l, prefix):
    for k, v in d.iteritems():
        k1= "_".join([prefix,k ])
        if issubclass(dict,OrderedDict) or isinstance(v, dict):
            l = iterme(v,l, k1)
        else:
            l.append([k1,v])
    return(l)



def process_post_sawdust(args):
    for c, f   in enumerate(args.infiles):
           
        a = yaml.load(f, OrderedDictYAMLLoader)

        m = iterme(a[0],[], 'in')
        m.extend(iterme(a[1],[], 'out'))

        if c ==0:
            out = [str(i[0]).rstrip('\n') for i in m]
            out.insert(0, "id") 
            out_line= "\t".join(out) 
            ae = re.sub(r'_sub', '', out_line)
            args.outfile.write(ae  +'\n')        

        out = [str(i[1]).rstrip('\n') for i in m]        
        out.insert(0, args.fnames[c]) 
        out_line= "\t".join(out) 
        args.outfile.write(out_line  +'\n')


def process_cutadapt(args):
#    print (yaml.dump(args))
    _pc(args, "R1",args.infiles,args.fnames )
    _pc(args, "R2", args.secondary_files, args.fnames_sec)

def _pc(args, read_type, mylist, fnames):

    p         = re.compile('^Total reads processed:')
    pclipped  = re.compile('^Reads with adapters:')
    pnumber   = re.compile('\d+')
    plength   = re.compile('^length')

    
    for c, f   in enumerate(mylist):
        process = None
        for line in f:
            if  not line.strip():
                continue
            line = line.rstrip('\n') 
            if  p.match(line):
                tmp = re.sub(r',', '', line)
                res = pnumber.search(tmp)
                total_reads = res.group()
                out_line= "\t".join([fnames[c], read_type, total_reads, args.step_id ]) 
                args.outfile.write(out_line  +'\n')
                
            if pclipped.match(line):
                tmp = re.sub(r',', '', line)
                res = pnumber.search(tmp)
                clipped_reads = res.group()
                not_clipped_reads =str(int(total_reads) - int(clipped_reads))  
                zero_line= "\t".join([fnames[c], read_type,  args.step_id, '0', not_clipped_reads  ]) 


            if plength.match(line):
                process = True
                args.outfile2.write(zero_line  +'\n')
                continue


            if process:
                tmp = line.split("\t")
                out_line= "\t".join([fnames[c], read_type,  args.step_id,tmp[0], tmp[1]  ]) 
                args.outfile2.write(out_line  +'\n')
                        


def process_htseq_count(args):
    logger = logging.getLogger()
    #d->id->gene->count
    genes = dict()

    gene_ids = dict()
    if args.htseq_strip_dot_gencode:
        logger.warning("--htseq-strip-dot-gencode == TRUE")

    p_underscore   = re.compile('^__*')

    for c, f   in enumerate(args.infiles):
        genes[args.fnames[c]] = dict()
        for line in f:
            if  not line.strip():
                continue
                
            if p_underscore.match(line):
                continue

            line = line.rstrip('\n') 
            if args.htseq_strip_dot_gencode:
                line = re.sub(r'\.\d+', '', line)

            tmp  = line.split("\t")
            ids   = tmp[0]
            count = tmp[1]
            gene_ids[ids] = 1
            genes[args.fnames[c]][ids] = count


    
    header = []        
    header.append("id")
    header.extend(args.fnames)
    header_line= "\t".join(str(x) for x in header)
    args.outfile.write(header_line  +'\n')


    for gene_id in gene_ids.keys():
        out = []
        out.append(gene_id)
        for fname in args.fnames:
            try: 
                res = genes[fname][gene_id]
            except:
                res = 0 
            out.append(res)

        out_line= "\t".join(str(x) for x in out)
        args.outfile.write(out_line  +'\n')

def main(args):
    if args.program == "post_sawdust":
        process_post_sawdust(args)

    elif args.program == "cutadapt":
        process_cutadapt(args)

    elif args.program == "htseq-count":
        process_htseq_count(args)

    else:
        print "no progL:"





__version__ = '0.001'
if __name__ == '__main__':
    args = read_arguments()
    args = eval_args(args)
    main(args)


