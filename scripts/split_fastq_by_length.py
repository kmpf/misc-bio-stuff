#!/usr/bin/python
#split_fastq_by_length.py

import sys
import argparse

  


def read_arguments():
    parser = argparse.ArgumentParser(description="Takes a fastq file and writes new ones for each nucleotide length")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin, help ="Infile: default=stdin")

    parser.add_argument('--out-prefix', default='Sample', help ="Outfile: default=stdout")
    parser.add_argument('--out-filename_extension', default='fastq', help ="common are fastq, fq")

   
    return parser.parse_args()


def eval_arguments(args):
        return(args)


## called from main                     
def open_fh_specific_length(outstreams, length):
    filename = "".join(map(str, [args.out_prefix,'-', length,'bp', '.', args.out_filename_extension]))
    outstreams[length]  = open (filename, 'w')
    return (outstreams)

def parse_fastq_record(temp):
    fastq_record = dict()
    fastq_record['id'] = temp[0]
    fastq_record['seq'] = temp[1]
    fastq_record['qual_id'] = temp[2]
    fastq_record['qual'] = temp[3]
    fastq_record['seq_length'] = len(temp[1])
    return fastq_record

def write_fastq_record(outstreams, fastq_record):
    fh = outstreams[fastq_record['seq_length']]
    fh.write(fastq_record['id'] + '\n')
    fh.write(fastq_record['seq']+ '\n')
    fh.write(fastq_record['qual_id']+ '\n')
    fh.write(fastq_record['qual']+ '\n')


def main(args):
  
    # print(args)
    line_count = 0
    temp = []
    outstreams = dict()

    for counter, line in enumerate(args.infile):

        line_count += 1
        line = line.rstrip('\n')
        temp.append(line)
        number_fastq_records = (counter +1 ) / 4
        if line_count == 4:
            fastq_record = parse_fastq_record(temp)
            line_count = 0
            temp = []
            if not fastq_record['seq_length'] in outstreams:
                outstreams = open_fh_specific_length(outstreams, fastq_record['seq_length'])
            
            write_fastq_record(outstreams,fastq_record)

            
            if number_fastq_records  % 1e6 == 0:
                sys.stderr.write( " ".join(['fastq records processed:', str(number_fastq_records), '\n']) )
            


        

    sys.stderr.write( " ".join(['fastq records processed:', str(number_fastq_records), 'done!', '\n']) )
    for fh in outstreams:
        outstreams[fh].close()



if __name__ == '__main__':
    args = read_arguments()
    args = eval_arguments(args)
    main(args)

__version__ = '0.00000001'
