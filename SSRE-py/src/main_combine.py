import os, argparse, getopt
from os.path import join, abspath, dirname
import numpy as np
# from visualization import tsne_bo
import sys

sys.path.append(dirname(abspath(__file__)))

from src.SSRE import SSRE_combine
from src.tsne_l import tsne_l


# def parse_args():
#     '''
#     Parses the SSRE arguments.
#     '''
#     parser = argparse.ArgumentParser(description="Run SSRE")
#
#     parser.add_argument('--path', nargs='?',
#                         help='Input path')
#
#     parser.add_argument('--cluster_num', type=int,
#                         help='Number of cluster.')
#     return parser.parse_args()


def usage():
    help_msg = "\n  usage: main.py -d <file_data> -l <file_label> -p <parameter ρ⁄λ> [options]\n\
  -d,--data\tinput file path and name of data\n\
  -l,--label\tinput file path and name of label\n\
  -p,--param\tparameter ρ⁄λ\n\
  -h,--help\tprint help message.\n"
    print(help_msg)
    sys.exit()


def main(argv):
    filedata = None
    filelabel = None
    param = 0
    try:
        opts, args = getopt.getopt(argv, "d:l:p:h",
                                   ["data=","label=", "param=","help"])
    except getopt.GetoptError as err:
            print(str(err))
            usage()
            sys.exit()

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-d', '--data'):
            filedata = arg
        elif opt in ('-l', '--label'):
            filelabel = arg
        elif opt in ('-p', '--param'):
            param = int(arg)
        else:
            print("Error in arguments: unknown option\n")
            usage()
            sys.exit()

    if (filedata is None):
        sys.stderr.write("Error in arguments: must specify -d\n")
        usage()
    else:
        if not os.path.exists(filedata):
            sys.stderr.write("Error: " + filedata + " does not exists.\n")
            sys.exit()
        if (filelabel is None):
            sys.stderr.write("Error in arguments: must specify -l\n")
            usage()
        else:
            if not os.path.exists(filelabel):
                sys.stderr.write("Error: " + filelabel + " does not exists.\n")
                sys.exit()
    if param == 0:
        sys.stderr.write("Error in arguments: must specify -p\n")
        usage()
        sys.exit()
    if filedata != None and filelabel != None and param != 0:
        data = np.loadtxt(filedata,delimiter='\t')
        label = np.loadtxt(filelabel,delimiter='\t',dtype=int)
        [NMI,ARI] = SSRE_combine(data, label, param)
        print('Results: SSR+pearson, SSR+spearman, SSR+cosine, SSR+pearson+spearman, SSR+pearson+cosine, SSR+spearman+cosine')
        print(NMI,ARI)

if __name__ == "__main__":
    main(sys.argv[1:])
