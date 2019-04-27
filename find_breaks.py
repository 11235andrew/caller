import vcf
#import numpy as np
#import matplotlib.pyplot as plt
from print_record import open_file
from print_record import print_to_file


#def func(x):
#    return x*x*x +3* x*x -2*x+1

def find_breaks(vcf_file_name):
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    all_sign = {}
    signature = {}
    chrm = None
    step = 500000
    for record in vcf_reader:
        if record.CHROM in ['chrM',  'chrX',  'chrY']:
            continue
        
        if chrm != record.CHROM:
            if chrm is not None:
                all_sign[chrm] = signature
            homo = {}
            for sample in record.samples:
                signature[sample.sample] = []
                homo[sample.sample] = 0
            chrm = record.CHROM
            segment = 0
        
        if int(record.POS / step) > segment:
            segment = int(record.POS / step)
            print('Chromosom: ' + chrm + ', Segment: ' + str(segment))
            for sample in record.samples:
                key = sample.sample
                signature[key].append(homo[key])
                homo[key] = 0
        
        for sample in record.samples:
            if sample.gt_type == 2:
                homo[sample.sample] += 1
    res_file_name = '/home/andrey/work/Caller/caller/case_187/signature.json'
    print_to_file(sign_all,  res_file_name)
    
#    X = np.linspace(-2.5,  2.5,  50)
#    Y = func(X)
#    print('X: ' + str(X))
#    print('Y: ' + str(Y))
#    plt.plot(X,  Y)
    
















if __name__ == '__main__':
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    find_breaks(vcf_file_name)
    print('Ok')
