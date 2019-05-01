import vcf
import sys
from print_record import open_file
from print_record import print_to_file



def homo(sample):
    if sample.gt_type == 2:
        return 1
    else:
        return 0

def noise(sample):
    if sample.gt_type is None:
        return 1
    else:
        return 0

def find_breaks(vcf_file_name, res_file_name, func):
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    all_sign = {}
    signature = {}
    chrm = None
    step = 10000
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
            homo[sample.sample] += func(sample)
    all_sign[chrm] = signature
    print_to_file(all_sign,  res_file_name)













if __name__ == '__main__':
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    if sys.argv[1] == '--homo':
        res_file_name = '/home/andrey/work/Caller/caller/case_187/signature.json'
        find_breaks(vcf_file_name,  res_file_name,  homo)
    elif sys.argv[1] == '--noise':
        res_file_name = '/home/andrey/work/Caller/caller/case_187/noise.json'
        find_breaks(vcf_file_name,  res_file_name,  noise)
    else:
        print('Unknown key: ' + sys.argv[1])
    print('Ok')
