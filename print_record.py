import sys
import vcf
import json

def print_record(record):
    rec_dict = record.__dict__
    rec_dict['aaf'] = record.aaf()
    samples = []
    for call in rec_dict['samples']:
        sample = call.__str__()
        samples.append(sample)
#        for key in rec_dict:
#            rec_dict[key] = rec_dict[key].__str__()
    rec_dict['samples'] = samples
    rec_dict['alleles'] = rec_dict['alleles'].__str__()
    rec_dict['ALT'] = rec_dict['ALT'].__str__()
    res = json.dumps(rec_dict,  indent=4)
    return res




if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Use following format: python print_record.py CHROM POS')
        sys.exit()
    chm = sys.argv[1]
    pos = sys.argv[2]
    
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    try:
        vcf_file = open(vcf_file_name, 'r')
    except IOError:
        print('File "' + vcf_file_name + '" not found.')
        sys.exit()
    
    vcf_reader = vcf.Reader(vcf_file)
    for record in vcf_reader:
        if record.CHROM != chm or str(record.POS) != pos:
            continue
        res = print_record(record)
        record_file_name = '/home/andrey/work/Caller/caller/case_187/record.json'
        record_file = open(record_file_name,  'w')
        record_file.write(res)
        record_file.close()
        vcf_file.close()
        print('Ok.')
        sys.exit()
    vcf_file.close()
    print('record CHROM=' + chm + ', POS=' + pos + 'not found.')
