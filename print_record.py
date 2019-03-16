import sys
import vcf
import json
import gzip


def print_gnomAD(chrm,  pos):
    gnomad = '/data/exp/trifon/vault/xl_BGM0187/fdata.json.gz'
    recs = []
    with gzip.open(gnomad, "rb") as inp:
        for line in inp:
            rec_data = json.loads(line)
            if str(rec_data['Start_Pos']) == pos and rec_data['Chromosome'] == chrm:
                recs.append(rec_data)
    inp.close()
    return recs

def print_record(record):
    rec_dict = record.__dict__
    rec_dict['aaf'] = record.aaf
    rec_dict['heterozygosity'] = record.heterozygosity
    samples = []
    hets = []
    for call in rec_dict['samples']:
        sample = call.__str__()
        hets.append(call.is_het)
        samples.append(sample)
    rec_dict['hets'] = hets
    rec_dict['samples'] = samples
    rec_dict['alleles'] = rec_dict['alleles'].__str__()
    rec_dict['ALT'] = rec_dict['ALT'].__str__()
    #res = json.dumps(rec_dict,  indent=4)
    return rec_dict


def print_to_file(data,  file_name):
    try:
        data_file = open(file_name, 'w')
    except IOError:
        print('File "' + file_name + '" not found.')
        return
    data_file.write(json.dumps(data,  indent=4))
    data_file.close()


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
        print_to_file(res,  record_file_name)
        vcf_file.close()
        
        freq = print_gnomAD(chm,  pos)
        if freq is None:
            print('Record in GnomAD for this variant is not found.')
        else:
            gnomAD_file_name = '/home/andrey/work/Caller/caller/case_187/gnomAD.json'
            print_to_file(freq,  gnomAD_file_name)
        
        print('Ok.')
        sys.exit()
    vcf_file.close()
    print('record CHROM=' + chm + ', POS=' + pos + 'not found.')
