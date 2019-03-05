import vcf
import json

def print_record(record):
    rec_dict = record.__dict__
    samples = []
    for call in rec_dict['samples']:
        sample = call.__dict__
        samples.append(sample)
#        for key in rec_dict:
#            rec_dict[key] = rec_dict[key].__str__()
    rec_dict['samples'] = samples
    rec_dict['alleles'] = rec_dict['alleles'].__str__()
    rec_dict['ALT'] = rec_dict['ALT'].__str__()
    res = json.dumps(rec_dict,  indent=4)
    return res
    
    
def intersection(list1,  list2):
    res = []
    for el in list1:
        if el in list2:
            res.append(el)
    return res

def find_false_negative(vcf_file_name,  cand_file_name):
    try:
        vcf_file = open(vcf_file_name, 'r')
    except IOError:
        print('File "' + vcf_file_name + '" not found.')
        return
    
    vcf_reader = vcf.Reader(vcf_file)
    pos_neg = []
    count = 0
    for record in vcf_reader:
        count += 1
        if count % 100 == 0:
            print('Record #' + str(count))
        if record.CHROM == 'chrM':
            continue
        for fr in record.AF:
            if record.INFO.AF > 0.01:
                continue
        rec = {}
        rec['CHROM'] = record.CHROM
        rec['POS'] = record.POS
        rec['AF'] = ','.join(record.INFO.AF)
        rec['QUAL'] = record.QUAL
        rec['owns'] = []
        for sample in record.samples:
            au = sample.sample[-2]
            if au == 'u':
                AD = sample.data.AD
                if AD[1] > 0:
                    rec['owns'].append(sample.sample)
            else:
                GT = sample.data.GT
                if GT in ['0/1',  '0/1', '1/1']:
                    rec['owns'].append(sample.sample)
        if rec['owns'] !=  []:
            pos_neg.append(str(rec))
#    print(str(len(f_neg)) + ' false negative records were found.')
#    print(str(len(f_pos)) + ' false positive records were found.')
    #print(str(len(intersection(f_neg,  f_pos))))
    print(str(len(pos_neg)) + ' varients were found.')
    
    vcf_file.close()
    
    
    try:
        cand_file = open(cand_file_name, 'w')
    except IOError:
        print('File "' + cand_file_name + '" not found.')
        return
    cand_file.write(json.dumps(pos_neg,  indent=4))
    cand_file.close()
    
#    try:
#        f_negative_file = open(f_negative_file_name, 'w')
#    except IOError:
#        print('File "' + f_negative_file_name + '" not found.')
#        return
#    f_negative_file.write(json.dumps(f_neg,  indent=4))
#    f_negative_file.close()
#    
#    try:
#        f_positive_file = open(f_positive_file_name, 'w')
#    except IOError:
#        print('File "' + f_positive_file_name + '" not found.')
#        return
#    f_positive_file.write(json.dumps(f_pos,  indent=4))
#    f_positive_file.close()
    





if __name__ == '__main__':
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    f_negative_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative.json'
    f_positive_file_name = '/home/andrey/work/Caller/caller/case_187/false_positive.json'
    candidats_file_name = '/home/andrey/work/Caller/caller/case_187/candidats.json'
    find_false_negative(vcf_file_name,  candidats_file_name)
    print('Ok.')
