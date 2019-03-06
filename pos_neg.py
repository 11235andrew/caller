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
    max_qual = 0
    filters = []
    gq = 0
    qd = 0
    fs = 0
    freq = 0
    pass_f = 0
    for record in vcf_reader:
        count += 1
        if count % 1000 == 0:
            print('Record #' + str(count))
        if record.CHROM == 'chrM':
            continue
        flag = False
        for fr in record.INFO['AF']:
            if fr < 0.01:
                flag = True
        if not flag:
            continue
        freq += 1
        for flt in record.FILTER:
            if flt not in filters:
                filters.append(flt)
        if record.FILTER != []:
            continue
        pass_f += 1
        #if record.QUAL < 500000:
        #    continue
        if record.QUAL > max_qual:
            max_qual = record.QUAL
        
#        if 'GQ_MEAN' in record.INFO and record.INFO['GQ_MEAN'] > 20:
#            gq += 1
#        else:
#            continue
        if 'QD' in record.INFO and record.INFO['QD'] > 4:
            qd += 1
        else:
            continue
        if 'FS' in record.INFO and record.INFO['FS'] < 30:
            fs += 1
        else:
            continue
        
        rec = {}
        rec['CHROM'] = record.CHROM
        rec['POS'] = record.POS
        rec['AF'] = str(record.INFO['AF'])
        rec['QUAL'] = record.QUAL
        rec['owns'] = []
        for sample in record.samples:
            au = sample.sample[-2]
            GQ = sample.data.GQ
            if au == 'u':
                AD = sample.data.AD
                if AD[1] > 0 and GQ > 20:
                    rec['owns'].append(sample.sample + '(' + str(AD[0]) + ',' + str(AD[1]) + ')')
            else:
                GT = sample.data.GT
                if GT in ['0/1',  '0/1', '1/1'] and GQ > 20:
                    rec['owns'].append(sample.sample  + '(' + GT + ')')
        if rec['owns'] !=  []:
            pos_neg.append(str(rec))
#    print(str(len(f_neg)) + ' false negative records were found.')
#    print(str(len(f_pos)) + ' false positive records were found.')
    #print(str(len(intersection(f_neg,  f_pos))))
    print('Maximal Quality is ' + str(max_qual))
    print('Allele freqency < 0.01: ' + str(freq) + ' variants.')
    print('Quality_PASS: ' + str(pass_f) + ' variants.')
    print('GQ_MEAN > 20: ' + str(gq) + ' variants.')
    print('QD > 4: ' + str(qd) + ' variants.')
    print('FS < 30: ' + str(fs) + ' variants.')
    print(str(len(pos_neg)) + ' variants were found.')
    
    vcf_file.close()
    
    
    print_to_file(pos_neg,  cand_file_name)
    print_to_file(filters,  'case_187/filters.json')
    
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
    

def print_to_file(data,  file_name):
    try:
        data_file = open(file_name, 'w')
    except IOError:
        print('File "' + file_name + '" not found.')
        return
    data_file.write(json.dumps(data,  indent=4))
    data_file.close()


if __name__ == '__main__':
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    f_negative_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative.json'
    f_positive_file_name = '/home/andrey/work/Caller/caller/case_187/false_positive.json'
    candidats_file_name = '/home/andrey/work/Caller/caller/case_187/candidats.json'
    find_false_negative(vcf_file_name,  candidats_file_name)
    print('Ok.')
