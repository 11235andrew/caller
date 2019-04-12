import vcf
import random
from random import WichmannHill
from print_record import get_json_from_file
from print_record import print_to_file
from print_record import json_to_csv
from print_record import open_file
from print_record import get_frequency
from print_record import json_to_csv_columns


def is_unaffected(name):
    au = name[-2]
    return au == 'u'

def drop_out():
    f_neg_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative_main.json'
    f_neg = get_json_from_file(f_neg_file_name)
    short_f_neg = []
    freq_count = 0
    for rec in f_neg:
        if rec['gnomAD_AF'] > 0.01:
            continue
        freq_count += 1
        unaff_flag = True
        ind = rec['ALT_index'] + 1
        doute = 0
        for sample in rec['samples']:
            if is_unaffected(sample['sample']):
                if sample['AD'][ind] != 0:
                    unaff_flag = False
            else:
                if sample['GT'] not in ['0/1',  '1/0',  '1/1']:
                    doute += 1
        if not unaff_flag or doute > 1:
            continue
        
        short_f_neg.append(rec)
    short_f_neg_file_name = '/home/andrey/work/Caller/caller/case_187/short_false_negative_main.json'
    print_to_file(short_f_neg,  short_f_neg_file_name)
    print('Rare: ' + str(freq_count) + ' variants')
    print('In short version of false negative ' + str(len(short_f_neg))) + ' variants.'



def samples(count,  vars_file_name):
    vars = get_json_from_file(vars_file_name)
    #res = []
    num_list = range(len(vars))
    r_h = WichmannHill(3948)
    r_h.shuffle(num_list)
    new_vars = []
    for k in num_list[:30]:
        new_vars.append(vars[k])
    
#    if len(vars)<count:
#        print('In this list there are ' + len(vars) + ' variants only.')
#        res = vars
#    else:
#        inds = []
#        for i in range(count):
#            ind = int(random.random() * len(vars))
#            while ind in inds:
#                ind = int(random.random() * len(vars))
#            inds.append(ind)
#        for i in range(len(inds)):
#            res.append(vars[inds[i]])

    new_file_name = vars_file_name[:-5] + '_' + str(count) +'_samples.json'
    print_to_file(new_vars,  new_file_name)
    csv_file_name = new_file_name[:-5] + '.csv'
    json_to_csv(new_file_name,  csv_file_name)
    return new_vars


def atlas(radius,  base_variants,  vcf_file_name):
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    res = []
    centers = []
    for var in base_variants:
        res.append([])
        centers.append(None)
    for record in vcf_reader:
        for k in range(len(base_variants)):
            var = base_variants[k]
            if record.CHROM != var['CHROM']:
                continue
            if abs(record.POS - var['POS']) > radius:
                continue
            if record.POS == var['POS']:
                centers[k] = len(res[k])
                print('Chart ' + var['CHROM'] + ':' + str(var['POS']))
            rec = {}
            rec['CHROM'] = record.CHROM
            rec['POS'] = record.POS
            rec['remoteness_(base_pairs)'] = record.POS - var['POS']
            rec['QUAL'] = record.QUAL
            rec['start'] = record.start
            rec['end'] = record.end
            rec['REF'] = str(record.REF)
            rec['ALT'] = ''
            for alt in record.ALT:
                rec['ALT'] += str(alt)
            quality = ''
            if record.FILTER == []:
                quality += 'YES/'
            else:
                quality += 'NO/'
            if 'FS' in record.INFO:
                quality += str(record.INFO['FS']) + '/'
            else:
                quality += 'None/'
            if 'QD' in record.INFO:
                quality += str(record.INFO['QD']) + '/'
            else:
                quality += 'None/'
            GQ_str = '('
            for sample in record.samples:
                GQ_str += str(sample.data.GQ) + '|'
            quality += GQ_str[:-1] + ')'
            rec['PASS/FS/QD/GQ'] = quality
            rec['ExAC_AF'] = ''
            alls = range(len(record.INFO['AF']))
            minors = ''
            for all in alls:
                frequency = get_frequency(record.INFO['CSQ'], str(record.ALT[all]))
                rec['ExAC_AF'] += str(frequency) + '/'
                minor = False
                for sample in record.samples:
                    if is_unaffected(sample.sample):
                        if sample.data.AD[all+1] > 0:
                            minor = False
                            break
                    else:
                        if not (sample.gt_type > 0 or sample.gt_type == 0 and sample.data.AD[all+1] > 0):
                            minor = False
                            break
                        if sample.gt_type == 0 and sample.data.AD[all+1]>0:
                            minor = True
                if minor:
                    minors += '+|'
                else:
                    minors += '-|'
            rec['minor_read'] = minors[:-1]
            rec['ExAC_AF'] = rec['ExAC_AF'][:-1]
            for sample in record.samples:
                rec[sample.sample + '_zygosity'] = sample.gt_type
            res[k].append(rec)
    
    new_res = []
    for k in range(len(res)):
        new_res.append([])
        for m in range(len(res[k])):
            res[k][m]['remoteness_(variants)'] = m - centers[k]
            if abs(m-centers[k]) <= 100:
                new_res[k].append(res[k][m])
    
    for k in range(len(res)):
        atlas_file_name = '/home/andrey/work/Caller/caller/atlas/chart_'+ str(k) + '.json'
        print_to_file(new_res[k],  atlas_file_name)
        json_to_csv_columns(atlas_file_name)






if __name__=='__main__':
    drop_out()
    vars_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative_main.json'
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    print('Generate samples...')
    base_variants = samples(30, vars_file_name)
    print('Compose atlas...')
    atlas(10000,  base_variants,  vcf_file_name)
    print('Ok.')
