import random
from print_record import get_json_from_file
from print_record import print_to_file
from print_record import json_to_csv


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
    res = []
    if len(vars)<count:
        print('In this list there are ' + len(vars) + ' variants only.')
        res = vars
    else:
        inds = []
        for i in range(count):
            ind = int(random.random() * count)
            while ind not in inds:
                ind = int(random.random() * count)
        for i in range(len(inds)):
            res.append(vars[inds[i]])
    new_file_name = vars_file_name[:-5] + '_' + str(count) +'_samples.json'
    print_to_file(res,  new_file_name)
    csv_file_name = new_file_name[:-5] + '.csv'
    json_to_csv(new_file_name,  csv_file_name)







if __name__=='__main__':
    drop_out()
    vars_file_name = '/home/andrey/work/Caller/caller/case_187/candidats_main.json'
    samples(30, vars_file_name)
    print('Ok.')
