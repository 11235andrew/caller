import vcf
import json
import gzip
from print_record import get_frequency
from print_record import print_to_file
from print_record import print_record
from print_record import print_short_record
from drop_out import is_unaffected
from print_record import open_file



def sample_to_dict(sample):
    res = {}
    res['sample'] = sample.sample
    res['AD'] = sample.data.AD
    res['GQ'] = sample.data.GQ
    res['GT'] = sample.data.GT
    return res

def get_work_version(record,  AFs):
    rec = {}
    rec['CHROM'] = record.CHROM
    rec['POS'] = record.POS
    rec['AF'] = str(record.INFO['AF'])
    rec['QUAL'] = record.QUAL
    rec['REF'] = record.REF
    rec['ALT'] = record.ALT
    if 'FS' in record.INFO:
        rec['FS'] = record.INFO['FS']
    else:
        rec['FS'] = None
    if 'QD' in record.INFO:
        rec['QD'] = record.INFO['QD']
    else:
        rec['QD'] = None
    rec['ExAC_AF'] = []
    if str(record.POS) in AFs:
        rec['gnomAD_AF'] = AFs[str(record.POS)]
    else:
        rec['gnomAD_AF'] = 'None'
    alls = range(len(record.INFO['AF']))
    for all in alls:
        frequency = get_frequency(record.INFO['CSQ'], str(record.ALT[all]))
        rec['ExAC_AF'].append(frequency)
    rec['samples'] = []
    for sample in record.samples:
        rec['samples'].append(sample_to_dict(sample))
    return rec

def get_readible_version(record,  AFs):
    rec = {}
    rec['CHROM'] = record.CHROM
    rec['POS'] = record.POS
    rec['AF'] = str(record.INFO['AF'])
    rec['QUAL'] = record.QUAL
    rec['REF'] = record.REF
    rec['ALT'] = record.ALT
    if 'FS' in record.INFO:
        rec['FS'] = record.INFO['FS']
    else:
        rec['FS'] = None
    if 'QD' in record.INFO:
        rec['QD'] = record.INFO['QD']
    else:
        rec['QD'] = None
    rec['ExAC_AF'] = []
    if str(record.POS) in AFs:
        rec['gnomAD_AF'] = AFs[str(record.POS)]
    else:
        rec['gnomAD_AF'] = 'None'
    alls = range(len(record.INFO['AF']))
    for all in alls:
        frequency = get_frequency(record.INFO['CSQ'], str(record.ALT[all]))
        rec['ExAC_AF'].append(frequency)
    rec['samples'] = []
    for sample in record.samples:
        rec['samples'].append(str(sample_to_dict(sample)))
    return rec

def get_gnomAD_frequency():
    print('Research frequencies in gnomAD...')
    gnomad = '/data/exp/trifon/vault/xl_BGM0187/fdata.json.gz'
    AFs = {}
    with gzip.open(gnomad, "rb") as inp:
        for line in inp:
            rec_data = json.loads(line)
            if rec_data['Chromosome'] in ['chrM',  'chrX',  'chrY']:
                continue
            AFs[str(rec_data['Start_Pos'])] = rec_data['gnomAD_AF']
#            if rec_data['Start_Pos'] == pos and rec_data['Chromosome'] == chrm:
#                inp.close()
#                return rec_data['gnomAD_AF']
    inp.close()
    return AFs
    
def intersection(list1,  list2):
    res = []
    for el in list1:
        if el in list2:
            res.append(el)
    return res

def is_f_neg(record):
    if record.CHROM == 'chrM' or record.CHROM == 'chrX' or record.CHROM == 'chrY':
        return
    alls = range(len(record.INFO['AF']))
#    if 'QD' not in record.INFO or record.INFO['QD'] < 4:
#        return
#    if 'FS' not in record.INFO or record.INFO['FS'] > 30:
#        return
    if record.FILTER != []:
        return
    for all in alls:
#        frequency = get_frequency(record.INFO['CSQ'], str(record.ALT[all]))
#        if frequency is None:
#            continue
        minor = False
        for sample in record.samples:
#            if sample.data.GQ is None or sample.data.GQ < 20:
#                minor = False
#                break
            if is_unaffected(sample.sample):
                if sample.gt_type != 0 or sample.data.AD[all+1] > 0:
                    minor = False
                    break
            else:
                if not (sample.gt_type > 0 or sample.gt_type == 0 and sample.data.AD[all+1] > 0):
                    minor = False
                    break
                if sample.gt_type == 0 and sample.data.AD[all+1]>0:
                    minor = True
        if minor:
            return all
    return
    
def all_records(vcf_file_name,  f_neg_file_name):
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    f_neg_main = []
    f_neg = []
    count = 0
    AFs = get_gnomAD_frequency()
    for record in vcf_reader:
        count += 1
        if count % 1000 == 0:
            print('Record #' + str(count))
        rec = get_work_version(record,  AFs)
        rec_read = get_readible_version(record,  AFs)
        is_neg = is_f_neg(record)
        if is_neg is not None:
            rec['ALT_index'] = is_neg
            rec_read['ALT_index'] = is_neg
            f_neg.append(rec_read)
            f_neg_main.append(rec)
    vcf_file.close()
    print(str(len(f_neg)) + ' false negative records were found.')
    print_to_file(f_neg,  f_neg_file_name)
    print_to_file(f_neg_main,  f_neg_file_name[:-5] + '_main.json')
    

def rude_classificator(vcf_file_name,  cand_file_name, f_pos_file_name,  f_neg_file_name):
    try:
        vcf_file = open(vcf_file_name, 'r')
    except IOError:
        print('File "' + vcf_file_name + '" not found.')
        return
    
    vcf_reader = vcf.Reader(vcf_file)
    pos_neg = []
    f_pos = []
    f_neg = []
    pos_neg_main = []
    f_pos_main = []
    f_neg_main = []
    count = 0
    max_qual = 0
    filters = []
    gq = 0
    qd = 0
    fs = 0
    freq = 0
    pass_f = 0
    in_gnom = 0
    AFs = get_gnomAD_frequency()
    print('In Gnom_AD ' + str(len(AFs)) + ' records.')
    for record in vcf_reader:
        count += 1
        if count % 1000 == 0:
            print('Record #' + str(count))
        if record.CHROM == 'chrM' or record.CHROM == 'chrX' or record.CHROM == 'chrY':
            continue
        
#        if str(record.POS) in AFs:
#            frequency = AFs[str(record.POS)]
#            if frequency is None or frequency>0.01:
#                continue
#        else:
#            continue
#        in_gnom += 1
        
        #flag = False
        alls = range(len(record.INFO['AF']))
#        min_fr = 3
#        min_all = None
#        for all in range(len(record.INFO['AF'])):
#            if record.INFO['AF'][all] < min_fr:
#                min_fr = record.INFO['AF'][all]
#                min_all = all
#        if min_all is not None:
#            alls.append(min_all)
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
        

        if 'QD' in record.INFO and record.INFO['QD'] > 4:
            qd += 1
        else:
            continue
        if 'FS' in record.INFO and record.INFO['FS'] < 30:
            fs += 1
        else:
            continue
        
        rec = get_work_version(record,  AFs)
        rec_read = get_readible_version(record,  AFs)
        for all in alls:
            rec['ALT_index'] = all
            rec_read['ALT_index'] = all
            frequency = get_frequency(record.INFO['CSQ'], str(record.ALT[all]))
            if frequency is not None and frequency > 0.01:
                continue
            owns = []
            pos_flag = False
            neg_flag = True
            neg_flag2 = False
            GQ_flag = True
            for sample in record.samples:
                au = sample.sample[-2]
                GT = sample.data.GT
                AD = sample.data.AD
                GQ = sample.data.GQ
                if GQ is None or GQ <= 20:
                    GQ_flag = False
                if au == 'u':
                    if GT == '0/0':
                        owns.append(sample.sample)
                    if AD[all+1] != 0:
                        pos_flag = True
                else:
                    if GT in ['0/1',  '0/1', '1/1']:
                        owns.append(sample.sample)
                        neg_flag2 = True
                    if AD[all+1] == 0:
                        neg_flag = False
            if not GQ_flag:
                continue
            if len(owns) ==  len(record.samples):
                pos_neg.append(rec_read)
                pos_neg_main.append(rec)
                if pos_flag:
                    f_pos.append(rec_read)
                    f_pos_main.append(rec)
                break
            elif neg_flag2 and neg_flag:
                f_neg.append(rec_read)
                f_neg_main.append(rec)
        

    print(str(len(f_neg)) + ' false negative records were found.')
    print(str(len(f_pos)) + ' false positive records were found.')
    #print(str(len(intersection(f_neg,  f_pos))))
    print('Need variants in gnomAD: ' + str(in_gnom))
    print('Maximal Quality is ' + str(max_qual))
    print('Allele freqency < 0.01: ' + str(freq) + ' variants.')
    print('Quality_PASS: ' + str(pass_f) + ' variants.')
    print('GQ_MEAN > 20: ' + str(gq) + ' variants.')
    print('QD > 4: ' + str(qd) + ' variants.')
    print('FS < 30: ' + str(fs) + ' variants.')
    print(str(len(pos_neg)) + ' variants were found.')
    
    vcf_file.close()
    
    
    print_to_file(pos_neg,  cand_file_name)
    print_to_file(f_pos,  f_pos_file_name)
    print_to_file(f_neg,  f_neg_file_name)
    print_to_file(pos_neg_main,  cand_file_name[:-5] + '_main.json')
    print_to_file(f_pos_main,  f_pos_file_name[:-5] + '_main.json')
    print_to_file(f_neg_main,  f_neg_file_name[:-5] + '_main.json')
    print_to_file(filters,  'case_187/filters.json')
    



if __name__ == '__main__':
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    f_negative_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative.json'
    f_positive_file_name = '/home/andrey/work/Caller/caller/case_187/false_positive.json'
    candidats_file_name = '/home/andrey/work/Caller/caller/case_187/candidats.json'
    #rude_classificator(vcf_file_name,  candidats_file_name,  f_positive_file_name,  f_negative_file_name)
    all_records(vcf_file_name,  f_negative_file_name)
    print('Ok.')
