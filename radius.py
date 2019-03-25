import vcf
from print_record import open_file
from print_record import print_to_file
from print_record import get_frequency
from print_record import get_json_from_file


def get_boundary(base_variants,  radius, measure, vcf_file_name):
    print('Find boundaries...')
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    all_vars = {}
    count = 0
    for record in vcf_reader:
        if record.CHROM not in all_vars:
            all_vars[record.CHROM] = []
        all_vars[record.CHROM].append(record.POS)
        count += 1
        if count % 1000 == 0:
            print('Record #' + str(count))
    
    res = []
    for variant in base_variants:
        if variant['POS'] in all_vars[variant['CHROM']]:
            rec = {}
            rec['CHROM'] = variant['CHROM']
            rec['POS'] = variant['POS']
            if measure == 'variants':
                ind = all_vars[variant['CHROM']].index(variant['POS'])
                if ind - radius < 0:
                    rec['left'] = all_vars[variant['CHROM']][0]
                else:
                    rec['left'] = all_vars[variant['CHROM']][ind-radius]
                if ind + radius >= len(all_vars[variant['CHROM']]):
                    rec['right'] = all_vars[variant['CHROM']][-1]
                else:
                    rec['right'] = all_vars[variant['CHROM']][ind+radius]
            else:
                for pos in all_vars[variant['CHROM']]:
                    if pos >= variant['CHROM'] - radius and 'left' not in rec:
                        rec['left'] = pos
                    if pos > variant['pos'] + radius and 'right' not in rec:
                        ind = all_vars.index(pos)
                        rec['right'] = all_vars[variant['CHROM']][ind - 1]
                if 'right' not in rec:
                    rec['right'] = all_vars[variant['CHROM']][-1]
            res.append(rec)
        else:
            print('Variant ' + variant['CHROM'] + ':' + str(variant['POS']) + ' not found in the vcf-file')
            continue
    vcf_file.close()
    return res

def neighbourhood(base_variants,  affected,  unaffected,  radius,
                                            measure,  frequency,  vcf_file_name):
    boundaries = get_boundary(base_variants,  radius,  measure,  vcf_file_name)
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    count = 0
    for bound in boundaries:
        bound['aproved_at_the_right'] = 0
        bound['aproved_at_the_left'] = 0
        bound['not_aproved_at_the_right'] = 0
        bound['not_aproved_at_the_left'] = 0
        bound['the_most_distant_at_the_right'] = bound['right'] - bound['POS']
        bound['the_most_distant_at_the_left'] = bound['POS'] - bound['left']
        bound['the_main_condition_at_the_left'] = 0
        bound['the_main_condition_at_the_right'] = 0
        bound['not_main_condition_right'] = None
        bound['not_main_condition_left'] = None
    for record in vcf_reader:
        count += 1
        if count % 1000 == 0:
            print('Record #' + str(count))
        exac = get_frequency(record.INFO['CSQ'], str(record.ALT[all]))
        aproved =  exac is not None and exac > frequency and record.QUAL == []
        main = True
        for sample in record.samples:
            if sample.sample in affected:
                if sample.GT not in ['1/0',  '0/1',  '1/1']:
                    main = False
                    break
            if sample.sample in unaffected:
                if sample.GT in ['1/0',  '0/1',  '1/1']:
                    main = False
                    break

        for bound in boundaries:
            if bound['CHROM'] != record.CHROM or bound['left'] > record.POS or bound['right'] < record.POS:
                continue
            if record.POS < bound['POS']:
                if aproved:
                    bound['aproved_at_the_left'] += 1
                else:
                    bound['not_aproved_at_the_left'] += 1
                if main:
                    bound['the_main_condition_at_the_left'] += 1
                else:
                    bound['not_main_condition_left'] = record.POS
            if record.POS > bound['POS']:
                if aproved:
                    bound['aproved_at_the_right'] += 1
                else:
                    bound['not_aproved_at_the_right'] += 1
                if main:
                    bound['the_main_condition_at_the_right'] += 1
                else:
                    if bound['not_main_condition_right'] is None:
                        bound['not_main_condition_right'] = record.POS
    vcf_file.close()
    
    file_name = '/home/andrey/work/Caller/caller/case_187/false_negative_neighbourhoods.json'
    print_to_file(file_name,  boundaries)




if __name__ == '__main__':
    f_neg_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative.json'
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    affected = ['bgm0187a1',  'bgm0187a2',  'bgm0187a3',  'bgm0187a4',  'bgm0187a5']
    unaffected = ['bgm0187u1']
    radius = 100
    frequency = 0.1
    f_neg = get_json_from_file(f_neg_file_name)
    neighbourhood(f_neg,  affected,  unaffected, 
                                radius,  'variants',  frequency,  vcf_file_name)
    print('Ok.')
