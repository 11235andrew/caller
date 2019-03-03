import vcf
import json

def find_false_negative(vcf_file_name,  f_negative_file_name):
    try:
        vcf_file = open(vcf_file_name, 'r')
    except IOError:
        print('File "' + vcf_file_name + '" not found.')
        return
    
    vcf_reader = vcf.Reader(vcf_file)
    calls = []
    count = 0
    for record in vcf_reader:
        calls.append(record.__dict__())
        count += 1
        if count > 10:
            break
    vcf_file.close()
    
    try:
        f_negative_file = open(f_negative_file_name, 'r')
    except IOError:
        print('File "' + f_negative_file_name + '" not found.')
        return
    f_negative_file.write(json.dumps(calls,  indent=4))
    f_negative_file.close()





if __name__ == '__main__':
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    f_negative_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative_bgm0187.vcf'
    find_false_negative(vcf_file_name,  f_negative_file_name)
    print('Ok.')
