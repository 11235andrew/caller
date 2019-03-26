import sys
import vcf
import json
import gzip


def get_chrms(vcf_file_name):
    vcf_file = open_file(vcf_file_name,  'r')
    vcf_reader = vcf.Reader(vcf_file)
    chrms = []
    for record in vcf_reader:
        if record.CHROM not in chrms:
            chrms.append(record.CHROM)
    chrms_file_name = '/home/andrey/work/Caller/caller/case_187/chrms.json'
    print_to_file(chrms,  chrms_file_name)

def get_frequency(csq, allele):
    infos_file_name = '/home/andrey/work/Caller/caller/case_187/infos.json'
    format = get_json_from_file(infos_file_name)
    format = format['Format']
    
    n = None
    freq = None
    for n in range(len(format)):
        if format[n] == 'ExAC_AF':
            break
    if n is None or csq == []:
        return
    for csq_el in csq:
        data = csq_el.split('|')
        if data[0] != allele:
            continue
        if data[n] == '':
            return
        freq = float(data[n])
        break
    return freq

def infos(vcf_reader):
    infos = vcf_reader.infos
    descr = infos['CSQ'][3]
    place = descr.find('Format')
    format = descr[place + 8:]
    csq = {}
    csq['Description'] = descr[:place]
    csq['Format'] = format.split('|')
    infos_file_name = '/home/andrey/work/Caller/caller/case_187/infos.json'
    print_to_file(csq,  infos_file_name)
    return format.split('|')

def print_gnomAD(chrm,  pos):
    gnomad = '/data/exp/trifon/vault/xl_BGM0187/fdata.json.gz'
    with gzip.open(gnomad, "rb") as inp:
        for line in inp:
            rec_data = json.loads(line)
            if str(rec_data['Start_Pos']) == pos and rec_data['Chromosome'] == chrm:
                return rec_data
    inp.close()
    return

def print_short_record(record, freq):
    res = {}
    res['CHROM'] = record.CHROM
    res['POS'] = record.POS
    res['AF'] = record.INFO['AF']
    res['gnomAD_AF'] = freq['gnomAD_AF']
    
   
    freqs = []
    for all in range(len(record.INFO['AF'])):
        frequency = get_frequency(record.INFO['CSQ'], format,  str(record.ALT[all]))
        freqs.append(frequency)
    res['ExAC_AF'] = freqs
    
    res['FS'] = record.INFO['FS']
    res['QD'] = record.INFO['QD']
    samples = []
    for sample in record.samples:
        samples.append(str(sample))
    res['samples'] = samples
    return res

def print_record(record, freq, format):
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
    rec_dict['gnomAD_AF'] = freq['gnomAD_AF']
    csq = []
    for data in rec_dict['INFO']['CSQ']:
        vert = data.split('|')
        csq_dict = {}
        if len(vert) != len(format):
            print('Defferent lengths: ' + len(format) + ' (format) and ' + len(vert) + ' (csq).')
        else:
            for n in range(len(vert)):
                csq_dict[format[n]] = vert[n]
            csq.append(csq_dict)
    rec_dict['INFO']['CSQ'] = csq
    return rec_dict


def open_file(file_name,  mode):
    try:
        file = open(file_name,  mode)
    except IOError:
        print('File "' + file_name + '" not found.')
        sys.exit()
    return file

def get_json_from_file(file_name):
    file = open_file(file_name,  'r')
    data = json.loads(file.read())
    file.close()
    return data

def print_to_file(data,  file_name):
    try:
        data_file = open(file_name, 'w')
    except IOError:
        print('File "' + file_name + '" not found.')
        return
    data_file.write(json.dumps(data,  indent=4))
    data_file.close()

def json_to_csv(json_file_name,  csv_file_name):
    data = get_json_from_file(json_file_name)
    if data == []:
        return
    keys = data[0].keys()
    csv = ''
    for el in data:
        line = ''
        for key in keys:
            line +=str(el[key]) + '\t'
        csv += line + '\n'
    csv_file = open_file(csv_file_name,  'w')
    csv_file.write(csv)
    csv_file.close()


if __name__ == '__main__':
    if len(sys.argv) not in [3,  4]:
        print('Use following format: python print_record.py CHROM POS')
        sys.exit()
    chm = sys.argv[1]
    pos = sys.argv[2]
    short = False
    if len(sys.argv) == 4:
        if sys.argv[3] != '-s':
            print('Unkown key ' + sys.argv[3])
            sys.exit()
        short = True
    
    vcf_file_name = '/data/bgm/cases/bgm0187/bgm0187_wes_run2_xbrowse.vep.vcf'
    vcf_file = open_file(vcf_file_name, 'r')
    
    vcf_reader = vcf.Reader(vcf_file)
    format = infos(vcf_reader)
    for record in vcf_reader:
        if record.CHROM != chm or str(record.POS) != pos:
            continue
        
        freq = print_gnomAD(chm,  pos)
        if freq is None:
            print('Record in GnomAD for this variant is not found.')
        else:
            gnomAD_file_name = '/home/andrey/work/Caller/caller/case_187/gnomAD.json'
            print_to_file(freq,  gnomAD_file_name)
        
        if short:
            res = print_short_record(record,  freq)
        else:
            res = print_record(record,  freq,  format)
        record_file_name = '/home/andrey/work/Caller/caller/case_187/record.json'
        record_file = open(record_file_name,  'w')
        print_to_file(res,  record_file_name)
        vcf_file.close()
        
        print('Ok.')
        sys.exit()
    vcf_file.close()
    print('record CHROM=' + chm + ', POS=' + pos + ' not found.')
