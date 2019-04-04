from print_record import get_json_from_file
from print_record import print_to_file


def drop_out():
    f_neg_file_name = '/home/andrey/work/Caller/caller/case_187/false_negative.json'
    f_neg = get_json_from_file(f_neg_file_name)
    short_f_neg = []
    for rec in f_neg:
        if rec['gnomAD_AF'] > 0.01:
            continue
        short_f_neg.append(rec)
    short_f_neg_file_name = '/home/andrey/work/Caller/caller/case_187/short_false_negative.json'
    print_to_file(short_f_neg,  short_f_neg_file_name)
    print('In short version of false negative ' + str(len(short_f_neg))) + ' variants.'












if __name__=='__main__':
    drop_out()
    print('Ok.')
