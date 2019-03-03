import os
import json
from case_utils import parse_fam_file
from sortedcontainers import SortedDict


def get_people_data(dir):
    relatives = []
    list = os.listdir(dir)
    for folder in list:
        #print('Dir: ' + folder)
        files = os.listdir(dir + '/' + folder)
        for file in files:
            #print('File: ' + file)
            if file[-4:] == '.fam':
                sample = parse_fam_file(dir + '/' + folder + '/' + file)
                relatives.append(sample)
    
    json_file_name = 'relatives/samples.json'
    json_file = open(json_file_name,  'w')
    json_file.write(json.dumps(relatives,  indent=4))
    json_file.close()
    return relatives


def find_relatives(dir):
    families = {}
    samples = get_people_data(dir)
    for sample in samples:
        for key in sample:
            fam_id = sample[key]['family']
            if fam_id not in families:
                families[fam_id] = {}
                families[fam_id]['people'] = {}
                families[fam_id]['count_of_affected'] = 0
            families[fam_id]['people'][key] = sample[key]
            if sample[key]['affected']:
                families[fam_id]['count_of_affected'] += 1
    
    for fam_id in families:
        families[fam_id] = dict(SortedDict(families[fam_id]))
    
    families_file_name = 'relatives/families.json'
    families_file = open(families_file_name,  'w')
    families_file.write(json.dumps(families,  indent=4))
    families_file.close()
    return families

def get_all_ancestors(name,  family, ancestors,  compare,  line):
    if name == '0':
        return False
    if len(ancestors) > len(family):
        return False
    if name not in family:
        return False
    if name in compare:
        return True
    if family[name]['mother'] != '0':
        ancestors.append(family[name]['mother'])
        new_line = line[:]
        new_line.append(family[name]['mother'])
        res = get_all_ancestors(family[name]['mother'],  family,  ancestors,  compare,  new_line)
        if res:
            line.append(family[name]['mother'])
            return res
        
    if family[name]['father'] != '0':
        ancestors.append(family[name]['father'])
        new_line = line[:]
        new_line.append(family[name]['father'])
        res = get_all_ancestors(family[name]['father'],  family,  ancestors,  compare,  new_line)
        if res:
            line.append(family[name]['father'])
            return res
    return False

def find_ancestors(families):
    ancestors = {}
    for fam_id in families:
        ancestors[fam_id] = {}
        for key in families[fam_id]['people']:
            anc = []
            get_all_ancestors(key,  families[fam_id]['people'],  anc, [], [])
            ancestors[fam_id][key] = anc
    
    anc_file_name = 'relatives/ancestors.json'
    anc_file = open(anc_file_name,  'w')
    anc_file.write(json.dumps(ancestors,  indent=4))
    anc_file.close()
    return ancestors

def intersection(list1,  list2):
    res = []
    for el in list1:
        if el in list2:
            res.append(el)
    return res

def common_ancestors(families,  ancestors):
    couples = []
    different = {}
    for fam_id in ancestors:
        people = families[fam_id]['people']
        for key in ancestors[fam_id]:
            print('Name: ' + key)
            if people[key]['affected']:
                for name in people:
                    if people[name]['affected'] and name != key:
                        inter = intersection(ancestors[fam_id][key],  ancestors[fam_id][name])
                        if inter != [] or name in ancestors[fam_id][key] or key in ancestors[fam_id][name]:
                            if key in different and name not in different[key]:
                                couple = {}
                                couple['family'] = fam_id
                                couple['first'] = key
                                couple['second'] = name
                                couple['common_ancestors'] = inter
                                couples.append(couple)
                    if name not in different:
                        different[name] = []
                    different[name].append(key)
    print('Couples: ' + str(len(couples)))
    couples_file_name = 'relatives/couples.json'
    couples_file = open(couples_file_name,  'w')
    couples_file.write(json.dumps(couples,  indent=4))
    couples_file.close()
    return couples

def draw_common_ancestors(couples,  families,  common_file_name):
    common_file = open(common_file_name,  'w')
    
    for couple in couples:
#        title = '\n'
#        for key in couple:
#            title = key + ': ' + couple[key] + '\t'
        title = '\n' + 'Family: ' + couple['family'] + ', affected: ' + couple['first'] + ', ' + couple['second']
        common_file.write(title)
        lines = draw_parents(couple['first'],  families[couple['family']]['people'],  [])
        common_file.write('\n')
        for line in lines:
            common_file.write('\n' + line)
        common_file.write('\n')
        lines = draw_parents(couple['second'],  families[couple['family']]['people'],  [])
        for line in lines:
            common_file.write('\n' + line)
        end = '\n'
        for n in range(50):
            end += '-'
        common_file.write(end)
    common_file.close()

def distant_relatives(couples,  families):
    distant = []
    for couple in couples:
        people = families[couple['family']]['people']
        if people[couple['first']]['mother'] == people[couple['second']]['mother']:
            continue
        if people[couple['first']]['father'] == people[couple['second']]['father']:
            continue
        if people[couple['first']]['mother'] == couple['second']:
            continue
        if people[couple['first']]['father'] == couple['second']:
            continue
        if people[couple['second']]['mother'] == couple['first']:
            continue
        if people[couple['second']]['father'] == couple['first']:
            continue
        distant.append(couple)
    print(str(len(distant)) + ' distant relatives were found.')
    return distant

#def find_common_ancestor(families):
#    couples = {}
#    for fam_id in families:
#        if families[fam_id]['count_of_affected'] < 2:
#            continue
#        couples[fam_id] = {}
#        people = families[fam_id]['people']
#        for key in people:
#            if people[key]['affected']:
#                ancestors = []
#                get_all_ancestors(key,  people,  ancestors,  [],  [])
#                print('Family: ' + fam_id + ', name: ' + key + ', ancestors: ' + str(ancestors))
#                for man in people:
#                    if man != key and people[man]['affected']:
#                        if man in couples[fam_id] and key in couples[fam_id][man]:
#                            continue
#                        line = []
#                        if get_all_ancestors(man,  people,  [],  ancestors,  line):
#                            other_line = []
#                            get_all_ancestors(key,  people,  [],  line,  other_line)
#                            couple = {}
#                            couple[key] = other_line
#                            couple[man] = line
#                            if key not in couples[fam_id]:
#                                couples[fam_id][key] = {}
#                            couples[fam_id][key][man] = couple
#    
#    couples_file_name = 'relatives/couples.json'
#    couples_file = open(couples_file_name,  'w')
#    couples_file.write(json.dumps(couples,  indent=4))
#    couples_file.close()


def get_count_of_ancestors(name,  family):
    if name == '0' or name not in family:
        return 0
    else:
        count = 1
        count += get_count_of_ancestors(family[name]['father'],  family)
        count += get_count_of_ancestors(family[name]['mother'],  family)
        return count

def draw_parents(name,  family,  lines):
    if name == '0' or name not in family:
        return lines
    father = family[name]['father']
    if father == '0' or father not in family:
        f_level = 0
    else:
        f_level = get_count_of_ancestors(family[father]['father'],  family)
    f_count = get_count_of_ancestors(father,  family)
    mother = family[name]['mother']
    if mother == '0' or mother not in family:
        m_level = f_count + 1
    else:
        m_level = f_count + 1 + get_count_of_ancestors(family[mother]['father'],  family)
    m_count = get_count_of_ancestors(mother,  family)
    #print('f_level: ' + str(f_level) + ', m_level: ' + str(m_level))
    
    if lines == []:
        lines = ['' for n in range(f_count + m_count + 1)]
#    print('Name: ' + name + ', f_count: ' + str(f_count) + ', ,_count: ' + str(m_count) +', lines: ')
#    for line in lines:
#        print('"' + line + '"')
#    print('Len: ' + str(len(lines)))
    if family[name]['affected']:
        aff = 'a'
    else:
        aff = 'u'
    key_line = name + '(' + aff + ')'
    indent = ''
    for k in range(len(key_line)):
        indent += ' '
    if f_count > 0:
        for n in range(f_level):
            lines[n] += indent + '     '
        lines[f_level] += indent + '  |--'
        for n in range(f_level + 1,  f_count):
            lines[n] += indent + '  |  '
    lines[f_count] += key_line + '--|  '
    if m_count > 0:
        for n in range(f_count+1,  m_level):
            lines[n] += indent + '  |  '
        lines[m_level] += indent + '  |--'
        for n in range(m_level+1,  f_count + m_count+1):
            lines[n] += indent + '     '
    
#    print('After: ')
#    for line in lines:
#        print('"' + line + '"')
    
    new_lines = draw_parents(father,  family,  lines[:f_count])
    new_lines.append(lines[f_count])
    new_lines.extend(draw_parents(mother,  family,  lines[f_count +1:]))
    return new_lines

def draw_all_parents(families):
    trees_file_name = 'relatives/trees.txt'
    trees_file = open(trees_file_name,  'w')
    for fam_id in families:
        for key in families[fam_id]['people']:
            trees_file.write('\nFamily: ' + fam_id + ', name: ' + key)
            lines = draw_parents(key,  families[fam_id]['people'],  [])
            for line in lines:
                trees_file.write('\n' + line)
            trees_file.write('\n')
    trees_file.close()

if __name__ == '__main__':
    dir = 'bgm'
    families = find_relatives(dir)
    #find_common_ancestor(families)
    ancestors = find_ancestors(families)
    couples = common_ancestors(families,  ancestors)
    common_file_name = 'relatives/common_ancestors.txt'
    draw_common_ancestors(couples, families, common_file_name)
    distant = distant_relatives(couples,  families)
    distant_file_name = 'relatives/distant_relatives.txt'
    draw_common_ancestors(distant, families, distant_file_name)
    lines = draw_parents('bgm0187a1',  families['1']['people'],  [])
    for line in lines:
        print(line)
    draw_all_parents(families)
    print('Ok.')
