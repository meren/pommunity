# -*- coding: utf-8 -*-

import os
import sys
import random
import ConfigParser

# everything in this code is based QC processes implemented in MBL. Once a 454
# run is done, quality control, chimera-checking and trimming steps take place
# to make sure that sequences with low quality is removed from the library. in
# normal conditions, expected error rate of 454 machine is 1%. however, once
# the pipeline in MBL done with QC/QA, expected error rate drops down to 0.25%.
# which is equivalent to the probability of 0.0025 to have an erronous base call:
expected_error_rate = 0.0025

# these substutition probabilities are coming from Sue Huse et al.'s
# 454 quality studies. they simply show what is the probability of a given base
# to be substituted with one of the other three possible bases, *if* there is
# a substutition error at the given location:
substutition_probabilities = {'T': 'C' * 23 + 'A' * 6  + 'G' * 6,
                              'C': 'A' * 7  + 'T' * 11 + 'G' * 1,
                              'G': 'A' * 7  + 'T' * 3  + 'C' * 4,
                              'A': 'T' * 7  + 'C' * 3  + 'G' * 22}

# following is important. even though expected_error_rate gives us the expected
# error rate, it doesn't reflect the fact that some bases more probable
# to have a substitution error than others. for instance, among all errors, 32% of errors 
# observed in A and %35 in T, according to the substutition_probabilities Huse computed. therefore
# expected_error_rate should be normalized in such a way that resulting dictionary
# should converge to expected_error_rate in the long run, but include base biases. following
# formula normalizes those results in such a way.
nucleotide_normalized_error_probabilities = {}
for n in ['A', 'T', 'C', 'G']:
    nucleotide_normalized_error_probabilities[n] = int(round(1 / ((len(substutition_probabilities[n]) / 100.0) * 4 * expected_error_rate)))

def pp(n):
    """Pretty print function for very big numbers.."""
    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()

    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

def has_error(location, base):
    if location < 400:
        return random.randint(1, 400) == 1
    if location > 400:
        return random.randint(1, 200) == 1

def update_nucleotide(nucleotide):
    # here it is an error that is not associated with homopolymer
    # regions. it may be a substitution, or indel. FIXME: these
    # probabilities must be updated.
    #
    # assumptions:
    #
    # %5  percent of errors are insertions.
    # %5  percent of errors are deletions.
    # %90 percent of errors are substitutions.

    r = random.randint(1, 100)
    if r <= 5:
        # this is a deletion
        return ''
    if r <= 10:
        # this is an insertion, we don't know what to insert,
        # ergo we will likelihood of the insertion of any base
        # is equal
        return nucleotide + ['A', 'T', 'C', 'G'][random.randint(0, 3)]
    else:
        # here it is a substitution.. probability of change is different
        # for every base. we will use probabilities provided by Susan Huse
        # coded in substutition_probabilities dict.
        return random.sample(substutition_probabilities[nucleotide], 1)[0]


class Member:
    def __init__(self):
        self.id = None
        self.sequence = None
        self.ratio = 0
        self.abundance = 0

    def __str__(self):
        return self.id

class CommunityConfiguration:
    def __init__(self):
        self.output_file = None
        self.platform = None
        self.total_sequences = 0

        self.members = []

    def compute_member_distribution(self):
        ratios = sum([member.ratio for member in self.members])
        self.member_distribution = []
        for member in self.members:
            abundance = int(round(((member.ratio + 0.0) / ratios) * self.total_sequences))
            member.abundance = abundance
            self.member_distribution += [member] * abundance
        
        random.shuffle(self.member_distribution)

def main(config):
    sys.stderr.write('\nCommunity membership distribution is being computed...')
    config.compute_member_distribution()
    sys.stderr.write('\rCommunity membership distribution for "%s" (a total of %s sequences):\n\n' % 
                                                       (config.community, pp(config.total_sequences)))
    sys.stderr.flush()

    sys.stderr.write('\t\t%20s %20s %20s\n' % ('member id', 'ratio', 'abundance'))
    sys.stderr.write('\t\t%s\n' % ('-' * 80))
    for member in config.members:
        sys.stderr.write('\t\t%20s %18sX %20s\n' % (member.id, member.ratio, pp(member.abundance)))
    sys.stderr.write('\t\t%s\n\n' % ('-' * 80))

   
    output = open(config.output_file, 'w')
    sys.stderr.write('Introducing random errors... ')
    for i in range(0, len(config.member_distribution)):
        if i % 1000:
            sys.stderr.write('\rIntroducing random erros -- %.2d%%' % (i * 100 / config.total_sequences))
            sys.stderr.flush()

        member = config.member_distribution[i]
        
        sequence_with_error = ''
        number_of_errors = 0
        for j in range(0, len(member.sequence)):
            nucleotide = member.sequence[j]
            # testing every nucleotide in the read
            if has_error(j, nucleotide):
                number_of_errors += 1
                # if we are here, it means we hit the prob of 0.0025.
                # it is time to find out what type error we have,
                # with what kind of outcome.
                nucleotide = update_nucleotide(nucleotide)
            
            sequence_with_error += nucleotide

        output.write('>%s | type: %s | errors: %s\n' % 
                              ('_'.join([config.community, i.__str__()]),
                               member.id,
                               number_of_errors))
        output.write(sequence_with_error + '\n')
               
    output.close()
    sys.stderr.write('\n')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Microbial Pseudo Community Sample Generation')
    parser.add_argument('-c', '--configuration', required=True, metavar = 'CONFIG_FILE', 
                                                 help = 'Configuration parameters of the community')

    args = parser.parse_args()
    user_config = ConfigParser.ConfigParser()
    user_config.read(args.configuration)

    config = CommunityConfiguration()

    config.output_file = user_config.get('general', 'output_file')
    config.platform = user_config.get('general', 'platform')
    config.community = user_config.get('general', 'community')
    config.total_sequences = int(user_config.get('general', 'total_sequences'))

    for section in user_config.sections():
        if section.startswith('member'):
            member = Member()
            member.id = '-'.join(section.replace('member', '').strip().split())
            member.sequence = user_config.get(section, 'sequence')
            member.ratio = int(user_config.get(section, 'ratio'))
            
            config.members.append(member)

    sys.exit(main(config))
