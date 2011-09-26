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


# the probabilty of having an insertion or deletion in a homopolymer region
# is correlated with the length of the sequences. following parameters
# of the empirical distributions in homopolymer errors are derived from
# a the paper titled "Characteristics of 454 pyrosequencing data —enabling
# realistic simulation with flowsim" by Suzanne Balzer et al. 
# 
# Paper can be obtained via this address:
#
# http://bioinformatics.oxfordjournals.org/content/26/18/i420.full
#
# Mean and standard deviation of normal distribution around homopolymer
# lengths of 6, 7, 8 and 9 are derived from the linear regression of
# this rule: {mean: n, std: 0.03494 + n * 0.06856} where n denotes the
# length of the homopolymer. 
homopolymer_error_probabilities = {0: {'mu': 0.1230, 'sigma': 0.0737},
                                   1: {'mu': 1.0193, 'sigma': 0.1230},
                                   2: {'mu': 2.0006, 'sigma': 0.1585},
                                   3: {'mu': 2.9934, 'sigma': 0.2188},
                                   4: {'mu': 3.9962, 'sigma': 0.3168},
                                   5: {'mu': 4.9550, 'sigma': 0.3863},
                                   6: {'mu': 6.0000, 'sigma': 0.4462},
                                   7: {'mu': 7.0000, 'sigma': 0.5148},
                                   8: {'mu': 8.0000, 'sigma': 0.5834},
                                   9: {'mu': 9.0000, 'sigma': 0.6519}}


# necessity for normalization: even though expected_error_rate gives us the expected
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

def has_in_del_sub_error(location, base):
    # this function basically doubles the base error probability
    # after position 400. Even though it makes sense when the resulting
    # quality scores from the 454 machine is analyzed, here I shall put
    # a FIXME, so at some point someone would take a better look.

    if location < 400:
        return random.randint(1, nucleotide_normalized_error_probabilities[base]) == 1
    if location > 400:
        return random.randint(1, nucleotide_normalized_error_probabilities[base] / 2) == 1

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
        # ergo we will assume the probability of any base to
        # be inserted here is equal
        return nucleotide + ['A', 'T', 'C', 'G'][random.randint(0, 3)]
    else:
        # here it is a substitution.. probability of change is different
        # for every base. we will use probabilities provided by Susan Huse
        # coded in substutition_probabilities dict.
        return random.sample(substutition_probabilities[nucleotide], 1)[0]

def update_homopolymer_region(homopolymer):
    # if, for some reason, 'homopolymer' is actually not a homopolymer,
    # return it ('not my problem' mode is [ON]).
    if not len(set(homopolymer)) == 1:
        return homopolymer

    # FIXME: I have to come up with an explanation for this one here:
    if random.randint(1, 5) != 1:
        return homopolymer

    e = homopolymer_error_probabilities[len(homopolymer)]

    return homopolymer[0] * int(round(random.gauss(e['mu'], e['sigma'])))

class Member:
    def __init__(self, sequence):
        self.sequence = sequence
        
        self.hp_regions = []
        self.hp_free_base_locs = []
        self.get_hp_free_base_locs()
      
        self.id = None
        self.ratio = 0
        self.abundance = 0

        self.error_distribution = {}

    def get_hp_free_base_locs(self):
        sequence = self.sequence + ' '
        regionIsHP = lambda x: len(set(x)) == 1
        baseInHP   = lambda l: filter(lambda tpl: l >= tpl[0] and l < tpl[1], self.hp_regions)

        index = 0
        while 1:
            if index == len(sequence):
                break

            j = 1
            while regionIsHP(sequence[index:index + j]):
                j += 1
                if j + index >= len(sequence):
                    break
            
            if j > 3:
                self.hp_regions.append((index, index + j - 1))
            index += j - 1

        for i in range(0, len(sequence.strip())):
            if not baseInHP(i):
                self.hp_free_base_locs.append(i)
 
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
            sys.stderr.write('\rIntroducing random erros -- %.2d%%' % (int(round(i * 100.0 / config.total_sequences))))
            sys.stderr.flush()

        member = config.member_distribution[i]
        
        sequence_with_errors = ''
        number_of_errors = 0
        number_of_hp_errors = 0
        
        j = 0
        while True:
            if j == len(member.sequence):
                break
            
            if j in member.hp_free_base_locs:
                nucleotide = member.sequence[j]
                # testing every nucleotide that is not in a homopolymer region in the
                # template sequence for insertions, deletions and substitutions:
                if has_in_del_sub_error(j, nucleotide):
                    # if we are here, it means we hit the prob of 0.0025.
                    # it is time to find out what type error we have,
                    # with what kind of outcome.
                    nucleotide = update_nucleotide(nucleotide)
                    
                    number_of_errors += 1
                
                sequence_with_errors += nucleotide
                
                j += 1
                     
            else:
                # potential FIXME: if we are in this else clause, we are about to process a homopolymer region.
                # at this point we ignore the fact that there may be substitution errors in homopolymer
                # regions, but I can't see why it would be the case. current assumption makes it easier to
                # implement this section due to the fact that we can rely on the template to find out where
                # homopolymer regions are (since they are not being changed by regular sub/in/del errors).
                # but this is something to think about. this approach might have to change.
                hp_tpl = [region for region in member.hp_regions if region[0] <= j and region[1] > j][0]
                
                original_hp_region = member.sequence[hp_tpl[0]:hp_tpl[1]]
                updated_hp_region = update_homopolymer_region(original_hp_region)
                
                sequence_with_errors += updated_hp_region
                
                if updated_hp_region != original_hp_region:
                    number_of_errors += 1
                    number_of_hp_errors += 1

                j += len(original_hp_region)
        

        member.error_distribution[number_of_errors] = member.error_distribution[number_of_errors] + 1 \
                                                         if member.error_distribution.has_key(number_of_errors) else 1
        
        output.write('>%s | type: %s | all_errors: %d | hp_errors: %d\n' % 
                              ('_'.join([config.community, i.__str__()]),
                               member.id,
                               number_of_errors,
                               number_of_hp_errors))
        
        output.write(sequence_with_errors + '\n')
    
    output.close()
    
    sys.stderr.write('\n\nBasic stats:\n')
    for member in config.members:
        sys.stderr.write('\n\tMember "%s"\n\t%s\n\n' % (member.id, '-'*50))
        sys.stderr.write('\t\t:: Abundance         : %s\n' % (pp(member.abundance)))
        total_bases = member.abundance * len(member.sequence)
        sys.stderr.write('\t\t:: Total bases       : %s\n' % (pp(total_bases)))
        total_errors = sum([x * member.error_distribution[x] for x in member.error_distribution])
        sys.stderr.write('\t\t:: Total errors      : %s\n' % (pp(total_errors)))
        error_base_ratio = total_errors / float(total_bases)
        sys.stderr.write('\t\t:: Error-base ratio  : 1 in %d\n' % (1 / error_base_ratio))
        sys.stderr.write('\t\t:: Error distribution:\n')
        for key in sorted(member.error_distribution.keys()):
            sys.stderr.write('\t\t   - %d Error%s: %s (%%%.4f)\n' % (key, 's' if key != 1 else ' ',
                                                                  pp(member.error_distribution[key]),
                                                                  member.error_distribution[key] * 100.0 / sum(member.error_distribution.values())))

    sys.stderr.write('\n\nFASTA File for community %s is ready: "%s"\n\n' % (config.community, config.output_file))


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
            member = Member(user_config.get(section, 'sequence'))
            member.id = '-'.join(section.replace('member', '').strip().split())
            member.ratio = int(user_config.get(section, 'ratio'))
            
            config.members.append(member)
            config.members.sort(reverse=True)

    sys.exit(main(config))
