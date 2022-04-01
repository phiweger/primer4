    def map_to_genomic_manually(self, db):
        # TODO: Not congruent w/ automatic for NM_015015.2:c.2441+1G>A
        name = self.tx
        c_pos = self.base
        offset = self.offset

        coding = []
        for i in db.children(f'rna-{name}', featuretype='CDS', order_by='start'):
            coding.append(i)
    
        strand = db[f'rna-{name}'].strand
        l = 0
        if strand == '-':
            coding = coding[::-1]  # reverse CDS order
            
            for i in coding:
                if l + len(i) >= c_pos:
                    break
                else:
                    l += len(i)
            # i.end holds genomic coordinate
            g_pos = i.end - (c_pos - l) + 1
        
        else:
            for i in coding:
                if l + len(i) >= c_pos:
                    break
                else:
                    l += len(i)
            g_pos = i.start + c_pos - 1

        if offset:
            # The offset is in relation to the coding sequence, eg "+3" means
            # three bases downstream of the previous coding sequence.
            if strand == '-':
                g_pos -= offset
            else:
                g_pos += offset

        # click.echo(log(f'Variant on chromosome {chromosome}, g. position {g_pos}'))
        return g_pos


    def load_variation2_(self, databases):
        for name, db in databases.items():
            variants = VariantFile(db)

            # dbSNP names chromosomes like "NC_000007.13", others like "7"
            if name != 'dbSNP':
                chrom = convert_chrom(self.feat.chrom)
            else:
                chrom = self.feat.chrom
            vv = variants.fetch(chrom, self.feat.start, self.feat.end)
            self.variation[name] = vv

            for i in vv:
                # .info.get(...) raises ValueError: Invalid header if not there
                info = dict(i.info)
                
                if name == 'dbSNP':
                    if info.get('COMMON'):
                        self.mask.add(i.pos - self.feat.start - 1)
                
                elif name == '1000Genomes':
                    if info['AF'][0] >= 0.01:
                        self.mask.add(i.pos - self.feat.start - 1)
                
                elif name == 'ESP':
                    if float(info['MAF'][0]) >= 1:
                        self.mask.add(i.pos - self.feat.start - 1)

                else:
                    print(f'"{name}" is not a valid variant database')
        return None










def get_mrna(transcript, feature_db):







# TODO: Have a universal translate fn that translates coords btw/ genomic, transcript, coding and relative -- these need to be object methods of the Template obj.


for i in list(db.children(f'rna-{v.tx}', featuretype='CDS', order_by='start')):
    print(i.frame)




def qpcr(template, feature_db, params):
    '''
    One primer inside the exon, one outside
    '''
    pass


'''
1 inside, 1 outside

should not be too difficult
'''


s = mrna(tmp, db)



mrna = db[f'rna-{v.tx}']

def get_exon(db, name, exon):
    try:
        return db[f'exon-{name}-{exon}']
    except gffutils.exceptions.FeatureNotFoundError:
        return None


# https://pythonhosted.org/gffutils/autodocs/gffutils.FeatureDB.create_introns.html
def get_intron(db, name, exon, relative=-1):



qry = 3
tx = v.tx

ex = retrieve_exon(db, tx, qry)



# https://github.com/seandavi/GFFutils#imputing-introns
exons = {}

for e in db.children('rna-NM_000546.6', featuretype='exon', order_by='start'):
    exons[int(e.id.split('-')[-1])] = e

mrna = ''
for k in sorted(exons.keys()):
    mrna += exons[k].sequence(genome)
    # Validated manually using screed that this considers strand, ie for "-"
    # we get the revcomp sequence.
'''
As a sanity check I can the composed sequence in Blastn and it returned:
"Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA" with 100%
identity and 100% query cover.
'''

# qrna
exons = list(db.children('rna-NM_000546.6', featuretype='exon', order_by='start'))
introns = list(db.interfeatures(exons))
# Naturally, we have more exons than introns
assert len(exons) == len(introns) + 1




l = twolists(exons, introns)
ix = [int(i.id.split('-')[-1]) if i.id else None for i in l]

ix.index(1)
# and then l[20+1] and l[20-1]









def qpcr(template, feature_db, params, n_exon):
    '''
    qpcr(tmp, db, params, 5)
    '''
    mn, mx = params['size_range_qPCR']

    exons = list(feature_db.children(
        template.feat.id, featuretype='exon', order_by='start'))
    introns = list(feature_db.interfeatures(exons))
    l = twolists(exons, introns)
    ix = [int(i.id.split('-')[-1]) if i.id else None for i in l]
    mid = ix.index(n_exon)
    left = mid - 1 
    right = mid + 1

    l[left].start
    rlb = template.relative_pos(l[left].start)  # rlb .. relative left bound
    rmb = template.relative_pos(l[mid].start)  # rmb .. mid

    return {
        'only_here': (
            (rlb, len(l[left])),
            (rmb, len(l[mid]))
            ),
        'size_range': (mn, mx)
        }


'''
db = gffutils.create_db(
    fp, 
    dbfn='hg19-p13_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True)
'''
# https://github.com/daler/gffutils/issues/111
introns = list(db.create_introns())
# 819,463
db.update(
    introns,
    disable_infer_transcripts=True,
    disable_infer_genes=True,
    verbose=True,
    merge_strategy='merge')
# backup=True


# f'exon-{name}-{exon}
# intron = db['exon-NR_026818.1-2', 'exon-NR_026818.1-3']


def mrna(template, feature_db, params):
    '''
    One primer in exon 1, the other in exon 2
    '''
    pass


def get_contraints(x):
    return x.start, x.end - x.start



# TODO: closest features
# https://github.com/seandavi/GFFutils#closest-features


'''
input: just name exon for now, later bold on the syntax to parse automatically

NM_000546.6, exon 6

class mRNA()

has a dict {1: (start, end)}

then select exon + and - (exon +- 1)

place primer there (constrains, see syntax of the primer3 param)

possibility to select only one splice site (eg +1 or -1 if no arg take both)
'''







# -----------------------------------------------------------------------------


'''
Problems encountered/ solved:

- DSL
- SNVs database preparations
- SNVs different chromosome namings and options fields (AF vs. MAF vs. COMMON)
- recursion necessary for multiple pairs (otherwise get 10000 takes forever but bc/ combinatorics still in the same place)
- automatic transcript conversion (no tests yet)
- test suite


also:

visit mibi



'''




# -----------------------------------------------------------------------------

'''
hdp = JSONDataProvider(['cdot-0.2.1.refseq.grch37_grch38.json.gz'])
db = gffutils.FeatureDB('/Users/phi/Dropbox/repos/primer/data/hg19-p13_annotation.db', keep_order=True)
genome = Fasta('/Users/phi/Dropbox/repos/primer/data/GRCh37_latest_genomic.fna')


import json
with open('/Users/phi/Dropbox/repos/primer/config.json', 'r') as file:
    params = json.load(file)


# v = Variant('NM_015015.2:c.2441+1G>A', hdp, db)
v = Variant('NM_000546.6:c.215C>G', hdp, db)
# ex = context(v, db, 'exon')
t = Template(v, db)

# TODO: Template() needs to work w/ ('NM_000546.6', '4') exon coords, too
# We have genomic coords in this case, so good.
# -- Create another class Exon() and do have the same interface, but do the
# manual coord conversion in there.

# Actually, we can parse this rather easily:
# 
# g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
# c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
# 
# Genomic we can look up, coding we'd have to translate using existing code.



constraints = mask_sanger(v, t, params)
# TODO: the mask_x fn should not need the variant
# ((0, 7544), (7882, 19070))
design_primers('PCR', params, t.get_sequence(genome), constraints)


# TODO design(search_space, template, params)



method = 'mRNA'

before = neighbor(ex, db, -1)
after = neighbor(ex, db, 1)


design(template, placements, mask, params)


placements = set(
    {'left': (10, 15), 'right': (45, 89)},
    {'left': (45, 89), 'right': (91, 98)},
    )

# then pass each to primer3



Target(v, ex, params, 'mRNA')

'''








##fileformat=VCFv4.2
##fileDate=20210513
##source=dbSNP
##dbSNP_BUILD_ID=155
##reference=GRCh37.p13
##phasing=partial

# def remove_overlapping(primers, length):
#     occupied = [0] * length
#     for k, v in primers.items():
#         for i in ['fwd', 'rev']:
#             start, end = v[i]['start'], v[i]['end']
#             overlap = sum(occupied[start:end])
#             if overlap < 10:
#                 print(k, i, overlap, v['penalty'])
#             # Not +1, primer3 py package has pythonic coords
#             for j in range(start, end + 1):
#                 occupied[j] = 1



















# Sanger
method = 'sanger'
code = 'NM_000546.6:c.215C>G'      # -- strand -, no offset
# code = 'NM_015015.3:c.2441+1G>A'   # -- strand +, offset 1
# code = 'NM_000546.6:c.672+3C>G'    # -- strand -, offset 3

v = Variant(code, hdp, db)
tmp = Template(v, db)
# TODO: Add this to tests
if tmp.feat.strand == '-':
    assert tmp.c_to_g[v.start] - v.start_offset == v.g_start
else:
    assert tmp.c_to_g[v.start] + v.start_offset == v.g_start
# ((0, 7544), (7882, 11188))
# (250, 600)
tmp.load_variation_(vardbs)
# masked = tmp.mask_sequence(genome)

masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)

# tmp.mask_sequence(genome, unmasked='-')[:1000]
constraints = tmp.apply(method, db, params)
primers = [p for p in next(design_primers(masked, constraints, params, []))]
# [7461-7481:8019-8039, loss: 0.1454,
#  7521-7541:8076-8096, loss: 0.2891,
#  7408-7429:7940-7960, loss: 1.7632,
#  7498-7520:7996-8017, loss: 4.1442]



'''
# qPCR, starting from HGVS code
method = 'qpcr'
code = 'NM_000546.6:c.(?_560-1)_(672+1_?)del'
v = ExonDelta(code, db)

if len(v.data) == 1:
    exon = next(iter(v.data))
    n_exon = int(exon.id.split('-')[-1])

# For qPCR, the fwd primer can be in the 5' intron or exon, the rev in exon or
# 3' intron, so we have two sets of constraints.
primers = []
for constraints in tmp.apply('qpcr', db, params, n_exon):
    x = [p for p in next(design_primers(masked, constraints, params, []))]
    primers.extend(x)
'''


# qPCR, starting from explicit exon annotation
method = 'qpcr'
# code = ('NM_000546.6', 5)
code = ('NM_001145408.2', 6)

v = SingleExon(*code)
tmp = Template(v, db)
tmp.load_variation_(vardbs)
# masked = tmp.mask_sequence(genome)
masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)

primers = []
for constraints in tmp.apply('qpcr', db, params):
    print(constraints)
    x = [p for p in next(design_primers(masked, constraints, params, []))]
    primers.extend(x)
# [10527-10548:10612-10632, loss: 1.1868,
#  10568-10592:10667-10686, loss: 5.3035,
#  10860-10880:10990-11010, loss: 0.4314,
#  10822-10842:10944-10964, loss: 2.2463,
#  10783-10803:10896-10918, loss: 2.619]


# mRNA

# reconstruct mRNA as template
# NM_000546.6:c.215C>G
code = ('NM_000546.6', 6, 7)
method = 'mrna'

v = ExonSpread(*code)
tmp = Template(v, db)

tmp.mrna = reconstruct_mrna(tmp.feat, db, genome, vardbs)
constraints = tmp.apply('mrna', db, params)
masked = tmp.mrna[0]
'''
As a sanity check I can the composed sequence in Blastn and it returned:
"Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA" with 100%
identity and 100% query cover.
'''
primers = [p for p in next(design_primers(masked, constraints, params, []))]
# [777-797:888-908, loss: 0.5132,
#  706-727:816-836, loss: 1.474,
#  735-755:861-879, loss: 6.427]



# primers.save('foo')







'''
makeblastdb -in GRCh37_latest_genomic.fna -dbtype nucl

blastn -dust no -word_size 7 -evalue 10 -outfmt 6 -query ../test.fna -db GRCh37_latest_genomic.fna -out result


import pandas as pd
names = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(' ')
df = pd.read_csv('result', sep='\t', names=names)

df[]

> I have logic that does round-robin pairwise comparisons (i.e. compares every hit to every other hit) to determine if each hit pair could result in an amplicon (on same chromosome, opposite-stranded, perfect match at 3'-ends, plus-strand hit is upstream of minus-strand hit, are separated by <= 1000 bp).


https://bioinformatics.stackexchange.com/questions/7262/hits-in-primer-blast-not-found-with-programmatic-blastn-query
https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
https://www.metagenomics.wiki/tools/blast/megablast

blastn -task megablast -dust no -word_size 15 -evalue 2 -outfmt "6 sam" -query ../test.fna -db GRCh37_latest_genomic.fna -out result -perc_identity 0.9 -strand both && cat result



'''




# TODO: now back translate primer coord from blast


results = check_for_multiple_amplicons(primers, fp_genome)

tmp.g_to_c[7578213]







'''
TODO: 

We now could use pybedtools 

https://daler.github.io/pybedtools/index.html

nearby = genes.closest(intergenic_snps, ...)

if primer coord (start, end) in exon, use g_to_c, else use pybedtools get_closest and then minus the corresponding start coord

nearby = genes.closest(intergenic
'''



'''
{'qseqid': '0605c189-3b95-49ac-a4d3-bc52947a8f0c.fwd',
 'sseqid': 'NC_000002.11',
 'pident': 100.0,
 'length': 16,
 'mismatch': 0,
 'gapopen': 0,
 'qstart': 1,
 'qend': 16,
 'sstart': 71186481,
...

'''



'''
    tmp = tempfile.TemporaryDirectory()
    p = tmp.name
    print(f'Aligning {Path(target).name} to {Path(query).name}')

    steps = [
        f'foldseek createdb {target} {p}/targetDB',
        f'foldseek createdb {query} {p}/queryDB',
        f'foldseek search {p}/queryDB {p}/targetDB {p}/aln {p}/tmp -a --cov-mode {mode} --tmscore-threshold {minscore}',
        f'foldseek aln2tmscore {p}/queryDB {p}/targetDB {p}/aln {p}/aln_tmscore',
        f'foldseek createtsv {p}/queryDB {p}/targetDB {p}/aln_tmscore {p}/aln_tmscore.tsv'
    ]

    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)
'''





import ast
import click

class PythonLiteralOption(click.Option):
    '''
    TODO: Subcommands vs weird command line:
    
    - https://click.palletsprojects.com/en/8.0.x/commands/
    - https://stackoverflow.com/questions/47631914/how-to-pass-several-list-of-arguments-to-click-option
    '''
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)