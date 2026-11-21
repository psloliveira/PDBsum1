/***********************************************************************

  dcommon.h - Include file for common parameters in DisaStr programs

***********************************************************************/

/* Number of amino acids */
#define NAMINO               20

/* Array parameters */
#define CATH_LEN             30
#define D_FILENAME_LEN      512
#define D_LINELEN          1024
#define D_LONG_LINELEN   100000
#define D_MAX_SEQLEN     100000
#define D_MAX_TOKEN      100000
#define ISOFORM_NAME_LEN     30
#define MAX_ISOFORM         100
#define MAX_PATTERNS          5
#define MAX_SITES             5
#define MAX_STORED           10
#define NAME_LEN            512
#define OMIM_LEN             15
#define UNIPROT_LEN          20

/* Experiment types */
#define XRAY                  0
#define NMR                   1
#define THEORETICAL_MODEL     2
#define OTHER_EXPT            3

static char *EXPERIMENT_TYPE[4] = { 
  "XRAY", "NMR", "THEORETICAL_MODEL", "OTHER"
};

/* Minimum interactions percentage */
#define MIN_PERCEN          ((float) 1.0)

/* Search comparison resuls */
#define D_LOWER              -1
#define D_EQUAL               0
#define D_HIGHER              1

/* Source of variant */
#define D_VARIANT             0
#define D_GNOMAD              1

/* Variant types */
#define NO_VARIANT           -9
#define NAT_VARIANT          -1
#define DISEASE               0
#define DISEASE_NOTE          1
#define MUTAGEN_EXPT          2
#define DISEASE_DDD           3
#define CLINVAR_DATA          4
#define ADDED_VAR             5

#define NVAR_TYPES            6

static int VAR_PRIORITY[NVAR_TYPES + 1] = { 
  0, 5, 2, 1, 3, 4, 6
};

static char *VARIANT_TYPE[NVAR_TYPES] = { 
  "Disease variant", "Disease note", "Mutagenesis experiment",
  "DDD variant", "ClinVar variant", "Report variant"
};
static char *VARIANT_SHORT[NVAR_TYPES] = { 
  "Disease", "Disease note", "Mutagen expt", "DDD", "ClinVar", "Report"
};
static char *VARIANT_DESC[NVAR_TYPES] = { 
  "disease variant", "disease note", "mutagen expt", "DDD variant",
  "clinvar variant", "added variant"
};
static char VARIANT_CHAR[NVAR_TYPES + 1]
= { 'd', 'n', 'm', 'D', 'C', 'A' };

static int VAR_TYPE_COLOUR[NVAR_TYPES + 1] = { 
  BLACK, RED, GREY3, GREEN, BLUE, ORANGE, PURPLE
};

/* Natural variant types */
#define BENIGN                0
#define DAMAGING              1

#define NGNOMAD_TYPES         2

static char *GNOMAD_DESC[NGNOMAD_TYPES] = { 
  "benign variant", "possibly harmful variant"
};
static int NATVAR_TYPE_COLOUR[NGNOMAD_TYPES + 1] = { 
  BLACK, CYAN, ORANGE
};

/* Cutoff for benign/deleterious score */
#define MAX_BENIGN           18

static char *PYMOL_COLOUR[NVAR_TYPES] = { 
  "red", "grey", "green", "blue", "orange", "purple"
};

static int VARIANT_SCORE[NVAR_TYPES] = { 20, 8, 5, 20, 20, 30 };

/* Image dimension subscripts */
#define COLUMN                0
#define ROW                   1

/* Max-min labels */
#define MIN                  0
#define MAX                  1
#define XMIN                 0
#define XMAX                 1
#define YMIN                 2
#define YMAX                 3

static char *MAP_TYPE[4] = { "XMIN", "XMAX", "YMIN", "YMAX" };

/* Letter definitions */
#define CAPITAL_I             9
#define CAPITAL_W            22
#define CAPITAL_X            23

/* Feature types */
#define ANNOT_NONE         -1
#define FT_ACT_SITE         0
#define FT_BINDING          1
#define FT_DISULFID         2
#define FT_DNA_BIND         3
#define FT_DOMAIN           4
#define FT_METAL            5
#define FT_MOTIF            6
#define FT_MUTAGEN          7
#define FT_REGION           8
#define FT_VARIANT          9
#define FT_ZN_FING         10
#define PFAMA              11
#define PFAMB              12
#define CATH               13
#define WHOLE_SEQUENCE     14
#define UNASSIGNED         15
#define SST_COIL           16
#define SST_HELIX          17
#define SST_STRAND         18
#define PDB_MODEL          19
#define PDB_ENTRY          20
#define CAT_RES            21
#define LIGAND_INT         22
#define DNA_INT            23
#define METAL_INT          24
#define ACT_SITE           25
#define DISULPHIDE         26
#define RES_CONS           27
#define RES_CONS_AREA      28
#define DDD_DATA           29
#define PROT_PROT_INT      30
#define GNOMAD             31
#define GNOMAD_AREA        32
#define ADDED_VARIANT      33
#define LOCAL_3D           34
#define LOCAL_3D_AREA      35
#define BINDING_SITE       36
#define COVID_DATA         37
#define AF_MODEL           38

#define NANNOT             39

/* Feature type names */
static char *ANNOT_TYPE[NANNOT] = {
  "ACT_SITE", "BINDING", "DISULFID", "DNA_BIND", "DOMAIN",
  "METAL", "MOTIF", "MUTAGEN", "REGION", "VARIANT", "ZN_FING",
  "Pfam-A", "Pfam-B", "CATH", "Whole", "Unassigned", "Coil",
  "Helix", "Strand", "Structure", "PDB entry", "Catalytic residue",
  "To ligand", "To DNA", "To metal", "Active site", "Disulphide",
  "Residue conservation", "Res cons area", "DDD variant",
  "Prot-prot", "Natural variant", "Natural area", "Local 3D",
  "Local 3D area", "Added variant", "Bind site", "Covid variant",
  "AlphaFold model"
 };

/* Feature description to be used if none supplied */
static char *ANNOT_DESC[NANNOT] = {
  "Active site", "Binding residue(s)", "Disulphide", "DNA binding", "\0",
  "Metal binding", "\0", "\0", "\0", "\0", "Zinc finger",
  "\0", "\0", "\0", "Whole", "Unassigned", "\0",
  "\0", "\0", "\0", "\0", "\0",
  "\0", "\0", "\0", "\0", "\0",
  "\0", "\0", "\0",
  "\0", "\0", "\0", "\0",
  "\0", "\0", "\0", "\0", "\0"
 };

/* Order of priority of domains */
static int DOMAIN_PRIORITY[NANNOT] = {
  0, 0, 0, 7, 4,
  0, 0, 0, 5, 0, 6,
  1, 2, 3, 9, 8, 0,
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0, 0
 };

/* Data in pdb.align file */
#define ALIGN_SEQ_NO         0
#define ALIGN_AA_CODE        1
#define ALIGN_PDB_CODE       2
#define ALIGN_CHAIN          3
#define ALIGN_ENTRY_NO       4
#define ALIGN_DOMAIN         5
#define ALIGN_ISEQ           6
#define ALIGN_PDB_AA         7
#define ALIGN_RES_NAME       8
#define ALIGN_RES_NUM        9
#define ALIGN_START         10
#define ALIGN_END           11

#define NALIGN_FIELDS       12

/* PDB data in pdb.align file */
#define PDB_CODE             0
#define PDB_CHAIN            1
#define PDB_NAME             2
#define PDB_EXPT_TYPE        3
#define PDB_IDENT            4
#define PDB_RESOL            5
#define PDB_RFACT            6
#define PDB_SW               7
#define PDB_E_VALUE_STRING   8
#define PDB_UNIPROT_MATCH    9
#define PDB_UNIPROT_ID      10

#define NPDB_FIELDS         11

/* Columns in covid_variants.tsv file */
#define VAR_NEW_NAME          0
#define VAR_ID                1
#define VAR_DESC              2
#define VAR_PROTEIN           3
#define VAR_UNIPROT_ACC       4
#define VAR_VARIANT           5
#define VAR_SEQ_FROM          6
#define VAR_SEQ_TO            7
#define VAR_AA_CODE           8
#define VAR_AA_MUT            9
#define VAR_REF              10
#define VAR_NOTE             11

#define NVAR_FIELDS          12

/* Contact types */
#define NOT_WANTED           -9
#define TO_PROTEIN            0
#define TO_DNA                1
#define TO_LIGAND             2
#define TO_METAL              3

#define CONTACT_TO            4

/* Contact descriptions */
static const int CONTACT_ANNOT[CONTACT_TO] =
  { PROT_PROT_INT, DNA_INT, LIGAND_INT, METAL_INT };
static char *CONTACT_DESC[CONTACT_TO] =
  { "PROTEIN", "DNA", "LIG", "MET" };
static const char *CONTACT_NAME[CONTACT_TO] =
  { "protein", "DNA/RNA", "ligand", "metal" };
static const char *CONTACT_NAME2[CONTACT_TO] =
  { "another protein chain", "DNA/RNA", "ligand", "metal" };
static const char *CONTACT_FNAME[CONTACT_TO] =
  { "dim", "nuc", "lig", "met" };
static const int START_POS[CONTACT_TO] = { 7, 3, 6, 5 };
static const int CONTACT_SCORE[CONTACT_TO] = { 1, 3, 2, 4 };

  
/* Starting colours for each annotation type */
static int ANNOT_START_COLOUR[NANNOT] = {
  0, 0, 0, 18, 30,
  19, 35, 0, 36, 0, 10,
  0, 10, 1, 0, 46, 0,
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  0, 0, 0,
  0, 0, 0, 0,
  0, 0, 5, 0
 };

/* Residue width for different image magnifications in DisaStr */
static const int RESIDUE_WIDTH[NIMG_MAG] = { 1, 3, 7, 35 };

/* List of valid chromosomes */
static char CHROMO[] = {
  "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M"};

/* List of valid bases */
static char BASES[] = {"A C G T"};

/* ClinVar fields */
#define CLIN_CHROMOSOME           0
#define CLIN_COORDS               1
#define CLIN_ID                   2
#define CLIN_ORIG_BASE            3
#define CLIN_VAR_BASE             4
#define CLIN_QUAL                 5
#define CLIN_FILTER               6
#define CLIN_INFO                 7

#define CLIN_NINPUT_FIELDS        8

/* ClinVar subfields */
#define CLINSUB_UNIPROT_ACC       0
#define CLINSUB_COORDS            1
#define CLINSUB_ID                2
#define CLINSUB_ORIG_BASE         3
#define CLINSUB_VAR_BASE          4
#define CLINSUB_QUAL              5
#define CLINSUB_FILTER            6
#define CLINSUB_INFO              7

#define CLINSUB_NFIELDS           8

#define NCLINSIG                 12
#define SIG_CUTOFF               10

/* ClinVar significance */
#define CLINSIG_OTHER             0
#define CLINSIG_UNKNOWN_SIGNIF    1
#define CLINSIG_UNCERTAIN_SIGNIF  2
#define CLINSIG_PROTECTIVE        3
#define CLINSIG_BENIGN            4
#define CLINSIG_LIKELY_BENIGN     5
#define CLINSIG_DRUG_RESPONSE     6
#define CLINSIG_ASSOCIATION       7
#define CLINSIG_AFFECTS           8
#define CLINSIG_RISK_FACTOR       9
#define CLINSIG_LIKELY_PATHOGENIC 10
#define CLINSIG_PATHOGENIC        11

/* Clinical significance terms in ClinVar data */
static char *CLNSIG[NCLINSIG] = {
  "other",
  "not_provided",
  "Uncertain_significance",
  "protective",
  "Benign",
  "Likely_benign",
  "drug_response",
  "association",
  "Affects",
  "risk_factor",
  "Likely_pathogenic",
  "Pathogenic"
};

/* Modified clinical significance terms */
static char *CLNSIG_MOD[NCLINSIG] = {
  "other",
  "unknown signif",
  "uncertain signif",
  "protective",
  "benign",
  "likely benign",
  "drug response",
  "association",
  "affects",
  "risk factor",
  "likely pathogenic",
  "pathogenic"
};

/* Fields output by API call to VEP */
#define VEP_IDENTIFIER       0
#define VEP_LOCATION         1
#define VEP_ALLELE           2
#define VEP_GENE_ID          3
#define VEP_FEATURE_ID       4
#define VEP_FEATURE_TYPE     5
#define VEP_CONSEQUENCE      6
#define VEP_CDNA_PSTN        7
#define VEP_CDS_PSTN         8
#define VEP_PROTEIN_PSTN     9
#define VEP_AA_CHANGE       10
#define VEP_CODON_CHANGE    11
#define VEP_VARIATION_ID    12
#define VEP_EXTRA_FIELDS    13

#define NVEP_FIELDS         14

static char *VEP_FIELD[NVEP_FIELDS] = {
  "VEP_IDENTIFIER      ",
  "VEP_LOCATION        ",
  "VEP_ALLELE          ",
  "VEP_GENE_ID         ",
  "VEP_FEATURE_ID      ",
  "VEP_FEATURE_TYPE    ",
  "VEP_CONSEQUENCE     ",
  "VEP_CDNA_PSTN       ",
  "VEP_CDS_PSTN        ",
  "VEP_PROTEIN_PSTN    ",
  "VEP_AA_CHANGE       ",
  "VEP_CODON_CHANGE    ",
  "VEP_VARIATION_ID    ",
  "VEP_EXTRA_FIELDS    ",
};

/* Fields read in from isoforms.txt file */
#define ISOFORM_UNIPROT_ACC  0
#define ISOFORM_COUNT        1
#define ISOFORM_ISOFORM      2
#define ISOFORM_RES_START    3
#define ISOFORM_SHIFT        4
#define ISOFORM_NOMAP_LEN    5
#define ISOFORM_ORDER        6

#define NISOFORM_FIELDS      7


/* Fields output by coorun.pl */
#define CHROMOSOME           0
#define COORDS               1
#define USER_BASE            2
#define USER_VARIANT         3
#define ENSEMBL_BASE         4
#define VEP_CODING_BASE      5
#define GENE                 6
#define GENE_ACC             7
#define REFSEQ_GENE_ACC      8
#define TRANSCRIPT           9
#define REFSEQ_TRANSCRIPT   10
#define HGVS_C              11
#define HGVS_P              12
#define STRAND_DIR          13
#define CODON_CHANGE        14
#define VEP_AA              15
#define UNIPROT_AA          16
#define AA_CHANGE           17
#define POLYPHEN_SCORE      18
#define SIFTS_SCORE         19
#define CADD_PHRED          20
#define CADD_RAW            21
#define CADD_MARK           22
#define UNIPROT_ACCESSION   23
#define PROTEIN_NAME        24
#define SEQ_NO              25
#define CHANGE_TYPE         26
#define ALL_TRANSCRIPTS     27
#define NOTE                28
#define GNOMAD_ALLELE_FREQ  29
#define NEGATIVE            30
#define USER_ID             31

#define NINPUT_FIELDS       32

static char *COORUN_FIELD[NINPUT_FIELDS] = {
  "CHROMOSOME          ",
  "COORDS              ",
  "USER_BASE           ",
  "USER_VARIANT        ",
  "ENSEMBL_BASE        ",
  "VEP_CODING_BASE     ",
  "GENE                ",
  "GENE_ACC            ",
  "REFSEQ_GENE_ACC     ",
  "TRANSCRIPT          ",
  "REFSEQ_TRANSCRIPT   ",
  "HGVS_C              ",
  "HGVS_P              ",
  "STRAND_DIR          ",
  "CODON_CHANGE        ",
  "VEP_AA              ",
  "UNIPROT_AA          ",
  "AA_CHANGE           ",
  "POLYPHEN_SCORE      ",
  "SIFTS_SCORE         ",
  "CADD_PHRED          ",
  "CADD_RAW            ",
  "CADD_MARK           ",
  "UNIPROT_ACCESSION   ",
  "PROTEIN_NAME        ",
  "SEQ_NO              ",
  "CHANGE_TYPE         ",
  "ALL_TRANSCRIPTS     ",
  "NOTE                ",
  "GNOMAD_ALLELE_FREQ  ",
  "NEGATIVE            ",
  "USER_ID             ",
};

static char *COORUN_HEADER[NINPUT_FIELDS] = {
  "Chromosome",
  "Chromosome coordinate",
  "User Ref allele",
  "User Alt allele",
  "Ensembl ref allele",
  "VEP coding allele",
  "Gene name",
  "Ensembl gene ID",
  "gene Refseq(HGNC)",
  "Ensembl spliced transcript ID(s) for Uniprot canonical isoform",
  "Transcript refseq(s)(Ensembl)",
  "HGVSc",
  "HGVSp",
  "Strand direction",
  "Codon change",
  "VEP original ref AA",
  "Swissprot ref AA",
  "Amino acid change",
  "Polyphen score",
  "SIFT score",
  "CADD PHRED",
  "CADD raw",
  "CADD mark",
  "Uniprot accession of canonical isoform",
  "Swissprot protein name",
  "Uniprot AA position",
  "VEP Variant consequence",
  "All isoforms",
  "Notes",
  "gnomAD allele freq",
  "Negative T-F",
  "User id"
};

/* Extra fields output to results.dbt file */
#define SYNONYM_FLAG         (NINPUT_FIELDS + 0)

#define NEXTRA_FIELDS        (NINPUT_FIELDS + 1)

/* All-transcript subfields, with added info */
#define SUB_ENST             0
#define SUB_REFSEQ           1
#define SUB_HGVS_C           2
#define SUB_UNIPROT_ACC      3
#define SUB_HGVS_P           4
#define SUB_SEQ_NO           5
#define SUB_AA_CHANGE        6
#define SUB_ASTERISK         7
#define SUB_GENE             8
#define SUB_POLYPHEN         9
#define SUB_SIFTS           10
#define SUB_CADD_PHRED      11
#define SUB_CADD_RAW        12
#define SUB_CONSEQUENCE     13
#define SUB_CODON_CHANGE    14
#define SUB_STRAND_DIR      15

#define SUB_NFIELDS         16

/* All-transcript subfields, without added info */
#define SUB2_ENST             0
#define SUB2_REFSEQ           1
#define SUB2_UNIPROT_ACC      2
#define SUB2_SEQ_NO           3
#define SUB2_AA_CHANGE        4
#define SUB2_CONSEQUENCE      5

#define SUB2_NFIELDS          6

/* gnomAD data fields */
#define GNOMAD_CHROM                 0
#define GNOMAD_POS                   1
#define GNOMAD_REF                   2
#define GNOMAD_ALT                   3
#define GNOMAD_AC                    4
#define GNOMAD_AN                    5
#define GNOMAD_AF                    6
#define GNOMAD_ALLELE                7
#define GNOMAD_CONSEQUENCE           8
#define GNOMAD_IMPACT                9
#define GNOMAD_SYMBOL               10
#define GNOMAD_GENE                 11
#define GNOMAD_FEATURE              12
#define GNOMAD_BIOTYPE              13
#define GNOMAD_HGVSC                14
#define GNOMAD_CDNA_POSITION        15
#define GNOMAD_CDS_POSITION         16
#define GNOMAD_PROTEIN_POSITION     17
#define GNOMAD_AMINO_ACIDS          18
#define GNOMAD_CODONS               19
#define GNOMAD_EXISTING_VARIATION   20
#define GNOMAD_ALLELE_NUM           21
#define GNOMAD_VARIANT_CLASS        22
#define GNOMAD_CANONICAL            23
#define GNOMAD_TSL                  24
#define GNOMAD_ENSP                 25
#define GNOMAD_UNIPROT              26
#define GNOMAD_SIFT                 27
#define GNOMAD_POLYPHEN             28
#define GNOMAD_DOMAINS              29

#define NGNOMAD_FIELDS              30

static char *GNOMAD_FIELD[NGNOMAD_FIELDS] =
  { 
    "CHROM", "POS", "REF", "ALT",
    "AC", "AN", "AF", "ALLELE",
    "CONSEQUENCE", "IMPACT", "SYMBOL", "GENE",
    "FEATURE", "BIOTYPE", "HGVSC", "CDNA_POSITION",
    "CDS_POSITION", "PROTEIN_POSITION", "AMINO_ACIDS", "CODONS",
    "EXISTING_VARIATION", "ALLELE_NUM", "VARIANT_CLASS", "CANONICAL",
    "TSL", "ENSP", "UNIPROT", "SIFT",
    "POLYPHEN", "DOMAINS"
 };

/* gnomAD populations */
#define GNOMAD_ALL_AF                   0
#define GNOMAD_AFR_AF                   1
#define GNOMAD_AMR_AF                   2
#define GNOMAD_ASJ_AF                   3
#define GNOMAD_EAS_AF                   4
#define GNOMAD_FIN_AF                   5
#define GNOMAD_NFE_AF                   6
#define GNOMAD_OTH_AF                   7
#define GNOMAD_SAS_AF                   8

#define NGNOMAD_POPS                    9

static char *GNOMAD_POPULATION[NGNOMAD_POPS] = {
  "gnomAD_AF=",
  "gnomAD_AFR_AF=",
  "gnomAD_AMR_AF=",
  "gnomAD_ASJ_AF=",
  "gnomAD_EAS_AF=",
  "gnomAD_FIN_AF=",
  "gnomAD_NFE_AF=",
  "gnomAD_OTH_AF=",
  "gnomAD_SAS_AF="
};

static char *GNOMAD_POP_ID[NGNOMAD_POPS] = {
  "ALL", "AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "SAS"
};

static char *ALLELE_POP[NGNOMAD_POPS] = {
  "_AF[", "_AFR[", "_AMR[", "_ASJ[", "_EAS[", "_FIN[", "_NFE[",
  "_OTH[", "_SAS["
};

static char *GNOMAD_POP_NAME[NGNOMAD_POPS] = {
  "Combined",
  "African",
  "American",
  "Ashkenazi Jewish",
  "East Asian",
  "Finnish",
  "Non-Finnish European",
  "Other combined",
  "South Asian"
};

/* Fields in the Ensembl to UniProt map file */
#define ENSMAP_UNIPROT_ACC              0
#define ENSMAP_GENE_NAME                1
#define ENSMAP_ENSG                     2
#define ENSMAP_ENST                     3
#define ENSMAP_ENSP                     4
#define ENSMAP_ISOFORM                  5
#define ENSMAP_MAPPING_TYPE             6
#define ENSMAP_SEQ_DIFF                 7
#define ENSMAP_SEQUENCE                 8

#define NENSMAP_FIELDS                  9

/* Fields in HGNC_refseq.tsv file */
#define HGNC_ID                         0
#define HGNC_SYMBOL                     1
#define HGNC_NAME                       2
#define HGNC_REFSEQ_ID                  3
#define HGNC_UNIPROT_ID                 4
#define HGNC_ENSEMBL_GENE_ID            5

/* SIFT scores for variants from gnomAD */
#define SIFT_NONE                      -1
#define SIFT_UNKNOWN                    0
#define SIFT_TOLERATED                  1
#define SIFT_POSS_TOLERATED             2
#define SIFT_POSS_DELETERIOUS           3
#define SIFT_DELETERIOUS                4

#define NSIFT_BANDS                     5

static char *SIFT_BAND[NSIFT_BANDS] = {
  "unknown",
  "tolerated",
  "tolerated_low_confidence",
  "deleterious_low_confidence",
  "deleterious",
};

static float SIFT_SCORE[NSIFT_BANDS] = { 
  1.1, 0.95, 0.90, 0.5, 0.0
};

/* Polyphen scores for variants from gnomAD */
#define POLYPHEN_NONE                  -1
#define POLYPHEN_UNKNOWN                0
#define POLYPHEN_BENIGN                 1
#define POLYPHEN_POSS_BENIGN            2
#define POLYPHEN_POSS_DAMAGING          3
#define POLYPHEN_PROB_DAMAGING          4

#define NPOLYPHEN_BANDS                 5

static char *POLYPHEN_BAND[NPOLYPHEN_BANDS] = {
  "unknown",
  "benign",
  "poss_benign",
  "possibly_damaging",
  "probably_damaging",
};

static float POLYPHEN_VALUE[NPOLYPHEN_BANDS] = { 
  0.0, 0.10, 0.15, 0.85, 1.0
};

/* CADD score markers */
#define NCADD_BANDS                     4
static float CADD_VALUE[NCADD_BANDS] = { 
  0.0, 20.0, 25.0, 30.0
};   // green, orange, purple, red
static int CADD_COLOUR[NCADD_BANDS] = { 
 GREEN, ORANGE, PURPLE, RED
};
#define TOP_CADD_SCORE                 50

/* Different variant types and their representations */
#define NSNP                   0
#define NS_SEQ                 1
#define MISSING                2
#define TERMINATION            3
#define DELETION               4
#define INSERTION              5
#define FRAME_SHIFT            6

#define NMUT_TYPES             7

static char MUT_CODE[NMUT_TYPES] = { ' ', '#', '%', 'X', '-', '+', '!' };
static char *MUT_TYPE_DESC[NMUT_TYPES] = {
  "nsSNP",
  "Variant seq segment",
  "Missing segment",
  "Termination",
  "Deletion",
  "Insertion",
  "Frame shift"
};
static char *MUT_TYPE_DESC_SHORT[NMUT_TYPES] = {
  "?",
  "Seg",
  "Mis",
  "Ter",
  "Del",
  "Ins",
  "Fshft"
};

static char *MUT_TYPE_MESSAGE[NMUT_TYPES] = {
  "?",
  "This variant corresponds to an inserted residue sequence at this position",
  "This variant corresponds to a missing residue sequence at this position",
  "This variant corresponds to a premature stop codon which results"
  " in a truncation of the protein at this position - with likely"
  " loss of the protein's function",
  "This variant corresponds to a residue deletion at this position",
  "This variant corresponds to a residue insertion at this position",
  "This variant corresponds to a frame shift at this position"
};

/* a.a. tables read in from .tsv files */
#define BASE_CHANGES            0
#define GNOMAD_VARIANTS         1
#define UNIPROT_VARIANTS        2

#define NTABLE_TYPES            3

/* Map types */
#define ANATOMAP             0
#define BRAINMAP             1

#define NMAP_TYPES           2

static char *MAP_NAME[NMAP_TYPES] = { "anatoman", "Brainmap" };
static char *MAP_OUT[NMAP_TYPES] = { "anatomap", "brainmap" };
static char MAP_CHAR[NMAP_TYPES] = { 'A', 'B' };
static int MAP_FONT[NMAP_TYPES] = { FONT_CALIBRI_18, FONT_CALIBRI_26 };

/* Tissue types */
#define EXPRESSED            0
#define IN_DISEASE           1

#define NTISSUE_TYPES        2

static char *TISSUE_TYPE[NTISSUE_TYPES] = { "expressed", "disease" };

static char *line_prefix[NMAP_TYPES] = { "AMAP", "BMAP" };

/* Fields in the results.tsv file */
#define RES_CHROMOSOME           0
#define RES_COORDS               1
#define RES_USER_BASE            2
#define RES_USER_VARIANT         3
#define RES_ENSEMBL_BASE         4
#define RES_VEP_CODING_BASE      5
#define RES_GENE                 6
#define RES_GENE_ACC             7
#define RES_REFSEQ_GENE_ACC      8
#define RES_TRANSCRIPT           9
#define RES_REFSEQ_TRANSCRIPT   10
#define RES_HGVS_C              11
#define RES_HGVS_P              12
#define RES_STRAND_DIR          13
#define RES_CODON_CHANGE        14
#define RES_VEP_AA              15
#define RES_UNIPROT_AA          16
#define RES_AA_CHANGE           17
#define RES_POLYPHEN_SCORE      18
#define RES_SIFTS_SCORE         19
#define RES_CADD_PHRED          20
#define RES_CADD_RAW            21
#define RES_CADD_MARK           22
#define RES_UNIPROT_ACCESSION   23
#define RES_PROTEIN_NAME        24
#define RES_SEQ_NO              25
#define RES_CHANGE_TYPE         26
#define RES_ALL_TRANSCRIPTS     27
#define RES_NOTE                28
#define RES_GNOMAD_ALLELE_FREQ  29
#define RES_NEGATIVE            30
#define RES_USER_ID             31
#define RES_SYNONYMOUS          32
#define RES_HAVE_PDB            33
#define RES_PDB_UNIPROT_MATCH   34
#define RES_CLOSEST_PDB_CODE    35
#define RES_PDB_CHAIN           36
#define RES_PDB_PROTEIN_NAME    37
#define RES_PDB_EXPT_TYPE       38
#define RES_PDB_RESOLUTION      39
#define RES_PDB_RFACT           40
#define RES_PDB_UNIPROT_ACC     41
#define RES_PDB_IDENTITY        42
#define RES_PDB_SW_SCORE        43
#define RES_PDB_E_VALUE         44
#define RES_RES_NAME            45
#define RES_RES_NUM             46
#define RES_SST                 47
#define RES_CAT_RES             48
#define RES_DISULPHIDE          49
#define RES_NTO_DNA             50
#define RES_NTO_LIGAND          51
#define RES_NTO_METAL           52
#define RES_NTO_PROTEIN         53
#define RES_NPDB_RES            54
#define RES_LIGANDS             55
#define RES_METALS              56
#define RES_PFAM_DOMAIN         57
#define RES_PFAM_NAME           58
#define RES_CATH_DOMAIN         59
#define RES_CATH_NAME           60
#define RES_RES_CONSERVATION    61
#define RES_NCONS_SEQS          62
#define RES_DISEASES            63
#define RES_DISEASE_VARIANTS    64
#define RES_NVARIANTS           65
#define RES_NAT_VARIANTS        66
#define RES_NOTE_TYPE           67
#define RES_DUPLIC              68

#define NRES_FIELDS             69

/* Amino acid substitution preferences (-4 to +3, as above), taken from:
   M.J. Betts, R.B. Russell. Amino acid properties and consequences of subsitutions.
   In Bioinformatics for Geneticists, M.R. Barnes, I.C. Gray eds, Wiley, 2003.. */
static int AA_SUBST[NAMINO + 1][NAMINO + 1] = {
   0,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2, -5,  // Ala
   0,  0, -3, -4, -1, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, -5,  // Cys
  -2, -3,  0,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -2, -4, -3, -5,  // Asp
  -1, -4,  2,  0, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2, -5,  // Glu
  -2, -2, -3, -3,  0, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3, -5,  // Phe
   0, -3, -1, -2, -3,  0, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3, -5,  // Gly
  -2, -3, -1,  0, -1, -2,  0, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2, -5,  // His
  -1, -1, -3, -3,  0, -4, -3,  0, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1, -5,  // Ile
  -1, -3, -1,  1, -3, -2, -1, -3,  0, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2, -5,  // Lys
  -1, -1, -4, -3,  0, -4, -3,  2, -2,  0,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1, -5,  // Leu
  -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  0, -2, -2,  0, -1, -1, -1,  1, -1, -1, -5,  // Met
  -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  0, -2,  0,  0,  1,  0, -3, -4, -2, -5,  // Asn
  -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  0, -1, -2, -1, -1, -2, -4, -3, -5,  // Pro
  -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  0,  1,  0, -1, -2, -2, -1, -5,  // Gln
  -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  0, -1, -1, -3, -3, -2, -5,  // Arg
   1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  0,  1, -2, -3, -2, -5,  // Ser
   0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  0,  0, -2, -2, -5,  // Thr
   0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  0, -3, -1, -5,  // Val
  -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3,  0,  2, -5,  // Trp
  -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  0, -5,  // Tyr
  -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5  // Ter
};

/* Covid-19 proteins */
#define NCOVID19          27

static char *COVID19_NAME[NCOVID19][3] = {
  {"P0DTC2", "S - Spike glycoprotein", "S"},
  {"P0DTC4", "E - Envelope small membrane protein", "E"},
  {"P0DTC5", "M - Membrane protei", "M"},
  {"P0DTC9", "N - Nucleocapsid phosphoprotein", "N"},
  {"P0DTD1_1", "NSP1 - Mediates RNA replication and processing", "nsp1"},
  {"P0DTD1_2", "NSP2 - Function unknown", "nsp2"},
  {"P0DTD1_3", "NSP3 - Involved in viral replication", "nsp3"},
  {"P0DTD1_4", "NSP4 - Involved in viral replication", "nsp4"},
  {"P0DTD1_5", "NSP5 - Main protease", "nsp5"},
  {"P0DTD1_6", "NSP6 - Involved in viral replication", "nsp6"},
  {"P0DTD1_7", "NSP7 - Part of RNA polymerase", "nsp7"},
  {"P0DTD1_8", "NSP8 - Part of RNA polymerase", "nsp8"},
  {"P0DTD1_9", "NSP9 - ssRNA binding", "nsp9"},
  {"P0DTD1_10", "NSP10 - Essential for NSP16 methyltransferase activity", "nsp10"},
  {"P0DTD1_11", "NSP11 - Short peptide", "nsp11"},
  {"P0DTD1_12", "NSP12 - RNA polymerase", "nsp12"},
  {"P0DTD1_13", "NSP13 - Helicase/triphosphatase", "nsp13"},
  {"P0DTD1_14", "NSP14 - 3'-5' exonuclease", "nsp14"},
  {"P0DTD1_15", "NSP15 - Uridine-specific endoribonuclease", "nsp15"},
  {"P0DTD1_16", "NSP16 - RNA-cap methyltransferase", "nsp16"},
  {"P0DTC3", "Protein 3a - Activates the NLRP3 inflammasome", "ORF 3a"},
  {"P0DTC6", "Protein 6 - Type 1 IFN antagonist", "ORF 6"},
  {"P0DTC7", "Protein 7a - Virus-induced apoptosis", "ORF 7a"},
  {"P0DTD8", "Protein 7b - Unknown function", "ORF 7b"},
  {"P0DTC8", "Protein 8 - Unknown function", "ORF 8"},
  {"P0DTD2", "Protein 9b - Type 1 IFN antagonist", "ORF 9b"},
  {"A0A663DJA2", "Protein 10 - Unknown function}", "ORF 10"}
};

/* Data files for indexing */
#define UNIPROT_ALL_FASTA    0
#define UNIPROT_ID_FASTA     1
#define UNIPROT_ENSEMBL_MAP  2
#define ENSEMBL_REFSEQ       3
#define HGNC_REFSEQ          4

#define NINDEX_FILES         5

static const char *INDEX_FILE_NAME[NINDEX_FILES] = {
  "uniprot_all.fasta",
  "uniprot_all.fasta",
  //  "uniprot_ensembl_mapping_2019_08v2.csv",
  //"human_sequences_mappings_Nov_20.csv",
  "human_sequences_mappingsV4.csv",
  "ensembl_refseq.tsv",
  //  "HGNC_refseq.tsv"
  "HGNC_refseq_Nov_2020.tsv"
};

static const char INDEX_TOKEN_SEP[NINDEX_FILES]
   = { '|', '|', ',', '\t', '\t' };

static const int INDEX_MIN_TOKENS[NINDEX_FILES] = { 3, 3, 4, 3, 6 };

static const int INDEX_KEY_TOKEN[NINDEX_FILES] = { 1, 2, 3, 1, 5 };

static const char INDEX_END_CHAR[NINDEX_FILES]
   = { '\0', ' ', '\0', '\0', '\0' };



/* Structure for each gnomad data item
   ---------------------------------- */
struct gnomad
{
  char               *field[NGNOMAD_FIELDS];
  struct gnomad      *next_gnomad_ptr;
};

/* Structure for each idx data item
   -------------------------------- */
struct idx
{
  char               *key;
  long               offset;
  int                nentry;
  struct idx         *next_idx_ptr;
};

/* Structure for variant for which CADD data have been written out
   --------------------------------------------------------------- */
struct cadd
{
  char                    *variant_id;
  float                   cadd_raw;
  float                   cadd_phred;
  struct cadd             *next_cadd_ptr;
};

/* Structure for isoform records
   ----------------------------- */
struct isoform
{
  char               *uniprot_acc;
  int                count;
  int                isoform;
  int                res_start;
  int                shift;
  int                nomap_len;
  int                order;
  struct isoform     *next_isoform_ptr;
};

