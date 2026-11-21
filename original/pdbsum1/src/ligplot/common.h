/**********************************************************************

  common.h - Include file for common global vraibles

***********************************************************************/

/* Standard system include files */
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

/* Include file for reading gzipped files */
#ifdef TEST
#else
#include <zlib.h>
#endif

/* Array parameters */
#define C_CATH_LEN           30
#define C_LINELEN          1024
#define C_FILENAME_LEN     1024
#define C_MAX_TOKEN        1000
#define MAX_PDB_CODE_LEN      4
#define C_STRING_LEN        512
#define C_MAXDOMCHAIN       200
#define C_DOMLEN             16

/* Flag-values */     
#define FALSE                0
#define TRUE                 1
#define UNSET                2

static char *T_F[3] = { "FALSE", "TRUE", "*UNSET*" };

typedef short  BOOL;

/* Directory separator */
static char DIR_SEP[] = "/";

/* Operating system */
#define LINUX                0
#define WINDOWS              1
#define MAC                  2

#define NOP_SYS              3

static char *OS[NOP_SYS] = { "linux", "win", "mac" };

/* Path format */
#define NOT_DOS              0
#define DOS                  1

#define NPATH_FORMATS        2

/* O/S commands */
static char *CP_COPY[NPATH_FORMATS] = { "cp", "copy" };
static char *MV_REN[NPATH_FORMATS] = { "mv", "ren" };
static char *RM_DEL[NPATH_FORMATS] = { "rm", "del" };
static char *SLEEP[NPATH_FORMATS] = { "sleep 1", "timeout /t 1" };

/* Special characters */
static char null = '\0';
static const char *TAB = "\t";

/* Array variables */
#define PPM_LINELEN         512
#define PPM_MAX_TOKEN       100

/* Field for storing current date and time 
   Access by:  time(&current_time);
               printf("Current time: %s\n",ctime(&current_time)); */
//time_t current_time;

/* Months of the year */
static char *MONTH[12] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
  "Nov", "Dec" };

/* Operating system we are running on */
#define LINUX_32              0
#define LINUX_64              1

#define DISULPHIDE_DIST      3.0
#define DISULPHIDE_DIST_SQRD (DISULPHIDE_DIST * DISULPHIDE_DIST)

/* PDB file types */
#define FIND_FIRST           -9
#define NOT_FOUND            -1
#define STANDARD_PDB          0
#define UPLOAD_PDB            1
#define MCSG_PDB              2
#define ALPHA_FOLD_PDB        3
#define PDBSUM1               4

/* PDB code length */
#define PDB_CODELEN          4

/* Template characters (length must match PDB code length above) */
static char TEMPLATE_CHARS[] = "abcd";

/* PDB file locations */
#define PDB_CURRENT          0
#define PDB_OBSOLETE         1
#define PDB_MODELS           2
#define PDB_MODELS_OBSOLETE  3
#define PDB_MISC             4
#define PDB_UPLOAD           5
#define PDB_MCSG             6
#define PDB_ALPHA_FOLD       7
#define PDB_PDBSUM_DIR       8

#define NPDB_DIR             9

/* Names of the PDB directories in CATHPARAM */
static char *PDB_DIR_TEMPLATE[NPDB_DIR] = { 
  "FTP_PDB_TEMPLATE", "FTP_OBSOLETE_TEMPLATE",
  "FTP_MODELS_TEMPLATE", "FTP_OBSOLETE_MODELS_TEMPLATE",
  "MISC_PDB_TEMPLATE", "UPLOAD_PDB_TEMPLATE", "MCSG_PDB_TEMPLATE",
  "ALPHAFOLD_PDB_TEMPLATE", "PDBSUM_PDB_TEMPLATE"
};

/* PDBsum file locations */
#define PDBSUM_CURRENT       0
#define PDBSUM_UPLOAD        1

#define NPDBSUM_DIR          2

/* Names of the PDBsum directory templates */
static char *PDBSUM_DIR_TEMPLATE[NPDBSUM_DIR] = { 
  "PDBSUM_REAL_TEMPLATE", "UPLOAD_PDBSUM_TEMPLATE"
};

/* Names of the PDBsum data directories */
static char *PDBSUM_DATA_DIR[NPDBSUM_DIR] = { 
  "PDBSUM_REAL_DATA_DIR", "UPLOAD_PDBSUM_DATA_DIR"
};

/* Executable file locations */
#define ACROREAD_EXE         0
#define BABEL_EXE            1
#define CONVERT_EXE          2
#define CRYPT_EXE            3
#define ESEARCH_PROG         4
#define FASTA                5
#define GHOSTSCRIPT_EXE      6
#define HBADD_EXE            7
#define HBPLUS_EXE           8
#define HSEARCH_PROG         9
#define JAVA_EXE            10
#define LIGPLOT_EXE         11
#define LIGSEARCH_PYMOL_EXE 12
#define NACCESS_EXE         13
#define NCBI_FETCH          14
#define NEW_PYMOL_EXE       15
#define OLD_CONVERT_EXE     16
#define PDFEXTRACT_EXE      17
#define PDFIMAGES_EXE       18
#define PPMQUANT_EXE        19
#define PPMTOGIF_EXE        20
#define PROSITE_SCAN        21
#define PSTOIMG_EXE         22
#define PSTOTEXT_EXE        23
#define PS2PDF              24
#define PYMOL_EXE           25
#define PYMOL_64_EXE        26
#define RASTER3D            27
#define RUNROMLAS_CMND      28
#define SPEEDFILL_EXE       29
#define SVMTEST_EXE         30
#define SVMTORCH_EXE        31
#define SW_JAVA             32
#define VIEWFIGS_PL         33
#define WGET                34

#define NEXE_DIR            35

/* Names of the PDBsum directories in CATHPARAM */
static char *EXE_NAME[NEXE_DIR] = { 
  "ACROREAD_EXE", "BABEL_EXE", "CONVERT_EXE", "CRYPT_EXE", "ESEARCH_PROG",
  "FASTA",
  "GHOSTSCRIPT_EXE", "HBADD_EXE", "HBPLUS_EXE", "HSEARCH_PROG", "JAVA_EXE",
  "LIGPLOT_EXE", "LIGSEARCH_PYMOL_EXE", "NACCESS_EXE", "NCBI_FETCH",
  "NEW_PYMOL_EXE", "OLD_CONVERT_EXE", "PDFEXTRACT_EXE",
  "PDFIMAGES_EXE", "PPMQUANT_EXE", "PPMTOGIF_EXE", "PROSITE_SCAN",
  "PSTOIMG_EXE", "PSTOTEXT_EXE", "PS2PDF", "PYMOL_EXE", "PYMOL_64_EXE",
  "RASTER3D", "RUNROMLAS_CMND", "SPEEDFILL_EXE", "SVMTEST_EXE",
  "SVMTORCH_EXE", "SW_JAVA", "VIEWFIGS_PL", "WGET"
};

/* Other file locations */
#define ALPHAFOLD_PDB_DIR                         0
#define ALPHAFOLD_PDBSUM_DATA_DIR                 1
#define ARCHCATH_OLD_REAL_DATA_DIR                2
#define ARCHCATH_REAL_DATA_DIR                    3
#define ARCHSCHEMA_REAL_DATA_DIR                  4
#define ATLAS_REAL_DIR                            5
#define BINDDB_DIR                                6
#define CADD_DATA_DIR                             7
#define CADD_DATA_DIR_37                          8
#define CADD_DATA_DIR_38                          9
#define CATH_ALIGNMENTS_DIR                      10
#define CATH_DATA_DIR                            11
#define CATH_LIST                                12
#define CATH_DOMAINS                             13
#define CATH_TITLES                              14
#define CITEXPLORE_CLASSPATH                     15
#define COORUN_DIR                               16
#define CULLED_DIR                               17
#define DIR_SEPARATOR                            18
#define DISASTR_REAL_DATA_DIR                    19
#define DISASTR_REAL_DIR                         20
#define DISASTR_REAL_DOMAINS_DIR                 21
#define DISASTR_REAL_PATHWAYS_DIR                22
#define DISASTR_REAL_UNIPROT_DIR                 23
#define DISASTR_REAL_TEMPLATES_DIR               24
#define DISASTR_REAL_TEMPLATES_PPM_DIR           25
#define DNA_REPORTS_DIR                          26
#define DNA_TEMPLATES_DIR                        27
#define DOMALL_FILE                              28
#define DRAT_REAL_DATA_DIR                       29
#define DRUGBANK_DIR                             30
#define DRUGBANK_XML                             31
#define DRUGCARDS_FILE                           32
#define DRUGPORT_REAL_DATA_DIR                   33
#define DRUGPORT_REAL_DIR                        34
#define EBI                                      35
#define EBI_HOME_PAGE                            36
#define EBI_INCLUDE_DIR                          37
#define EHEADERS_FILE                            38
#define ENSEMBL_DIR                              39
#define ENTRIES_DIR                              40
#define ENTRIES_REAL_DIR                         41   
#define ENZ_REPORTS_DIR                          42
#define ENZ_TEMPLATES_DIR                        43
#define ENZYMES_GIF_DIR                          44
#define ENZYMES_HOME_PAGE                        45
#define ENZYMES_MOLGIF_DIR                       46
#define ENZYMES_REAL_DIR                         47
#define ENZYMES_REAL_DATA_DIR                    48   
#define ENZYMES_TEMPLATES_DIR                    49
#define EXPASY_ENZYME_DIR                        50
#define FAKE_PDB_ENTRY                           51
#define FTP_PDB_MMCIF_TEMPLATE                   52
#define FTP_PFAM_DIR                             53
#define GENE3D_DATA_DIR                          54
#define GENOMES_LOCAL_DIR                        55
#define GENOMES_LOCAL_DAT_DIR                    56
#define GET_ECPAGE_PL                            57
#define GET_PDBSUM_ENTRY                         58
#define GET_PYMOL_FILE                           59
#define GET_SASPAGE_PL                           60
#define GET_UNIPROT_DATA                         61
#define GNOMAD_DIR                               62
#define GNOMAD_FILE                              63
#define GOA_DIR                                  64
#define GOA_LINK                                 65
#define GO_LOCAL_DIR                             66
#define GO_NAMES_FILE                            67
#define GO_ONTOLOGY_DIR                          68
#define HETGIF_DIR                               69
#define HETGIF_REAL_DIR                          70
#define HETGROUP_DICTIONARY                      71
#define HETGROUPS_DIR                            72
#define HETGROUPS_REAL_DIR                       73
#define HINX_DNA_TEMPLATES_DIR                   74
#define HINX_ENZ_TEMPLATES_DIR                   75
#define HINX_LIG_TEMPLATES_DIR                   76
#define HINX_TEMPLATES_DIR                       77
#define HTTP_BASE                                78
#define JAX_LIB                                  79
#define KEGG_COMPOUND_FILE                       80
#define KEGG_DATA_DIR                            81
#define KEGG_DRUG_FILE                           82
#define KEGG_GIF_DIR                             83
#define KEGG_MOL_DIR                             84
#define LIG_REPORTS_DIR                          85
#define LIG_TEMPLATES_DIR                        86
#define LIGANDS_DIR                              87
#define LIGANDS_REAL_DIR                         88
#define LIG_CLUSTERS_REAL_DIR                    89
#define LIGPLOT_BIN_DIR                          90
#define LIGSEARCH_REAL_RESULTS                   91
#define LOCAL_PATTERNS_FILE                      92
#define NEW_CONSURF_DIR                          93
#define NISTHET_NAME                             94
#define NISTLIST_NAME                            95
#define OPERATING_SYSTEM                         96
#define PATHLIST_TXT                             97
#define PDB_XML_DIR                              98
#define PDBLIB_FILE                              99
#define PDBSUM_CGI_DIR                          100
#define PDBSUM1_DATA_DIR                        101
#define PDBSUM_DIR                              102
#define PDBSUM_DOCS_DIR                         103
#define PDBSUM_EXE_DIR                          104
#define PDBSUM_GENDATA_DIR                      105
#define PDBSUM_GIF_DIR                          106
#define PDBSUM_HOME_PAGE                        107
#define PDBSUM_PARAMS_DIR                       108
#define PDBSUM_PROFUNC_DIR                      109
#define PDBSUM_REAL_DATA_DIR                    110
#define PDBSUM_REAL_TEMPLATES_DIR               111
#define PDBSUM_TEMPLATE                         112
#define PDBSUM_TEMPLATES_DIR                    113
#define PDBSUM_TEMPLATES_GIF_DIR                114
#define PDBSUM1_RESULTS_DIR                     115
#define PFAM_DATA_DIR                           116
#define PQS_DIR                                 117
#define PQS_DATA_DIR                            118
#define PROCHECK_EXE_64_DIR                     119
#define PROFUNC_CGI_DIR                         120
#define PROFUNC_DATA_DIR                        121
#define PROFUNC_PARAMS_DIR                      122
#define PROFUNC_TEMPLATES_DIR                   123
#define PROSITE_FILE                            124
#define PROSITE_PATTERNS_FILE                   125
#define PS2PDF_OUT                              126
#define PUBMED_REAL_DIR                         127
#define RESULTS_DIR                             128
#define S20REP_DOMAINS                          129
#define SAS_CGI_DIR                             130
#define SAS_CSS_DIR                             131
#define SAS_EXE_DIR                             132
#define SAS_GIF_DIR                             133
#define SAS_HOME_PAGE                           134
#define SAS_REAL_TEMPLATES_DIR                  135
#define SPECIES_DIR                             136
#define SPECIES_HTML_FILE                       137
#define SWISS_CODE                              138
#define SWISS_NAME                              139
#define SWISSPROT_DIR                           140
#define SWISSPROT_LOCAL_DIR                     141
#define TAXONOMY_DIR                            142
#define TMP_DIR                                 143
#define TOF_MATCHES_FILE                        144
#define UNIPROT_DIR                             145
#define UPLOAD_PDBLIB_FILE                      146
#define UPLOAD_PROSITE_FILE                     147
#define UPLOAD_TOF_MATCHES_FILE                 148
#define URL_EBI                                 149
#define URL_GNOMAD                              150
#define URL_GXA_HOME                            151
#define URL_GXA_SEARCH                          152
#define URL_GXA_UNIPROT                         153
#define URL_NIST_DIRECT                         154
#define URL_PFAM                                155
#define URL_PFAM_CODE                           156
#define URL_UNIPROT                             157
#define VARIANTS_REAL_TEMPLATE                  158
#define VARIANTS_TEMPLATE                       159
#define VARSITE_REAL_DATA_DIR                   160
#define VARSITE_DATA_DIR                        161
#define VARSITE_REAL_DIR                        162
#define VARSITE_REAL_DOMAINS_DIR                163
#define VARSITE_REAL_PATHWAYS_DIR               164
#define VARSITE_REAL_UNIPROT_DIR                165
#define VARSITE_REAL_TEMPLATES_DIR              166
#define VARSITE_REAL_TEMPLATES_PPM_DIR          167
#define VARSITE_UPLOAD_DIR                      168
#define ZPDB_TO_GO_FILE                         169
#define ZSPROT_DAT                              170
#define ZSPROT_FASTA                            171
#define ZTREMBL_FASTA                           172
#define ZUNIPROT_TO_GO_FILE                     173

#define NOTHER_DIR                              174

/* Names of the PDBsum directories in CATHPARAM */
static char *OTHER_DIR_NAME[NOTHER_DIR] = { 
  "ALPHAFOLD_PDB_DIR", "ALPHAFOLD_PDBSUM_DATA_DIR",
  "ARCHCATH_OLD_REAL_DATA_DIR", "ARCHCATH_REAL_DATA_DIR",
  "ARCHSCHEMA_REAL_DATA_DIR", "ATLAS_REAL_DIR", "BINDDB_DIR",
  "CADD_DATA_DIR", "CADD_DATA_DIR_37", "CADD_DATA_DIR_38",
  "CATH_ALIGNMENTS_DIR", "CATH_DATA_DIR",
  "CATH_LIST", "CATH_DOMAINS", "CATH_TITLES",
  "CITEXPLORE_CLASSPATH", "COORUN_DIR", "CULLED_DIR",
  "DIR_SEP", "DISASTR_REAL_DATA_DIR", "DISASTR_REAL_DIR",
  "DISASTR_REAL_DOMAINS_DIR", "DISASTR_REAL_PATHWAYS_DIR",
  "DISASTR_REAL_UNIPROT_DIR",
  "DISASTR_REAL_TEMPLATES_DIR",  "DISASTR_REAL_TEMPLATES_PPM_DIR",
  "DNA_REPORTS_DIR", "DNA_TEMPLATES_DIR", "DOMALL_FILE",
  "DRAT_REAL_DATA_DIR", "DRUGBANK_DIR", "DRUGBANK_XML", "DRUGCARDS_FILE",
  "DRUGPORT_REAL_DATA_DIR", "DRUGPORT_REAL_DIR", "EBI", "EBI_HOME_PAGE",
  "EBI_INCLUDE_DIR", "EHEADERS_FILE", "ENSEMBL_DIR", "ENTRIES_DIR",
  "ENTRIES_REAL_DIR",
  "ENZ_REPORTS_DIR", "ENZ_TEMPLATES_DIR", "ENZYMES_GIF_DIR",
  "ENZYMES_HOME_PAGE", "ENZYMES_MOLGIF_DIR", "ENZYMES_REAL_DIR",
  "ENZYMES_REAL_DATA_DIR", "ENZYMES_TEMPLATES_DIR", "EXPASY_ENZYME_DIR",
  "FAKE_PDB_ENTRY", "FTP_PDB_MMCIF_TEMPLATE", "FTP_PFAM_DIR",
  "GENE3D_DATA_DIR", "GENOMES_LOCAL_DIR",
  "GENOMES_LOCAL_DAT_DIR", "GET_ECPAGE_PL", "GET_PDBSUM_ENTRY",
  "GET_PYMOL_FILE", "GET_SASPAGE_PL", "GET_UNIPROT_DATA",
  "GNOMAD_DIR", "GNOMAD_FILE",
  "GOA_DIR", "GOA_LINK", "GO_LOCAL_DIR", "GO_NAMES_FILE",
  "GO_ONTOLOGY_DIR", "HETGIF_DIR", "HETGIF_REAL_DIR",
  "HETGROUP_DICTIONARY", "HETGROUPS_DIR", "HETGROUPS_REAL_DIR",
  "HINX_DNA_TEMPLATES_DIR", "HINX_ENZ_TEMPLATES_DIR",
  "HINX_LIG_TEMPLATES_DIR", "HINX_TEMPLATES_DIR",
  "HTTP_BASE", "JAX_LIB", "KEGG_COMPOUND_FILE",
  "KEGG_DATA_DIR", "KEGG_DRUG_FILE", "KEGG_GIF_DIR", "KEGG_MOL_DIR",
  "LIG_REPORTS_DIR", "LIG_TEMPLATES_DIR",
  "LIGANDS_DIR", "LIGANDS_REAL_DIR", "LIG_CLUSTERS_REAL_DIR",
  "LIGPLOT_BIN_DIR", "LIGSEARCH_REAL_RESULTS", "LOCAL_PATTERNS_FILE",
  "NEW_CONSURF_DIR", "NISTHET_NAME", "NISTLIST_NAME", 
  "OPERATING_SYSTEM", "PATHLIST_TXT",
  "PDB_XML_DIR", "PDBLIB_FILE", "PDBSUM_CGI_DIR", "PDBSUM1_DATA_DIR",
  "PDBSUM_DIR", "PDBSUM_DOCS_DIR", "PDBSUM_EXE_DIR", "PDBSUM_GENDATA_DIR",
  "PDBSUM_GIF_DIR", "PDBSUM_HOME_PAGE", "PDBSUM_PARAMS_DIR",
  "PDBSUM_PROFUNC_DIR", "PDBSUM_REAL_DATA_DIR",
  "PDBSUM_REAL_TEMPLATES_DIR", "PDBSUM_TEMPLATE", "PDBSUM_TEMPLATES_DIR",
  "PDBSUM_TEMPLATES_GIF_DIR", "PDBSUM1_RESULTS_DIR", "PFAM_DATA_DIR",
  "PQS_DIR", "PQS_DATA_DIR", "PROCHECK_EXE_64_DIR", "PROFUNC_CGI_DIR",
  "PROFUNC_DATA_DIR", "PROFUNC_PARAMS_DIR", "PROFUNC_TEMPLATES_DIR",
  "PROSITE_FILE", "PROSITE_PATTERNS_FILE", "PS2PDF_OUT",
  "PUBMED_REAL_DIR", "RESULTS_DIR", "S20REP_DOMAINS",
  "SAS_CGI_DIR", "SAS_CSS_DIR", "SAS_EXE_DIR", "SAS_GIF_DIR",
  "SAS_HOME_PAGE", "SAS_REAL_TEMPLATES_DIR", "SPECIES_DIR",
  "SPECIES_HTML_FILE", "SWISS_CODE", "SWISS_NAME", "SWISSPROT_DIR",
  "SWISSPROT_LOCAL_DIR", "TAXONOMY_DIR", "TMP_DIR", "TOF_MATCHES_FILE",
  "UNIPROT_DIR", "UPLOAD_PDBLIB_FILE", "UPLOAD_PROSITE_FILE",
  "UPLOAD_TOF_MATCHES_FILE", "URL_EBI", "URL_GNOMAD", "URL_GXA_HOME",
  "URL_GXA_SEARCH", "URL_GXA_UNIPROT", "URL_NIST_DIRECT",
  "URL_PFAM", "URL_PFAM_CODE",
  "URL_UNIPROT", "VARIANTS_REAL_TEMPLATE", "VARIANTS_TEMPLATE",
  "VARSITE_REAL_DATA_DIR", "VARSITE_DATA_DIR", "VARSITE_REAL_DIR",
  "VARSITE_REAL_DOMAINS_DIR", "VARSITE_REAL_PATHWAYS_DIR",
  "VARSITE_REAL_UNIPROT_DIR",
  "VARSITE_REAL_TEMPLATES_DIR",  "VARSITE_REAL_TEMPLATES_PPM_DIR",
  "VARSITE_UPLOAD_DIR", "ZPDB_TO_GO_FILE",
  "ZSPROT_DAT", "ZSPROT_FASTA", "ZTREMBL_FASTA", "ZUNIPROT_TO_GO_FILE"
};

/* Coordinate subscripts */
#define X                    0
#define Y                    1
#define Z                    2

static char COORD[3] = { 'x', 'y', 'z' };

/* Large and small constants */
#ifndef INFINITY
#define INFINITY   (9999999999.9)
#endif
#ifndef INT_INFINITY
#define INT_INFINITY   INT_MAX
#endif
#define TINY       (0.0000000001)
#define ZERO                 0

/* Angle constants */
#define PI             ((double)  3.141592654)
#define RADDEG         ((double)  (180.0 / PI))
#define PI_BY_2        ((double)  (PI / 2.0))
#define TWO_PI         ((double)  (2.0 * PI))
#define SQRT_TWO_PI    ((double)  (sqrt(TWO_PI)))
#define ROOT2          ((double)  sqrt(2.0))
#define ONE_OVER_ROOT2 ((double)  (1.0 / ROOT2))

/* Other constants */
#define E_VALUE        ( exp(1) )
#define E_VALUE_SQRD   ( E_VALUE * E_VALUE )

/* Open output file flags */
#define TERMINATE      TRUE
#define WARN_ONLY      FALSE

/* RGB colour codes */
#define BLACK                 0
#define WHITE             ( 256 * 256 * 256 - 1 )
#define OFF_WHITE         ( WHITE - 1 )
#define RED               ( 256 * 256 * 255 )
#define GREEN             ( 256 * 255 )
#define BLUE              ( 255 )
#define BROWN             ( 256 * 256 * 204 + 256 * 153 +   0 )
#define CORNFLOWER_BLUE   ( 256 * 256 * 100 + 256 * 149 + 237 )
#define CREAM             ( 256 * 256 * 255 + 256 * 255 + 179 )
#define CYAN              ( 256 * 256 *   0 + 256 * 255 + 255 )
#define DARK_SLATE_BLUE   ( 256 * 256 *  72 + 256 *  61 + 139 )
#define DEEP_PINK         ( 256 * 256 * 255 + 256 *  20 + 147 )
#define DEEP_SKY_BLUE     ( 256 * 256 *   0 + 256 * 191 + 255 )
#define OLIVE_DRAB        ( 256 * 256 * 107 + 256 * 142 +  35 )
#define GOLD3             ( 256 * 256 * 205 + 256 * 173 +   0 )
#define GREY1             ( 256 * 256 * 128 + 256 * 128 + 128 )
#define GREY2             ( 256 * 256 * 168 + 256 * 168 + 168 )
#define GREY3             ( 256 * 256 * 188 + 256 * 188 + 188 )
#define GREY4             ( 256 * 256 * 208 + 256 * 208 + 208 )
#define GREY5             ( 256 * 256 * 228 + 256 * 228 + 228 )
#define LAWN_GREEN        ( 256 * 256 * 124 + 256 * 252 +   0 )
#define MAGENTA           ( 256 * 256 * 255 + 256 * 0   + 255 )
#define MIDNIGHT_BLUE     ( 256 * 256 *  25 + 256 *  25 + 112 )
#define ORANGE            ( 256 * 256 * 255 + 256 * 165 +   0 )
#define PINK              ( 256 * 256 * 255 + 256 * 105 + 180 )
#define PURPLE            ( 256 * 256 * 204 + 256 * 102 + 255 )
#define SANDY_BROWN       ( 256 * 256 * 244 + 256 * 164 +  96 )
#define SKY_BLUE          ( 256 * 256 * 135 + 256 * 206 + 235 )
#define SLATE_BLUE        ( 256 * 256 * 106 + 256 *  90 + 205 )
#define STEEL_BLUE        ( 256 * 256 *  70 + 256 * 130 + 180 )
#define YELLOW            ( 256 * 256 * 255 + 256 * 255 +   0 )


/* Cutoff for deciding whether a colour is a shade of grey */
#define GREY_CUTOFF                 100
#define STDEV_CUTOFF                 20
#define SPAN_CUTOFF                 20

/* Shadow parameters */
#define NSHADOW_COLOURS              5

/* Image sizes */
#define SMALL_IMAGE                  0
#define MEDIUM_IMAGE                 1
#define LARGE_IMAGE                  2
#define EXTRA_LARGE_IMAGE            3

#define NIMG_MAG                     4

static const int SHADOW_RANGE[NIMG_MAG] = { 2, 3, 4, 8 };

#define MAX_RANGE                    8

static const int SHADOW_COLOUR[NSHADOW_COLOURS]
= { GREY1, GREY2, GREY3, GREY4, GREY5 };

static const int NEMPTY[NIMG_MAG][NSHADOW_COLOURS] =
  { {0, 1, 2, 3, 4},
    {0, 1, 2, 3, 4},
    {0, 1, 2, 3, 4},
    {0, 3, 5, 6, 7 }};

static const int APPLY_SHADE[NIMG_MAG][MAX_RANGE] =
  { {0, 2, 4, 4, 4, 4, 4, 4},
    {0, 2, 3, 4, 4, 4, 4, 4},
    {0, 1, 3, 4, 4, 4, 4, 4},
    {0, 0, 0, 1, 1, 3, 3, 4 }};

/* Chain colours */
#define NCHAIN_COLOURS        8
static char *CHAIN_LIST
= "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890!#_-=:<>|abcdefghijklmnopqrstuvwxyz";

static char *CHAIN_COLOUR[NCHAIN_COLOURS]
= {
  "purple", "red", "brown", "pink",
  "blue", "green", "cyan", "yellow" 
};

/* Defined colours */
#define NCOLOURS             24
static char *COLOUR_LIST[] = {
  "black", "red", "green", "blue", "yellow", "purple", "brown",
  "sky_blue", "skyblue", "orange", "cyan", "grey", "white", "magenta",
  "pink", "dark_grey", "cream", "paleblue","greenblue", "redorange",
  "violet", "peach", "limegreen", "olivegreen"
};
static double RGB_TRIPLET[NCOLOURS][3] = {
  {0.000,  0.000,  0.000},  /* black */
  {1.000,  0.000,  0.000},  /* red */
  {0.000,  1.000,  0.000},  /* green */
  {0.000,  0.000,  1.000},  /* blue */
  {1.000,  1.000,  0.000},  /* yellow */
  {0.800,  0.400,  1.000},  /* purple */
  {0.800,  0.600,  0.000},  /* brown */
  {0.000,  0.600,  1.000},  /* sky_blue */
  {0.000,  0.600,  1.000},  /* skyblue */
  {1.000,  0.400,  0.000},  /* orange */
  {0.200,  0.800,  1.000},  /* cyan */
  {0.500,  0.500,  0.500},  /* grey */
  {1.000,  1.000,  1.000},  /* white */
  {0.938,  0.000,  0.500},  /* magenta */
  {1.000,  0.412,  0.706},  /* pink */
  {0.200,  0.200,  0.200},  /* dark_grey */
  {1.000,  1.000,  0.700},  /* cream */
  {0.680,  0.850,  0.900},  /* paleblue */
  {0.200,  0.600,  0.400},  /* greenblue */
  {1.000,  0.400,  0.000},  /* redorange */
  {0.930,  0.510,  0.930},  /* violet */
  {1.000,  0.800,  0.600},  /* peach */
  {0.500,  1.000,  0.000},  /* limegreen */
  {0.412,  0.545,  0.133}   /* olivegreen */
};

/* Colours for CATH domains */
#define NCATH_COLOURS              8

static int CATH_COLOUR[NCATH_COLOURS] = {
  45,  // Black
  1,   // Red
  2,   // Blue
  0,   // Green
  37,  // SaddleBrown
  23,  // Orange
  5,   // Cyan
  13  // HotPink
};

/* Amino acid colours by property */
static char AMINO[]={"ACDEFGHIKLMNPQRSTVWYX"};

static char *AMINO_CODE[] =
  {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "UNK"
  };

static char *AMINO_NAME[] =
  {
    "Alanine", "Cysteine", "Aspartic acid", "Glutamic acid", "Phenylalanine",
    "Glycine", "Histidine", "Isoleucine", "Lysine", "Leucine",
    "Methionine", "Asparagine", "Proline", "Glutamine", "Arginine",
    "Serine", "Threonine", "Valine", "Tryptophan", "Tyrosine",
    "Unknown"
  };

static char *AMINO_COLOUR_NAME[] = { 
  "DimGrey", "Orange", "Red", "Red", "Purple", "SaddleBrown", "Blue", "DimGrey",
  "Blue", "DimGrey", "DimGrey", "Green", "SaddleBrown", "Green", "Blue",
  "Green", "Green", "DimGrey", "Purple", "Purple", "Black" };

static char *AMINO_HTML_COLOUR_NAME[] = { 
  "gray", "orange", "red", "red", "purple", "brown", "blue", "gray",
  "blue", "gray", "gray", "green", "brown", "green", "blue",
  "green", "green", "gray", "purple", "purple", "black" };

/* Colours for Pfam domains */
#define NPFAM_COLOURS             46
static char *COLOUR_DEFN[NPFAM_COLOURS + 1][3] = {
  {"Green",           "#00FF00", "0 255 0"},     //  0 - A
  {"Red",             "#FF0000", "255 0 0"},     //  1 - B
  {"Blue",            "#0000FF", "0 0 255"},     //  2 - C
  {"Gold",            "#FFD700", "255 215 0"},   //  3 - D
  {"Purple",          "#A020F0", "160 32 240"},  //  4 - E
  {"Cyan",            "#00FFFF", "0 255 255"},   //  5 - F
  {"Firebrick",       "#B22222", "178 34 34"},   //  6 - G
  {"LimeGreen",       "#32CD32", "50 205 50"},   //  7 - H
  {"DeepSkyBlue",     "#00BFFF", "0 191 255"},   //  8 - I
  {"LightGoldenrod1", "#FFEC8B", "255 236 139"}, //  9 - J
  {"OrangeRed",       "#FF4500", "255 69 0"},    // 10 - K
  {"LightGreen",      "#90EE90", "144 238 144"}, // 11 - L
  {"SkyBlue",         "#87CEEB", "135 206 235"}, // 12 - M
  {"HotPink",         "#FF69B4", "255 105 180"}, // 13 - N
  {"SeaGreen",        "#2E8B57", "46 139 87"},   // 14 - O
  {"MidnightBlue",    "#191970", "25 25 112"},   // 15 - P
  {"DarkRed",         "#8B0000", "139 0 0"},     // 16 - Q
  {"Peru",            "#CD853F", "205 133 63"},  // 17 - R
  {"SpringGreen",     "#00FF7F", "0 255 127"},   // 18 - S
  {"Magenta",         "#FF00FF", "255 0 255"},   // 19 - T
  {"DimGrey",         "#696969", "105 105 105"}, // 20 - U
  {"DarkGoldenrod",   "#B8860B", "184 134 11"},  // 21 - V
  {"Maroon",          "#B03060", "176 48 96"},   // 22 - W
  {"Orange",          "#FFA500", "255 165 0"},   // 23 - X
  {"OliveDrab",       "#6B8E23", "107 142 35"},  // 24 - Y
  {"CornflowerBlue",  "#6495ED", "100 149 237"}, // 25 - Z
  {"DeepPink",        "#FF1493", "255 20 147"},  // 26 - a
  {"DarkSlateBlue",   "#483D8B", "72 61 139"},   // 27 - b
  {"SandyBrown",      "#F4A460", "244 164 96"},  // 28 - c
  {"DarkKhaki",       "#BDB76B", "189 183 107"}, // 29 - d
  {"Tomato",          "#FF6347", "255 99 71"},   // 30 - e
  {"SteelBlue",       "#4682B4", "70 130 180"},  // 31 - f
  {"Chocolate",       "#D2691E", "210 105 30"},  // 32 - g
  {"PaleGreen",       "#98FB98", "152 251 152"}, // 33 - h
  {"MediumOrchid",    "#BA55D3", "186 85 211"},  // 34 - i
  {"DarkSlateGrey",   "#2F4F4F", "47 79 79"},    // 35 - j
  {"Salmon",          "#FA8072", "250 128 114"}, // 36 - k
  {"SaddleBrown",     "#8B4513", "139 69 19"},   // 37 - l
  {"SlateGrey",       "#708090", "112 128 144"}, // 38 - m
  {"LawnGreen",       "#7CFC00", "124 252 0"},   // 39 - n
  {"Burlywood",       "#DEB887", "222 184 135"}, // 40 - o
  {"DarkGreen",       "#006400", "0 100 0"},     // 41 - p
  {"SlateBlue",       "#6A5ACD", "106 90 205"},  // 42 - q
  {"DarkOliveGreen",  "#556B2F", "85 107 47"},   // 43 - r
  {"Grey",            "#BEBEBE", "190 190 190"}, // 44 - s
  {"Black",           "#000000", "0 0 0"},       // 45 - t
  {"White",           "#FFFFFF", "255 255 255"}  // 46 - u
};
static char COLOUR_CODE[]
= "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstu";
#define BLACK_COL            45
#define WHITE_COL            46
#define PFAMA_STARTCOL        0
#define PFAMB_STARTCOL       10
#define HALF_TEXT           300.0

/* Colour of dots indicating break in sequence */
static int RGB_BREAK[3] = { 64, 64, 64 };

/* Colour of labels showing names of Pfam domains */
static int RGB_TEXT[3] = { 170, 170, 170 };

/* Het Groups that aren't drugs, despite being in DrugBank */
static char NOT_DRUGS_LIST[] = "|010|ACY|CIT|DMS|GLY|GOL|NAG|SAM|TRS|";


/* Fonts used when generating .ppm image files */
#define FONT_TIMES_NEW_ROMAN_8    0
#define FONT_TIMES_NEW_ROMAN_40   1
#define FONT_CALIBRI_7            2
#define FONT_CALIBRI_8            3
#define FONT_CALIBRI_9            4
#define FONT_CALIBRI_10           5
#define FONT_CALIBRI_11           6
#define FONT_CALIBRI_12           7
#define FONT_CALIBRI_13           8
#define FONT_CALIBRI_14           9
#define FONT_CALIBRI_16          10
#define FONT_CALIBRI_18          11
#define FONT_CALIBRI_20          12
#define FONT_CALIBRI_24          13
#define FONT_CALIBRI_25          14
#define FONT_CALIBRI_26          15


#define FONT_CALIBRI_40          16
#define FONT_CALIBRI_50          17
#define FONT_CALIBRI_90          18
#define FONT_COURIER_13          19
#define FONT_NU_7                20
#define FONT_COURIER_NUMBERS_7   21
#define FONT_COURIER_NUMBERS_20  22
#define FONT_COURIER_NUMBERS_35  23
#define FONT_SCHOOLBOOK_11       24
#define FONT_GREEK_10            25
#define FONT_GREEK_50            26

#define NFONT_TYPES              27

static char *FONT_FILE[NFONT_TYPES]
= {
  "font8", "font40", "font_calibri7", "font_calibri8", "font_calibri9", 
  "font_calibri10", "font_calibri11", "font_calibri12", "font_calibri13", 
  "font_calibri14", "font_calibri16", "font_calibri18",
  "font_calibri20", "font_calibri24", "font_calibri25", "font_calibri26",
  "font_calibri40", "font_calibri50", "font_calibri90",
  "font_courier13",
  "font_nu7", "fontnum7",  "fontnum20", "fontnum35", "font_school11",
  "greek10", "greek50"
};

/* Letters */
static char LETTER[]
= "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

/* Number of letter shades */
#define NSHADES              20

/* Minimum allele frequency */
#define MIN_ALLELE_FREQ       10

/* Covid protein definitions */
#define NCOVID_PROTEINS      28

static char *COVID_PROTEIN[NCOVID_PROTEINS]
= {
  "P0DTC2",
  "P0DTC3",
  "P0DTC4",
  "P0DTC5",
  "P0DTC6",
  "P0DTC7",
  "P0DTC8",
  "P0DTC9",
  "P0DTD1_1",
  "P0DTD1_2",
  "P0DTD1_3",
  "P0DTD1_4",
  "P0DTD1_5",
  "P0DTD1_6",
  "P0DTD1_7",
  "P0DTD1_8",
  "P0DTD1_9",
  "P0DTD1_10",
  "P0DTD1_11",
  "P0DTD1_12",
  "P0DTD1_13",
  "P0DTD1_14",
  "P0DTD1_15",
  "P0DTD1_16",
  "A0A663DJA2",
  "P0DTD2",
  "P0DTD",
  "P0DTD8"
 };

/* Alpha Fold confidence scores */
#define NAF_RANGES             4

/* Upper limits of each score range */
static int AF_SCORE[NAF_RANGES] = { 50, 70, 90, 100 };

/* Colours for each score range */
static char *AF_COLOUR[NAF_RANGES] = {
  "#FF7D45", "#FFDB13", "#65CBF3", "#0053D3" };
static char *AF_COLOUR_RGB[NAF_RANGES] = {
  "1.000, 0.490, 0.270",
  "1.000, 0.860, 0.080",
  "0.400, 0.800, 0.950",
  "0.000, 0.330, 0.840"
};
static int AF_COLOUR_IRGB[NAF_RANGES][3] = {
  { 255, 125,  69 },
  { 255, 219,  20 },
  { 102, 204, 242 },
  {   0,  84, 214 }
};
static char *AF_COLOUR_NAME[NAF_RANGES] = {
  "Tomato", "Gold", "SkyBlue", "MidnightBlue" };
static char *AF_COLOUR_CONF[NAF_RANGES] = {
  "VeryLowConf", "LowConf", "Confident", "VeryHigh" };

static int AF_COLOUR_NUM[NAF_RANGES] = { 30, 3, 12, 2 };

/* Structure for common directory names
   ------------------------------------ */
struct c_file_data
{
  char                    *pdb_template[NPDB_DIR];
  char                    *pdbsum_exe_dir;
  char                    *pdbsum_gen_dir;
  char                    *pdbsum_template[NPDBSUM_DIR];
  char                    *pdbsum_data_dir[NPDBSUM_DIR];
  char                    *exe_name[NEXE_DIR];
  char                    *other_dir[NOTHER_DIR];
  int                     dos;
};

/* Structure for the given font
   ---------------------------- */
struct font
{
  char                    *name;
  int                     ncharacters;
  int                     height;
  char                    *character;
  struct letter           **letter_index_ptr;
};

/* Structure for letters
   --------------------- */
struct letter
{
  char                    character;
  int                     height;
  int                     width;
  int                     miny;
  int                     maxy;
  char                    *array;
};

/* Structure for CATH domain name and description
   ---------------------------------------------- */
struct c_cath
{
  char                    *code;
  char                    *name;
  struct c_cath           *next_c_cath_ptr;
};

/* Structure for domain record
   --------------------------- */
struct c_domain
{
  int                     domain_no;
  struct c_cath           *c_cath_ptr;
  struct c_domlist        *first_c_domlist_ptr;
  struct c_domain         *next_c_domain_ptr;
};

/* Structure for domain segment
   ---------------------------- */
struct c_domlist
{
  int                     segment_no;
  struct c_domain         *c_domain_ptr;
  struct c_domlist        *next_c_domlist_ptr;
};

/* Structure for Alpha Fold scores
   ------------------------------- */
struct c_score
{
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  double                  score;
  int                     range;
  struct c_score          *next_c_score_ptr;
};

