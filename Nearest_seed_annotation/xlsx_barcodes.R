# install.packages("readxl")
# install.packages('seqinr')
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(version = '3.15')
  

# BiocManager::install("Biostrings")

library(readxl)
library(seqinr)
library(dplyr)
library(Biostrings)
library(stringr)
library(data.table)
library(tidyr)
library(utils)

# test reading ----
test <- read_excel("tag_all/Tag list FL_amplicon_1_LM1_correct.xlsx",
                   sheet = 1)

# test["1", 1]

# for LM
test <- test[complete.cases(test), ]

# for 16s
test <- test[!(is.na(test$`Oligo Name`)), ]

nchar(test[1, 5])


# rownames(test)
seq_vec <- c()
seq_names <- c()


for (x in rownames(test)) {
  seq_names <- c(seq_names, test[x, 1])
  seq_vec <- c(seq_vec, paste0("", test[x, 5])) #5
}

length(seq_vec)
substring(seq_vec, 9, 28)

a <- table(substring(seq_vec, 8, 28))
as.vector(a)
str(as.vector(a))
1 %in% a

names(a)
a[names(a) == "ACTCAGCAGCCGCGGTAATAC"] ==1

# https://statisticsglobe.com/extract-values-and-names-from-table-object-in-r
# https://stackoverflow.com/questions/10104400/extract-a-row-from-a-table-object



length(unique(substring(seq_vec, 7, 10)))


# seq_vec
# write.fasta(sequences = as.list(seq_vec),# need as.list to correctly give multiple fasta
#             names = seq_names,
#             file.out = "barcodes/test_barcodes.fasta",
#             open = "w",
#             nbchar = 60,
#             )

# on all ----

# unlist(gregexpr("_correct.xlsx", 
#                 "Tag list FL_amplicon_1_LM1_correct.xlsx")
#        ) # 27
# 
# substr("Tag list FL_amplicon_1_LM1_correct.xlsx",
#        10, 26)
# 
# paste0("1", "2", "fasta")

for (files in list.files("tag_all/")) {
  # print(files)
  string_ind <- unlist(gregexpr("_correct.xlsx",
                                files
                                )
                       ) - 1
  
  seq_vec <- c()
  seq_names <- c()
  
  f_togo <- read_excel(paste0("tag_all/",files),
                       sheet = 1)
  
  if("_LM" %in% files) {
    f_togo <- f_togo[complete.cases(f_togo), ]
  } else {
    f_togo <- f_togo[!(is.na(f_togo$`Oligo Name`)), ]
  }
  
  for (x in rownames(f_togo)) {
    seq_names <- c(seq_names, f_togo[x, 1])
    seq_vec <- c(seq_vec, paste0("", f_togo[x, 5]))
  }
  
  write.fasta(sequences = as.list(seq_vec),# need as.list to correctly give multiple fasta
              names = seq_names,
              file.out = paste0("barcodes/",
                               substr(files, 10, string_ind),
                               ".fasta"),
              open = "w",
              nbchar = 60,
  )
}

# min overalp length test ----

# for forward ----
for (files in list.files("tag_all/")) {
  seq_vec <- c()
  
  f_togo <- read_excel(paste0("tag_all/",files),
                       sheet = 1)
  
  if("_LM" %in% files) {
    f_togo <- f_togo[complete.cases(f_togo), ]
  } else{
    f_togo <- f_togo[!(is.na(f_togo$`Oligo Name`)), ]
  }
  
  for (x in rownames(f_togo)) {
    seq_vec <- c(seq_vec, paste0("", f_togo[x, 3]))
  }
  
  start_ <- 1
  resolution <- length(seq_vec)
  while (resolution == length(seq_vec)) {
    resolution <- length(unique(substring(seq_vec, start_, 10)))
    start_ <- start_ + 1
  }
  print(files)
  print(start_ - 1)
}

# generally you need at least 5 bp barcodes to differntiate all

# however more specifically for every prefix ----

for (files in list.files("tag_all/")) {
  string_ind <- unlist(gregexpr("_correct.xlsx",
                                files
                                )
                       ) - 1
  
  seq_names <- c()
  seq_vec <- c()
  seq_final <- c()
  
  f_togo <- read_excel(paste0("tag_all/",files),
                       sheet = 1)
  
  # "_LM" %in% files
  if(grepl( "_LM", files, fixed = TRUE)) {
    f_togo <- f_togo[complete.cases(f_togo), ]
  } else {
    f_togo <- f_togo[!(is.na(f_togo$`Oligo Name`)), ]
  }

  for (x in rownames(f_togo)) {
    seq_names <- c(seq_names, f_togo[x, 1])
    seq_vec <- c(seq_vec, paste0("", f_togo[x, 5]))
  }
  # print(length(seq_names))
  # print(length(seq_vec))=
  
  if(grepl( "_LM", files, fixed = TRUE)) {
    for (y in seq(1, length(seq_vec))) {
      for (count_ in seq(10, 1)) {
        a <- table(substring(seq_vec, count_, 30))
        
        if (a[names(a) == substr(seq_vec[y], count_, 30)] == 1) {
          seq_final <- c(seq_final,
                         paste0(seq_vec[y],
                                ";min_overlap=",
                                as.character((11 - count_)
                                )
                         )
          )
          break
        }
        
      }
    }
  } else {
    for (y in seq(1, length(seq_vec))) {
      for (count_ in seq(10, 1)) {
        a <- table(substring(seq_vec, count_, 28))
        
        if (a[names(a) == substr(seq_vec[y], count_, 28)] == 1) {
          seq_final <- c(seq_final,
                         paste0(seq_vec[y],
                                ";min_overlap=",
                                as.character((11 - count_)
                                )
                         )
          )
          break
        }
        
      }
    }
  }
  
  
  
  # print(length(seq_vec))
  
  write.fasta(sequences = as.list(seq_final),# need as.list to correctly give multiple fasta
              names = seq_names,
              file.out = paste0("barcodes/",
                                substr(files, 10, string_ind),
                                ".fasta"
                                ),
              open = "w",
              nbchar = 60
              )
}
  


# for reverse ----
# for 19 ----
# GACTACCAGGGTATCTAATCCTGT
nchar("GACTACCAGGGTATCTAATCCTGT") # 24

rev_comp <-  DNAStringSet(c("GACTACCAGGGTATCTAATCCTGT"))

rev_comp <- reverseComplement(rev_comp) %>% as.character()
rev_comp

refs <- read.fasta("E:/adcadmic/bioinfo/docker_qiime2_share_win/FL2/all/19_ref.fasta",
                   seqtype = c("DNA"),
                   as.string = TRUE,
                   forceDNAtolower = FALSE)
# str(refs)
# refs[[1]][1]
# names(refs)[1]
for (x in names(refs)) {
  print(as.character(nchar(refs[[x]][1]
                           )
                     )
        )
}


302 - 5 - 5 - 244

seq_ref <- c()

for (x in seq(1, length(refs))) {
  seq_ref <- c(seq_ref, refs[[x]][1])
}

seq_ref

for (x in seq(23, 1)) {
  tt <- substr(rev_comp, x, 24)
  
  if (!(TRUE %in% str_detect(seq_ref, tt))) {
    print(x)
    break
  }
}#20

substr(rev_comp, 20, 24) # 5


# for LM ----
# frank e-mail ##
# The primers for the LM samples are:
#   
# Lm_amplicon_FL_F            -         Forward         tag-TGAACAATACCAATTTGCCC
# 
# Lm_amplicon_FL_R            -         Reverse         GAGCATTATGCCAAGGTTTG

#
# cat LM_ori_ref.txt | tr -d "[:blank:]" > LM_ori_ref.fasta


LM_refs <- read.fasta("E:/adcadmic/bioinfo/docker_qiime2_share_win/FL2/all/LM_ori_ref.fasta",
                     seqtype = c("DNA"),
                     as.string = TRUE,
                     forceDNAtolower = FALSE,
                     whole.header = TRUE)

for (x in names(LM_refs)) {
  print(as.character(nchar(LM_refs[[x]][1]
                           )
                     )
        )
}

# names(LM_refs)[1]
# nchar(LM_refs[[4]][1])
151*2 - (28 + 24) - 193

seq_to_save <- c()
name_to_save <- c()

for (x in seq(1, length(LM_refs))) {
  seq_to_save <- c(seq_to_save, LM_refs[[x]][1])
  name_to_save <- c(name_to_save, names(LM_refs)[x])
}

# seq_to_save
# name_to_save

write.fasta(sequences = as.list(seq_to_save),# need as.list to correctly give multiple fasta
            names = name_to_save,
            file.out = "LM_ref.fasta",
            open = "w",
            nbchar = 60,
            as.string = TRUE
            )

LM_refs <- read.fasta("LM_ref.fasta",
                      seqtype = c("DNA"),
                      as.string = TRUE,
                      forceDNAtolower = FALSE,
                      whole.header = TRUE)


to_look_for1 <- c("GAGCATTATGCCAAGGTTTG")
nchar(to_look_for1[1])
to_look_for <- DNAStringSet(to_look_for1)
to_look_for <- reverseComplement(to_look_for) %>% as.character()
print(to_look_for)



seq_ref <- c()

for (x in seq(1, length(LM_refs))) {
  seq_ref <- c(seq_ref, LM_refs[[x]][1])
}

# seq_ref

for (x in seq(20, 1)) {
  tt <- substr(to_look_for, x, 20)
  
  if (!(TRUE %in% str_detect(seq_ref, tt))) {
    print(x)
    break
  }
}#16

substr(to_look_for, 16, 20)

# finder ----
x <- DNAString("GAGCATTATGCCAAGGTTTG") %>% reverseComplement %>% as.character()
x
# x <- DNAString("") %>% as.character()
TRUE %in% str_detect(seq_ref, x)

# both refs rewrite with short name ----
# 19 ----
seq_names <- c()
seq_vec <- c()

# a <- names(refs)[17] %>% strsplit(";")
# a[[1]][6]
# refs[[3]][1]

for (x in seq(1, length(names(refs)))) {
  # print(x)
  a <- names(refs)[x] %>% strsplit(";")
  # print(length(a))
  
  seq_names <- c(seq_names, 
                 # paste0(";", a[[1]][7])
                 a[[1]][7])
  seq_vec <- c(seq_vec, refs[[1]][1])
}

write.fasta(sequences = as.list(seq_vec),# need as.list to correctly give multiple fasta
            names = seq_names,
            file.out = "19_sname_ref.fasta",
            open = "w",
            nbchar = 60,
            as.string = TRUE
            )

# LM ----
seq_names <- c()
seq_vec <- c()

# a <- names(refs)[2] %>% strsplit(";")
# a[[1]][7]
# refs[[3]][1]

for (x in seq(1, length(names(LM_refs)))) {
  # print(x)
  a <- names(LM_refs)[x] %>% strsplit(";")
  # print(length(a))
  
  seq_names <- c(seq_names, paste0(";", a[[1]][7]))
  seq_vec <- c(seq_vec, LM_refs[[1]][1])
}

write.fasta(sequences = as.list(seq_vec),# need as.list to correctly give multiple fasta
            names = seq_names,
            file.out = "LM_sname_ref.fasta",
            open = "w",
            nbchar = 60,
            as.string = TRUE
            )

# q2 taxonomy and fasta making ----
eg_taxonomy <- read.csv("qiime2/fileformat_eg/tax.txt",
                        sep = "\t",
                        header = FALSE)
# so it's tab delimited

eg_taxonomy[1, 2]
unlist(gregexpr(" ", eg_taxonomy[1, 2]))
# so there are spaces

# LM ----
LM_refs <- read.fasta("E:/adcadmic/bioinfo/docker_qiime2_share_win/FL2/all/LM_ori_ref.fasta",
                      seqtype = c("DNA"),
                      as.string = TRUE,
                      forceDNAtolower = FALSE,
                      whole.header = TRUE)
# LM_refs[[1]][1]
# gsub(";", "; ",names(LM_refs)[1])
# __Firmicutes; c__Bacilli; o__Bac
# a <- strsplit(names(LM_refs)[1], split = ";")[[1]][7]
# 
# unlist(gregexpr("\\(", a)
#        )
# 
# substr(a, 19, nchar(a) - 1 )


seq_name_1 <- c()
Seqs_1 <- c()

for (x in seq(1, length(names(LM_refs)))) {
  a <- strsplit(names(LM_refs)[x], split = ";")[[1]][7]
  pos <- unlist(gregexpr("\\(", a)
              ) + 1
  a <- substr(a, pos, nchar(a) - 1 )
  seq_name_1 <- c(seq_name_1, a)
  Seqs_1 <- c(Seqs_1, LM_refs[[x]][1])
}

write.fasta(sequences = as.list(Seqs_1),# need as.list to correctly give multiple fasta
            names = seq_name_1,
            file.out = "qiime2/LM_q2_seqs.fasta",
            open = "w",
            nbchar = 60,
            as.string = TRUE
)

# taxonomy
df1 <- data.frame(matrix("", nrow = 6, ncol = 2))

for (x in seq(1, length(names(LM_refs)))) {
  a <- strsplit(names(LM_refs)[x], split = ";")[[1]][7]
  pos <- unlist(gregexpr("\\(", a)
                ) + 1
  a <- substr(a, pos, nchar(a) - 1 )
  df1[x, "X1"] <- a
  df1[x, "X2"] <- gsub(";", "; ",names(LM_refs)[x])
}

write.table(x = df1, "qiime2/LM_q2_tax.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
# 19 ----
refs <- read.fasta("E:/adcadmic/bioinfo/docker_qiime2_share_win/FL2/all/raw data file V4 region extraction 18 Lm 244&245 length.txt",
                   seqtype = c("DNA"),
                   as.string = TRUE,
                   forceDNAtolower = FALSE
                   )

citro <- read.fasta("E:/adcadmic/bioinfo/docker_qiime2_share_win/FL2/all/16S_Citrobacterspp.fasta",
                    seqtype = c("DNA"),
                    as.string = TRUE,
                    forceDNAtolower = FALSE
                    )# 37

raou <- read.fasta("E:/adcadmic/bioinfo/docker_qiime2_share_win/FL2/all/16S_Raoultellaspp.fasta",
                    seqtype = c("DNA"),
                    as.string = TRUE,
                    forceDNAtolower = FALSE
                    )# 35

corrected_37 <- c("GGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAA")
# nchar(corrected_37)

refs[["37"]][1] == corrected_37
refs[["35"]][1] == corrected_37

citro_seqs <- c()
for (x in names(citro)) {
  citro_seqs <- c(citro_seqs, citro[[x]][1])
}

raou_seqs <- c()
for (x in names(raou)) {
  raou_seqs <- c(raou_seqs, raou[[x]][1])
}


str_detect(citro_seqs, refs[["35"]][1])
str_detect(citro_seqs, refs[["37"]][1])

str_detect(raou_seqs, refs[["35"]][1])
str_detect(raou_seqs, refs[["37"]][1])

# nchar(refs[["35"]][1])
# nchar(refs[["37"]][1])

refs[["35"]][1] == refs[["37"]][1]

seq_name_2 <- c()
Seqs_2 <- c()

names(refs)

for (x in seq(1, length(names(refs)))) {
  # a <- strsplit(names(refs)[x], split = ";")[[1]][7]
  # pos <- unlist(gregexpr("\\(", a)
  # ) + 1
  # a <- substr(a, pos, nchar(a) - 1 )
  a <- names(refs)[x]
  if (a == "Listeria") {
    seq_name_2 <- c(seq_name_2, "00")
  } else {
    seq_name_2 <- c(seq_name_2, a)
  }
  
  Seqs_2 <- c(Seqs_2, refs[[x]][1])
}

write.fasta(sequences = as.list(Seqs_2),# need as.list to correctly give multiple fasta
            names = seq_name_2,
            file.out = "qiime2/19_q2_seqs.fasta",
            open = "w",
            nbchar = 60,
            as.string = TRUE
            )

# taxonomy
df2 <- data.frame(matrix("", nrow = 19, ncol = 2))

for (x in seq(1, length(names(refs)))) {
  a <- strsplit(names(refs)[x], split = ";")[[1]][7]
  pos <- unlist(gregexpr("\\(", a)
  ) + 1
  a <- substr(a, pos, nchar(a) - 1 )
  
  if (a == "LM") {
    df2[x, "X1"] <- "00"
  } else {
    df2[x, "X1"] <- a
  }
  
  df2[x, "X2"] <- gsub(";", "; ",names(refs)[x])
}

write.table(x = df2, "qiime2/19_q2_tax.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# sequence similarity ----
LM_refs <- read.fasta("qiime2/LM_q2_seqs.fasta",
                      seqtype = c("DNA"),
                      as.string = TRUE,
                      forceDNAtolower = FALSE,
                      whole.header = TRUE)

refs <- read.fasta("qiime2/19_q2_seqs.fasta",
                   seqtype = c("DNA"),
                   as.string = TRUE,
                   forceDNAtolower = FALSE)


# CJ(names(LM_refs), names(LM_refs), unique = TRUE)
# expand.grid(names(LM_refs),names(LM_refs)) #%>% rownames()
# 
# expand.grid(names(LM_refs),names(LM_refs))[1,1] %>% typeof()
# 
# # crossing(names(LM_refs),names(LM_refs)) %>% rownames()
# 
# LM_refs[["544"]][1]
# 
# palign1 <- pairwiseAlignment(LM_refs[[544]][1] %>% DNAString(), 
#                              LM_refs[[2]][1] %>% DNAString())
# 
# palign1
# pid(palign1)

# min(c(1,2,3))
# LM
coll1 <- c()

df1 <- expand.grid(names(LM_refs),names(LM_refs))

for (x in rownames(df1)) {
  if (df1[x, 1] == df1[x, 2]) {
    next
  } else {
    palign <- pairwiseAlignment(LM_refs[[as.character(df1[x, 1])]][1] %>% DNAString(),
                                LM_refs[[as.character(df1[x, 2])]][1] %>% DNAString())
    coll1 <- c(coll1, pid(palign))
  }
}

min(coll1)# 92.74611
max(coll1)# 97.40933

# 19
comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}

comb(18, 2)

combn(c(1,2,3), 2)

coll2 <- c()

# df2 <- expand.grid(names(refs),names(refs))

df2 <- combn(names(refs), 2) %>% as.data.frame()

for (x in colnames(df2)) {
  if (df2[1, x] == df2[2, x]) {
    next
  } else {
    palign <- pairwiseAlignment(refs[[as.character(df2[1, x])]][1] %>% DNAString(),
                                refs[[as.character(df2[2, x])]][1] %>% DNAString())
    if (pid(palign) == 100) {
      print(x) # 22, 40
    }
    coll2 <- c(coll2, pid(palign))
  }
  
}

df2 <- rbind(df2, coll2 %>% as.data.frame() %>% t())

min(coll2)# 64.89796
max(coll2)# 99.59184


palign <- pairwiseAlignment(LM_refs[[as.character(df1[22, 1])]][1] %>% DNAString(),
                            LM_refs[[as.character(df1[22, 2])]][1] %>% DNAString())

pid(palign,type="PID4")
LM_refs[[as.character(df1[22, 1])]][1] == LM_refs[[as.character(df1[22, 2])]][1]


# df2[40, ]

# 1 to 18 ----
# refs[["152"]][1]

for (x in names(refs)) {
  print(x)
  palign <- pairwiseAlignment(DNAString("AGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGAATGTGAAATCCCTGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTATGGTAGAGGGTAGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACTACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAA"),
                              refs[[x]][1] %>% DNAString())
  print(pid(palign))
}

# meta making ----

# grepl( "_LM", "Metadata_FL_3_LM1.xlsx", fixed = TRUE)

# try remove columns with NA and communicate with Frank later

# strsplit("Metadata_FL_4_16S3.xlsx", "\\.")[[1]][1] %>% substr(10,nchar("Metadata_FL_4_16S3.xlsx"))

for (files in list.files("meta/")) {
  if (grepl( "zip", files, fixed = TRUE)){
    next
  } else {
    meta <- read_excel(paste0("meta/", files),
                       sheet = 1,
                       # col_names = FALSE
                       )
    meta <- meta[ , colSums(is.na(meta)) == 0]
    
    write.table(x = meta, 
                file = paste0("qiime2/", 
                              strsplit(files, "\\.")[[1]][1],
                              "/",
                              strsplit(files, "\\.")[[1]][1],
                              "-metadata.txt"),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)
  }
  
}

# small alignment ----
refs[["37"]][1]

# 5a1b2771718f934391c0b2dd9cadf66f	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
palign1 <- pairwiseAlignment(DNAString("GGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAA"),
                             refs[["35"]][1] %>% DNAString())
pid(palign1)#100

palign2 <- pairwiseAlignment(DNAString("GGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAA"),
                             refs[["35"]][1] %>% DNAString())
pid(palign2)

