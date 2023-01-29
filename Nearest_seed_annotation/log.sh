# gunzip notmatching and concat ----
    # /data/docker_qiime2_share_container_E/FL2/all/not_matching
    parallel -j 10 gunzip {} \
    ::: `ls /data/docker_qiime2_share_container_E/FL2/all/need_cat`

    ######### cat1
    cat_list1=""
    cat_list2=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/need_cat`
    do
    if grep -q "amplicon_2" <<< $x && grep -q "1.fastq" <<< $x
    then
    cat_list1+="/data/docker_qiime2_share_container_E/FL2/all/need_cat/$x "
    elif grep -q "amplicon_2" <<< $x && grep -q "2.fastq" <<< $x
    then
    cat_list2+="/data/docker_qiime2_share_container_E/FL2/all/need_cat/$x "
    fi
    done

    cat_list1=${cat_list1::-1}
    cat_list2=${cat_list2::-1}

    echo $cat_list1
    echo $cat_list2

    cat `echo $cat_list1` > /data/docker_qiime2_share_container_E/FL2/all/seq_all/NG-27123_FL_amplicon_2_lib470794_1111_1_1.fastq
    cat `echo $cat_list2` > /data/docker_qiime2_share_container_E/FL2/all/seq_all/NG-27123_FL_amplicon_2_lib470794_1111_1_2.fastq

    ######### cat2
    cat_list1=""
    cat_list2=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/need_cat`
    do
    if grep -q "amplicon_3" <<< $x && grep -q "1.fastq" <<< $x
    then
    cat_list1+="/data/docker_qiime2_share_container_E/FL2/all/need_cat/$x "
    elif grep -q "amplicon_3" <<< $x && grep -q "2.fastq" <<< $x
    then
    cat_list2+="/data/docker_qiime2_share_container_E/FL2/all/need_cat/$x "
    fi
    done

    cat_list1=${cat_list1::-1}
    cat_list2=${cat_list2::-1}

    echo $cat_list1
    echo $cat_list2

    cat `echo $cat_list1` > /data/docker_qiime2_share_container_E/FL2/all/seq_all/NG-28397_FL_amplicon_3_16S1_NEW_lib570224_1111_3_1.fastq
    cat `echo $cat_list2` > /data/docker_qiime2_share_container_E/FL2/all/seq_all/NG-28397_FL_amplicon_3_16S1_NEW_lib570224_1111_3_2.fastq

    ######### compress and mv
    # press_list=""
    # for x in `ls /data/docker_qiime2_share_container_E/FL2/all/not_matching`
    # do
    # if grep -q "1111" <<< $x
    # then
    # press_list+="/data/docker_qiime2_share_container_E/FL2/all/not_matching/$x "
    # fi
    # done

    # press_list=${press_list::-1}

    # echo $press_list

    # parallel -j 4 gzip --best {} ::: `echo $press_list`


    parallel -j 4 gzip /data/docker_qiime2_share_container_E/FL2/all/seq_all/{} \
    ::: `ls /data/docker_qiime2_share_container_E/FL2/all/seq_all`

# copy left to seq_all ----
# fastqc on all ----
    go_list=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/seq_all/`
    do
    go_list+="/data/docker_qiime2_share_container_E/FL2/all/seq_all/$x "
    done

    go_list=${go_list::-1}

    echo $go_list

    time /data/docker_qiime2_share_container_E/bioinformatics-tools/fastqc/FastQC/fastqc \
    -o /data/docker_qiime2_share_container_E/FL2/all/Fastqc/op \
    -t 16 \
    `echo $go_list` > /data/docker_qiime2_share_container_E/FL2/all/fastaqc-log.txt 2>&1

    real    7m54.167s
    user    39m4.704s
    sys     4m44.149s

# multiqc ----
    # cd /data/docker_qiime2_share_container_E/bioinformatics-tools/multiqc
    # wget https://github.com/ewels/MultiQC/archive/refs/tags/v1.12.tar.gz

    pip install multiqc -i https://pypi.tuna.tsinghua.edu.cn/simple

    time multiqc /data/docker_qiime2_share_container_E/FL2/all/Fastqc/op/*.zip \
    -o /data/docker_qiime2_share_container_E/FL2/all/Multiqc -v > /data/docker_qiime2_share_container_E/FL2/all/multiqc-log.txt 2>&1


# fastp ----
    # apt-get install libisal2
    # /data/docker_qiime2_share_container_E/bioinformatics-tools/fastp/fastp-0.23.2/fastp --help > /data/docker_qiime2_share_container_E/FL2/all/fastp_help.sh 2>&1

    # rename ----
    # 1  NG-26931_FL_amplicon_1_lib460790_7284_1_1.fastq.gz
    # 2  NG-26931_FL_amplicon_1_lib460790_7284_1_2.fastq.gz
    # 3  NG-27123_FL_amplicon_2_lib470794_1111_1_1.fastq.gz
    # 4  NG-27123_FL_amplicon_2_lib470794_1111_1_2.fastq.gz

    # to 
    # 1  NG-26931_FL_amplicon_1_LM1_lib460790_7284_1_1.fastq.gz
    # 2  NG-26931_FL_amplicon_1_LM1_lib460790_7284_1_2.fastq.gz
    # 3  NG-27123_FL_amplicon_2_LM1_lib470794_1111_1_1.fastq.gz
    # 4  NG-27123_FL_amplicon_2_LM1_lib470794_1111_1_2.fastq.gz

    go_list=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/seq_all/`
    do
    go_list+="/data/docker_qiime2_share_container_E/FL2/all/seq_all/${x::-10} "
    done

    go_list=${go_list::-1}

    go_list=`echo $go_list | tr ' ' '\n' | sort | uniq | tr '\n' ' '`

    echo $go_list

    for x in `echo $go_list`
    do
    echo `echo $x`1.fastq.gz
    echo $x | cut -d "/" -f 7
    done

    for x in `echo $go_list`
    do
    time /data/docker_qiime2_share_container_E/bioinformatics-tools/fastp/fastp-0.23.2/fastp \
    -i `echo $x`1.fastq.gz \
    -I `echo $x`2.fastq.gz \
    -o /data/docker_qiime2_share_container_E/FL2/all/Fastp/`echo $x | cut -d "/" -f 7`1.fastq \
    -O /data/docker_qiime2_share_container_E/FL2/all/Fastp/`echo $x | cut -d "/" -f 7`2.fastq \
    --verbose --disable_adapter_trimming --dont_eval_duplication --disable_trim_poly_g --cut_front --cut_tail \
    --cut_right --cut_mean_quality 20 --length_required 145 -w 16 \
    -j /data/docker_qiime2_share_container_E/FL2/all/Fastp/`echo $x | cut -d "/" -f 7`.json \
    -h /data/docker_qiime2_share_container_E/FL2/all/Fastp/`echo $x | cut -d "/" -f 7`.html
    done >  /data/docker_qiime2_share_container_E/FL2/all/fastp-log.txt 2>&1

# barcodes generation using R from excel ----
# cut adapt ----

    # renaming amplicon 1 and amplicon 2 to match bracode file name
    # apt-get update
    # apt-get install rename

    # cd /data/docker_qiime2_share_container_E/FL2/all/Fastp

    # rename 's/amplicon_1/amplicon_1_LM1/g' *
    # rename 's/amplicon_2/amplicon_2_LM1/g' *

    # dir making
    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/barcodes/ | cut -d "." -f 1`
    do
    mkdir /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x
    done

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`
    do
    for y in `seq 1 4`
    do
    mkdir /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op$y
    done
    done

    # 1 for first round
    # "/data/docker_qiime2_share_container_E/FL2/all/Fastp/${x::-7} "

    go_list=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/Fastp/*.fastq`
    do
    go_list+="${x::-7} "
    done

    go_list=${go_list::-1}

    go_list=`echo $go_list | tr ' ' '\n' | sort | uniq | tr '\n' ' '`

    echo $go_list

    for x in $go_list
    do
    echo $x
    done

    #
    # go_list2=""

    # for x in `ls /data/docker_qiime2_share_container_E/FL2/all/Fastp/ | cut -d "." -f 1`
    # do
    # go_list2+="${x:9:-19} "
    # done
    # echo $go_list2

    # for x in `ls /data/docker_qiime2_share_container_E/FL2/all/barcodes/ | cut -d "." -f 1`
    # do
    # for y in `echo $go_list`
    # do
    # if grep -q $x <<< $y
    # then
    # echo $x
    # echo $y
    # fi
    # done
    # done

    # need lower -e for some too similar adapter
    # no need -O cause already set within fasta file
    # 151 - 30 - 1 = 120
    # 151 -28 -1 = 122
    


    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/barcodes/ | cut -d "." -f 1`
    do
    for y in `echo $go_list`
    do
    if grep -q $x <<< $y
    then
    time cutadapt -j 16 -e 0 --no-indels --max-n 2 --action trim --report full -m 120 \
    -g file:/data/docker_qiime2_share_container_E/FL2/all/barcodes/$x.fasta \
    -o /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op1/{name}.R1.fastq \
    -p /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op1/{name}.R2.fastq \
    `echo $y`1.fastq \
    `echo $y`2.fastq
    fi
    done
    done > /data/docker_qiime2_share_container_E/FL2/all/cut_log.txt 2>&1

    # 2 for second round
    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/barcodes/ | cut -d "." -f 1`
    do
    time cutadapt -j 16 -e 0 --max-n 2 --no-indels --action trim --report full -m 120 \
    -g file:/data/docker_qiime2_share_container_E/FL2/all/barcodes/$x.fasta \
    -o /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op2/{name}.R1.fastq \
    -p /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op2/{name}.R2.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op1/unknown.R2.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op1/unknown.R1.fastq
    done > /data/docker_qiime2_share_container_E/FL2/all/cut2_log.txt 2>&1

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/barcodes/ | cut -d "." -f 1`
    do
    # echo `ls -lh /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op2/unknown*`
    ls -lh /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op2/unknown*
    done


    # 3 for cating
    go_list3=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`
    do
    for y in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op1/`
    do
    if grep -q "unknown" <<< $y
    then 
    continue
    else
    go_list3+="`echo $x`/op1/`echo $y` "
    fi
    done
    done

    go_list3=${go_list3::-1}
    echo $go_list3

    for x in `echo $go_list3`
    do
    echo $x| sed -e "s/op1/op2/g"
    done

    # parallel -j 16 --header echo {folder1} {folder2} \
    # ::: folder1 `echo $go_list3` \
    # ::: folder2 `echo $go_list3 | sed -e "s/op1/op2/g"`

    parallel -j 16 echo {} '{= s/op1/op2/g =}' '{= s/op1/op3/g =}' ::: `echo $go_list3`

    # more than 16 cause io mission
    parallel -j 100 "cat /data/docker_qiime2_share_container_E/FL2/all/cutadapt/{} /data/docker_qiime2_share_container_E/FL2/all/cutadapt/'{= s/op1/op2/g =}' > /data/docker_qiime2_share_container_E/FL2/all/cutadapt/'{= s/op1/op3/g =}'" \
    ::: `echo $go_list3` > /data/docker_qiime2_share_container_E/FL2/all/cut3_log.txt 2>&1

    # 4 for tail trimming
    go_list4=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`
    do
    for y in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op3/`
    do
    if grep -q "R2" <<< $y
    then
    continue
    else
    go_list4+="`echo $x`/op3/`echo $y` "
    fi
    done
    done

    go_list4=${go_list4::-1}

    for x in `echo $go_list4`
    do
    echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"
    done

    # do not need $ sign in adapter since it require the adapter to present in full length

    # The minimum overlap length cannot be set for anchored adapters as these always need to occur at full length.

    time for x in `echo $go_list4`
    do
    if grep -q "_LM" <<< $x
    then
    echo "<<<" $x
    time cutadapt -j 16 -e 0 -O 5 --max-n 2 --no-indels --action trim --report full \
    -G GAGCATTATGCCAAGGTTTG \
    -o /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R1.fastq \
    -p /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R2.fastq \
    --untrimmed-output /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R1.untrimed.fastq \
    --untrimmed-paired-output /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R2.untrimed.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R1.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R2.fastq
    else
    echo "<<<" $x
    time cutadapt -j 16 -e 0 -O 5 --max-n 2 --no-indels --action trim --report full \
    -G GACTACCAGGGTATCTAATCCTGT \
    -o /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R1.fastq \
    -p /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R2.fastq \
    --untrimmed-output /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R1.untrimed.fastq \
    --untrimmed-paired-output /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1 | sed -e "s/op3/op4/g"`.R2.untrimed.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R1.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R2.fastq
    fi
    done > /data/docker_qiime2_share_container_E/FL2/all/cut4_log.txt 2>&1

# FLASH ----
    # FLASH help
    /data/docker_qiime2_share_container_E/bioinformatics-tools/FLASH/FLASH-1.2.11-Linux-x86_64/flash --help >  /data/docker_qiime2_share_container_E/FL2/all/FLASH_help.txt 2>&1

    # mkdir
    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt`
    do
    mkdir /data/docker_qiime2_share_container_E/FL2/all/FLASH/$x
    done

    # go list
    go_list5=""

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`
    do
    for y in `ls /data/docker_qiime2_share_container_E/FL2/all/cutadapt/$x/op4/`
    do
    if grep -q "untrimed" <<< $y || grep -q "R2" <<< $y
    then
    continue
    else
    go_list5+="`echo $x`/op4/`echo $y` "
    fi
    done
    done

    go_list5=${go_list5::-1}

    for x in `echo $go_list5`
    do
    # echo $x | cut -d "/" -f 3 | cut -d "." -f 1
    # echo $x | cut -d "/" -f 1
    echo $x
    done

    date ; time for x in `echo $go_list5`
    do
    echo "<<<" $x
    if grep -q "_LM" <<< $x
    then 
    time /data/docker_qiime2_share_container_E/bioinformatics-tools/FLASH/FLASH-1.2.11-Linux-x86_64/flash \
    -m 50 -M 100 -x 0 -t 16 -o `echo $x | cut -d "/" -f 3 | cut -d "." -f 1` \
    -d /data/docker_qiime2_share_container_E/FL2/all/FLASH/`echo $x | cut -d "/" -f 1`/ \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R1.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R2.fastq
    else
    time /data/docker_qiime2_share_container_E/bioinformatics-tools/FLASH/FLASH-1.2.11-Linux-x86_64/flash \
    -m 4 -M 48 -x 0 -t 16 -o `echo $x | cut -d "/" -f 3 | cut -d "." -f 1` \
    -d /data/docker_qiime2_share_container_E/FL2/all/FLASH/`echo $x | cut -d "/" -f 1`/ \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R1.fastq \
    /data/docker_qiime2_share_container_E/FL2/all/cutadapt/`echo $x | cut -d "." -f 1`.R2.fastq
    fi
    done > /data/docker_qiime2_share_container_E/FL2/all/FLASH_log.txt 2>&1
# fastqc multiqc after FLASH ----
    # fastqc
    # dir making
    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH`
    do
    mkdir /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH/$x
    done

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH`
    do
    echo "<<<"  $x
    a_go_list=""
    for y in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH/$x`
    do
    if grep -q extendedFrags <<< $y
    then
    a_go_list+="/data/docker_qiime2_share_container_E/FL2/all/FLASH/$x/$y "
    fi
    done
    a_go_list=${a_go_list::-1}
    for z in `echo $a_go_list`
    do
    echo $z
    done
    done > /data/docker_qiime2_share_container_E/FL2/all/fastaqc-post-FLASH-log.txt 2>&1

    # for x in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH`
    # do
    # echo "<<<"  $x
    # a_go_list=""
    # for y in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH/$x`
    # do
    # if grep -q "extendedFrags" <<< $y
    # then
    # a_go_list+="/data/docker_qiime2_share_container_E/FL2/all/FLASH/$x/$y "
    # fi
    # done
    # a_go_list=${a_go_list::-1}
    # /data/docker_qiime2_share_container_E/bioinformatics-tools/fastqc/FastQC/fastqc \
    # -o /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH/$x \
    # `echo $a_go_list` -t 16
    # done > /data/docker_qiime2_share_container_E/FL2/all/fastaqc-post-FLASH-log.txt 2>&1

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH`
    do
    echo "<<<"  $x
    a_go_list=""
    for y in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH/$x`
    do
    if grep -q "extendedFrags" <<< $y
    then
    a_go_list+="/data/docker_qiime2_share_container_E/FL2/all/FLASH/$x/$y "
    fi
    done
    a_go_list=${a_go_list::-1}
    time parallel -j 5 /data/docker_qiime2_share_container_E/bioinformatics-tools/fastqc/FastQC/fastqc \
    -o /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH/$x/ -t 2 {} \
    ::: `echo $a_go_list`
    done > /data/docker_qiime2_share_container_E/FL2/all/fastaqc-post-FLASH-log.txt 2>&1

    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH`
    # do
    # echo "<<<"  $x
    # # a_go_list=""
    # for y in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH/$x`
    # do
    # if grep -q "extendedFrags" <<< $y
    # then
    # time /data/docker_qiime2_share_container_E/bioinformatics-tools/fastqc/FastQC/fastqc \
    # -o /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH/$x/ -t 16 \
    # /data/docker_qiime2_share_container_E/FL2/all/FLASH/$x/$y
    # fi
    # done
    # done > /data/docker_qiime2_share_container_E/FL2/all/fastaqc-post-FLASH-log.txt 2>&1







    # multiqc
    for x in `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH`
    do
    mkdir /data/docker_qiime2_share_container_E/FL2/all/Multiqc/post_FLASH/$x
    done

    # for x in `ls /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH`
    # do
    # echo "<<<" $x
    # time multiqc /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH/$x \
    # -o /data/docker_qiime2_share_container_E/FL2/all/Multiqc/post_FLASH/$x/*.zip
    # done > /data/docker_qiime2_share_container_E/FL2/all/multiqc-post-FLASH-log.txt 2>&1

    time parallel -j 16 multiqc /data/docker_qiime2_share_container_E/FL2/all/Fastqc/After_FLASH/{}/*.zip \
    -o /data/docker_qiime2_share_container_E/FL2/all/Multiqc/post_FLASH/{} \
    ::: `ls /data/docker_qiime2_share_container_E/FL2/all/FLASH` > /data/docker_qiime2_share_container_E/FL2/all/multiqc-post-FLASH-log.txt 2>&1



# uchime ----
# cd /data/docker_qiime2_share_container_E/bioinformatics-tools/usearch
# wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz

# chmod +x /data/docker_qiime2_share_container_E/bioinformatics-tools/usearch/usearch11.0.667_i86linux32

# test
# /data/docker_qiime2_share_container_E/FL2/all/FLASH/FL_amplicon_1_LM1/Lm_amplicon_FL_F1.extendedFrags.fastq \
# a3_s1.fasta
# 19_ref.fasta
# LM_ref.fasta

/data/docker_qiime2_share_container_E/bioinformatics-tools/usearch/usearch11.0.667_i86linux32 -uchime2_ref \
/data/docker_qiime2_share_container_E/FL2/all/a3_s4.fasta \
-db /data/docker_qiime2_share_container_E/FL2/all/19_ref.fasta \
-uchimeout /data/docker_qiime2_share_container_E/FL2/all/uchime/test_out.txt \
-uchimealnout /data/docker_qiime2_share_container_E/FL2/all/uchime/test_aln.txt \
-notmatched /data/docker_qiime2_share_container_E/FL2/all/uchime/test_notmatched.fasta \
-strand plus -mode balanced -threads 16

/data/docker_qiime2_share_container_E/bioinformatics-tools/usearch/usearch11.0.667_i86linux32 -annot \
/data/docker_qiime2_share_container_E/FL2/all/a3_s4.fasta \
-db /data/docker_qiime2_share_container_E/FL2/all/19_ref.fasta \
-tabbedout /data/docker_qiime2_share_container_E/FL2/all/uchime/test_annote.txt -threads 16

Mugelia

#### LM

/data/docker_qiime2_share_container_E/bioinformatics-tools/usearch/usearch11.0.667_i86linux32 -uchime2_ref \
/data/docker_qiime2_share_container_E/FL2/all/a1_LM1_repseq.fasta \
-db /data/docker_qiime2_share_container_E/FL2/all/LM.udb \
-uchimeout /data/docker_qiime2_share_container_E/FL2/all/uchime/test_out.txt \
-uchimealnout /data/docker_qiime2_share_container_E/FL2/all/uchime/test_aln.txt \
-strand plus -mode balanced -threads 16

/data/docker_qiime2_share_container_E/bioinformatics-tools/usearch/usearch11.0.667_i86linux32 -annot \
/data/docker_qiime2_share_container_E/FL2/all/a1_LM1_repseq.fasta \
-db /data/docker_qiime2_share_container_E/FL2/all/LM_ref.fasta \
-tabbedout /data/docker_qiime2_share_container_E/FL2/all/uchime/test_annote.txt -threads 16

# qiime2 ----
deblur with positive filter

# muscle ----

# chmod +x /data/docker_qiime2_share_container_E/bioinformatics-tools/muscle/muscle5.1.linux_intel64

/data/docker_qiime2_share_container_E/bioinformatics-tools/muscle/muscle5.1.linux_intel64 \
-align /data/docker_qiime2_share_container_E/FL2/all/uchime/2_fragi.fasta \
-output /data/docker_qiime2_share_container_E/FL2/all/uchime/2_fragi_muscled.fasta -threads 16

/data/docker_qiime2_share_container_E/bioinformatics-tools/muscle/muscle5.1.linux_intel64 \
-align /data/docker_qiime2_share_container_E/FL2/all/qiime2/19_q2_seqs.fasta \
-output /data/docker_qiime2_share_container_E/FL2/all/qiime2/19_q2_seqs_muscled.fasta -threads 16