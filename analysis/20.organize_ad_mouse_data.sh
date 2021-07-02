#!/bin/bash

# pull summarized tables together - leave images alone
datadir="/datasets/merfish_vizgen_ADmouse_XiangminXu/analyzed_data/"
files=$(find $datadir -maxdepth 3 -name *.csv)

for file in $files; do
    echo $file
    newfile=${file/$datadir/}
    # newfile=${newfile/"/"/"_"}
    newfile=$(echo ${newfile} | tr "/" "_")
    newfile="../data/raw_admouse/"$newfile
    echo $newfile
    ln -s $file $newfile
done

# simplify sample name
sample1="VS8_MsBrain_Xulab-2-5_Abstain_VS8_V3_BW_05-13-2021"
simp1="Xulab_2_5"
sample2="VS8_MsBrain_Xulab-2-6_Abstain_VS8_V3_BW_05-14-2021"
simp2="Xulab_2_6"
sample3="VS8_VizMsBrain_VS8_V3_BW_05-09-2021"
simp3="Vizgen"

for file in `ls ../data/raw_admouse/*.csv`; do
    newfile=${file/$sample1/$simp1}
    newfile=${newfile/$sample2/$simp2}
    newfile=${newfile/$sample3/$simp3}
    mv $file $newfile
done

