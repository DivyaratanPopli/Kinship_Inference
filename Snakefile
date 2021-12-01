import pysam
import pybedtools
import pandas as pd
import numpy as np

from input_preparation_func import *
from hbd_hmm_func_new import *

interval=1e7
#split_bams = "./bamfiles/stephane"

#HN_diff0=0.36 chagyrskaya
#HN_diff1=0.55 chagyrskaya
#HN_diff0=0.35 vindija
#HN_diff1=0.475 vindija

#HN_diff=[HN_diff0,HN_diff1]
#target_ind=config["target_ind"]
#contam_ind=config["contam_ind"]

is_contam=int(config["is_contam"])

if is_contam==1:
    tc_diff=int(config["tc_diff"])
    if tc_diff==0:
        tar_ind1=config["tar_ind"]
        contam_ind1=config["contam_ind"]

        contam_diff_vcf=config["contam_diff.vcf.gz"]
        contam_diff_ind=contam_diff_vcf + '.tbi'
    elif tc_diff==1:
        tar_ind1=0
        contam_ind1=0
        contam_diff_vcf="contam_diff.vcf.gz"
        with open(contam_diff_vcf,'w') as f:
            print("NA")
        contam_diff_ind=contam_diff_vcf+'.tbi'
        with open(contam_diff_ind,'w') as f:
            print("NA")
        with open("test_fil0.vcf",'w') as f:
            print("NA")
        with open("test_fil1.vcf",'w') as f:
            print("NA")
        with open("contam_diff_fil0.txt",'w') as f:
            print("NA")
        with open("contam_diff_fil1.txt",'w') as f:
            print("NA")


elif is_contam==0:
    tar_ind1='none'
    contam_ind1='none'
    phased1='none'
    contam_diff_vcf="contam_diff.vcf.gz"
    with open(contam_diff_vcf,'w') as f:
        print("NA")
    contam_diff_ind=contam_diff_vcf+'.tbi'
    with open(contam_diff_ind,'w') as f:
        print("NA")
    contam_estF="contam_est"
    with open(contam_estF,'w') as f:
        print("NA")



with open("targets.txt") as f:
    libraries = [line.strip() for line in f]

totalch=list(range(1,23))

listf=[]
for i, l1 in enumerate(libraries):
    for l2 in libraries[(i+1):]:
        s = '%s_%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
        listf.append(s)

with open("identical.txt") as f:
    listidentical = [line.strip() for line in f]


with open("unrelated.txt") as f:
    lib_un = [line.strip() for line in f]


listun=[]
for i, l1 in enumerate(lib_un):
    for l2 in lib_un[(i+1):]:
        s = '%s_%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
        listun.append(s)


Mfilist=[0,1]


#run name_split.sh to first create right names in bamfiles/stephane

rule index_bam:
    input: "{name}.bam"
    output: "{name}.bam.bai"
    shell: "samtools index {input}"


rule sort:
    input: "bamfiles/stephane/{name}.bam"
    output: "bamfiles/sorted/{name}.bam"
    shell: "samtools sort {input} | bam-rmdup -r -o {output}"

rule splitbams:
    input:
        bam="bamfiles/sorted/{bam}.bam",
        bai="bamfiles/sorted/{bam}.bam.bai"
    output:
        "bamfiles/split_bams/chrm{ch}/{bam}.bam"
    shell: "samtools view -b {input.bam} {wildcards.ch} > {output}"



rule split_bed:
    input:
        bed="allsites.bed"
    output:
        split="chrm{ch}_allsites.bed"
    shell: "awk '{{if ($1 == {wildcards.ch}) print $0}}' {input.bed} > {output.split}"


rule make_hapProbs:
    input:
        bam="bamfiles/split_bams/chrm{ch}/{lib}.bam",
        bai="bamfiles/split_bams/chrm{ch}/{lib}.bam.bai",
        bed="chrm{ch}_allsites.bed"
    output:
        ps="hapProbs/chrm{ch}/{lib}.csv",
        ds="identicaldiffs/chrm{ch}/{lib}.csv"
    run:

        makeHapProbs(inbam=input.bam, inbai=input.bai, inbed=input.bed, outprob=output.ps, outdiff_id=output.ds)


rule hapProbs_all:
    input:
        pslist=expand("hapProbs/chrm{{ch}}/{lib}.csv", lib=libraries),
    output:
        probs="hapProbs/all_libs/chrm{ch}/probs.csv",
        psf="hapProbs/all_libs/chrm{ch}/ps.csv",
        chrmf="hapProbs/all_libs/chrm{ch}/chrm.csv",
        posf="hapProbs/all_libs/chrm{ch}/pos.csv"
    run:
        hapProbsAll(haplist=input.pslist, outf=output.probs, hap='noid', posf=output.posf, chrmf=output.chrmf, psf=output.psf)


rule idhapProbs_all:
    input:
        dslist=expand("identicaldiffs/chrm{{ch}}/{lib}.csv", lib=libraries)
    output:
        id_diffs="identicaldiffs/all_libs/chrm{ch}/idDiff.csv"
    run:
        hapProbsAll(haplist=input.dslist, outf=output.id_diffs, hap='id', posf='none', chrmf='none', psf='none')

rule basic:
    input:
        h=expand("identicaldiffs/all_libs/chrm{ch}/idDiff.csv", ch=2),
        i=expand("hapProbs/all_libs/chrm{ch}/probs.csv", ch=2)


rule Mfilter_tped:      #makes monomorphic filter and READ input file
    input:                      #each chromosome
        infiles="hapProbs/all_libs/chrm{ch}/ps.csv",
        posf="hapProbs/all_libs/chrm{ch}/pos.csv"
    output:
        fil="pshap/chrm{ch}/pos_filtered{Mfil}",
        read="readfiles/chrm{ch}/allsamples_fil{Mfil}.tped"
    run:
        filterPseudohap(filist=input.infiles, posfile=input.posf, outfilter=output.fil, readfile=output.read,
            fil=wildcards.Mfil, fil_version='hap')


rule Mfilter_prob:      #makes monomorphic filter and READ input file
    input:                      #each chromosome
        infiles="hapProbs/all_libs/chrm{ch}/probs.csv",
        posf="hapProbs/all_libs/chrm{ch}/pos.csv"
    output:
        fil="hapProbs/chrm{ch}/pos_filtered{Mfil}",
    run:
        filterPseudohap(filist=input.infiles, posfile=input.posf, outfilter=output.fil, readfile='none',
            fil=wildcards.Mfil, fil_version='prob')


rule READ_input:
    input:
        tped=expand("readfiles/chrm{ch}/allsamples_fil{{Mfil}}.tped", ch=totalch),
        posf=expand("pshap/chrm{ch}/pos_filtered{{Mfil}}", ch=totalch)
    output:
        read="readfiles_fil{Mfil}/allsamples.tped",
        tfam="readfiles_fil{Mfil}/allsamples.tfam"

    run:
        READInput(filist=input.tped,read=output.read, posfile=input.posf,
                tfam=output.tfam, libraries=libraries)


rule find_diff:
    input:
        indf="hapProbs/all_libs/chrm{ch}/probs.csv",
        infilter="hapProbs/chrm{ch}/pos_filtered{Mfil}",
        posallf="hapProbs/all_libs/chrm{ch}/pos.csv",
    output:
        diff="diffs_fil{Mfil}/chrm{ch}/alldiffs.csv",

    run:
        findDiff(indfile=input.indf,infilter=input.infilter, posall=input.posallf,
            difffile=output.diff,fil=wildcards.Mfil)



rule identical_diff:
    input:
        ind="identicaldiffs/all_libs/chrm{ch}/idDiff.csv",
        infilter="hapProbs/chrm{ch}/pos_filtered{Mfil}",
        posallf="hapProbs/all_libs/chrm{ch}/pos.csv",

    output:
        diff="identicaldiffs_fil{Mfil}/chrm{ch}/alldiffs.csv",

    run:
        identicalDiff(indfile=input.ind,posall=input.posallf,
                infilter=input.infilter, difffile=output.diff,fil=wildcards.Mfil)



rule get_win:
    input:
        pwfile="diffs_fil{Mfil}/chrm{ch}/alldiffs.csv",
        posf="hapProbs/chrm{ch}/pos_filtered{Mfil}",

    output:
        df="windows_fil{Mfil}/chrm{ch}/allwind.csv",
        tf="windows_fil{Mfil}/chrm{ch}/allwint.csv",

    run:
        getWin(diff=input.pwfile,
            dfile=output.df, tfile=output.tf, posf=input.posf,
            interval=interval)

rule identical_win:
    input:
        pwfile="identicaldiffs_fil{Mfil}/chrm{ch}/alldiffs.csv",
        posf="hapProbs/chrm{ch}/pos_filtered{Mfil}",

    output:
        df="identicalwindows_fil{Mfil}/chrm{ch}/idwind.csv",
        tf="identicalwindows_fil{Mfil}/chrm{ch}/idwint.csv",

    run:
        getWin(diff=input.pwfile,
            dfile=output.df, tfile=output.tf, posf=input.posf,
            interval=interval)


rule merge_chrm:
    input:
        dinfiles= expand("windows_fil{{Mfil}}/chrm{ch}/allwind.csv",ch=totalch),
        tinfiles= expand("windows_fil{{Mfil}}/chrm{ch}/allwint.csv",ch=totalch)

    output:
        dfile= "mergedwin_fil{Mfil}/merged_wind.csv",
        tfile= "mergedwin_fil{Mfil}/merged_wint.csv",
        chfile= "mergedwin_fil{Mfil}/merged_chrm.csv",
    run:

        mergeChrm(dinfiles=input.dinfiles,
                    tinfiles=input.tinfiles,
                    dfile=output.dfile,
                    tfile=output.tfile,
                    chfile=output.chfile)


rule merge_identicalchrm:
    input:
        dinfiles= expand("identicalwindows_fil{{Mfil}}/chrm{ch}/idwind.csv",ch=totalch),
        tinfiles= expand("identicalwindows_fil{{Mfil}}/chrm{ch}/idwint.csv",ch=totalch)
    output:
        dfile= "identicalmergedwin_fil{Mfil}/id_wind.csv",
        tfile= "identicalmergedwin_fil{Mfil}/id_wint.csv",
        chfile= "identicalmergedwin_fil{Mfil}/id_chrm.csv",

    run:

        mergeChrm(dinfiles=input.dinfiles,
                            tinfiles=input.tinfiles,
                            dfile=output.dfile,
                            tfile=output.tfile,
                            chfile=output.chfile)



rule contam_file:
    input:
        cfile="contam_est"
    output:
        pairf="contam_est_pairwise",
        idfile="contam_est_2"
    run:
        contamFile(infile=input.cfile, outfile=output.pairf, targets=libraries, idfile=output.idfile, iscnt=is_contam)


rule merge_pos:
    input:
        poslist=expand("hapProbs/chrm{ch}/pos_filtered{{Mfil}}", ch=totalch)
    output:
        filbed="filtered_bed{Mfil}.bed"
    run:
        mergePos(poslist=input.poslist, filbed=output.filbed)



rule nh_inputFile:
    input:
        genf=contam_diff_vcf,
        index=contam_diff_ind,
        bed="filtered_bed{Mfil}.bed"
    output:
        test="test_fil{Mfil}.vcf",
        diff="contam_diff_fil{Mfil}.txt",
    params:
        cnt=is_contam,
        tc=tc_diff,
        tar_ind=tar_ind1,
        contam_ind=contam_ind1
    shell:
        """
        (
        if [[ {params.cnt} -eq 1 ]] && [[ {params.tc} -eq 0 ]]
        then
            bcftools view {input.genf} -R {input.bed} -s {params.tar_ind},{params.contam_ind} -e 'GT="mis"'>{output.test}
            bcftools query -f '%CHROM\t%POS[\t%GT]\n' {output.test}>{output.diff}
        elif [[ {params.cnt} -eq 0 ]] || [[ {params.tc} -eq 1 ]]
        then
            echo "NA">{output.test}
            echo "NA">{output.diff}
        fi
        )
        """


def second_input(wildcards):
    xxstr="contam_diff_fil%s.txt" %(str(wildcards.Mfil))
    if os.path.exists(xxstr):
        return ""
    else:
        return xxstr

rule nhfile:
    input:
        alldiff=second_input,
    output:
        avgdiff="nhfile_fil{Mfil}.txt"
    params:
        is_contam=is_contam,
    run:
        nhFile(alldiff=input.alldiff, avgdiff=output.avgdiff, is_contam=params.is_contam)





rule contam_all:
    input:
        dfile="mergedwin_fil{Mfil}/merged_wind.csv",
        tfile= "mergedwin_fil{Mfil}/merged_wint.csv",
        contam_est="contam_est_pairwise",
        nh_file="nhfile_fil{Mfil}.txt"


    output:
        difffile="mergedwin_contam_fil{Mfil}/contam_diff.csv",
        totalfile="mergedwin_contam_fil{Mfil}/contam_total.csv"
    run:

        contamAll(dfile=input.dfile, tfile=input.tfile, cfile=input.contam_est,
        Dnhfile=input.nh_file, difffile=output.difffile,
        totalfile=output.totalfile, iscnt=is_contam)


rule contam_id:
    input:
        dfile="identicalmergedwin_fil{Mfil}/id_wind.csv",
        tfile="identicalmergedwin_fil{Mfil}/id_wint.csv",
        contam_est="contam_est_2",
        nh_file="nhfile_fil{Mfil}.txt"

    output:
        difffile="identicalmergedwin_contam_fil{Mfil}/id_diff.csv",
        totalfile="identicalmergedwin_contam_fil{Mfil}/id_total.csv",
    run:
        contamAll(dfile=input.dfile, tfile=input.tfile, cfile=input.contam_est,
        Dnhfile=input.nh_file[0], difffile=output.difffile,
        totalfile=output.totalfile, iscnt=is_contam)


rule get_highDiv:
    input:
        difffile = "mergedwin_contam_fil{Mfil}/contam_diff.csv",
        totalfile = "mergedwin_contam_fil{Mfil}/contam_total.csv",
    output:
        highdiv="windows_with_diffProp_higher_than3sd{Mfil}.txt"
    run:
        getHighDiv(diff=input.difffile, total=input.totalfile, fout=output.highdiv)



rule rem_highDiv:
    input:
        diff="mergedwin_contam_fil{Mfil}/contam_diff.csv",
        total="mergedwin_contam_fil{Mfil}/contam_total.csv",
        winfile="windows_with_diffProp_higher_than3sd{Mfil}.txt"
    output:
        outd="mergedwin_remHighdiv_fil{Mfil}/pw_diff.csv",
        outt="mergedwin_remHighdiv_fil{Mfil}/pw_total.csv",
    run:
        remHighDiv(diffname=input.diff, totalname=input.total, winfile=input.winfile, outd=output.outd, outt=output.outt)




rule rem_highDiv_id:
    input:
        diff="identicalmergedwin_fil{Mfil}/id_wind.csv",
        total="identicalmergedwin_fil{Mfil}/id_wint.csv",
        winfile="windows_with_diffProp_higher_than3sd{Mfil}.txt"
    output:
        outd="identicalmergedwin_remHighdiv_fil{Mfil}/id_diff.csv",
        outt="identicalmergedwin_remHighdiv_fil{Mfil}/id_total.csv",
    run:
        remHighDiv(diffname=input.diff, totalname=input.total, winfile=input.winfile, outd=output.outd, outt=output.outt)




rule get_p1:
    input:
        difffile="mergedwin_remHighdiv_fil{Mfil}/pw_diff.csv",
        totalfile="mergedwin_remHighdiv_fil{Mfil}/pw_total.csv",
        id_difffile="identicalmergedwin_remHighdiv_fil{Mfil}/id_diff.csv",
        id_totalfile="identicalmergedwin_remHighdiv_fil{Mfil}/id_total.csv",

    output:
        pfile="hmm_parameters_fil{Mfil}/p1_file",
        p_all="hmm_parameters_fil{Mfil}/p_all",
        phalf="hmm_parameters_fil{Mfil}/phalf_file",
        phalf_all="hmm_parameters_fil{Mfil}/phalf_all",
        newpairs="goodpairs_fil{Mfil}.txt",
        newid="identical_fil{Mfil}.txt",
        over="overlap_fil{Mfil}",
        id_over="overlap_fil{Mfil}_id"

    run:
        getP(dfile=input.difffile, tfile=input.totalfile, targets=listf,
        pfile=output.pfile,goodpairs=output.newpairs,
        allp=output.p_all, overf=output.over)

        getP(dfile=input.id_difffile, tfile=input.id_totalfile, targets=libraries,
        pfile=output.phalf,goodpairs=output.newid,
        allp=output.phalf_all, overf=output.id_over)



rule hmm:
    input:
        p1file="hmm_parameters_fil{Mfil}/p1_file",
        dfile="identicalmergedwin_remHighdiv_fil{Mfil}/id_diff.csv",
        tfile="identicalmergedwin_remHighdiv_fil{Mfil}/id_total.csv",
        chfile="identicalmergedwin_fil{Mfil}/id_chrm.csv",

    output:
        resfile="hbd_win_fil{Mfil}/pw_{listind}.csv",
        likfile="hbd_win_fil{Mfil}/lik/pw_{listind}.csv",


    run:
        l=wildcards.listind

        hmm(dfile=input.dfile, tfile=input.tfile, chrmf=input.chfile, p1file=input.p1file, ind=l,
                        resfile=output.resfile, likfile=output.likfile,
                        upA='initial', targets=libraries)


rule all:
    input:
        over=expand("overlap_fil{Mfil}", Mfil=Mfilist),
        iover=expand("overlap_fil{Mfil}", Mfil=Mfilist),
        pfile=expand("hmm_parameters_fil{Mfil}/p1_file", Mfil=Mfilist),
        phalf=expand("hmm_parameters_fil{Mfil}/phalf_file", Mfil=Mfilist),
        tped=expand("readfiles_fil{Mfil}/allsamples.tped", Mfil=Mfilist),
        hbd=expand("hbd_win_fil{Mfil}/pw_{listid}.csv", Mfil=Mfilist, listid=libraries)
