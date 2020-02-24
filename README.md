- [Clean GWAS](#clean-gwas)


The scripts depends on the [Anaconda](https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-2019.10-Linux-x86_64.sh) development environment:

Anaconda3-2019.10-Linux-x86_64.sh.

# Clean GWAS

Summary:
    This a script that conducts basic and advanced quality control for a gwas summary data.

Dependencies: Anaconda

Usage:
    1) Configure basic variables in the head of this script.
    2) using following command: python3 cleangwas.py -h for help.

Details:
    1) SNPs meeting following conditions will be removed:
       beta<0,
       se<0,
       P<0,
       info<0,
       freq<0 or freq>1,
       or < 0,
    2) SNPs meeting will be removed if you set additional parameters:
       p>threshold,
       info<threshold,
       freq<maf,
       N < threshold,
       duplicated SNPs,
       palindrome SNPs,
       indels,
       SNPs with missing values,
    3) Alleles of SNPs will be uppercased.
    4) Uppercase RS number
    5) Eligible SNPs will be written into file you specified.
