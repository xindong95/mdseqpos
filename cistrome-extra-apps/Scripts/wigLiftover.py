#!/usr/bin/env python

import sys
import os
import subprocess

def main():
    inputWig = sys.argv[1]
    overChainProgram = sys.argv[2]
    outputWig = sys.argv[3]

    #converting the input wig file to a bed file
    subprocess.call(['wigToBed.py', inputWig, 'tmpBed1.bed'])

    #performing liftover on bed file
    #supposing liftOver has been installed to the system PATH.
    subprocess.call(['liftOver','tmpBed1.bed', overChainProgram, 'tmpBed2.bed', 'UnMapped'])

    #converting the new bed file to wig file
    subprocess.call(['bedToWig2.py', inputWig, 'tmpBed2.bed', outputWig])

    #removing the temporary bed files
#    os.remove('./tmpBed1.bed')
 #   os.remove('./tmpBed2.bed')
             


if __name__ =='__main__': main()
