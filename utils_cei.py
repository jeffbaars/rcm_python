#!/usr/bin/python
import os, sys, glob, re, math
import numpy as np
from utils_date import *
from utils_ghcnd_obs import *

def get_daily_perc(mod_all, mod_dts_all):

    bdt = '1970010100'
    tx10 = {}
    tn10 = {}
    tx90 = {}
    tn90 = {}
    for d in range(0,364):
        fdate = ct10(bdt, (24 * d))
    $basemmddhh = substr($fdate,4,6);
    @txunsorted=();
    @tnunsorted=();
    $t=0;
    for ($y=$basefirst;$y<=$baselast;$y++) {
	$basedate = "$y"."$basemmddhh";
	for ($l=-2;$l<=2;$l++) {
	    $date = &ct10($basedate,(24*$l));
	    $txunsorted[$t] = $tx{$date};
	    $tnunsorted[$t] = $tn{$date};
	    $t++;
	}
    }
    @txsorted=();
    @txsorted = sort { $a <=> $b } @txunsorted;
    $num = $#txsorted;
    $n10 = ($num+1)*0.1;
    $n90 = ($num+1)*0.9;
    $tx10{$fdate} = $txsorted[$n10];
    $tx90{$fdate} = $txsorted[$n90];
    @tnsorted=();
    @tnsorted = sort { $a <=> $b } @tnunsorted;
    $num = $#tnsorted;
    $n10 = ($num+1)*0.1;
    $n90 = ($num+1)*0.9;
    $tn10{$fdate} = $tnsorted[$n10];
    $tn90{$fdate} = $tnsorted[$n90];
    printf TX10 ("%10.10d,%4.1f\n",$fdate,$tx10{$fdate});
    printf TX90 ("%10.10d,%4.1f\n",$fdate,$tx90{$fdate});
    printf TN10 ("%10.10d,%4.1f\n",$fdate,$tn10{$fdate});
    printf TN90 ("%10.10d,%4.1f\n",$fdate,$tn90{$fdate});
}
close TX10;
