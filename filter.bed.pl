#!/usr/bin/perl -w
#we remove all lines with fields ~=5 p,value (field 5) >1000
while(<>){
	@a=split;
	print if (4==$#a && $a[4]<=1000);
}
