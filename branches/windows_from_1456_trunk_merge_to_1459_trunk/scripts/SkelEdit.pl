#!/bin/sh
# The following line is for platform independence.
PERL=`type perl | cut -f3 -d\ `
#
# see whether perl was found
RESP=`echo $PERL | grep perl` 
if [ $? -eq 1 ] ; then
#
#   perl wasn't found
  echo "Couldn't find perl in your path. Please add perl to your path"
  echo "and try again."
  exit 1;
fi
$PERL -x $0 "$@" ; exit $?
#!perl
#------------------------------------------------------------------------------
#
# Author Ed.Rice & J.Gerschwitz
# COPYRIGHT Ed Rice & J.Gerschwitz  2004- All Rights Reserved.
# This software is the property of Ed.Rice & J.Gerschwitz.
# No part of this code may be copied in any form or by any means without the
# expressed written consent of Ed.Rice & J.Gerschwitz.
#
# Perl script to perform skeleton driven file creation
# Based on a skeleton configuration file this script will 
# edit a skeleton file replacing "placekeeper" text with user specified
# text generating updated files
#
# This script expects the name of the skeleton file as its first
# argument - additional arguments depend on the skeleton
#-----------------------------------------------------------------------------
use strict 'vars';

#
# some global variables because I don't know how to make perl pass 
# array variables around to subroutines
#
my $nReplacements = 0;
my @originalText;
my @replacementText;

{
   # There must be at least one command-line argument was provided - the
   # name of the skeleton file
   # Additional arguments may be required by the Skeleton file
   if ($#ARGV < 0) {
      print("Program Usage: SkelEdit <SkeletonFile> <possible skeleton ",
            "file arguments>\n\n",
	    "  where SkeletonFile - is the name of the skeleton file to ",
	    "process\n\n");

      exit(-1);
   }

   my $skeletonFile = $ARGV[0];

   #
   # Open the skeleton file and process it
   #
   if (!open(SKELH, "$skeletonFile")) {
       print("SkelEdit: Unable to open skeleton file ($skeletonFile).\n");
       exit(-1);
   }

   my $currentBlock = "none";
   my $nArgs = 0;
   my $maxArgs = $#ARGV; 
   my $nLineUsage = 0;
   my @usageText;
   my $nFiles;
   my @filesToCheck;
   while(<SKELH>) {
       if($currentBlock eq "none") {
	   if(/^\s*BeginMessage/) {
	       $currentBlock="message";
	   }
	   elsif(/^\s*BeginUsage/) {
	       $currentBlock="usage";
               $nLineUsage = 0;
	       if(/^\s*BeginUsage\s+(\d*)/) {
		   $nArgs = $1;
	       }
	   }
	   elsif(/^\s*BeginReplacements/) {
	       $currentBlock="replacements";
	       $nReplacements = 0;
	   }
	   elsif(/^\s*BeginFileCheck/) {
	       $currentBlock="filecheck";
	       $nFiles = 0;
	   }
	   elsif(/^\s*BeginFile\s*/) {
	       $currentBlock="file";
	       if(/^\s*BeginFile\s+(\S*)/) {
		   my $file = doReplacements($1);
		   if(-f $file) {
		       print "Error in File block - output file exists\n";
		       print "Attempting to output file $file\n";
		       exit(-1);
		   }
		   if (!open(OUTPUTH, ">$file")) {
		       print("Unable to open output file ($file)\n");
		       exit(-1);
		   }
		   print "Processing $file";
	       }
	       else {
		   print "Error in File block - no output file specified\n";
		   print "$_\n";
		   exit(-1);
	       }
	   }
       }
       elsif($currentBlock eq "message") {
	   if(/^\s*EndMessage/) {
	       $currentBlock="none";
	   }
	   else {
	       print doReplacements($_);
	   }
       }
      elsif($currentBlock eq "usage") {
	   if(/^\s*EndUsage/) {
               if($nArgs>$maxArgs) {
		   print "Error - too few arguments specified\n";
		   for(my $i = 0; $i<$nLineUsage; $i++) {
		       print $usageText[$i];
		   }
                   exit(-1);
	       }
	       $currentBlock="none";
	   }
	   else {
	       $usageText[$nLineUsage] = $_;
               $nLineUsage++;
	   }
       }
      elsif($currentBlock eq "replacements") {
	  if(/^\s*EndReplacements/) {
	      $currentBlock="none";
	  }
	  else {
	      if(/^\s*(<[^>]*>)\s*(\S.*)/) {
		  my $orig = $1;
		  my $repl = $2;
		  my $newText;
                  if($repl =~ /arg:(\d*)/) {
		      my $argNum = $1;
                      if($argNum >= $nArgs) {
			  print("Problem with replacement block\n",
				$_,
				"Argument number exceeds maximum specified",
				"in usage block\n",
				"Required argument index : $argNum\n",
				"Number of arguments in usage : $nArgs\n",
				"Remember argument indices are zero-based\n");
			  exit(-1);
		      }
		      $newText = $ARGV[$argNum+1];
		  }
		  elsif($repl =~ /date:(\S*)/) {
		      my $dateFormat = $1;
		      if($dateFormat eq "yyyymmdd") {
			  my ($sec,$min,$hour,$mday,$mon,
			      $year,$wday,$yday,$isdst) = localtime(time);
			  $mon++;
			  if($mon<10) {
			      $mon = "0".$mon;
			  }
			  $year+=1900;
			  if($mday<10) {
			      $mday = "0".$mday;
			  }
			  $newText = $year.$mon.$mday;
		      }
		      else {
			  print("Problem with replacement block\n",
				$_,
				"date format ($dateFormat) unknown\n");
			  exit(-1);
		      }
		  }
		  else {
		      print("Problem with replacement block\n",
			    $_,
			    "unknown replacement format ($repl)\n");
		      exit(-1);
		  }
		  $originalText[$nReplacements] = $orig;
		  $replacementText[$nReplacements] = $newText;
                  $nReplacements++;
	      }
	  }
      }
      elsif($currentBlock eq "filecheck") {
	   if(/^\s*EndFileCheck/) {
	       my $nFailed=0;
	       my @failedFiles;
	       for(my $k=0; $k<$nFiles; $k++) {
		   if(-f $filesToCheck[$k]) {
		       $failedFiles[$nFailed]=$filesToCheck[$k];
		       $nFailed++;
		   }
	       }
	       if($nFailed>0) {
		   print("Error processing FileCheck block\n",
			 "The following files already exist - \n");
		   for(my $i=0; $i<$nFailed; $i++) {
		       print "$failedFiles[$i]\n";
		   }
		   exit(-1);
	       }
	       $currentBlock="none";
	   }
	   else {
	       # chomp not working cleanly, remove \n and \r
	       s/\n//;
	       s/\r//;
	       $filesToCheck[$nFiles] = doReplacements($_);
               $nFiles++;
	   }
       }
      elsif($currentBlock eq "file") {
	   if(/^\s*EndFile\s+$/) {
	       close(OUTPUTH);
	       $currentBlock="none";
	       print " .....ok\n";
	   }
	   else {
	       print OUTPUTH doReplacements($_);
	   }
       }

   }
   if($currentBlock ne "none") {
       print "Error - end of skeleton file <$skeletonFile> reached\n";
       print "But still in $currentBlock block\n";
       print "Missing End keyword ???\n";
       exit(-1);
   }
   close(SKELH);
}


sub doReplacements {
  my $text=$_[0];
  for(my $j=0; $j<$nReplacements; $j++) {
      $text =~ s/$originalText[$j]/$replacementText[$j]/g;
  }
  return $text;
}
