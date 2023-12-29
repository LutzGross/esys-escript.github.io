#!/usr/bin/perl
#

# environment variables supported:
# EXACT_DRIVER - mpirun command.  Defaults to mpirun -np
#   For solaris, it should be
#   "/Net/local/mpi/build/solaris/ch_p4/bin/mpirun -np"
# ZOLTAN_ARCH - architechture of the machine.  Choices are generally:
#   generic, solaris, sun, linux, qed.  Defaults to "sun"
# ZOLTAN_ZDRIVE - location of zdrive, defaults to
#   "../../src/Obj_$ZOLTAN_ARCH/zdrive".  Usually not set if ZOLTAN_ARCH
#   has already been set, except for special purposes.

use Sys::Hostname;
my $host = hostname();

use Data::Dumper;

# defaults, can be overridden in experiment.xml
my ($mpi,$zarch,$zdrive);	# filled in below based on env
my $zinp_def = "zdrive.inp";	# static, used below for checks
my $zinp = $zinp_def;		# zdrive.inp used
my $np=0;

# values to track and print by parse() routine
my (%bal, %cutl, %cutn);
my $iters;
my $ckok = '';	# output string for pass/fail of checksums
my $erreason = "";  # if set, report in exit_comment


# use specified driver if environment var is set.  Set up defaults here
if ($ENV{EXACT_DRIVER}) {
  $mpi = $ENV{EXACT_DRIVER};
} else {
  $mpi = "mpirun -np";
}

if ($ENV{ZOLTAN_ARCH}) {
  $zarch = $ENV{ZOLTAN_ARCH};
} else {
  $zarch = "sun";
}

if ($ENV{ZOLTAN_ZDRIVE}) {
  $zdrive = $ENV{ZOLTAN_ZDRIVE};
} else {
  $zdrive = "../../src/Obj_$zarch/zdrive";
}

# contents of $infile, after parseinfile()
my %args;

# checksums of all the answers files
my %cksum;

# list of output filenames, to be checksummed
my @checkfiles=();

# mesh/model name
my $name;

my $delimiter='"""';
my $cmdline;	# constructed later

my $outfile = $ARGV[1];
my $logfile = $ARGV[2];


sub parseinfile {
  # Process infile (variables, levels, etc)
  # argument: opened file descriptor

  my $infile = @_[0];

  while (<$infile>) {
    my ($key,$a) = split /\t/, $_, 2;
    chomp $a;
    if ($a eq $delimiter) {
      $a = '';
      while (<$infile>) {
	if (m/$delimiter/) {		# single line with delimiter
	  chomp $a;
	  break;
	} else {
	  $a .= $_;
	}
      }
    }
    $args{$key} = $a;
  }
# loop through factor levels for vars
  for (my $i=1;; $i++) {
    my $x= "_factor_${i}_name";
    if (defined $args{$x}) {
      my $y= "_factor_${i}_value";
      if ($args{$x} eq "numproc") {		# possible factors
	$np = $args{$y};
      } elsif ($args{$x} eq "zdrive") {
	$zdrive = $args{$y};
      } elsif ($args{$x} eq "zinput") {
	$zinp = $args{$y};
      }
    } else {
      last;
    }
  }
# set to what was found in tests filename for unset variables and overrides
# (you can have a factor numproc->np=4, and override in file with np=6)
  $zdrive = $args{zdrive} if (defined $args{zdrive});
  $zinp = $args{zinp} if (defined $args{zinp});
  $np = $args{np} if (defined $args{np});

}

sub parseoutfile {
  # parse zoltan output file
  # arguments: INFD LOGFD
  # INFD is a pipe to running command
  # LOGFD is the logfile, to dump output into

# example output
#| 0/ 4|: ( 0, 0)/[ 2, 2] ->STATS Runs 10  bal  CURRENT 1.099466  MAX 1.099466  MIN 1.099466  AVG 1.099466
#| 0/ 4|: ( 0, 0)/[ 2, 2] ->STATS Runs 10  cutl CURRENT 28482.000000  MAX 29833.000000  MIN 26321.000000  AVG 28602.500000
#| 0/ 4|: ( 0, 0)/[ 2, 2] ->STATS Runs 10  cutn CURRENT 25237.000000  MAX 26335.000000  MIN 23499.000000  AVG 25365.500000
#| 0/ 4|: ( 0, 0)/[ 2, 2] ->STATS Runs 10  moveVol CURRENT 17839440.000000  MAX 17943520.000000  MIN 17619200.000000  AVG 17723204.000000
#| 0/ 4|: ( 0, 0)/[ 2, 2] ->STATS Runs 10  repart CURRENT 20687640.000000  MAX 20888820.000000  MIN 20283204.000000  AVG 20583454.000000


  $iters = 0;

  ($bal{MAX},$bal{MIN},$bal{AVG})		= ("ERROR","ERROR","ERROR");
  ($cutl{MAX},$cutl{MIN},$cutl{AVG})		= ("ERROR","ERROR","ERROR");
  ($cutn{MAX},$cutn{MIN},$cutn{AVG})		= ("ERROR","ERROR","ERROR");
  ($movevol{MAX},$movevol{MIN},$movevol{AVG})	= ("ERROR","ERROR","ERROR");
  ($repart{MAX},$repart{MIN},$repart{AVG})	= ("ERROR","ERROR","ERROR");
  $timeavg = -1;

  my $infd = @_[0];
  my $logfd = @_[1];

  while (<$infd>) {
    print $logfd $_;
    if (/.*STATS Runs ([0-9]+)/) {
      my $vname="";
      $iters = "$1";
      if (/ bal (.*)/) {
	%bal = (split ' ', $1);
      } elsif (/ cutl (.*)/) {
	%cutl = (split ' ', $1);
      } elsif (/ cutn (.*)/) {
	%cutn = (split ' ', $1);
      } elsif (/ moveVol (.*)/) {
	%movevol = (split ' ', $1);
      } elsif (/ repart (.*)/) {
	%repart = (split ' ', $1);
      }
    } elsif (/^FILE .*Average:  ([\d.e+]+) seconds per Iteration/) {
      $timeavg = $1 + 0;
    }
  }

#  print "iters: $iters\n";
#  print "bal: $bal{MAX} $bal{MIN} $bal{AVG}\n";
#  print "cutl: $cutl{MAX} $cutl{MIN} $cutl{AVG}\n";
#  print "cutn: $cutn{MAX} $cutn{MIN} $cutn{AVG}\n";
#  print "time: $timeavg\n";

}


sub printvars {
  my $filename = @_[0];
  open FD, ">$filename" or die "printout open failed on $filename";

  print FD "BalMax	numeric/double	$bal{MAX}\n";
  print FD "BalMin	numeric/double	$bal{MIN}\n";
  print FD "BalAvg	numeric/double	$bal{AVG}\n";
  print FD "CutlMax	numeric/double	$cutl{MAX}\n";
  print FD "CutlMin	numeric/double	$cutl{MIN}\n";
  print FD "CutlAvg	numeric/double	$cutl{AVG}\n";
  print FD "CutnMax	numeric/double	$cutn{MAX}\n";
  print FD "CutnMin	numeric/double	$cutn{MIN}\n";
  print FD "CutnAvg	numeric/double	$cutn{AVG}\n";
  print FD "MVolMax	numeric/double	$movevol{MAX}\n";
  print FD "MVolMin	numeric/double	$movevol{MIN}\n";
  print FD "MVolAvg	numeric/double	$movevol{AVG}\n";
  print FD "RepartMax	numeric/double	$repart{MAX}\n";
  print FD "RepartMin	numeric/double	$repart{MIN}\n";
  print FD "RepartAvg	numeric/double	$repart{AVG}\n";
  print FD "SolverTime	numeric/double	$timeavg\n";
  print FD "Iterations	numeric/integer	$iters\n";
  print FD $ckok;
  print FD "exit_comment text/string \"$erreason\"\n" if ($erreason);
  print FD "exit_status	numeric/integer	0\n";
  close FD;
}

##
## MAIN ROUTINE
##
if (! @ARGV) {
   print "$0 - Zoltan experiment driver\n";
   print "\n";
   print "   This script is launched by 'exact' to execute and process\n";
   print "   tests from an experiment file.\n";
   exit;
}

open LOGFILE, ">$logfile" or die "can't open $logfile";
open SAVEERR, ">&STDERR";
open STDERR, ">&LOGFILE" or die "can't dup stderr to logfile";

select STDERR;	# all output is now going to logfile, as STDERR
$| = 1;		# unbuffered writes

open INFILE, "<$ARGV[0]" or die "can't open $infile";
parseinfile(INFILE);		# populate %args
close INFILE;

if ($args{_exact_debug}) {
  print "DEBUG INFO ON\n";
  print "Input file: ", Dumper(\%args);
}

# can configure a short script to create a custom zdrive here, if desired
# but it uses the commandline to specify the name


#
# Setup command line
#

if (! $np) {
  die "number of processors undefined ($np)";
}

$cmdline = "$mpi $np $zdrive $zinp";
#print "CMDLINE: $cmdline\n" if ($args{_exact_debug});
print "CMDLINE: $cmdline\n";

open ZIN, $zinp or die "couldn't open $zinp";
my $temprs = $/;
undef $/;		# turn off input record separator
my $ztext = <ZIN>;	# slurp whole zdrive.inp
close ZIN;
$/ = $temprs;

unless ($ztext =~ /text output *= *0+\b/) {
  if ($zinp eq $zinp_def) {
    print "Running with default $zinp ... skipping diff check\n";
    $erreason .= "Using default $zinp, skipping diffs - ";
    break;
  }

  # construct list of output files
  my $l = length (sprintf "%d", $np);	# length of nproc
  my $ll = sprintf "%0${l}d", 0;	# 5 becomes "05" when 9 < np < 100
  if ($ztext =~ /File Name\s*=\s*(\w+)/) {
    $name = $1;
  } else {
    print "Couldn't determine file name ... skipping diff check\n";
    $erreason .= "File Name in zinp missing - ";
    break;
  }
  my $fullname = "$name.out.$np"; # needs the node number at end
  for (my $i = 0; $i < $np; $i++) {
    push @checkfiles,"$fullname.$ll";
    $ll++;				# string increment
  }

  # also need to read in checksum file
  if (!defined (open ZIN, "checksums") ) {
    $erreason .= "checksums missing (run mktests), skipping diff check - ";
    @checkfiles = ();
  } else {
    while (<ZIN>) {
      /(\d+)\s+(\d+)\s+(\S+)/;
      $cksum{"$3"} = [$1, $2];
    }
    close ZIN;
  }
}

# Launch the test, finally, dumping output to LOGFILE
# Before running, copy $zinp to zdrive.inp in case we're running zfdrive 
# (which doesn't take command-line arguments).

open ZOLTAN, "/bin/cp -f $zinp zdrive.inp; $cmdline |" or die "Can't run $cmdline";
parseoutfile(ZOLTAN, LOGFILE);
close ZOLTAN;

# check output file checksums, if @checkfiles non-empty
my $pass = 0;
my $fail = 0;
mkdir "output", 0755;	# place to store the output

foreach $f (@checkfiles) {
  my $sum = `csh -f ../bin/cksum-w $f`;
  $sum =~ /(\d+)\s+(\d+)/;
  my $sum0 = $1;
  my $sum1 = $2;
  my $partname = $zinp;
  my $fullname = $f;
  $partname =~ s/^.*${zinp_def}\.//;
  $fullname =~ s/^($name)\.out\./$1.$partname./;
  print "CRC: $fullname ?= $f\n" if ($args{_exact_debug});
  if (! defined($cksum{$fullname})) {
    $erreason .= "answer file $fullname not in checksums - ";
  }
  if ($sum eq "") {
    $erreason .= "output file $f not found! - ";
  }
  print "  cmp $sum0/$sum1 to $cksum{$fullname}[0]/$cksum{$fullname}[1]\n"
  	if ($args{_exact_debug});
  if ($cksum{$fullname}[0] == $sum0 && $cksum{$fullname}[1] == $sum1 ) {
    $pass++;
  } else {
    $fail++;
    $erreason .= "bad checksum for $f - ";
  }
  rename $f, "output/$fullname";		
}
$ckok = "OutputPass\tnumeric/integer\t$pass\n";
$ckok .= "OutputFail\tnumeric/integer\t$fail\n";
$ckok .= "AllPass\tnumeric/boolean\t" . (($pass > 0 && $pass == $np) ? "1\n" : "0\n");

printvars($outfile);

# fix things back
open STDERR, ">&SAVEERR";
select STDOUT;

