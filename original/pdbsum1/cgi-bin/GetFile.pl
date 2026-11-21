#!/usr/bin/perl -sw

#------------------------------------------------------------------------
#
# GetFile.pl - PerlScript to download PDBsum1 zip file and log user
#              details (where entered)
#
# Roman Laskowski, 22 Jul 2022
#
#------------------------------------------------------------------------

use lib('./');
use CGI qw/:standard/;

#-- Required libraries
$dir = `pwd`;
$dir =~ m/(.*)\/.*/;
$cathparam = $1.'/param/CATHPARAM';
$ENV{'CATHPARAM'} = "$cathparam";
require "ReadCATHPARAM.pl";
require "dumpvar.pl";

$SIG{'TERM'} = 'SigHand';
$SIG{CHLD} = 'IGNORE';

#
#-- Initialise
#
$date = `date`;
chop($date);

$outcome = "OK";
$option = " ";
$zip_file = "NONE";

#
#-- Get the parameters (paths to directories, etc)
#
%param = ReadCATHPARAM("$cathparam");

#
#-- Set all the required variables using the parameter data
#
SetVars();

#
#-- Get URL
#
$address = $ENV{'REMOTE_ADDR'};

#
#-- Get the input parameters
#
$q = new CGI;

#
#-- Retrieve the input parameters from the form
#
$fname         = $q->param('fname');
$lname         = $q->param('lname');
$email         = $q->param('email');

$debug         = $q->param('debug');

#
#-- Initialise uninitialised variables
#
if (!$debug)
  {
     $debug = "FALSE";
  }

#
#-- If any details not entered, set to none
#
if (!$fname)
  {
     $fname = "-*-";
  }
if (!$lname)
  {
     $lname = "-*-";
  }
if (!$email)
  {
     $email = "=+=";
  }


#
#-- Form the file names
#
$dir = $pdbsum1_distrib_dir;
$version_txt = $dir."/version.txt";
$zip_file = $dir."/pdbsum1.zip";
$log_file = $logs_dir."/PDBsum1.log";

#
#-- Get the program version from the version.txt fil
#
$version = `cat $version_txt`;
chomp($version);

#
#-- If all OK, then return the zip file
#
if ($debug eq "FALSE")
  {
    #
    #-- Prepare log details
    #
    $log_details = $date." : v.".$version." Name: ".$fname." ".$lname;
    $log_details = $log_details." From: ".$address."  ".$email;

    #
    #-- Write the details to the log file
    #
    open LOGGER, ">>$log_file";
    print LOGGER $log_details."\n";
    close(LOGGER);

    #
    #-- Define output type as download
    #
    print "Content-disposition: attachment; filename=pdbsum1.zip\n\n";

    #
    #-- Return the licence file
    #
    open(F,"/bin/cat $zip_file |");
    while(<F>)
      {
        print;
      }
    close(F);
  }

#
#-- Otherwise, display debug message
#
else
  {
    #
    #-- Define MIME type as HTML page
    #
    print "Content-type: text/html\n\n";

    #
    #-- Show error message
    #
    print "<H1>Debugging</H1>\n";
    print "<P>\n";
    print "<H2>Outcome</H2>\n";
    print $outcome."\n";
    print "<H2>Option</H2>\n";
    print $option."\n";
    print "<H2>CATHPARAM</H2>\n";
    print $cathparam."\n";
    print "<H2>First name:</H2>\n";
    print $fname."\n";
    print "<H2>Last name:</H2>\n";
    print $lname."\n";
    print "<H2>E-mail:</H2>\n";
    print $email."\n";
    print "<H2>Zip file:</H2>\n";
    print $zip_file."\n";
    print "<H2>Version:</H2>\n";
    print $version."\n";
  }

exit(0);

#-----------------------------------------------------------------------------

#
# Various subroutines
#

sub SetVars {
    $logs_dir                        = $param{'LOGS_DIR'};
    $pdbsum1_distrib_dir             = $param{'PDBSUM1_DISTRIB_DIR'};
}

sub SigHand {
  exit(1);
}

