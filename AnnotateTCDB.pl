#use warnings;
$indb = $ARGV[0];
if($indb !~/\w/ || !-s "UR100vsTCDB.m8"){
        print "must specify the database file, ex: perl AnnotateTCDB.pl UNIREF100_INFO_OCT_2021.txt\n";
        print "and must have aligned the UR100 proteins to the TCDB.dmnd. See my GitHub: https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database/wiki/Home/_edit\n";
        die;}
$outdb=$indb;
$outdb=~s/\..*//;
$outdb.="wtcdbs.txt";
open(OUTINFO, ">", $outdb)||die;

#input transported compounds
qx{wget -N  https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py};
open(INSUB, "getSubstrates.py")||die;
while(<INSUB>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $tcdb, my $cpds)=split("\t",$_);
        @CHEBI = ($cpds=~/(CHEBI.\d+)/g);
        @CHEBI = sort(@CHEBI);
        $cpd=join(";",@CHEBI);
        if($tcdb=~/\w/){$TCDB_CPDS{$tcdb}=$cpd;}
}
$tc=keys %TCDB_CPDS;
print "$tc tcdbs with cpds\n";


#input the matchs
open(INTVU, "UR100vsTCDB.m8")||die;
while(<INTVU>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_,-1);
        $ur100 = $stuff[0];
        if($stuff[2]=~/([^\|]+$)/){$tcdb = $1;}
        else{$tcdb='';}
        $pid  = $stuff[9];
        if($pid < 60){next;}
        $cov  = $stuff[11];
        if($cov < 50){next;}
        $sco  = $pid*$cov;

        #get best hit
        $FOUND{$tcdb}=1;
        if($TCDB_CPDS{$tcdb}=~/CHEBI/){$HASCPD{$tcdb}=1;}
        $UR_TCDB{$ur100}{$tcdb}=$sco;
}
$um = keys %UR_TCDB;
$hc = keys %HASCPD;
$ft = keys %FOUND;
undef(%FOUND);
undef(%HASCPD);
print "have $um ur100s and $ft found tcdbs $hc have cpds\n";


#get top tcdb and top tcdb with compounds (for some reason tcdb has other proteins ::shrugs::
foreach my $ur100 (keys %UR_TCDB){
        $max=0;
        $top='';
        $twf='';
        $cpds='';
        foreach my $tcdb (sort{$UR_TCDB{$ur100}{$b}<=>$UR_TCDB{$ur100}{$a}} keys %{$UR_TCDB{$ur100}}){
                if($UR_TCDB{$ur100}{$tcdb}>$max){$max=$UR_TCDB{$ur100}{$tcdb}; $top=$tcdb;}
                if($TCDB_CPDS{$tcdb}=~/\w/ && $twf eq ''){ $twf=$tcdb; $cpds=$TCDB_CPDS{$tcdb}; last;}
        }
        delete($UR_TCDB{$ur100});
        if($top eq $twf){$UR_TCDB{$ur100}=$top;}
        else{$UR_TCDB{$ur100}=$top.";".$twf;}
        if($cpds=~/\w/){$UR_CPDS{$ur100}=$cpds;}
}
$tc=keys %UR_TCDB;
$uc=keys %UR_CPDS;
print "have $tc total ur100 with tcdb and $uc with compounds\n";

#in/output annotatd UNIREF100
open(ININFO, $indb)||die;
while(<ININFO>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_,-1);
        $ur100 = $stuff[0];
        if(exists($UR_TCDB{$ur100})){$stuff[11]=$UR_TCDB{$ur100};}
        if(exists($UR_CPDS{$ur100})){$stuff[24]=$UR_CPDS{$ur100}; $added++;}
        $out=join("\t",@stuff);
        print OUTINFO "$out\n";
        $on++;
        if($on%1000000==0){print "on $on added $added\n";}
}
print "added $added out of total prots $on\n";
close(ININFO);
close(OUTINFO);

if(-s $outdb > -s $indb){ qx{mv $outdb $indb}; }
