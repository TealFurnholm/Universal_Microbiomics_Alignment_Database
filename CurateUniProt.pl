#use warnings;

#GET UPDATE INFO
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$month = $1; $year = $2;
$version=$1."_".$2;

#CHECK INPUTS
if(!-s "UniRef100_Len_AA.txt" && !-s "uniref100.fasta.gz"){
        print "missing uniref100.fasta.gz or UniRef100_Len_AA.txt, downloading...\n";
        qx{wget -N https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz};
}
if(!-s "uniprot_to_uniref.txt.gz"){
        if(!-s "idmapping.dat.gz"){
                print "missing idmapping.dat.gz, downloading...\n";
                qx{wget -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz};
        }
        print "missing uniprot_to_uniref.txt.gz, extracting data...\n";
        qx{zcat idmapping.dat.gz | grep -P "(UniRef100|UniRef90)" | gzip > uniprot_to_uniref.txt.gz};
}
$inarx = 'all_reaction_info_'.$month.'_'.$year.'.txt';  if(!-s $inarx ){print "missing $inarx - please place into this folder and restart\n"; $fail=1;}
qx{wget -N https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py};
$inupr = 'all_upid_rxns_'.$month.'_'.$year.'.txt';      if(!-s $inupr){ print "missing $inupr - please place into this folder and restart\n"; $fail=1;}
$inkeg = 'all_kegg_rxns_'.$month.'_'.$year.'.txt';      if(!-s $inkeg ){print "missing $inkeg - please place into this folder and restart\n"; $fail=1;}
$inecs = 'all_ec_rxns_'.$month.'_'.$year.'.txt';        if(!-s $inecs ){print "missing $inecs - please place into this folder and restart\n"; $fail=1;}
$intax = 'TAXONOMY_DB_'.$year.'.txt';                   if(!-s $intax ){print "missing $intax - please place into this folder and restart\n"; $fail=1;}
$inmon = 'all_mono_rxns_'.$month.'_'.$year.'.txt';      if(!-s $inmon ){print "missing $inmon - please place into this folder and restart\n"; $fail=1;}
$inup  = 'uniprot-all.tab.gz';                          if(!-s $inup ){ print "missing $inup - please place into this folder and restart\n";  $fail=1;}
if($fail==1){
        print "missing needed input files.\n";
        print "see my GitHub repositories: https://github.com/TealFurnholm/ for further guidance\n";
        die;
}
open(DEBUG, ">", "curate_uniprot_".$month.'_'.$year.".debug")||die;



#GET PROT LENGTH AND ODD AA DISTs
$on=0;
print "INPUT FASTA\n";
if(-s "UniRef100_Len_AA.txt"){
        print "Input Existing Prot Lengths\n";
        open(INODD, "UniRef100_Len_AA.txt")||die;
        while(<INODD>){
                if($_!~/\w/){next;}
                $_=uc($_);
                $_=~s/[\n\r]+//;
                (my $ur100, my $len, my $odd)=split(";",$_);
                $PROT_LEN{$ur100}=$len;
                $ODD_AA{$ur100}=$odd;
                if($on%1000000==0){print "on $on len $len odd $odd\n";} $on++;
        }
}
else{   #GET ODD AA DIST PER PROTEIN
        open(OUTODD,">>","UniRef100_Len_AA.txt")||die;
        open(INURFA, "gunzip -c uniref100.fasta.gz |")||die;
        while(<INURFA>){
                if($_!~/\w/){next;}
                $_=uc($_);
                $_=~s/[\n\r]+//;

                if($_=~/(UNIREF\S+)/){
                        $np=$1;
                        $seq=join("",@SEQ);
                        if($seq=~/^[A-Z]{10,}$/ && !exists($PROT_LEN{$prot})){
                                $seq=~s/[^A-Z]+//g;
                                $len = length($seq);
                                $PROT_LEN{$prot}=$len;
                                @AA=();
                                for my $i (A..Z){
                                        my $count = () = $seq =~ /$i/g;
                                        $per=$count*100/$len;
                                        $per=~s/\..*//;
                                        if($per >= 10){$aa=$i.":".$per; push(@AA,$aa);}
                                }
                                $odd=join("|",@AA);
                                if($odd=~/[A-Z]/){ $ODD_AA{$prot}=$odd; }
                                print OUTODD "$prot\t$PROT_LEN{$prot}\t$ODD_AA{$prot}\n";
                        }
                        if($on%1000000==0){$time=localtime; print "on $on time $time prot $prot odd $ODD_AA{$prot} len $PROT_LEN{$prot}\n";} $on++;
                        $prot=$np; @SEQ=();
                }
                if($_=~/^[A-Z]+$/){push(@SEQ,$_);}
}       }
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;



#GET UNIPROT TO UNIREF IDS
$on=0; $time=localtime;
print "INPUT UNIREF TO UNIPROT MAPPING $time\n";
open(INMAP, "gunzip -c uniprot_to_uniref.txt.gz |")||die;
while(<INMAP>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $prot, my $type, my $id)=split("\t",$_);
        if($type eq "UNIREF100"){$PID_UR100{$prot}=$id;}
        if($type eq "UNIREF90"){ $PID_UR90{$prot}=$id;}
        if($on%1000000==0){
                $time=localtime; $tpo = keys %PID_UR100; $tpn=keys %PID_UR90;
                print "on $on time $time prot $prot tpo $tpo tpn $tpn\n";
        }
        $on++;
}
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;



$time=localtime;
print "INPUT REACTION INFO $time\n";
open(INARX,$inarx)||die;
while(<INARX>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
        $rx =$stuff[0];
        $lc =$stuff[2];
        $lc.=";".$stuff[3];
        $lc.=";".$stuff[4];
        $rc =$stuff[5];
        $rc.=";".$stuff[6];
        $rc.=";".$stuff[7];
        $rx.=";".$stuff[8];
        $rx.=";".$stuff[9];
        $rx.=";".$stuff[10];
        @rx=split(";",$rx); %seen=(); @rx = grep { !$seen{$_}++ } @rx;
        foreach my $x (@rx){
                if($x !~ /\w/){next;}
                $LEFT_CPDS{$x}=$lc;
                $RITE_CPDS{$x}=$rc;
                foreach my $y (@rx){
                        if($y !~ /\w/){next;}
                        $RXN_ALTS{$x}{$y}=1;
}       }       }
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;




$time=localtime;
print "INPUT TRANSPORTER SUBSTRATES $time\n";
open(INTRCH, "getSubstrates.py")||die "unable to open getSubstrates.py: $!\n";
while(<INTRCH>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\n\r]+//;
        (my $tcdb, my $chebs)=split("\t", $_);
        $tcdb="TCDB:".$tcdb;
        @CHEBS=();
        @CHEBS = ( $chebs =~ /(CHEBI.\d+)/g );
        $chebi=join(";",@CHEBS);
        $TRCH{$tcdb}=$chebi;
}
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$TRCH{$tcdb}=$chebs;




$time=localtime;
print "INPUT UNIPROT REACTIONS $time\n";
open(INUPR, $inupr)||die;
while(<INUPR>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
        $pid=shift(@stuff);
        foreach my $x (@stuff){ if($x=~/\w/){$UPID_RXNS{$pid}{$x}=1;} }
}
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$TRCH{$tcdb}=$chebs;
#$UPID_RXNS{$pid}{$x}=1;





$time=localtime;
print "INPUT KEGG REACTIONS $time\n";
open(INKEG, $inkeg)||die;
while(<INKEG>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $gene, my $krxn)=split("\t",$_);
        if($gene!~/\w/ || $krxn!~/\w/){next;}
        $KGEN_KRXN{$gene}=$krxn;
}
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$TRCH{$tcdb}=$chebs;
#$UPID_RXNS{$pid}{$x}=1;
#$KGEN_KRXN{$gene}=$krxn;




$time=localtime;
print "INPUT EC REACTIONS $time\n";
open(INECS, $inecs)||die;
while(<INECS>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $ec, my $rxns)=split("\t",$_);
        if($ec !~/\w/ || $rxns !~/\w/){next;}
        $EC_RXNS{$ec}=$rxns;
}
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$TRCH{$tcdb}=$chebs;
#$UPID_RXNS{$pid}{$x}=1;
#$KGEN_KRXN{$gene}=$krxn;
#$EC_RXNS{$ec}=$rxns;




$time=localtime;
print "INPUT TAXONOMY $time\n";
open(INTAX, $intax)||die;
while(<INTAX>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
        $tid=shift(@stuff);
        $lin=join(";",@stuff);
        if($stuff[6]=~/\w/){$SPECIES{$stuff[6]}=join(";",@stuff[0..6]);}
        if($stuff[5]=~/\w/){$GENUS{$stuff[5]}=join(";",@stuff[0..5]);}
        $PHY{$tid}=$lin;
}
close(INTAX);
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$TRCH{$tcdb}=$chebs;
#$UPID_RXNS{$pid}{$x}=1;
#$KGEN_KRXN{$gene}=$krxn;
#$EC_RXNS{$ec}=$rxns;
#$SPECIES{$stuff[6]}=lin0..6
#$GENUS{$stuff[5]}=lin0..5
#$PHY{$tid}=$lin;




$time=localtime;
print "INPUT MONOMERS $time\n";
open(INMON, $inmon)||die;
while(<INMON>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $mono, my $rxns)=split("\t",$_);
        $MONO_RXNS{$mono}=$rxns;
}
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
#$PID_UR100{$prot}=$id;
#$PID_UR90{$prot}=$id;
#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$TRCH{$tcdb}=$chebs;
#$UPID_RXNS{$pid}{$x}=1;
#$KGEN_KRXN{$gene}=$krxn;
#$EC_RXNS{$ec}=$rxns;
#$SPECIES{$stuff[6]}=lin0..6
#$GENUS{$stuff[5]}=lin0..5
#$PHY{$tid}=$lin;
#$MONO_RXNS{$mono}=$rxns;





$on=0;
print "INPUT UNIPROT\n";
open(INFO, "gunzip -c $inup |")||die;
open(OUTPINT, ">", "UNIPROT_INFO_".$version.".txt")||die;
print OUTPINT "PID\tName\tLength\tUR100\tUR90\tTaxonId\tLineage\t";
print OUTPINT "SigPep\tTMS\tDNA\tMetal\tTCDB\tLoc\tCOG\tPfam\tTigr\tGene_Ont\tInterPro\tECs\tkegg\trhea\tbiocyc\tLeftCPDs\tRightCPDs\tTransCPDs\n";
while(<INFO>){
        if($_!~/\w/){next;}
        $on++;
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
        if($stuff[3]!~/\d/){next;}
        $pid   = $stuff[0];
        $plen  = $stuff[2];
        $ur100=$PID_UR100{$pid};
        $ur90 =$PID_UR90{$pid};
        if($PROT_LEN{$ur100}!~/\d/){$PROT_LEN{$ur100}=$plen;}
        if($plen!~/\d/){$plen=$PROT_LEN{$ur100};}

        $first="$pid\t$name\t$plen\t$ur100\t$ur90";

        #GET OR FIX TAXONOMY
        $tid=$stuff[3];
        if(!exists($PHY{$tid})){
                $genus = $stuff[4];
                $species = $stuff[5]; $species =~ s/\s*\(.*//; $species=fix_names($species);
                 $strain = $stuff[6];  $strain =~ s/\s\/\s.*//; $strain=fix_names($strain);
                   if(exists($GENUS{$genus})){ @LIN=split(";",$GENUS{$genus}); push(@LIN,$species); if($strain ne $species){ push(@LIN, $strain); }}
                elsif(exists($SPECIES{$species})){ @LIN=split(";",$SPECIES{$species}); if($strain ne $species){ push(@LIN, $strain); }}
                else{ print DEBUG "missing tid\t$tid\n"; next;}
                $lin=join(";",@LIN);
                $PHY{$tid}=$lin;
                $lin=join("\t",@LIN);
                open(OUTTAX, ">>", $intax)||die;
                print OUTTAX "$tid\t$lin\n";
        }
        $lin=$PHY{$tid};



        #GET FUNCTIONS
        @TMS=(); @COG=(); @GOS=(); @ECS=(); @MON=(); @METS=(); @LOCS=(); @KEGS=(); @RHEA=();
        $sig='';        if($stuff[7]=~/SIGNAL\W+(\d+)\.\.(\d+)/){                                               $sig = "SIGNAL:".$1."..".$2;}
        $tm='';         if($stuff[8]=~/TRANSMEM/){      @TMS = ($stuff[8]=~/TRANSMEM\W(\d+\.\.\d+)/g);          $tm  = join(";",@TMS); $tm="TMS:".$tm;}
        $trc='';        if($stuff[9]=~/([\w\.]+)/){                                                             $trc = "TCDB:".$1;}
        $cog='';        if($stuff[10]=~/[CK]OG\d+/){    @COG = ($stuff[10]=~/\b([CK]OG\d\d\d\d)\b/g);           $cog = join(";",@COG);}
        $pfa='';        if($stuff[11]=~/(PF.*\d+)/){                                                            $pfa = $1;}
        $tig='';        if($stuff[12]=~/(TIGR.*\d+)/){                                                          $tig = $1;}
        $gos='';        if($stuff[13]=~/GO.\d+/){       @GOS = ($stuff[13]=~/(GO.\d+)/g);                       $gos = join(";",@GOS);}
        $ipr='';        if($stuff[14]=~/(IPR.*\d+)/){                                                           $ipr = $1;}
        $ecs='';        if($stuff[15]=~/[\d\.\-]+/){    @ECS = ($stuff[15]=~/([\d\.\-]+)/g);                    $ecs = join(";",@ECS);}
        $bioc='';       if($stuff[16]=~/^(\w.*\w)/){    @MON = split(";",$stuff[16]);                           }
        $dna='';        if($stuff[17]=~/DNA.BIND\s(\d+\.\.\d+).*?NOTE\=\"([^\"]+)/){$cln=CleanNames($2);        $dna = "DNA:".$1."|".$cln;}
        $met='';        if($stuff[18]=~/NOTE/){         @METS = ($stuff[18]=~/NOTE\W+(\w\S+?)[\"\;]+/g );
                                                        %seen=(); @METS = grep{ !$seen{$_}++ } @METS;           $met = join(";", @METS);}
        $loc='';        if($stuff[19]=~/LOCATION/){     @LOCS = ($stuff[19]=~/([^\.\:\;]+)\{|([^\.\:\;]+)NOTE/g);
                                                        for my $i (0..$#LOCS){$LOCS[$i]=CleanNames($LOCS[$i]);}
                                                        %seen=(); @LOCS = grep{ !$seen{$_}++ } @LOCS;           $loc = join("|", @LOCS);}
        $kegg='';       if($stuff[20]=~/\w+/){          $stuff[20]=~s/\s//g; @KEGS=split(";",$stuff[20]);}
        $rhea='';       if($stuff[21]=~/\w+/){          $stuff[21]=~s/\s//g; @RHEA=split(";",$stuff[21]);}

        #GET REACTIONS
        %RXNS=();
        #rhearxns
        foreach my $rh (@RHEA){ if($rh=~/\w/){$RXNS{$rh}++; }}
        #uprxns
        if(keys %{$UPID_RXNS{$pid}}>0){ foreach my $rx (keys %{$UPID_RXNS{$pid}}){ if($rx=~/\w/){$RXNS{$rx}++; }}}
        #ecrxns
        foreach my $ec (@ECS){  if($EC_RXNS{$ec}=~/\w/){   @RXNS=split(";",$EC_RXNS{$ec});   foreach my $rx (@RXNS){ if($rx=~/\w/){$RXNS{$rx}++; }}}}
        #keggrxns
        foreach my $kg (@KEGS){ if($KGEN_KRXN{$kg}=~/\w/){ @RXNS=split(";",$KGEN_KRXN{$kg}); foreach my $rx (@RXNS){ if($rx=~/\w/){$RXNS{$rx}++; }}}}
        #monomers
        foreach my $mn (@MON){  if($MONO_RXNS{$mono}=~/\w/){ @RXNS=split(";",$MONO_RXNS{$mono}); foreach my $rx (@RXNS){ if($rx=~/\w/){$RXNS{$rx}++; }}}}
        #rxnalts
        foreach my $rx (keys %RXNS){ if($rx=~/\w/){foreach my $alt (keys %{$RXN_ALTS{$rx}}){ if($alt=~/\w/){$RXNS{$alt}++; }}}}

        #GET TOP OF EACH RXN TYPE
        $rm=0; $km=0; $bm=0;
        @RHE=(); @KEG=(); @BIO=(); %ALTS=();
        foreach my $rxn (sort{$RXNS{$b}<=>$RXNS{$a}} keys %RXNS){
                if($rxn !~/\w/){next;}
                   if($rxn=~/^\d+$/){   if($RXNS{$rxn}>$rm){$rm=$RXNS{$rxn};} if($RHE[3]!~/\w/ && $RXNS{$rxn} > $rm/2){ push(@RHE,$rxn); $ALTS{$rxn}=1; }}
                elsif($rxn=~/^R\d+$/){  if($RXNS{$rxn}>$km){$km=$RXNS{$rxn};} if($KEG[3]!~/\w/ && $RXNS{$rxn} > $km/2){ push(@KEG,$rxn); $ALTS{$rxn}=1; }}
                elsif($rxn=~/RXN/){     if($RXNS{$rxn}>$bm){$bm=$RXNS{$rxn};} if($BIO[3]!~/\w/ && $RXNS{$rxn} > $bm/2){ push(@BIO,$rxn); $ALTS{$rxn}=1; }}
                else{next;}
        }
        $kegg=join(";",@KEG);
        $rhea=join(";",@RHE);
        $bioc=join(";",@BIO);


        #GET L/R/T COMPOUNDS
        %LEFTS=(); %RITES=();
        foreach my $rxn (keys %ALTS){
                if($rxn!~/\w/){next;}
                @LCPDS=split(";",$LEFT_CPDS{$rxn});
                foreach my $cpd (@LCPDS){ if($cpd=~/\w/){ $LEFTS{$cpd}=1; }}
                @RCPDS=split(";",$RITE_CPDS{$rxn});
                foreach my $cpd (@RCPDS){ if($cpd=~/\w/){ $RITES{$cpd}=1; }}
        }

        @CL=(); @KL=(); @BL=();
        foreach my $cpd (sort(keys %LEFTS)){
                if($cpd !~/\w/){next;}
                if($cpd=~/CHEBI/){      push(@CL,$cpd);}
                elsif($cpd=~/^C\d+$/){  push(@KL,$cpd);}
                else{                   push(@BL,$cpd);}
        }
        @LEFTS=(); $lcpd='';
        $cl=join(";",@CL); if($cl=~/\w/){push(@LEFTS,$cl);}
        $kl=join(";",@KL); if($kl=~/\w/){push(@LEFTS,$kl);}
        $bl=join(";",@BL); if($bl=~/\w/){push(@LEFTS,$bl);}
        $lcpd=join(";", @LEFTS);

        @CR=(); @KR=(); @BR=();
        foreach my $cpd (sort(keys %RITES)){
                if($cpd !~/\w/){next;}
                if($cpd=~/CHEBI/){      push(@CR,$cpd);}
                elsif($cpd=~/^C\d+$/){  push(@KR,$cpd);}
                else{                   push(@BR,$cpd);}
        }
        @RITES=(); $rcpd='';
        $cr=join(";",@CR); if($cr=~/\w/){push(@RITES,$cr);}
        $kr=join(";",@KR); if($kr=~/\w/){push(@RITES,$kr);}
        $br=join(";",@BR); if($br=~/\w/){push(@RITES,$br);}
        $rcpd=join(";",@RITES);
        $trcpd=''; if(exists($TRCH{$trc}) && $trc=~/\w/){ $trcpd=$TRCH{$trc}; }

        # FIX NAMES
        @NAMES=split('\(',$stuff[1]);
        @GN=();
        foreach my $n (@NAMES){ $n=CleanNames($n); if($n!~/[A-Z]{4}/){next;} push(@GN,$n);}
        @NAMES=();
        while($GN[0]=~/\w/){$n=shift(@GN); if(!grep /\Q$n\E/i, @NAMES && !grep /\Q$n\E/i, @GN){ push(@NAMES,$n); }}
        $name=join(";",@NAMES);
        $odd=$ODD_AA{$ur100};
        if($odd=~/[A-Z]/){$name.=";".$odd;}


        #PREP FUNC OUTPUT
        $first="$pid\t$name\t$plen\t$ur100\t$ur90\t$tid\t$lin";
        $funcs="$sig\t$tm\t$dna\t$met\t$trc\t$loc\t$cog\t$pfa\t$tig\t$gos\t$ipr\t$ecs\t$kegg\t$rhea\t$bioc\t$lcpd\t$rcpd\t$trcpd";

        if($on%1000000==0){$time=localtime; print "on $on time $time $first\n$funcs\n";}

        print OUTPINT "$first\t";
        print OUTPINT "$funcs\n";
}
print "on $on uniprots\n";

sub fix_names{
        my $name = $_[0];

        #fix the species/strain name
        $name =~ s/([A-Z]+)\s(PROTEOBACTER)(IA|IUM)/$1$2$3/;
        $name =~ s/\bPROPIONIBACTERIUM/CUTIBACTERIUM/g;
        $name =~ s/\bLEPIDOSAURIA/SAURIA/g;
        $name =~ s/ENDOSYMBIONT.OF\s+/ENDOSYMBIONT-/;
        $name =~ s/COMPOSITE.GENOME.*//;
        $name =~ s/MARINE.GROUP.(\w+)/MARINE-GROUP-$1/;
        $name =~ s/\s+METAGENOME//;
        $name =~ s/OOMYCETES/OOMYCOTA/;
        $name =~ s/LILIOPSIDA/MAGNOLIOPSIDA/;
        $name =~ s/^NR[^A-Z]//;
        $name =~ s/.*INCERTAE.SEDIS.*//;
        $name =~ s/\_(PHYLUM|CLASS|ORDER|FAMILY|GENUS)[\b\_].*/\_/;
        $name =~ s/ENRICHMENT.CULTURE.CLONES*|ENRICHMENT.CULTURES*//;
        $name =~ s/\_(SENSU\_LATO|AFF|GEN|CF)\_/\_/g;
        $name =~ s/^(SENSU\_LATO|AFF|GEN|CF)\_//g;
        $name =~ s/\b(SENSU\_LATO|AFF|GEN|CF)\_//g;

        #remove ambiguous junk
        $name =~ s/(CANDIDATUS|CANDIDATUAS|CANDIDATE.\S+|VOUCHERED|UNDESCRIBED|FRAGMENT|UNSCREENED|UNKNOWN|UNCULTIVATED|UNCULTURED)\s*/\_/g;
        $name =~ s/(UNIDENTIFIED|UNCLASSIFIED|CONTAMINATION|SCREENED|UNASSIGNED|PUTATIVE|\-*LIKE)\s*/\_/g;

        #remove junk punctuation/standardize
        $name =~ s/\s+/_/g;
        $name =~ s/[^\w]+/_/g;
        $name =~ s/\_+/\_/g;
        $name =~ s/(^\_+|\_+$)//g;
        $name =~ s/^(X|CF)\_//;

        return($name);
}


sub CleanNames{
        $name = $_[0];

        #remove pointless ambiguators
        $name =~ s/(CANDIDATUS|CANDIDATUAS|CANDIDATE|VOUCHERED|UNDESCRIBED|UNSCREENED|UNKNOWN|UNCULTIVATED|UNCULTURED)\s*/UNCHARACTERIZED\_/g;
        $name =~ s/(UNIDENTIFIED|UNCLASSIFIED|CONTAMINATION|SCREENED|UNASSIGNED|PUTATIVE|HYPOTHETICAL)\s/UNCHARACTERIZED\_/g;
        $name =~ s/(UNCHARACTERIZED\_)+/UNCHARACTERIZED\_/g;
        $name =~ s/\-*LIKE\s*/\_/g;

        #remove junk punctuation/standardize
        $name =~ s/\s+/_/g;
        $name =~ s/[^\w]+/_/g;
        $name =~ s/\_+/\_/g;
        $name =~ s/(^\_+|\_+$)//g;

        return($name);
}

