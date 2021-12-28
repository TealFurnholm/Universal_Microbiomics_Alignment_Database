#use warnings;
use Sort::Naturally 'nsort';
#die;


##############################################################
##############################################################
#GET INITIAL FILES
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$version=$1."_".$2;

#CHECK FILE STATUS
$inup  	= 'uniprot-all.tab.gz';			if(!-s $inup  ){print "missing $inup - please place into this folder and restart\n";  	$fail=1;}
$inupur	= 'uniprot_to_uniref.txt.gz';		if(!-s $inupur){print "missing $inupur - please place into this folder and restart\n"; 	$fail=1;}
$intc	= 'UR100vsTCDB.m8';			if(!-s $intc  ){print "missing $intc - please place into this folder and restart\n";  	$fail=1;}
$insu	= 'getSubstrates.py';			if(!-s $insu  ){print "missing $insu - please place into this folder and restart\n";    $fail=1;}
$inarx 	= 'all_reaction_info_'.$version.'.txt';	if(!-s $inarx ){print "missing $inarx - please place into this folder and restart\n"; 	$fail=1;}
$inupr 	= 'all_upid_rxns_'.$version.'.txt';	if(!-s $inupr ){print "missing $inupr - please place into this folder and restart\n"; 	$fail=1;}
$inkeg 	= 'all_kegg_rxns_'.$version.'.txt';	if(!-s $inkeg ){print "missing $inkeg - please place into this folder and restart\n"; 	$fail=1;}
$inecs 	= 'all_ec_rxns_'.$version.'.txt';	if(!-s $inecs ){print "missing $inecs - please place into this folder and restart\n"; 	$fail=1;}
$inmon 	= 'all_mono_rxns_'.$version.'.txt';	if(!-s $inmon ){print "missing $inmon - please place into this folder and restart\n"; 	$fail=1;}
$intax 	= 'TAXONOMY_DB_'.$version.'.txt';	if(!-s $intax ){print "missing $intax - please place into this folder and restart\n"; 	$fail=1;}
if($fail==1){
	print "missing needed input files.\n"; 
	print "see my GitHub repositories: https://github.com/TealFurnholm/ for further guidance\n";
	die;
}

#OPEN INPUT FILES
open(INFO, "gunzip -c $inup |")||die;
open(INMAP, "gunzip -c $inupur |")||die;
open(INTVU, $intc )||die;
open(INTRCH,$insu )||die;
open(INARX, $inarx)||die;
open(INUPR, $inupr)||die;
open(INKEG, $inkeg)||die;
open(INECS, $inecs)||die;
open(INMON, $inmon)||die;
open(INTAX, $intax)||die;

#OPEN OUTPUT FILES
open(OUTPINT, ">", "XUNIPROT_INFO_".$version."int.txt")||die;
open(OUTCOMBO,">", "XFUNCTION_COMBOS.txt")||die;
open(LOG, ">", "XCURATE_UNIPROT_".$version.".log")||die;
print LOG "Curatin UniProt Reference Database Version $version\n";
##############################################################
##############################################################





##############################################################
##############################################################
#GET PROT LENGTH AND ODD AA DISTs
$on=0;
$time = localtime;
$/=">";
print "GET PROT LEN AND ODD AA STATS\n";
open(INURFA, "gunzip -c uniref100.fasta.gz |")||die;
while(<INURFA>){
	if($_!~/\w/){next;}
	$_=uc($_);
	@stuff=split("\n",$_);
	$header = shift(@stuff);
	$seq = join("",@stuff);
	$header =~ /(UNIREF\S+)/;
	$ur100=$1;
	$seq =~ s/[^A-Z]//g;
	$len = length($seq);
	if($len < 10){next;}
	$PROT_LEN{$ur100}=$len;
	@AA=();
	for my $i (A..Z){ 
		my $count = () = $seq =~ /$i/g; 
		$per=$count*100/$len;	
		$per=~s/\..*//;
		if($per >= 10){$aa=$i.":".$per; push(@AA,$aa);}
	}
	$odd=join("|",@AA);
	if($odd=~/[A-Z]/){ $ODD_AA{$ur100}=$odd; }
	if($on%1000000==0){$time=localtime; 
		print "on $on time $time ur100 $ur100 odd $odd len $len\n";
	} $on++;
}
$/="\n";

#GET UNIPROT TO UNIREF IDS
$on=0; $time=localtime;
print "INPUT UNIREF TO UNIPROT MAPPING $time\n";
while(<INMAP>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $prot, my $type, my $id)=split("\t",$_);
	if($type eq "UNIREF100"){$PID_UR100{$prot}=$id;}
	if($type eq "UNIREF90"){ $PID_UR90{$prot}=$id;}
	if($on%1000000==0){
		$time=localtime; 
		$tpo = keys %PID_UR100; 
		$tpn = keys %PID_UR90; 
		print "on $on time $time prot $prot tpo $tpo tpn $tpn id $id\n";
	}
	$on++; 
}

#GET TCDB SUBSTRATES
$time=localtime;
print "INPUT TRANSPORTER SUBSTRATES $time\n";
while(<INTRCH>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\n\r]+//;
        (my $tcdb, my $chebs)=split("\t", $_);
        if($tcdb!~/TCDB/){$tcdb="TCDB:".$tcdb;}
        @CHEBS=();
        @CHEBS = ( $chebs =~ /(CHEBI.\d+)/g );
        $chebi=join(";",@CHEBS);
        $TCDB_CPDS{$tcdb}=$chebi;
}

#INPUT TCDB ALIGNMENTS
$time=localtime;
print "INPUT UNIREF100 vs TCDB MATCHES $time\n";
while(<INTVU>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_,-1);
        $ur100 = $stuff[0];
        if($stuff[2]=~/([^\|]+$)/){$tcdb = "TCDB:".$1;}
        else{$tcdb='';}
        $pid  = $stuff[9];
        if($pid < 60){next;}
        $cov  = $stuff[11];
        if($cov < 50){next;}
        $sco  = $pid*$cov;
 
        #get best hit, weight to tcdbID w/cpds
        $FOUND{$tcdb}=1;
        if(!exists($UR100_TCDB{$ur100})){
		$UR100_TCDB{$ur100}=$tcdb; $UR_TCDB_SCO{$ur100}=$sco; 
		if($TCDB_CPDS{$tcdb}=~/CHEBI/){$HASCPD{$ur100}=1;}
	}
        elsif($TCDB_CPDS{$tcdb}=~/CHEBI/ && !exists($HASCPD{$ur100})){    $UR100_TCDB{$ur100}=$tcdb; $UR_TCDB_SCO{$ur100}=$sco; $HASCPD{$ur100}=1; }
        elsif($TCDB_CPDS{$tcdb}=~/CHEBI/ && $sco > $UR_TCDB_SCO{$ur100}){ $UR100_TCDB{$ur100}=$tcdb; $UR_TCDB_SCO{$ur100}=$sco; $HASCPD{$ur100}=1;}
        elsif($sco > $UR_TCDB_SCO{$ur100}){ $UR100_TCDB{$ur100}=$tcdb; $UR_TCDB_SCO{$ur100}=$sco; }
        else{}
}
undef(%FOUND);
undef(%HASCPD);
undef(%UR_TCDB_SCO);

#$PROT_LEN{$ur100}	= $len;
#$ODD_AA{$ur100}	= $odd;
#$PID_UR100{$prot}	= $id;
#$PID_UR90{$prot}	= $id;
#$TCDB_CPDS{$tcdb}	= $chebs;
#$UR100_TCDB{$ur100}	= $tcdb;
##############################################################
##############################################################






##############################################################
######	    INPUT UMRAD COMPOUND AND REACTION DATA     #######
$time=localtime;
print "INPUT REACTION INFO $time\n";
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
}	}	}

$time=localtime;
print "INPUT UNIPROT REACTIONS $time\n";
while(<INUPR>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
	$pid=shift(@stuff);
	foreach my $x (@stuff){ if($x=~/\w/){$UPID_RXNS{$pid}{$x}=1;} }
}

$time=localtime;
print "INPUT KEGG REACTIONS $time\n";
while(<INKEG>){
	if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $gene, my $krxn)=split("\t",$_);
	if($gene!~/\w/ || $krxn!~/\w/){next;}
        $KGEN_KRXN{$gene}=$krxn;
}

$time=localtime;
print "INPUT EC REACTIONS $time\n";
while(<INECS>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $ec, my $rxns)=split("\t",$_);
	if($ec !~/\w/ || $rxns !~/\w/){next;}
	$EC_RXNS{$ec}=$rxns;
}

$time=localtime;
print "INPUT MONOMERS $time\n";
while(<INMON>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $mono, my $rxns)=split("\t",$_);
	$MONO_RXNS{$mono}=$rxns;
}

$time=localtime;
print "INPUT TAXONOMY $time\n";
while(<INTAX>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
	$tid=shift(@stuff);
	$lin=join(";",@stuff);
	$PHY{$tid}=$lin;
}
close(INTAX);

#$LEFT_CPDS{$x}=$lc;
#$RITE_CPDS{$x}=$rc;
#$RXN_ALTS{$x}{$y}=1;
#$UPID_RXNS{$pid}{$x}=1;
#$KGEN_KRXN{$gene}=$krxn;
#$EC_RXNS{$ec}=$rxns;
#$MONO_RXNS{$mono}=$rxns;
#$PHY{$tid}=$lin;
#$GENUS{$stuff[5]}=lin0..5
#$SPECIES{$stuff[6]}=lin0..6
##############################################################
##############################################################




 

##############################################################
##### 		INPUT UNIPROT DOWNLOAD			######
$on=0;
$time = localtime;
print "INPUT UNIPROT $time\n";
print OUTPINT "PID\tName\tLength\tUR100\tUR90\tTaxonId\tLineage\t";
print OUTPINT "SigPep\tTMS\tDNA\tMetal\tTCDB\tLoc\tCOG\tPfam\tTigr\tGene_Ont\tInterPro\tECs\tkegg\trhea\tbiocyc\tLeftCPDs\tRightCPDs\tTransCPDs\n";
while(<INFO>){
	if($_=~/^(entry|PID)/i){next;}
	if($_!~/\w/){next;}
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


	#GET TAXONOMY
	$tid=$stuff[3]; 
	$lin=$PHY{$tid};
	if(!exists($PHY{$tid})){ print LOG "missing tid\t$tid\n"; next;}
 

	#COMPILE/CLEAN UNIPROT ANNOTATIONS
	@FIN=();

	$sig=''; if($stuff[7]=~/(\d+)\.\.(\d+)/){$sig = "SIGNAL:".$1."..".$2;}

	$trc='';  $trcpd=''; @TCDB=(); @TCPD=();
	if(exists($UR100_TCDB{$ur100})){$stuff[9].=";".$UR100_TCDB{$ur100};}
	if($stuff[9]=~/([\w\.]+)/){
		@TX=split(";",$stuff[9]);
		@TCPD=split(";",$stuff[24]);
		foreach my $trc (@TX){
			if($trc!~/\d/){next;}
		 	if($trc !~ /TCDB/){$trc="TCDB:".$trc;}
			push(@TCDB,$trc); 
			if($TCDB_CPDS{$trc}=~/\w/){push(@TCPD,$TCDB_CPDS{$trc});}
	}	}
	%seen=(); @TCDB = grep{ !$seen{$_}++ } @TCDB;
	%seen=(); @TCPD = grep{ !$seen{$_}++ } @TCPD;
	$trc=join(";",@TCDB);
	$trcpd=join(";",@TCPD);
	
	$tm=''; @TMS=();
	if($stuff[8]=~/TRANSMEM|TMS/){
		@TMS = ($stuff[8]=~/(\d+\.\.\d+)/g); 
		$tm  = join(";",@TMS); 
		$tm="TMS:".$tm;
	}

	$cog=''; @COG=();	if($stuff[10]=~/[CK]OG\d+/){@COG  = ($stuff[10]=~/\b([CK]OG\d{4})\b/g);	$cog = join(";",@COG);}

	$pfa=''; @PFAM=();	if($stuff[11]=~/(PF\d+)/){  @PFAM = ($stuff[11]=~/(PF\d+)/g); 		$pfa = join(";",@PFAM);}

	$tig=''; @TIGR=();	if($stuff[12]=~/(TIGR\d+)/){@TIGR = ($stuff[12]=~/(TIGR\d+)/g); 	$tig = join(";",@TIGR);}

	$gos=''; @GOS=();	if($stuff[13]=~/GO.\d+/){   @GOS  = ($stuff[13]=~/(GO.\d+)/g); 		$gos = join(";",@GOS);}

	$ipr=''; @IPR=();	if($stuff[14]=~/(IPR\d+)/){ @IPR  = ($stuff[14]=~/(IPR\d+)/g);		$ipr = join(";",@IPR);}

	$ecs=''; @ECS=();	if($stuff[15]=~/[\d\.\-]+/){@ECS  = ($stuff[15]=~/([\d\.\-]+)/g); 	$ecs = join(";",@ECS);}

	$dna=''; if($stuff[17]=~/(\d+\.\.\d+).*?NOTE\=\"([^\"]+)/){ $cln=CleanNames($2); $dna = "DNA:".$1."|".$cln;}

	$met=''; @METS=();
	if($stuff[18]=~/NOTE/){
		@METS = ($stuff[18]=~/NOTE\W+(\w\S+?)[\"\;]+/g ); 
		$met = join(";", @METS);
	}

	$loc=''; @LOCS=(); @LCG=();
	if($stuff[19]=~/LOCATION/){ 
		$stuff[19]=~s/NOTE\=.*//g;
		if($stuff[19]=~/\{/){$stuff[19]=~s/\,/\{/g; @LOCS = ($stuff[19]=~/([^\.\:\;]+)\s*\{/g); }
		else{ @LOCS = ($stuff[19]=~/([^\.\,\:\;]+)\s*[\;\,\.]/g); }
	}
	for my $i (0..$#LOCS){$LOCS[$i]=CleanNames($LOCS[$i]);}
	while($LOCS[0]=~/./){ 
		$cl = shift(@LOCS); 
		if($cl=~/\d/){next;} 
		   if(grep /$cl/, @LOCS){}
		elsif(grep /$cl/, @LCG){}
		else{push(@LCG,$cl);}
	}
	$loc = join(";", @LCG);


	### COMPILE REACTION DATABASE INFO ###
	#uprxns
	if(keys %{$UPID_RXNS{$pid}}>0){ foreach my $rx (keys %{$UPID_RXNS{$pid}}){if($rx=~/\w/){$RXNS{$rx}++;}}}
	%RXNS=();

	#ecrxns
	@ECS=();  if($stuff[15]=~/[\d\.\-]+/){	@ECS  = ($stuff[15]=~/([\d\.\-]+)/g);} 
	foreach my $ec (@ECS){  
		if(  $EC_RXNS{$ec}=~/\w/){ @RXNS=split(";",$EC_RXNS{$ec});   
			foreach my $rx (@RXNS){ 
				if($rx=~/\w/){$RXNS{$rx}++; 
	}	}	}	}

	#monomers
	@MON=();  if($stuff[16]=~/MONOMER/){	@MON  = ($stuff[16] =~ /([\-\w]*MONOMER[\-\w]*)/g);} 
	foreach my $mn (@MON){	
		if($MONO_RXNS{$mn}=~/\w/){ @RXNS=split(";",$MONO_RXNS{$mn}); 
			foreach my $rx (@RXNS){ 
				if($rx=~/\w/){$RXNS{$rx}++;}	
	}	}	}

	#keggrxns
	@KEGS=(); if($stuff[20]=~/\w+/){	@KEGS = ($stuff[20]=~/([^\;]+)/g);}
	foreach my $kg (@KEGS){ 
		if($KGEN_KRXN{$kg}=~/\w/){ @RXNS=split(";",$KGEN_KRXN{$kg}); 
			foreach my $rx (@RXNS){ 
				if($rx=~/\w/){$RXNS{$rx}++;
	}	}	}	}

	#rhearxns
	@RHEA=(); if($stuff[21]=~/\w+/){	@RHEA = ($stuff[21]=~/(\d+)/g);}
	foreach my $rh (@RHEA){ if($rh=~/\w/){$RXNS{$rh}++; }}

	#GET MISSING RXNS USING %RXN_ALTS
	foreach my $rx (keys %RXNS){ 
		if($rx=~/(^\d{4,}$|^R\d+$|RXN)/){foreach my $alt (keys %{$RXN_ALTS{$rx}}){if($alt=~/\w/){$RXNS{$alt}++; }}}
		else{delete($RXNS{$rx});}
	}

	#GET TOP OF EACH RXN TYPE
	$rm=0; $km=0; $bm=0;
	@RHE=(); @KEG=(); @BIO=(); %ALTS=();
	foreach my $rxn (sort{$RXNS{$b}<=>$RXNS{$a}} keys %RXNS){
		if($rxn !~/\w/){next;}
		   if($rxn=~/^\d{4,}$/){if($RXNS{$rxn}>$rm){$rm=$RXNS{$rxn};} if($RHE[3]!~/\w/ && $RXNS{$rxn} > $rm/2){ push(@RHE,$rxn); $ALTS{$rxn}=1; }}
		elsif($rxn=~/^R\d+$/){ 	if($RXNS{$rxn}>$km){$km=$RXNS{$rxn};} if($KEG[3]!~/\w/ && $RXNS{$rxn} > $km/2){ push(@KEG,$rxn); $ALTS{$rxn}=1; }}
		elsif($rxn=~/RXN/){ 	if($RXNS{$rxn}>$bm){$bm=$RXNS{$rxn};} if($BIO[3]!~/\w/ && $RXNS{$rxn} > $bm/2){ push(@BIO,$rxn); $ALTS{$rxn}=1; }}
		else{next;}
	}
	@KEG=nsort(@KEG);
	@RHE=nsort(@RHE);
	@BIO=nsort(@BIO);
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
		if($cpd=~/CHEBI/){ 	push(@CL,$cpd);}
		elsif($cpd=~/^C\d+$/){ 	push(@KL,$cpd);}
		else{			push(@BL,$cpd);}
	}
	@LEFTS=(); $lcpd='';
	$cl=join(";",@CL); if($cl=~/\w/){push(@LEFTS,$cl);}
	$kl=join(";",@KL); if($kl=~/\w/){push(@LEFTS,$kl);}
	$bl=join(";",@BL); if($bl=~/\w/){push(@LEFTS,$bl);}
	$lcpd=join(";", @LEFTS);
	
	@CR=(); @KR=(); @BR=();
	foreach my $cpd (sort(keys %RITES)){
		if($cpd !~/\w/){next;} 
		if($cpd=~/CHEBI/){ 	push(@CR,$cpd);}
		elsif($cpd=~/^C\d+$/){ 	push(@KR,$cpd);}
		else{			push(@BR,$cpd);}
	}
	@RITES=(); $rcpd='';
	$cr=join(";",@CR); if($cr=~/\w/){push(@RITES,$cr);}
	$kr=join(";",@KR); if($kr=~/\w/){push(@RITES,$kr);}
	$br=join(";",@BR); if($br=~/\w/){push(@RITES,$br);}
	$rcpd=join(";",@RITES);


	# FIX NAMES
	@NAMES=split('\(',$stuff[1]);
	@GN=();
	foreach my $n (@NAMES){ $n=CleanNames($n); if($n!~/[A-Z]{4}/){next;} push(@GN,$n);}
	@NAMES=();
	while($GN[0]=~/\w/){ $n=shift(@GN); if(grep /$n/i, @NAMES || grep /$n/i, @GN){} else{ push(@NAMES,$n); }}
	@NAMES=nsort(@NAMES);
	$name=join(";",@NAMES);
	$odd=$ODD_AA{$ur100};
	if($odd=~/[A-Z]/){$name.=";".$odd;}


	#COMBINE AND CLEAN ALL FUNCS
	@FIN=();
	push(@FIN, $pid);
	push(@FIN, $name);
	push(@FIN, $plen);
	push(@FIN, $ur100);
	push(@FIN, $ur90);
	push(@FIN, $tid);
	push(@FIN, $lin);
	push(@FIN, $sig);
	push(@FIN, $tm);
	push(@FIN, $dna);
	push(@FIN, $met);
	push(@FIN, $trc);
	push(@FIN, $loc);
	push(@FIN, $cog);
	push(@FIN, $pfa);
	push(@FIN, $tig);
	push(@FIN, $gos);
	push(@FIN, $ipr);
	push(@FIN, $ecs);
	push(@FIN, $kegg);
	push(@FIN, $rhea);
	push(@FIN, $bioc);
	push(@FIN, $lcpd);
	push(@FIN, $rcpd);
	push(@FIN, $trcpd);


	#IF PROTEIN HAD A POOR NAME
	if($name !~/\w/){
		@NAMES=(); 
		@GN=();	
		for my $i (11,13..20){
			@FNX=split(";",$FIN[$i]);
			foreach my $id (@FNX){ 
				if($FUNC_NAME{$id}=~/\w/){ $n=CleanNames($FUNC_NAME{$id}); push(@GN,$n); last;}
		}	}
		while($GN[0]=~/\w/){ $n=shift(@GN); if(grep /$n/i, @NAMES || grep /$n/i, @GN){} else{ push(@NAMES,$n); }}
		$name=join(";",@NAMES);
		if($odd=~/[A-Z]/){$name.=";".$odd;}
		$FIN[1]=$name;
	}


	#GET 4 func combos - to fill in function blanks
	@FF=();
        @FF=@FIN[11,13..20];
        for my $q (0..$#FF-3){
                if($FF[$q]!~/\w/){next;}
                $j=$q+1;
                for my $r ($j..$#FF-2){
                        if($FF[$r]!~/\w/){next;}
                        $k=$r+1;
                        for my $s ($k..$#FF-1){
                                if($FF[$s]!~/\w/){next;}
                                $l=$s+1;
                                for my $t ($l..$#FF){
                                        if($FF[$t]!~/\w/){next;}
                                        $combo=$FF[$q]."&".$FF[$r]."&".$FF[$s]."&".$FF[$t];
                                        for my $x (0..8){
                                                if($x!=$q && $x!=$r && $x!=$s && $x!=$t){
                                                        @FS=split(";",$FF[$x]);
                                                        foreach my $id (@FS){ $FUNC_AC{$combo}{$id}++; }
                                        }       }
                                        if(keys %{$FUNC_AC{$combo}} < 1){next;} #if no additional functions, skip
                                        $FUNC_CC{$combo}++;
                                        if($FIN[10]=~/\w/){ @FS=split(";",$FIN[10]); foreach my $id (@FS){ $FUNC_MET{$combo}{$id}++; }}
                                        if($FIN[12]=~/\w/){ @FS=split(";",$FIN[12]); foreach my $id (@FS){ $FUNC_LOC{$combo}{$id}++; }}
                                        if($FIN[21]=~/\w/){ @FS=split(";",$FIN[21]); foreach my $id (@FS){ $FUNC_AC{$combo}{$id}++;  }}
                                        if($FIN[22]=~/\w/){ @FS=split(";",$FIN[22]); foreach my $id (@FS){ $FUNC_LEFT{$combo}{$id}++;}}
                                        if($FIN[23]=~/\w/){ @FS=split(";",$FIN[23]); foreach my $id (@FS){ $FUNC_RITE{$combo}{$id}++;}}
        }       }       }       }


	#PREP FUNC OUTPUT
        for my $i (11,13..20){
                @NF=split(";",$FIN[$i]);
                %seen=(); @NF = grep{ !$seen{$_}++ } @NF;
                @KEEP=(); foreach my $x (@NF){if($x=~/\w/){push(@KEEP,$x);}}
                @KEEP=nsort(@KEEP);
                $FIN[$i]=join(";",@KEEP);
        }
	$out=join("\t",@FIN);
	print OUTPINT "$out\n";

	@OUTFUNCS=@FIN[7..24];
	$outfuncs=join("\t",@OUTFUNCS);
	if($on%100000==0){$time=localtime; print "on $on time $time funcs $outfuncs\n";} $on++;

	delete($PID_UR100{$pid});
	delete($PID_UR90{$pid});
	delete($UPID_RXNS{$pid});
}
close(OUTPINT);
print "final on $on uniprots\n";
undef(%PID_UR100);
undef(%PID_UR90);
undef(%UPID_RXNS);
undef(%PROT_LEN);
undef(%ODD_AA);
undef(%UR100_TCDB);
undef(%LEFT_CPDS);
undef(%RITE_CPDS);
undef(%RXN_ALTS);
undef(%KGEN_KRXN);
undef(%EC_RXNS);
undef(%MONO_RXNS);







### OUTPUT FUNC COMBOS ###
##########################

# CIRCULATE THROUGH THE FUNCTIONS
$kc1=keys %FUNC_AC;
$kc2=keys %FUNC_LOC;
$kc3=keys %FUNC_MET;
$kc4=keys %FUNC_RITE;
$kc5=keys %FUNC_LEFT;
print "combofuncs b4screen ac $kc1 loc $kc2 met $kc3 rite $kc4 left $kc5\n";
$on=0;
foreach my $combo (keys %FUNC_CC){
        if($on%10000==0){print "on $on done cc $FUNC_CC{$combo} combo $combo\n";}
        if($FUNC_CC{$combo}<5 || keys %{$FUNC_AC{$combo}}<1){
                delete($FUNC_CC{$combo});
                delete($FUNC_AC{$combo});
                next;
        }

        #MUST be in at least 3 genes to keep the alt id
        foreach my $id (keys %{$FUNC_MET{$combo}}){ if($FUNC_MET{$combo}{$id} < 3 ){delete($FUNC_MET{$combo}{$id});}  else{print OUTCOMBO "$combo\tMET\t$id\t$FUNC_MET{$combo}{$id}\n";}}
        foreach my $id (keys %{$FUNC_AC{$combo}}){  if($FUNC_AC{$combo}{$id}  < 3 ){delete($FUNC_AC{$combo}{$id});}   else{print OUTCOMBO "$combo\tAC\t$id\t$FUNC_AC{$combo}{$id}\n";}}
        foreach my $id (keys %{$FUNC_LOC{$combo}}){ if($FUNC_LOC{$combo}{$id} < 3 ){delete($FUNC_LOC{$combo}{$id});}  else{print OUTCOMBO "$combo\tLOC\t$id\t$FUNC_LOC{$combo}{$id}\n";}}
        foreach my $id (keys %{$FUNC_LEFT{$combo}}){if($FUNC_LEFT{$combo}{$id}< 3 ){delete($FUNC_LEFT{$combo}{$id});} else{print OUTCOMBO "$combo\tLEFT\t$id\t$FUNC_LEFT{$combo}{$id}\n";}}
        foreach my $id (keys %{$FUNC_RITE{$combo}}){if($FUNC_RITE{$combo}{$id}< 3 ){delete($FUNC_RITE{$combo}{$id});} else{print OUTCOMBO "$combo\tRITE\t$id\t$FUNC_RITE{$combo}{$id}\n";}}
        if(keys %{$FUNC_MET{$combo}}<1){delete($FUNC_MET{$combo});}
        if(keys %{$FUNC_AC{$combo}}<1){delete($FUNC_AC{$combo});}
        if(keys %{$FUNC_LOC{$combo}}<1){delete($FUNC_LOC{$combo});}
        if(keys %{$FUNC_LEFT{$combo}}<1){delete($FUNC_LEFT{$combo});}
        if(keys %{$FUNC_RITE{$combo}}<1){delete($FUNC_RITE{$combo});}
        $on++;
}
undef(%FUNC_CC);
$kc1=keys %FUNC_AC;
$kc2=keys %FUNC_LOC;
$kc3=keys %FUNC_MET;
$kc4=keys %FUNC_RITE;
$kc5=keys %FUNC_LEFT;
print "combofuncs afterscrn ac $kc1 loc $kc2 met $kc3 rite $kc4 left $kc5\n";


 





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


































__END__

# CIRCULATE THROUGH THE FUNCTIONS
$kc1=keys %FUNC_AC;
$kc2=keys %FUNC_LOC;
$kc3=keys %FUNC_MET;
$kc4=keys %FUNC_RITE;
$kc5=keys %FUNC_LEFT;
print "combofuncs b4screen ac $kc1 loc $kc2 met $kc3 rite $kc4 left $kc5\n";
$on=0;
foreach my $combo (keys %FUNC_CC){
        if($on%10000==0){print "on $on done cc $FUNC_CC{$combo} combo $combo\n";}
        if($FUNC_CC{$combo}<5){
                delete($FUNC_CC{$combo});
                delete($FUNC_AC{$combo});
                delete($FUNC_LOC{$combo});
                delete($FUNC_MET{$combo});
                delete($FUNC_RITE{$combo});
                delete($FUNC_LEFT{$combo});
                next;
        }

        #MUST be in at least 3 genes to keep the alt id
        foreach my $id (keys %{$FUNC_AC{$combo}}){  if($FUNC_AC{$combo}{$id}   < 3 ){delete($FUNC_AC{$combo}{$id});}}
        foreach my $id (keys %{$FUNC_LOC{$combo}}){ if($FUNC_LOC{$combo}{$id}  < 3 ){delete($FUNC_LOC{$combo}{$id});}}
        foreach my $id (keys %{$FUNC_MET{$combo}}){ if($FUNC_MET{$combo}{$id}  < 3 ){delete($FUNC_MET{$combo}{$id});}}
        foreach my $id (keys %{$FUNC_RITE{$combo}}){if($FUNC_RITE{$combo}{$id} < 3 ){delete($FUNC_RITE{$combo}{$id});}}
        foreach my $id (keys %{$FUNC_LEFT{$combo}}){if($FUNC_LEFT{$combo}{$id} < 3 ){delete($FUNC_LEFT{$combo}{$id});}}
        $on++;
}
$kc1=keys %FUNC_AC;
$kc2=keys %FUNC_LOC;
$kc3=keys %FUNC_MET;
$kc4=keys %FUNC_RITE;
$kc5=keys %FUNC_LEFT;
print "combofuncs afterscrn ac $kc1 loc $kc2 met $kc3 rite $kc4 left $kc5\n";



$on=0;
open(OUTFIN, ">", "UNIPROT_INFO_".$version."finX.txt")||die;
print OUTFIN "PID\tName\tLength\tUR100\tUR90\tTaxonId\tLineage\t";
print OUTFIN "SigPep\tTMS\tDNA\tMetal\tTCDB\tLoc\tCOG\tPfam\tTigr\tGene_Ont\tInterPro\tECs\tkegg\trhea\tbiocyc\tLeftCPDs\tRightCPDs\tTransCPDs\n";
open(OUTPINT, "UNIPROT_INFO_".$version."intC.txt")||die;
while(<OUTPINT>){
        if($_=~/^PID/){next;}
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
        $begin = $_;

        #CHECK FOR TCDB ANNOTATION BY TMS
        @TMS = (); %POT=();
        if($stuff[8]=~/(\d+\.\.\d+)/){
                @TMS = ($stuff[8]=~/(\d+\.\.\d+)/g);
                $tc=@TMS;
                if($tc>=5){
                        foreach my $tcdb (keys %{$TMS_MAX_ST{$tc}}){
                                $gt=1;
                                for my $i (0..$#TMS){
                                        $TMS[$i]=~/(\d+)\.\.(\d+)/;
                                        $start=$1; $end=$2;
                                        if($start >= $TMS_MIN_ST{$tc}{$tcdb}{$i}-1 && $start <=$TMS_MAX_ST{$tc}{$tcdb}{$i}+1){$gs=1;}else{$gt=0;}
                                        if(  $end >= $TMS_MIN_EN{$tc}{$tcdb}{$i}-1 &&   $end <=$TMS_MAX_EN{$tc}{$tcdb}{$i}+1){$gs=1;}else{$gt=0;}
                                }
                                if($gt==1){$POT{$tcdb}=1;}
                                else{next;}
        }       }       }

        #GET 4 FUNC COMBOS _ TO FILL MISSING FUNCS
        @FF=();
        if($stuff[11]=~/\d/){push(@FF,$stuff[11]);}
        if($stuff[13]=~/\d/){push(@FF,$stuff[13]);}
        if($stuff[14]=~/\d/){push(@FF,$stuff[14]);}
        if($stuff[15]=~/\d/){push(@FF,$stuff[15]);}
        if($stuff[16]=~/\d/){push(@FF,$stuff[16]);}
        if($stuff[17]=~/\d/){push(@FF,$stuff[17]);}
        if($stuff[18]=~/\d/){push(@FF,$stuff[18]);}
        if($stuff[20]=~/\d/){push(@FF,$stuff[20]);}

        #IF HAS 4+ FUNCS, CHECK COMBOS, SCREEN POTENTIAL MISSING FUNCS
        %ALTF=(); %ALTM=(); %ALTR=(); %ALTE=(); %ALTL=();$ff=@FF;
        if($ff>3){ #has 4+ identifying functions
                for my $i (0..$#FF-3){ #loop thru funcs
                        $combo=join("&",@FF[$i..$i+3]);
                        if(exists($FUNC_AC{$combo})){  foreach my $id (keys %{$FUNC_AC{$combo}}){  $ALTF{$id} += $FUNC_AC{$combo}{$id};   }}
                        if(exists($FUNC_LOC{$combo})){ foreach my $id (keys %{$FUNC_LOC{$combo}}){ $ALTL{$id} += $FUNC_LOC{$combo}{$id};  }}
                        if(exists($FUNC_MET{$combo})){ foreach my $id (keys %{$FUNC_MET{$combo}}){ $ALTM{$id} += $FUNC_MET{$combo}{$id};  }}
                        if(exists($FUNC_RITE{$combo})){foreach my $id (keys %{$FUNC_RITE{$combo}}){$ALTR{$id} += $FUNC_RITE{$combo}{$id}; }}
                        if(exists($FUNC_LEFT{$combo})){foreach my $id (keys %{$FUNC_LEFT{$combo}}){$ALTE{$id} += $FUNC_LEFT{$combo}{$id}; }}
                }

                #get best funcs of each type
                $cm=0; $pm=0; $tm=0; $gm=0; $im=0; $em=0; $km=0; $hm=0; $xm=0; $nm=0;
                foreach my $id (sort{$ALTF{$b}<=>$ALTF{$a}} keys %ALTF){
                           if($id=~/^COG\d+$/       ){  if($ALTF{$id}>=$cm){$cm=$ALTF{$id};} if($ALTF{$id}<$cm-1){delete($ALTF{$id});}}
                        elsif($id=~/^PF\d+$/        ){  if($ALTF{$id}>=$pm){$pm=$ALTF{$id};} if($ALTF{$id}<$pm-1){delete($ALTF{$id});}}
                        elsif($id=~/^TIGR\d+$/      ){  if($ALTF{$id}>=$tm){$tm=$ALTF{$id};} if($ALTF{$id}<$tm-1){delete($ALTF{$id});}}
                        elsif($id=~/^GO.\d+$/       ){  if($ALTF{$id}>=$gm){$gm=$ALTF{$id};} if($ALTF{$id}<$gm-1){delete($ALTF{$id});}}
                        elsif($id=~/^IPR\d+$/       ){  if($ALTF{$id}>=$im){$im=$ALTF{$id};} if($ALTF{$id}<$im-1){delete($ALTF{$id});}}
                        elsif($id=~/^\d+$/          ){  if($ALTF{$id}>=$hm){$hm=$ALTF{$id};} if($ALTF{$id}<$hm-1){delete($ALTF{$id});}}
                        elsif($id=~/^[\-\.\d]+$/    ){  if($ALTF{$id}>=$em){$em=$ALTF{$id};} if($ALTF{$id}<$em-1){delete($ALTF{$id});}}
                        elsif($id=~/^R\d+$/         ){  if($ALTF{$id}>=$km){$km=$ALTF{$id};} if($ALTF{$id}<$km-1){delete($ALTF{$id});}}
                        elsif($id=~/RXN/            ){  if($ALTF{$id}>=$xm){$xm=$ALTF{$id};} if($ALTF{$id}<$xm-1){delete($ALTF{$id});}}
                        elsif($id=~/TCDB.[\-\.\dA-Z]+$/ ){ if($ALTF{$id}>=$nm){$nm=$ALTF{$id};} if($ALTF{$id}<$nm-1){delete($ALTF{$id});}}
                        else{delete($ALTF{$id});}
                }
                $xm=0; foreach my $id (sort{$ALTL{$b}<=>$ALTL{$a}} keys %ALTL){ if($ALTL{$id}>=$xm){$xm=$ALTL{$id};} if($ALTL{$id}<$xm-1){delete($ALTL{$id});}}
                $xm=0; foreach my $id (sort{$ALTM{$b}<=>$ALTM{$a}} keys %ALTM){ if($ALTM{$id}>=$xm){$xm=$ALTM{$id};} if($ALTM{$id}<$xm-1){delete($ALTM{$id});}}
                $xm=0; foreach my $id (sort{$ALTR{$b}<=>$ALTR{$a}} keys %ALTR){ if($ALTR{$id}>=$xm){$xm=$ALTR{$id};} if($ALTR{$id}<$xm-1){delete($ALTR{$id});}}
                $xm=0; foreach my $id (sort{$ALTE{$b}<=>$ALTE{$a}} keys %ALTE){ if($ALTE{$id}>=$xm){$xm=$ALTE{$id};} if($ALTE{$id}<$xm-1){delete($ALTE{$id});}}


                #incorporate missing functions
                foreach my $id (keys %ALTF){
                           if($id=~/^COG\d+$/       ){$stuff[13].=";".$id;}
                        elsif($id=~/^PF\d+$/        ){$stuff[14].=";".$id;}
                        elsif($id=~/^TIGR\d+$/      ){$stuff[15].=";".$id;}
                        elsif($id=~/^GO.\d+$/       ){$stuff[16].=";".$id;}
                        elsif($id=~/^IPR\d+$/       ){$stuff[17].=";".$id;}
                        elsif($id=~/^\d+$/          ){$stuff[20].=";".$id;}
                        elsif($id=~/^[\-\.\d]+$/    ){$stuff[18].=";".$id;}
                        elsif($id=~/^R\d+$/         ){$stuff[19].=";".$id;}
                        elsif($id=~/RXN/            ){$stuff[21].=";".$id;}
                        elsif($id=~/TCDB.[\-\.\dA-Z]+$/ && exists($POT{$id})){
                                $stuff[11].=";".$id; $stuff[24].=";".$TRCH{$id};
                        }
                        else{ next;}
                }
                foreach my $id (keys %ALTL){$stuff[12].=";".$id;}
                foreach my $id (keys %ALTM){$stuff[10].=";".$id;}
                foreach my $id (keys %ALTR){$stuff[23].=";".$id;}
                foreach my $id (keys %ALTE){$stuff[22].=";".$id;}
        }

        #sort and clean new ids
        for my $i (10..24){
                @NF=split(";",$stuff[$i]);
                %seen=(); @NF = grep{ !$seen{$_}++ } @NF;
                @KEEP=(); foreach my $x (@NF){if($x=~/\w/){push(@KEEP,$x);}}
                @KEEP=nsort(@KEEP);
                $stuff[$i]=join(";",@KEEP);
        }

        $out=join("\t",@stuff);
        if($begin ne $out){ print "beg $begin\nend $out\n";}
        print OUTFIN "$out\n";
        $on++;
}


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

