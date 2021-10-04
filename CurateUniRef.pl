#use warnings;


#GET UPDATE INFO
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$month = $1; $year = $2;
$version=$1."_".$2;

#CHECK INPUTS
qx{wget -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_all.xml.gz};
$inup  = 'UNIPROT_INFO_'.$version.'.txt.gz';
if(!-s $inup ){ print "missing $inup - please place into this folder and restart\n";  $fail=1;}
$intax = 'TAXONOMY_DB_'.$version.'.txt';
if(!-s $intax ){print "missing $intax - please place into this folder and restart\n"; $fail=1;}
$infnam= 'Function_Names_'.$version.'.txt';
if(!-s $infnam){print "missing $infnam - please run Get_Function_Names.pl\n"; $fail=1;}
if(!-s "uniref100.fasta.gz"){
        print "missing uniref100.fasta.gz downloading...\n";
        qx{wget -N https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz};
}
if($fail==1){ print "missing needed input files.\n";
        print "see my GitHub repositories: https://github.com/TealFurnholm/ for further guidance\n"; die;}
open(INTAX, $intax)||die;
open(INNAM, $infnam)||die;
open(INFO, "gunzip -c $inup |")||die;
open(DEBUG, ">", "curate_uniref_".$month.'_'.$year.".debug")||die;
open(UNIPARC, "gunzip -c uniparc_all.xml.gz |")||die;



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
#$PHY{$tid}=$lin;



$time=localtime;
print "INPUT FUNCNAMES $time\n";
while(<INNAM>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $func, my $name)=split("\t",$_);
        $FUNC_NAME{$func}=$name;
}
#$PHY{$tid}=$lin;
#$FUNC_NAME{$func}=$name;


#GET PROT LENGTH, NAME AND ODD AA DISTs
print "GET ODD AAS, NAMES, and LENGTHS\n";
if(-s "UNIREF100_AANSTAT".$month.'_'.$year.".txt"){
        print "USING EXISTING AA ANNOTATIONS\n";
        open(INODD, "UNIREF100_AANSTAT".$month.'_'.$year.".txt")||die;
        while(<INODD>){
                if($_!~/\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                (my $ur100, my $len, my $odd, my $name, my $tid)=split("\t",$_);
                $ODD_AA{$ur100}=$odd;
                $UR100_TID{$ur100}{$tid}=1;
                $UR100_NAME{$ur100}=$name;
                $PROT_LEN{$ur100}=$len;
        }
}
else{   print "CREATING AA ANNOTATIONS FILE\n";
        open(OUTODD,">","UNIREF100_AANSTAT".$month.'_'.$year.".txt")||die;
        open(INURFA, "gunzip -c uniref100.fasta.gz |")||die;
        while(<INURFA>){
                if($_!~/\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                if($_=~ /(UNIREF100\S+)\s(.*)N\=\d.*TAXID\=(\d+)/){
                        $line=$_;
                        if($ur100=~/\w/ && $SEQ[0]=~/[A-Z]/){
                                $seq=join("",@SEQ);
                                $seq=~s/[^A-Z]+//g;
                                $len = length($seq);
                                $PROT_LEN{$ur100}=$len;
                                for my $i (A..Z){ $aa='';
                                        my $count = () = $seq =~ /$i/g; $per=$count*100/$len;
                                        $per=~s/\..*//; if($per >= 15){$aa=$i.":".$per; push(@AA,$aa);}
                                }
                                $odd=join("|",@AA);
                                if($odd=~/[A-Z]/){ $ODD_AA{$ur100}=$odd; }
                                print OUTODD "$ur100\t$len\t$odd\t$name\t$tid\n";
                                if($on%1000000==0){$time=localtime; print "on $on time $time ur $ur100 len $len odd $odd tid $tid name $name seq $seq\n";} $on++;
                        }
                        $ur100=''; $name=''; $tid=''; @SEQ=(); @AA=();
                        $line=~ /(UNIREF100\S+)\s(.*)N\=\d.*TAXID\=(\d+)/;
                        $ur100=$1; $name=CleanNames($2); $tid=$3;
                        $UR100_TID{$ur100}{$tid}=1;
                        $UR100_NAME{$ur100}=$name;
                }
                if($_=~/^[A-Z]+$/){push(@SEQ,$_);}
}       }
#$PHY{$tid}=$lin;
#$FUNC_NAME{$func}=$name;
#$UR100_TID{$ur100}{$tid}=1;
#$UR100_NAME{$ur100}=$name;
#$PROT_LEN{$ur100}=$len;
#$ODD_AA{$ur100}=$odd;
$phykc=keys %PHY;
$fname=keys %FUNC_NAME;
$urtc=keys %UR100_TID;
$urnm=keys %UR100_NAME;
$prtl=keys %PROT_LEN;
$odda=keys %ODD_AA;
print DEBUG "1phykc $phykc fname $fname urtc $urtc urnm $urnm prtl $prtl odda $odda\n";


#LOAD VARIABLES
$on=0;
$/="\n";
print "INPUT UNIPROT\n";
while(<INFO>){
        if($_!~/\w/){next;}
        $on++;
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
        if($stuff[3]!~/\d/){next;}
        $pid=$stuff[0];
        $ur100=$stuff[3];
        $ur90=$stuff[4];
        if($ur100 !~ /\w/){next;}
        if($UR100_NAME{$ur100} !~/\w/){$UR100_NAME{$ur100}=$stuff[1];}
        if($PROT_LEN{$ur100}!~/\d/){$PROT_LEN{$ur100}=$stuff[2];}
        $UR100_PID{$ur100}{$pid}=1;
        if($ur90 =~ /\w/){
                $UR100_UR90{$ur100} = $ur90;
                $UR90_UR100{$ur90}{$ur100}=1;
        }
        @IDS=(); if($stuff[5] =~/^\d+$/){       @IDS=split(";",$stuff[5]);      foreach my $id (@IDS){  if($id=~/\w/){  $UR100_TID{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[7] =~/\w/){                                                                                  $UR100_SIG{$ur100}{$stuff[7]}++;}
        @IDS=(); if($stuff[8] =~/\w/){                                                                                  $UR100_TMS{$ur100}{$stuff[8]}++;}
        @IDS=(); if($stuff[9] =~/\w/){          @IDS=split(";",$stuff[9]);      foreach my $id (@IDS){  if($id=~/\w/){  $UR100_DNA{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[10]=~/\w/){          @IDS=split(";",$stuff[10]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_MET{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[11]=~/\w/){          @IDS=split(";",$stuff[11]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_TCDB{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[12]=~/\w+/){         @IDS=split(";",$stuff[12]);     foreach my $id (@IDS){  if($id=~/\w/ && $id !~ /NOTE/){ $UR100_LOC{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[13]=~/OG\d+/){       @IDS=split(";",$stuff[13]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_COG{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[14]=~/^PF\d+/){      @IDS=split(";",$stuff[14]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_PFAM{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[15]=~/^TIGR\d+/){    @IDS=split(";",$stuff[15]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_TIGR{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[16]=~/^GO.\d+/){     @IDS=split(";",$stuff[16]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_GO{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[17]=~/^IPR\d+/){     @IDS=split(";",$stuff[17]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_IPR{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[18]=~/\d/){          @IDS=split(";",$stuff[18]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_EC{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[19]=~/^R\d+/){       @IDS=split(";",$stuff[19]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_KEGG{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[20]=~/^\d+$/){       @IDS=split(";",$stuff[20]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_RHEA{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[21]=~/\w+/){         @IDS=split(";",$stuff[21]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_BIOC{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[22]=~/\w+/){         @IDS=split(";",$stuff[22]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_LEFT{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[23]=~/\w+/){         @IDS=split(";",$stuff[23]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_RITE{$ur100}{$id}++;}}}
        @IDS=(); if($stuff[24]=~/\w+/){         @IDS=split(";",$stuff[24]);     foreach my $id (@IDS){  if($id=~/\w/){  $UR100_TRC{$ur100}{$id}++;}}}

        if($on%1000000==0){ print "on $on loading uniprot\n";} $on++;
}
$phykc=keys %PHY;
$fname=keys %FUNC_NAME;
$urtc=keys %UR100_TID;
$urnm=keys %UR100_NAME;
$prtl=keys %PROT_LEN;
$odda=keys %ODD_AA;
print DEBUG "2phykc $phykc fname $fname urtc $urtc urnm $urnm prtl $prtl odda $odda\n";



#GO THROUGH UNIREF 90 - GET BEST FUNCTIONS - DISTRIBUT MISSING TO ITS UNIREF100s - OUTPUT THE UNIREF100 INFO DB
print "OUTPUT UNIREF\n";
open(OUTUR100,">","UNIREF100_INFO_".$version."X.txt")||die;
$first="UNIREF100\tNAME\tLENGTH\tUNIPROT_IDS\tUNIREF90\tTAXON_IDS\tLINEAGE";
$funcs="SIGALPEP\tTMS\tDNA\tMETAL\tTCDB\tLOCATION\tCOG\tPFAM\tTIGR\tGO\tIPR\tEC\tKEGG\tRHEA\tBIOCYC\tREACTANTS\tPRODUCTS\tTRANS_CPD";
print OUTUR100 "$first\t$funcs\n";
$on=0;
foreach my $ur90 (keys %UR90_UR100){
        if($ur90!~/\d/){next;}
        $n=$ur90;
        foreach my $ur100 (keys %{$UR90_UR100{$n}}){
                if($ur100!~/\d/){next;}
                $u=$ur100;
                %NAMES=(); @NAMES=();   @PIDS=(); @TIDS=(); @PHYL=();
                $name=''; $odd='';      $pids=''; $tids=''; $lin=''; $first=''; $funcs='';
                $sig=''; $tms=''; $dna=''; $met=''; $tcdb=''; $loc=''; $cog=''; $pfam=''; $tigr='';
                $go=''; $ipr=''; $ec=''; $kegg=''; $rhea=''; $bioc=''; $left=''; $rite=''; $trc='';

                #COMPILE BASIC INFO
                $name=$UR100_NAME{$u};
                $len=$PROT_LEN{$u};
                foreach my $x (  sort(keys %{$UR100_PID{$ur100}})){if($x=~/\w/){ push(@PIDS,$x); }}
                $pids=join(";",@PIDS);
                foreach my $x (  sort(keys %{$UR100_TID{$ur100}})){if($x=~/\w/){ push(@TIDS,$x); }}
                $tids=join(";",@TIDS);
                foreach my $tid  (@TIDS){     if($PHY{$tid}=~/^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)\;.*\;.*\;.*\;.*\;.*\;.*/){ push(@PHYL,$PHY{$tid}); }} $lin=MakeLCA(@PHYL);
                        if($PHYL[0]!~/\w/){ foreach my $tid  (@TIDS){ if($PHY{$tid}=~/^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/){ push(@PHYL,$PHY{$tid}); }} $lin=MakeLCA(@PHYL);
                }
                if($on%1000000==0){$time=localtime; print DEBUG "on $on ur90 $n ur100 $u len $len pids $pids tids $tids name $name lin $lin\n";}

                #COMPILE FUNCTION INFO
                $mx=0; foreach my $x ( sort{$UR100_SIG{$ur100}{$b}<=>$UR100_SIG{$ur100}{$a}} keys %{$UR100_SIG{$ur100}}){       $sig=$x; last;}
                $mx=0; foreach my $x ( sort{$UR100_TMS{$ur100}{$b}<=>$UR100_TMS{$ur100}{$a}} keys %{$UR100_TMS{$ur100}}){       $tms=$x; last;}
                $mx=0; foreach my $x ( sort{$UR100_DNA{$ur100}{$b}<=>$UR100_DNA{$ur100}{$a}} keys %{$UR100_DNA{$ur100}}){       $dna=$x; last;}
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_MET{$u}}){       $IDS{$id}+=$UR100_MET{$u}{$id};}  delete($UR100_MET{$u});       if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_MET{$ux}}){      $IDS{$id}+=$UR100_MET{$ux}{$id};}}}  foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $met=join(";",@IDS);     foreach my $id (@IDS){ $UR100_MET{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_TCDB{$u}}){      $IDS{$id}+=$UR100_TCDB{$u}{$id};} delete($UR100_TCDB{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_TCDB{$ux}}){     $IDS{$id}+=$UR100_TCDB{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $tcdb=join(";",@IDS);    foreach my $id (@IDS){ $UR100_TCDB{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_LOC{$u}}){       $IDS{$id}+=$UR100_LOC{$u}{$id};}  delete($UR100_LOC{$u});       if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_LOC{$ux}}){      $IDS{$id}+=$UR100_LOC{$ux}{$id};}}}  foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $loc=join(";",@IDS);     foreach my $id (@IDS){ $UR100_LOC{$u}{$id} +=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_COG{$u}}){       $IDS{$id}+=$UR100_COG{$u}{$id};}  delete($UR100_COG{$u});       if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_COG{$ux}}){      $IDS{$id}+=$UR100_COG{$ux}{$id};}}}  foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $cog=join(";",@IDS);     foreach my $id (@IDS){ $UR100_COG{$u}{$id} +=$IDS{$id};
                if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_PFAM{$u}}){      $IDS{$id}+=$UR100_PFAM{$u}{$id};} delete($UR100_PFAM{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_PFAM{$ux}}){     $IDS{$id}+=$UR100_PFAM{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $pfam=join(";",@IDS);    foreach my $id (@IDS){ $UR100_PFAM{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_TIGR{$u}}){      $IDS{$id}+=$UR100_TIGR{$u}{$id};} delete($UR100_TIGR{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_TIGR{$ux}}){     $IDS{$id}+=$UR100_TIGR{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $tigr=join(";",@IDS);    foreach my $id (@IDS){ $UR100_TIGR{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_GO{$u}}){        $IDS{$id}+=$UR100_GO{$u}{$id};}   delete($UR100_GO{$u});        if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_GO{$ux}}){       $IDS{$id}+=$UR100_GO{$ux}{$id};}}}   foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $go=join(";",@IDS);      foreach my $id (@IDS){ $UR100_GO{$u}{$id}  +=$IDS{$id};
                if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_IPR{$u}}){       $IDS{$id}+=$UR100_IPR{$u}{$id};}  delete($UR100_IPR{$u});       if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_IPR{$ux}}){      $IDS{$id}+=$UR100_IPR{$ux}{$id};}}}  foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $ipr=join(";",@IDS);     foreach my $id (@IDS){ $UR100_IPR{$u}{$id} +=$IDS{$id};
                if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_EC{$u}}){        $IDS{$id}+=$UR100_EC{$u}{$id};}   delete($UR100_EC{$u});        if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_EC{$ux}}){       $IDS{$id}+=$UR100_EC{$ux}{$id};}}}   foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $ec=join(";",@IDS);      foreach my $id (@IDS){ $UR100_EC{$u}{$id}  +=$IDS{$id};
                if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_KEGG{$u}}){      $IDS{$id}+=$UR100_KEGG{$u}{$id};} delete($UR100_KEGG{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_KEGG{$ux}}){     $IDS{$id}+=$UR100_KEGG{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $kegg=join(";",@IDS);    foreach my $id (@IDS){ $UR100_KEGG{$u}{$id}+=$IDS{$id};
                if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_RHEA{$u}}){      $IDS{$id}+=$UR100_RHEA{$u}{$id};} delete($UR100_RHEA{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_RHEA{$ux}}){     $IDS{$id}+=$UR100_RHEA{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $rhea=join(";",@IDS);    foreach my $id (@IDS){ $UR100_RHEA{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_BIOC{$u}}){      $IDS{$id}+=$UR100_BIOC{$u}{$id};} delete($UR100_BIOC{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_BIOC{$ux}}){     $IDS{$id}+=$UR100_BIOC{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $bioc=join(";",@IDS);    foreach my $id (@IDS){ $UR100_BIOC{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_LEFT{$u}}){      $IDS{$id}+=$UR100_LEFT{$u}{$id};} delete($UR100_LEFT{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_LEFT{$ux}}){     $IDS{$id}+=$UR100_LEFT{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $left=join(";",@IDS);    foreach my $id (@IDS){ $UR100_LEFT{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_RITE{$u}}){      $IDS{$id}+=$UR100_RITE{$u}{$id};} delete($UR100_RITE{$u});      if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_RITE{$ux}}){     $IDS{$id}+=$UR100_RITE{$ux}{$id};}}} foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $rite=join(";",@IDS);    foreach my $id (@IDS){ $UR100_RITE{$u}{$id}+=$IDS{$id}; }
                %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_TRC{$u}}){       $IDS{$id}+=$UR100_TRC{$u}{$id};}  delete($UR100_TRC{$u});       if(keys %IDS < 1 && $n=~/\w/){ foreach my $ux (keys %{$UR90_UR100{$n}}){
                                         foreach my $id (keys %{$UR100_TRC{$ux}}){      $IDS{$id}+=$UR100_TRC{$ux}{$id};}}}  foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
                $trc=join(";",@IDS);     foreach my $id (@IDS){ $UR100_TRC{$u}{$id}+=$IDS{$id}; }

                #COMPLETE NAME
                $name=$UR100_NAME{$u};
                $len=$PROT_LEN{$u};
                if($name=~/UNCHARACTERIZED/){
                        foreach my $nx (sort{$NAMES{$b}<=>$NAMES{$a}} keys %NAMES){
                                if($nx !~ /UNCHARACTERIZED/ && $NAMES[2]!~/\w/){push(@NAMES,$nx);}
                        }
                }
                if($NAMES[0]=~/\w/){ $name=join(";",@NAMES); }
                if($name !~ /\w/){$name=$UR100_NAME{$u};}
                if($ODD_AA{$u}=~/\w/){$name.=";".$ODD_AA{$u};}
                if(keys %NAMES > 4 && $on=~/7$/){$nkc = keys %NAMES; print DEBUG "on $on ur90 $n ur100 $u nkc $nkc odd $ODD_AA{$u} name $name atnames @NAMES\n";}
                if($on%1000000==0){$time=localtime; $nkc=keys %UR90_UR100; $ukc=keys %PROT_LEN; print DEBUG "on $on ur90 $n has $nkc remaining ur100 $u has $ukc remaining len $len pids $pids tids $tids name $name lin $lin\n";}

                $first="$ur100\t$name\t$len\t$pids\t$ur90\t$tids\t$lin";
                $funcs="$sig\t$tms\t$dna\t$met\t$tcdb\t$loc\t$cog\t$pfam\t$tigr\t$go\t$ipr\t$ec\t$kegg\t$rhea\t$bioc\t$left\t$rite\t$trc";
                print OUTUR100  "$first\t$funcs\n";

        }
        $on++;
        foreach my $ur100 (keys %{$UR90_UR100{$n}}){
                delete($UR90_UR100{$n});
                delete($PROT_LEN{$ur100});
                delete($ODD_AA{$ur100});
                delete($UR100_PID{$ur100});
                delete($UR100_TID{$ur100});
                delete($UR100_NAME{$ur100});
                delete($UR100_SIG{$ur100});
                delete($UR100_TMS{$ur100});
                delete($UR100_DNA{$ur100});
                delete($UR100_MET{$ur100});
                delete($UR100_TCDB{$ur100});
                delete($UR100_LOC{$ur100});
                delete($UR100_COG{$ur100});
                delete($UR100_PFAM{$ur100});
                delete($UR100_TIGR{$ur100});
                delete($UR100_GO{$ur100});
                delete($UR100_IPR{$ur100});
                delete($UR100_EC{$ur100});
                delete($UR100_KEGG{$ur100});
                delete($UR100_TRC{$ur100});
                delete($UR100_RITE{$ur100});
                delete($UR100_LEFT{$ur100});
                delete($UR100_BIOC{$ur100});
                delete($UR100_RHEA{$ur100});
}       }


$on=0;
print "INPUT UNIPARC DATA\n";
while(<UNIPARC>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        if($_=~/\s+.ACCESSION.(UPI\w+)/){$inac=1; $ur100="UNIREF100_".$1; %PIDS=(); %TIDS=(); %TIGRS=(); %PFAMS=(); %IPRS=(); @NAMES=(); $length=0; next;}
        if($inac==1){
                if($_=~/UNIPROTKB.*ID\=.(\w+)/){$PIDS{$1}=1;}
                if($_=~/NCBI.TAXONOMY.ID.*VALUE\=.(\d+)/){$TIDS{$1}=1;}
                if($_=~/TIGRFAMS.*ID\=.(TIGR\d+)/){$TIGRS{$1}=1;}
                if($_=~/PFAM.*ID\=.(PF\d+)/){$PFAMS{$1}=1;}
                if($_=~/IPR.NAME.*ID\=.(IPR\d+)/){$IPRS{$1}=1;}
                if($_=~/PROTEIN.NAME.*VALUE\=.([^\"]+)/){ if($NAMES[2]!~/\w/){ $name=CleanNames($1); push(@NAMES,$name); }}
                if($_=~/SEQUENCE.LENGTH\=.(\d+)/){$length=$1;}
        }
        if($_ =~/\<\/ENTRY\>/ && $inac==1){
                $inac=0;
                if(!exists($UR100_TID{$ur100})){next;}
                if($on%1000000==0){$pkc = keys %UR100_PFAM; print "on $on ur100 $ur100 name $name have pkc $pkc pfam ur100s\n";} $on++;
                $ur90='';
                if($UR100_NAME{$ur100}=~/UNKNOWN|HYPOTHETICAL|UNCLASSIFIED|UNCHARACTERIZED/ && $name !~/UNKNOWN|HYPOTHETICAL|UNCLASSIFIED|UNCHARACTERIZED/){$UR100_NAME{$ur100}=$name;}
                if($UR100_NAME{$ur100}!~/\w/){$UR100_NAME{$ur100}=$name;}
                if($length =~ /^\d{2,4}$/ && $PROT_LEN{$ur100} !~ /^\d{2,4}$/){$PROT_LEN{$ur100}=$length;}
                foreach my $id (keys %TIDS){  $UR100_TID{$ur100}{$id}++; }
                foreach my $id (keys %TIGRS){ $UR100_TIGR{$ur100}{$id}++;}
                foreach my $id (keys %PFAMS){ $UR100_PFAM{$ur100}{$id}++;}
                foreach my $id (keys %IPRS){  $UR100_IPR{$ur100}{$id}++; }
                foreach my $id (keys %PIDS){  $UR100_PID{$ur100}{$id}++; }
                foreach my $id (keys %PIDS){
                        if(exists($UR100_UR90{$u})){$ur90=$UR100_UR90{$u}; $UR100_UR90{$ur100}=$ur90; last;}
}       }       }





open(OUTUR100,">>","UNIREF100_INFO_".$version."X.txt")||die;
foreach my $ur100 (keys %PROT_LEN){
        if($ur100!~/\d/){next;}
        $u=$ur100;
        $ur90='';
        if(exists($UR100_UR90{$ur100})){$ur90=$UR100_UR90{$ur100};}
        %NAMES=(); @NAMES=();   @PIDS=(); @TIDS=(); @PHYL=();
        $name=''; $odd='';      $pids=''; $tids=''; $lin='';
        $sig=''; $tms=''; $dna=''; $met=''; $tcdb=''; $loc=''; $cog=''; $pfam=''; $tigr='';
        $go=''; $ipr=''; $ec=''; $kegg=''; $rhea=''; $bioc=''; $left=''; $rite=''; $trc='';

        #COMPILE BASIC INFO
        $name=$UR100_NAME{$ur100};
        $len=$PROT_LEN{$ur100};
        foreach my $x (  sort(keys %{$UR100_PID{$ur100}})){if($x=~/\w/){ push(@PIDS,$x); }}
        $pids=join(";",@PIDS);
        foreach my $x (  sort(keys %{$UR100_TID{$ur100}})){if($x=~/\w/){ push(@TIDS,$x); }}
        $tids=join(";",@TIDS);
        foreach my $tid  (@TIDS){     if($PHY{$tid}=~/^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)\;.*\;.*\;.*\;.*\;.*\;.*/){ push(@PHYL,$PHY{$tid}); }} $lin=MakeLCA(@PHYL);
                if($PHYL[0]!~/\w/){ foreach my $tid  (@TIDS){ if($PHY{$tid}=~/^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/){ push(@PHYL,$PHY{$tid}); }} $lin=MakeLCA(@PHYL);
        }

        #COMPILE FUNCTIONS
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_MET{$u}}){       $IDS{$id}+=$UR100_MET{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $met=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_TCDB{$u}}){      $IDS{$id}+=$UR100_TCDB{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $tcdb=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_LOC{$u}}){       $IDS{$id}+=$UR100_LOC{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $loc=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_COG{$u}}){       $IDS{$id}+=$UR100_COG{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);
        if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}}
        $cog=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_PFAM{$u}}){      $IDS{$id}+=$UR100_PFAM{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $pfam=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_TIGR{$u}}){      $IDS{$id}+=$UR100_TIGR{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $tigr=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_GO{$u}}){        $IDS{$id}+=$UR100_GO{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);
        if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}}
        $go=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_IPR{$u}}){       $IDS{$id}+=$UR100_IPR{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);
        if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}}
        $ipr=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_EC{$u}}){        $IDS{$id}+=$UR100_EC{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);
        if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}}
        $ec=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_KEGG{$u}}){      $IDS{$id}+=$UR100_KEGG{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);
        if($FUNC_NAME{$id}=~/\w/){$NAMES{$FUNC_NAME{$id}}+=$IDS{$id};}}}
        $kegg=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_RHEA{$u}}){      $IDS{$id}+=$UR100_RHEA{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $rhea=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_BIOC{$u}}){      $IDS{$id}+=$UR100_BIOC{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $bioc=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_LEFT{$u}}){      $IDS{$id}+=$UR100_LEFT{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $left=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_RITE{$u}}){      $IDS{$id}+=$UR100_RITE{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $rite=join(";",@IDS);
        %IDS=(); @IDS=(); $mx=0; foreach my $id (keys %{$UR100_TRC{$u}}){       $IDS{$id}+=$UR100_TRC{$u}{$id};}
        foreach my $id (sort{$IDS{$b}<=>$IDS{$a}} keys %IDS){ if($IDS{$id}>$mx){$mx=$IDS{$id};} if($IDS{$id}>=$mx/2){push(@IDS,$id);}}
        $trc=join(";",@IDS);

        #COMPLETE NAME
        if($name=~/UNCHARACTERIZED/){
                foreach my $nx (sort{$NAMES{$b}<=>$NAMES{$a}} keys %NAMES){
                        if($nx !~ /UNCHARACTERIZED/ && $NAMES[2]!~/\w/){push(@NAMES,$nx);}
                }
        }
        if($NAMES[0]=~/\w/){ $name=join(";",@NAMES); }
        if($ODD_AA{$ur100}=~/\w/){$name.=";".$ODD_AA{$ur100};}

        $first="$ur100\t$name\t$len\t$pids\t$ur90\t$tids\t$lin";
        $funcs="$sig\t$tms\t$dna\t$met\t$tcdb\t$loc\t$cog\t$pfam\t$tigr\t$go\t$ipr\t$ec\t$kegg\t$rhea\t$bioc\t$left\t$rite\t$trc";
        print OUTUR100  "$first\t$funcs\n";

}



sub MakeLCA{
        my @ARR;
        %seen=();
        @array1 = @_;
        @array1 = grep { !$seen{$_}++ } @array1;

        #get the kingdoms, JIC lineage is NCA
        %LET=();
        foreach my $lin (@array1){
                if($lin !~ /^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/i){next;}
                $lin =~ /^(.)/; $LET{$1}=1;
        }
        @LET=();
        foreach my $let (sort(keys %LET)){push(@LET,$let);}
        $let=join("",@LET);

        $len1 = @array1;
        if($len1 == 1){$LCA = $array1[0]; }
        elsif($len1 > 1){
                $first = $array1[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @array1);
                        $len2 = @matched;
                        if($len2 == $len1){push(@ARR, $alevel);}
                        else{last;}
                }
                $len3 = @ARR;
                if($len3 > 1){$LCA = join(";", @ARR);}
                elsif($len3==1){$LCA = $ARR[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }

        #add kingdoms to NCA
        if($LCA eq "NCA"){ $LCA.="-".$let; }
        return($LCA);
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
