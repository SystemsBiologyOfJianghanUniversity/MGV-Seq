#! /usr/bin/perl -w

use Bio::DB::Sam;
use Getopt::Long;
die "USAGES perl $0 -bam <bam> -ref <reference fasta> -bed <bed> -ssr <ssrFile>  \n
	<bam>: input alignment
	<reference fasta>: the reference genome in fasta format
	<bed>: the target regions and the SNPs within\n
	<ssr>: the tandem repeats among the target region\n
	
	" unless @ARGV >0;

my ($bam,$refFasta,$bed) = ();

GetOptions(
	'bam=s' => \$bam,
	'ref=s' => \$refFasta,
	'bed=s' => \$bed,
	'ssr=s' => \$ssrFile,
	
	
);



#------------------------------------#
#----- Preparation
#------------------------------------#

unless(defined $bam)
{
	die "Error! A bam file is required\n";
}
unless(defined $refFasta)
{
	die "Error! Rerence genome is required\n";
}
unless(defined $bed)
{
	die "Error! a bed file for target regions is required\n";
}
unless(defined $ssrFile)
{
	die "Error! a ssr file is required\n";
}




my $conseBP = 2;
my $indelBP = 1;



my %typrep  = ('1',3,'2',3,'3',3,'4',3,'5',3,'6',3,'7',3); 
my @typ = sort { $a <=> $b } keys %typrep;


#-------------------------------------#
#---------Process SSRs on reference
#-------------------------------------#
my %SSR_END;
my %SSR_MOTIF;
my %SSR_REPEAT;


open IN_SSR, "$ssrFile" or die "cannot open $ssrFile\n";
<IN_SSR>;
while(my $line=<IN_SSR>)
{
	my ($amp_id,$ssr_index,$ssr_p1,$ssr_p2,$ssr_motif,$ssr_repeat)=(split /\s+/,$line);
	$SSR_END{$amp_id}{$ssr_p1}=$ssr_p2;
	$SSR_MOTIF{$amp_id}{$ssr_p1}=$ssr_motif;
	$SSR_REPEAT{$amp_id}{$ssr_p1}=$ssr_repeat;	
}
close(IN_SSR);


#-------------------------------------#
#------process bam based on the region in bed
#-------------------------------------#

my $sam = Bio::DB::Sam->new(-bam=>$bam,-fasta=>$refFasta,-autoindex=>1);

#No need to scan every parts of the chromosomes, so segment the chromosome, and scan each segment
my $lastLine = "";
my $flankNum = 2;

open IN_BED, "$bed" or die "cannot open BED file $bed\n";
while(my $line=<IN_BED>)
{
    my %h=();
	my %s=(); #strand
	my ($chr,$segment_start,$segment_end,$segment_name)=(split /\s+/,$line)[0,1,2,3];
	my $segment_len = $segment_end-$segment_start+1;
	###process each segment/window
	#print "$chr $segment_start $segment_end: ";
	my @alignments = $sam->get_features_by_location(-seq_id => $chr,-start=>$segment_start,-end =>$segment_end);#!!
	my $temNum = scalar(@alignments);
	#print "$temNum Reads\n";
	
	my @markPos1=();
	my @markPos2=();
	
	
	#retrieve marked regions
	foreach my $ssr_p1 (sort{$a<=>$b}keys %{$SSR_END{$segment_name}})
	{
		if($SSR_REPEAT{$segment_name}{$ssr_p1}>=3) #all ssrs 
		{
			if($ssr_p1-$flankNum<=0)
			{
				push @markPos1,1;
			}else
			{
				push @markPos1,$ssr_p1-$flankNum;
			}
			my $ssr_p2=$SSR_END{$segment_name}{$ssr_p1};
			if($ssr_p2+$flankNum>$segment_len)
			{
				push @markPos2,$segment_len;
			}else
			{
				push @markPos2,$ssr_p2+$flankNum;
			}
		}
	}
		
		
	#merge overlapping region
	for(my $i=0;$i<scalar(@markPos1)-1;$i++)
	{
		for(my $j=$i+1;$j<@markPos1;$j++)
		{
			if($markPos2[$i]>=$markPos1[$j] and $markPos1[$j]!=-1)
			{
				$markPos2[$i] = $markPos2[$j];
				$markPos1[$j]=-1;
				$markPos2[$j]=-1;
			}
		}	
	}
	@markPos1 = grep{$_>0}@markPos1;
	@markPos2 = grep{$_>0}@markPos2;
		
	

	#go through each read covering the segment
	foreach my $align(@alignments)
	{
		my $read_start  = $align->start;#the start position of the alignment in the genome
		my $read_end    = $align->end;  #the end position of the alignment in the genome
		#print "$chr $read_start $read_end\n";
        next unless $read_start <= $segment_start-10 and $read_end >= $segment_end+10; #
		
		my $read_strand = $align->strand;
		#my $read_strand = $align->strand;
		my $read_cigar  = $align->cigar_str;
		my $read_name = $align->query->name;
		my $str = $align->aux;
		if($read_cigar=~m/^(\d+)S/)
        {
            next if $1>=50;       
        }
        if($read_cigar=~m/(\d+)S$/)
        {
            next if $1>=50;
        }
        my ($ref,$matches,$query) = $align->padded_alignment;
		
        	
		#--- mask the SSR in reads 
		
		$query = maskReadSSR($query);
			

		#--- mask snps near indels
		
		$query=maskSNPNearIndel($ref,$matches,$query);
		
		#--- mask consecutive snps
		
		$query=maskConsecutiveSNP($ref,$matches,$query);
		
		
		my $allele =  getAllele($ref,$matches,$query,$read_start,$read_cigar,$segment_start,$segment_end);
        #print "$allele\n";
		
		## mark the repeats (1).poly A/T/C/G (length >=5),and 5bp up/downstream
		$allele=getMarkedAllele($allele,\@markPos1,\@markPos2);
		
		if(exists $h{$allele})
		{
			$h{$allele}++;
			if($read_strand==1)
			{
				$s{$allele}++;
			}
			
		}else
		{
			$h{$allele}=1; 
			if($read_strand==1)
			{
				$s{$allele}=1;
			}else
			{
				$s{$allele}=0;
			}
		}

	}
    #output 
	
    if(@alignments>0)
    {
        my $i = 1;
        foreach my $allele(sort{$h{$b} <=> $h{$a}} keys %h)
        {
            print "$segment_name\t$i\t$allele\t$h{$allele}\t$s{$allele}\n";
			$i++;
        }
    }

}

close(IN_BED);


#---mask snps near indel 
sub maskSNPNearIndel
{
	my ($ref,$matches,$query)=@_;
	
	my @REF = (split //,$ref);
	my @MATCHES = (split //,$matches);
	my @QUERY = (split //,$query);
	for(my $i=0;$i<@MATCHES;$i++)
	{
		if($MATCHES[$i] eq " ")
		{
			if($QUERY[$i] eq "-" or $REF[$i] eq "-")#SNPs
			{
				my $p1= $i-$indelBP;
				my $p2= $i+$indelBP;
				if($p1<0){$p1=0}
				if($p2>=@MATCHES){$p2=scalar(@MATCHES)-1}
				for(my $j=$p1;$j<=$p2;$j++)
				{
					if($QUERY[$j] ne "-")
					{
						$QUERY[$j] = "N";
					}
				}
			}
		}
	}
	$query = join("",@QUERY);
	return($query);
}

#---mask consecutive mismathcs
sub maskConsecutiveSNP
{
	my ($ref,$matches,$query)=@_;
	
	my @REF = (split //,$ref);
	my @MATCHES = (split //,$matches);
	my @QUERY = (split //,$query);
	my @MISMATCH = ();
	for(my $i=0;$i<@MATCHES;$i++)
	{
		if($MATCHES[$i] eq " ")
		{
			push @MISMATCH,$i;
		}
	}
	@MISMATCH = sort{$a <=>$b}@MISMATCH;
	
	if(@MISMATCH > 0)
	{
		my @toRemove;
		for($i=1;$i<@MISMATCH;$i++)
		{
			if($MISMATCH[$i]-$MISMATCH[$i-1]<=$conseBP)
			{
				push @toRemove,$MISMATCH[$i];
				push @toRemove,$MISMATCH[$i-1];
			}
		}
		if(@toRemove>0)
		{
			for($i=0;$i<@toRemove;$i++)
			{
				if($QUERY[$toRemove[$i]] ne "-")
				{
					$QUERY[$toRemove[$i]] = "N";
				}
			}
		}
	}
	$query = join("",@QUERY);
	return($query);
}
#----remove Insert in target region
sub removeInsertInRegion
{
	my ($str_insert,$pos1,$pos2)=@_;
	while($str_insert=~m/(\d+)\|([ATCG]+)/ig)
	{
		if($1>=$pos1 and $1<=$pos2)
		{
			#print "$str_insert\n";
			$str_insert=~s/$1\|$2//;
			#print "$str_insert\n";
		}
	}
	$str_insert=~s/\_{2,}/\_/;
	$str_insert=~s/\_$//;
	$str_insert=~s/:\_+/:/;
	return($str_insert);
}


#-----mark positions within SSRs
sub getMarkedAllele
{
	my $allele=$_[0];
	my @markPos1=@{$_[1]};
	my @markPos2=@{$_[2]};
	
	if(@markPos1>0)
	{
		my ($mnp,$str_insert)=(split /,/,$allele);
		my $tmp_mnp="";

		
		for(my $i=0;$i<@markPos1;$i++)
		{
			#process mnp
			if($i==0)
			{
				$tmp_mnp = substr($mnp,0,$markPos1[$i]-1);
				$tmp_mnp = $tmp_mnp."N"x($markPos2[$i]-$markPos1[$i]+1)
			}else
			{
				$tmp_mnp = $tmp_mnp.substr($mnp,length($tmp_mnp),$markPos1[$i]-length($tmp_mnp)-1);
				$tmp_mnp = $tmp_mnp."N"x($markPos2[$i]-$markPos1[$i]+1)
			}
			#process insert
			$str_insert = removeInsertInRegion($str_insert,$markPos1[$i],$markPos2[$i]);
		}
		$tmp_mnp = $tmp_mnp.substr($mnp,length($tmp_mnp),length($mnp)-length($tmp_mnp));
		#print "$allele\n";
		$allele=$tmp_mnp.",".$str_insert;
		#print "$allele\n";
	}
	

	return($allele);
}
#-----get read genotype at a specific position

sub getAllele
{
  my ($ref,$matches,$query,$read_start,$read_cigar,$segment_start,$segment_end)=@_;
  #$segment_start = $segment_start-$read_start;
  #$segment_end = $segment_end-$read_start;
  
  my %insert ;
  my $str_insert="";
  my $res;
  #process soft clip
  if($read_cigar=~m/^(\d+)S/)
  {
    $ref=substr($ref,$1);
    $matches=substr($matches,$1);
    $query=substr($query,$1);


#    print "$ref\n$matches\n$query\n";

  }
  

  if($read_cigar=~m/(\d+)S$/)
  {
      $ref=substr($ref,0,length($ref)-$1);
      $matches=substr($matches,0,length($matches)-$1);
      $query=substr($query,0,length($query)-$1);
#      print "$ref\n$matches\n$query\n";
  }

  #process insertion 
  my @i_start;
  my @i_end; #the position
  while($ref=~m/(-+)/ig)
  {
	my $tmp_i_end = pos($ref)-1;
	my $tmp_i_start = $tmp_i_end-length($1)+1;
	push @i_start,$tmp_i_start;
	push @i_end,$tmp_i_end;
  }
  if(@i_start >0)
  {
	
	my $tmp_ref;
	my $tmp_matches;
	my $tmp_query;
	for(my $i=0;$i<@i_start;$i++)
	{
		if($i==0)
		{
			$tmp_ref=substr($ref,0,$i_start[$i]);
			$tmp_matches=substr($matches,0,$i_start[$i]);
			$tmp_query=substr($query,0,$i_start[$i]);
			
		}else
		{
			$tmp_ref=$tmp_ref.substr($ref,$i_end[$i-1]+1,$i_start[$i]-$i_end[$i-1]-1);
			$tmp_matches=$tmp_matches.substr($matches,$i_end[$i-1]+1,$i_start[$i]-$i_end[$i-1]-1);
			$tmp_query=$tmp_query.substr($query,$i_end[$i-1]+1,$i_start[$i]-$i_end[$i-1]-1);
			
		}
		my $i_pos = length($tmp_ref);
        
		$insert{$i_pos}=substr($query,$i_start[$i],$i_end[$i]-$i_start[$i]+1);
		if($i_pos+$read_start-1>=$segment_start and $i_pos+$read_start-1<=$segment_end)
		{
			#print "$tmp_ref\n";
			my $i_pos_relative=$i_pos+$read_start-$segment_start;
			$str_insert=$str_insert."$i_pos_relative|$insert{$i_pos}_";
		}

	}
	$tmp_ref=$tmp_ref.substr($ref,$i_end[-1]+1);
	$tmp_matches=$tmp_matches.substr($matches,$i_end[-1]+1);
	$tmp_query=$tmp_query.substr($query,$i_end[-1]+1);
	
	$ref = $tmp_ref;
	$matches = $tmp_matches;
	$query = $tmp_query;
	$str_insert=~s/_$//;

  }
  
  ## get consensus sequence
  
  @MATCHES=(split //,$matches);
  @QUERY=(split //,$query);
  for(my $i=0;$i<@MATCHES;$i++)
  {
	if($MATCHES[$i] eq "|")
	{
		$MATCHES[$i]="R";
	}else
	{
		$MATCHES[$i]=$QUERY[$i];
	}
  }
  
  $str_insert="";
  #$res = join("",@MATCHES);
  $res = join("",@QUERY);
  $res = substr($res,$segment_start-$read_start,$segment_end-$segment_start+1);
  $res=$res.",I:".$str_insert;
 
  
  #print "$res\n";
  return($res);
}

#-------detect ssr in reads
sub str_detection
{
    my $seq = shift @_;
	my (@start,@end,@motif,@repeats,@strs);
    ###detect STRs modified based on misa.pl
    ### each repeat has following features: 1) start postion; 2)end position; 3)motif; 4) repeated sequence; 5) number of repeats(optional)
    for (my $i=0; $i < scalar(@typ); $i++) #check each motif class
    {
        my $motiflen = $typ[$i];
        my $minreps = $typrep{$typ[$i]} - 1;
        
        my $search = "(([acgt]{$motiflen})\\2{$minreps,})";
		#print "Search string: $search\n";
        while ( $seq =~ /$search/ig ) #detect the target motif among the sequence
        {
            
			my $temMotif = uc $2;
			#print "temMotif: $temMotif\n";
            my $redundant = 0;; #reject false type motifs [e.g. (TT)6 or (ACAC)5]
            for (my $j = $motiflen - 1; $j > 0; $j--)
            {
				next if $motiflen/$j-1 < 1;
                my $redmotif = "([ACGT]{$j})\\1{".($motiflen/$j-1)."}";
				#print  "$1 $2\n";
				#print "redmotif: $redmotif\n";
				#print "motiflen: $motiflen\n";
                if(($motiflen/$j-1)!=int(($motiflen/$j-1)))
                {
					$redundant = 1;
                }else
				{	
	
                	$redundant = 1 if ( $temMotif =~ /$redmotif/ )
				}
				#print "$redundant\n";		
            }
            next if $redundant;
 
            my $ssr = uc $1;
            my $temRepeats= length($ssr) / $motiflen;
            my $temEnd = pos($seq)-1; 
            my $temStart = $temEnd - length($ssr)+1;
            
        #   print "$temStart $temEnd $temMotif $ssr\n";
      
            push @start,$temStart;
            push @end,$temEnd;
            push @motif,$temMotif;
            push @repeats,$temRepeats;
            push @strs,$ssr;
      
        }
    }
	#---reorder
	my @refs=&reOrderSTRs(\@start,\@end,\@motif,\@repeats);
	@start =@{$refs[0]};
	@end = @{$refs[1]};
	@motif = @{$refs[2]};
	@repeats = @{$refs[3]};
    return(\@start,\@end,\@motif,\@repeats);
}
#-----mark ssrs in reads
sub maskReadSSR
{
	my $query = shift @_;
	my $nt2mask = 2;# to be update
	
	my $newQuery = "";
	
	my $seq =$query;
	#print ("SEQ: $seq\n");
	#---locate the deletions (
	my @POS_STARTS = ();
	my @MINUS_STRS = ();
	while($seq=~m/(\-+)/g)
	{
		my $minusSTR = $1;
		my $posEnd = pos($seq)-1;
		my $posStart = $posEnd-length($minusSTR)+1;
		push @POS_STARTS,$posStart;
		push @MINUS_STRS,$minusSTR;
		#print "$minusSTR\t$posStart\t$posEnd\n";
	}

	
	
	#----locate the SSR
	$seq =~s/-//g;
	my @refs=&str_detection($seq);
    my @start =@{$refs[0]};
	my @end = @{$refs[1]};
    my @motif = @{$refs[2]};
    my @repeats = @{$refs[3]};
	
	
	#----mask STRs 
	my $newSeq = $seq;
	if(scalar(@start)>0)
	{
		for(my $i=0;$i<scalar(@start);$i++)
		{
			my $p1 = $start[$i]-$nt2mask;
			my $p2 = $end[$i]+$nt2mask;
			if($p1<0){$p1=0;}
			if($p2>=length($seq)){$p2 = length($seq)-1;}
			
			
			#process mnp
			if($p1==0)
			{
				
				$newSeq = "N"x($p2-$p1+1).substr($newSeq,$p2+1,length($newSeq)-$p2);
			}else
			{
				if($p2 != length($seq)-1)
				{
					$newSeq = substr($newSeq,0,$p1)."N"x($p2-$p1+1).substr($newSeq,$p2+1,length($newSeq)-$p2+1);
				}else
				{
					$newSeq = substr($newSeq,0,$p1)."N"x($p2-$p1+1);
				}
			}
			
		}
		#print "$newSeq\n";
	}
	
	#---put minus back
	if(scalar(@POS_STARTS)>0)
	{
		for(my $i=0;$i<scalar(@POS_STARTS);$i++)
		{
			my $p1 = $POS_STARTS[$i];
			if($p1 ==0)
			{
				$newSeq = $MINUS_STRS[$i].$newSeq;
			}else
			{
				my $tmpSTR1 = substr($newSeq,0,$p1);
				my $tmpSTR2 = substr($newSeq,$p1,length($newSeq)-$p1+1);
				$newSeq = $tmpSTR1.$MINUS_STRS[$i].$tmpSTR2;
			}
		}
	}
	
	$newQuery = $newSeq;
	
	return($newQuery);
}
#-----re-order the ssrs
sub reOrderSTRs
{

    my @start = @{$_[0]};
    my @end   = @{$_[1]};
    my @motif = @{$_[2]};
    my @repeats = @{$_[3]};
    
    my @index = (0..(scalar(@start)-1));
    my @order = sort { $start[$a] <=> $start[$b] } @index;
    my (@newStart2,@newEnd2,@newMotif2,@newRepeat2);
    for(my $j=0;$j<@order;$j++)
    {
        my $i=$order[$j];
        push @newStart2,$start[$i];
        push @newEnd2,$end[$i];
        push @newMotif2,$motif[$i];
        push @newRepeat2,$repeats[$i];
    }

    return(\@newStart2,\@newEnd2,\@newMotif2,\@newRepeat2);
}