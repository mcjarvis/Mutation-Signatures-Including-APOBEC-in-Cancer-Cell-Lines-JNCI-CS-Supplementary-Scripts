#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;

my $type = shift;
my $fasta = shift;
my $variants = shift;
my $stdout = 0;

my $time = time;

unless (defined $type && defined $fasta && defined $variants) {
	my $usage="Error: Not enough inputs\n\n\tUsage: perl count_trinuc_muts <maf, vcf, or pvcf> <reference> <variants>\\nn";
	die $usage;
}

my $name = "$variants";
my $trinucs = "#chr\tpos\t5'tetranuc\t3'tetranuc\ttrinuc\tmut\ttrinuc_mut\tstrand\tflank41bp\tCcount\tTCcount\tTCAcount\tTCTcount\tYTCAcount\tRTCAcount\tsample\n";


my %trinucHash;

my %compDict = (
	A => "T",
	T => "A",
	C => "G",
	G => "C"
);

my %normalization = (
	ACA =>	0.07,
	ACC =>	0.06,
	ACG =>	0.02,
	ACT =>	0.06,
	CCA =>	0.10,
	CCC =>	0.09,
	CCG =>	0.02,
	CCT =>	0.09,
	GCA =>	0.07,
	GCC =>	0.08,
	GCG =>	0.02,
	GCT =>	0.07,
	TCA =>	0.07,
	TCC =>	0.08,
	TCG =>	0.01,
	TCT =>	0.08
);

my @norm = (
	0.07,
	0.06,
	0.02,
	0.06,
	0.10,
	0.09,
	0.02,
	0.09,
	0.07,
	0.08,
	0.02,
	0.07,
	0.07,
	0.08,
	0.01,
	0.08
);

#system("awk 'BEGIN {FS="\t"}; {print "chr"$5"\t"$6"\t"$18"\t"$12"\t"$16}' $filein");

print "creating fasta db... ";
my $db = Bio::DB::Fasta->new($fasta);
print "complete!\n";

my %repeatHash;

open (VAR, "< $variants");
while (my $line = <VAR>) {
	chomp($line);
	next if $line =~ /^$/;
	my ($chr,$pos,$ref,$alt,$strand);
	my $sample = "";
	my @fields = split("\t",$line);
	unless ($type eq "maf") {
		$chr = shift @fields;
		$pos = shift @fields;
		$ref = shift @fields;
		$alt = shift @fields;
		$strand = shift @fields;
		$sample = shift @fields if $type eq "pvcf"; # if sample name is available
	} else {
		next unless $fields[9] eq "SNP";
		$chr = $fields[4];
		$pos = $fields[5];
		$ref = $fields[10];
		$alt = $fields[12];
		$strand = $fields[16];
		$sample = substr($fields[15],0,12);
	}
	next unless length($alt) == 1;
	#print "$chr $pos $ref $alt	$strand\n";
	if ($strand =~ /^\+$/) {
		my $TCcxt=0;
		my $GAcxt=0;
		my $CGcxt=0;
		my $TCAcxt=0;
		my $TGAcxt=0;
		my $TCTcxt=0;
		my $AGAcxt=0;
		my $RTCAcxt=0;
		my $TGAYcxt=0;
		my $YTCAcxt=0;
		my $TGARcxt=0;

		my $flank = uc($db->seq($chr, $pos-10 => $pos+10));
		my $flank40 = uc($db->seq($chr, $pos-20 => $pos+20));
		while ($flank40 =~ /TC/g) { $TCcxt++ }
		while ($flank40 =~ /GA/g) { $GAcxt++ }
		while ($flank40 =~ /TC[AT]/g) { $TCAcxt++ }
		while ($flank40 =~ /[AT]GA/g) { $TGAcxt++ }
		while ($flank40 =~ /TC[AT]/g) { $TCTcxt++ }
		while ($flank40 =~ /[AT]GA/g) { $AGAcxt++ }

		while ($flank40 =~ /[GA]TCA/g) { $RTCAcxt++ }
		while ($flank40 =~ /TGA[CT]/g) { $TGAYcxt++ }
		while ($flank40 =~ /[CT]TCA/g) { $YTCAcxt++ }
		while ($flank40 =~ /TGA[GA]/g) { $TGARcxt++ }

		my $TCNcxt = $TCcxt + $GAcxt;
		my $A3cxt1 = $TGAcxt + $TCAcxt;
		my $A3cxt2 = $AGAcxt + $TCTcxt;
		my $YtetraCxt = $YTCAcxt + $TGARcxt;
		my $RtetraCxt = $RTCAcxt + $TGAYcxt;

		while ($flank40 =~ /[CG]/g) { $CGcxt++ }
		if (exists $repeatHash{$flank}) {
			next;
		} else {
			$repeatHash{$flank} = 1;
		}
		my $seq = uc($db->seq($chr, $pos-2 => $pos+2));
		my @bases = split('', $seq);
		$seq = $bases[1] . $bases[2] . $bases[3];
		my $tetraseq5 = $bases[0] . $bases[1] . $bases[2] . $bases[3];
		my $tetraseq3 = $bases[1] . $bases[2] . $bases[3] . $bases[4];
		my $fullstring = "$bases[1]\[$ref>$alt\]$bases[3]";
		$trinucs .= "$chr\t$pos\t$tetraseq5\t$tetraseq3\t$seq\t$ref>$alt\t$fullstring\t$strand\t$flank40\t$CGcxt\t$TCNcxt\t$A3cxt1\t$A3cxt2\t$YtetraCxt\t$RtetraCxt\t$sample\n";

	} elsif ($strand =~ /^\-$/ ) {
		my $TCcxt=0;
		my $GAcxt=0;
		my $CGcxt=0;
		my $TCAcxt=0;
		my $TGAcxt=0;
		my $TCTcxt=0;
		my $AGAcxt=0;
		my $RTCAcxt=0;
		my $TGAYcxt=0;
		my $YTCAcxt=0;
		my $TGARcxt=0;

		my $flank = uc($db->seq($chr, $pos-10 => $pos+10));
		my $flank40 = uc($db->seq($chr, $pos-20 => $pos+20));
		while ($flank40 =~ /TC/g) { $TCcxt++ }
		while ($flank40 =~ /GA/g) { $GAcxt++ }
		while ($flank40 =~ /TCA/g) { $TCAcxt++ }
		while ($flank40 =~ /TGA/g) { $TGAcxt++ }
		while ($flank40 =~ /TCT/g) { $TCTcxt++ }
		while ($flank40 =~ /AGA/g) { $AGAcxt++ }

		while ($flank40 =~ /[GA]TCA/g) { $RTCAcxt++ }
		while ($flank40 =~ /TGA[CT]/g) { $TGAYcxt++ }
		while ($flank40 =~ /[CT]TCA/g) { $YTCAcxt++ }
		while ($flank40 =~ /TGA[GA]/g) { $TGARcxt++ }

		my $TCNcxt = $TCcxt + $GAcxt;
		my $A3cxt1 = $TGAcxt + $TCAcxt;
		my $A3cxt2 = $AGAcxt + $TCTcxt;
		my $YtetraCxt = $YTCAcxt + $TGARcxt;
		my $RtetraCxt = $RTCAcxt + $TGAYcxt;

		while ($flank40 =~ /[CG]/g) { $CGcxt++ }
		my @fbases = split('', $flank40);
		my @revcomp;
		foreach(@fbases) {
			unshift(@revcomp, $compDict{$_});
		}
		my $rcflank = join('', @revcomp);
		if (exists $repeatHash{$rcflank}) {
			next;
		} else {
			$repeatHash{$rcflank} = 1;
		}
		my $seq = uc($db->seq($chr, $pos-2 => $pos+3));
		my @bases = split('', $seq);
		$seq = $compDict{$bases[3]} . $compDict{$bases[2]} . $compDict{$bases[1]};
		my $tetraseq3 = $compDict{$bases[3]} . $compDict{$bases[2]} . $compDict{$bases[1]} . $compDict{$bases[0]};
		my $tetraseq5 = $compDict{$bases[4]} . $compDict{$bases[3]} . $compDict{$bases[2]} . $compDict{$bases[1]};
		my $fullstring = "$compDict{$bases[3]}\[$ref>$alt\]$compDict{$bases[1]}";
		$trinucs .= "$chr\t$pos\t$tetraseq5\t$tetraseq3\t$seq\t$ref>$alt\t$fullstring\t$strand\t$rcflank\t$CGcxt\t$TCNcxt\t$A3cxt1\t$A3cxt2\t$YtetraCxt\t$RtetraCxt\t$sample\n";

	} else { next }
}

if ($stdout == 0) {
	open (OUT, "> $name.$time.count.txt");
	print OUT $trinucs;
	close(OUT);
} else {
	print $trinucs;
}
# count all TCW:WGA and sum flank C:Gs

# count all C:Gs and sum flank TCW:WGA

### Process maf file
# data.w <- mapply("*", data.c, data.radj, SIMPLIFY=F)

exit;
