#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use Data::Dumper;
use Genome::Ranges;
use Genome::Ranges::IntSet;
use Data::Dumper;
use HOP::Parser qw(:all);
use HOP::Lexer  qw(string_lexer);
use HOP::Stream qw(:all);

use Getopt::Lucid qw(:all);

my @specs =
  (	Switch("help|h")->anycase,
	Param("expr|e"),    # required expression
	Keypair("define|d"),# defined set variables associated files
  );

my $opt = Getopt::Lucid->getopt(\@specs);

if ($opt->get_help || !$opt->get_define || !$opt->get_expr) {
  print STDERR<<'EOT';

beds_operations.pl v0.1 Copyright(C) 2009 Xueya Zhou <xueyazhou@gmail.com>

Usage:
  beds_operations.pl --define|-d <SYMBOL=filename> --expr|-e <EXPRESSION>

  Supported operations include: + (union), - (diff), * (intersect)
  and / (or ^, for exclusive or). The precedence is the same as conventional
  arithematic operations.

  For example:
     perl -I lib bin/beds_operations.pl --expr 'A * B + C / (A - B)' \
       --define A=A.bed --define B=B.bed --define C=C.bed

EOT
exit -1;
}

my %defs = $opt->get_define;

my %VAR;
foreach my $k (keys %defs) {
	$VAR{$k} = Genome::Ranges::IntSet->new($defs{$k}, { bed => 1 });
}

my $input = 'print '. $opt->get_expr . ';';

# Parse the input and calculate the results
my @input_tokens = (
					['TERMINATOR', qr/;\n*|\n+/                 ],
					['PRINT',      qr/\bprint\b/                ],
					['IDENTIFIER', qr|[A-Za-z]\w*|              ],
					['OP',         qr#[-=+*/^()]#                ],
					['WHITESPACE', qr/\s+/,          sub { "" } ],
				   );

my $lexer = iterator_to_stream(string_lexer($input,@input_tokens));

my ($base, $expression, $factor, $program, $statement, $term);
my $Expression = parser { $expression->(@_) };
my $Factor     = parser { $factor->(@_) };
my $Program    = parser { $program->(@_) };
my $Statement  = parser { $statement->(@_) };
my $Term       = parser { $term->(@_) };

$program   = concatenate(star($Statement), \&End_of_Input);

$statement = alternate(T(concatenate(lookfor('PRINT'),
                                     $Expression,
                                     lookfor('TERMINATOR')
									),
                         sub {
						   $_[1]->write(\*STDOUT, {bed => 1});
						 }),
                       T(concatenate(lookfor('IDENTIFIER'),
                                     lookfor(['OP', '=']),
                                     $Expression,
                                     lookfor('TERMINATOR')
                                    ),
                         sub {
						   $VAR{$_[0]} = $_[2];
						 }),
                      );

$expression =  operator($Term,
						[lookfor(['OP', '+']),
						 sub { $_[0] + $_[1] }],
						[lookfor(['OP', '-']),
						 sub { $_[0] - $_[1] }]
					   );

$term       =  operator($Factor,
						[lookfor(['OP', '*']),
						 sub { $_[0] * $_[1] }],
						[lookfor(['OP', '^']),
						 sub { $_[0] ^ $_[1] }],
						[lookfor(['OP', '/']),
						 sub { $_[0] / $_[1] }]
					   );


$factor     = alternate(lookfor('IDENTIFIER',
								sub { $VAR{$_[0][1]} or die "undefined $_[0][1]" }),
						T(concatenate(lookfor(['OP', '(']),
									  $Expression,
									  lookfor(['OP', ')'])),
						  sub { $_[1] })
					   );

$program->($lexer);


__END__

=head1 DESCRIPTION

The Genome::Ranges::IntSet Calculator.

It does not support self-complement which will yield infinite sets.

+ : union, - : diff, * : intersection, ^ : exclusive or
/ is not used here because it will confuse the parser with file path.


The complete grammer rule:

program   -> star(statement) 'End_of_Input'

statement -> 'PRINT' expression 'TERMINATOR'
           | 'IDENTIFIER' '=' expression 'TERMINATOR'

expression-> term star('+' term | '-' term)

term      -> factor star('*' factor | '/' factor)

factor    -> 'IDENTIFIER' | '(' expression ')'


=cut

