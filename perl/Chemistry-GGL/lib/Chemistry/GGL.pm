package Chemistry::GGL;

=head1

Chemistry::GGL - Perl interface to the Graph Grammar Library

=head1 SYNOPSIS

   use Chemistry::GGL;

   my $smiles = Chemistry::GGL::apply_rules(\@substrates, \@rules, 2);
   # or...
   my $canonical_smiles = Chemistry::GGL::canonical_smiles($smiles_string);

=head1 DESCRIPTION

This package is a Perl interface to the Graph Grammar Library. It currently
provides only two routines which may be useful for an artifical chemist.

=cut

our $VERSION='0.01';

use warnings;
use strict;

use Chemistry::Mol;
use Chemistry::File::SMILES;

use Carp;

BEGIN {
  require DynaLoader;
  @Chemistry::GGL::ISA = qw/DynaLoader/;

  eval {
    bootstrap Chemistry::GGL '0.01';
  };
  if ($@) {
    die "Couldn't load GGL binding\n\n$@\n";
  }
}

$Chemistry::GGL::mol_mode='s';

eval {
  require Chemistry::ReactionGrammar::File::GML;
};
if ($@) {
  $Chemistry::GGL::no_gml = 1;
}

=head1 FUNCTIONS

=over 4

=item apply_rules(\@substrates, \@rules, $outmode)

Take an array of substrate molecules, containing either SMILES strings or
Chemistry::Mol-objects (or even both!), apply the rules given as a reference to an array of GML strings and apply those rules (one iteration). $outmode specifies the return value. It can be either

=over 8

=item 1 (RXNS)

returns an array of hashes specifying reactions. The hash structure is as
follows:

    {
      substrates => [ $smiles1, $smiles2, ... ],
      products   => [ $smiles3, ... ],
    }

=item 2 (SMILES)

simply returns an array of SMILES strings of the products of the rule
applications, already canonicalized.

=back

=cut

sub apply_rules {
  my ($substrates, $rules, $outmode) = @_;

  my @substrates = map { _convert_mol($_) } @$substrates;
  my @rules = map { ref $_ ? sprintf("%s", $_) : $_ } @$rules;

  toyChem(\@substrates, \@rules, $outmode);
}

=item check_rules($rule1, $rule2, ...)

Takes a list of rule objects or strings (in GML format) and feeds it to the rule check of the GGL. Returns an array of hashes, one for each input rule, of the form:

=over 8

=item res

1 (true) if ok, 0 otherwise

=item msg

If not ok, contains a brief error string.

=back

=cut

sub check_rules {
  my @rules = map { ref $_ ? sprintf("%s", $_) : $_ } @_;
  my $ret = checkRules(\@rules);
  wantarray ? @$ret : $ret;
}


sub _convert_mol {
  my $mol = shift;
  my $outstr;
  if (!ref $mol) {
    $outstr = $mol;
  }
  elsif ($Chemistry::GGL::mol_mode eq 'g') {
    croak("GML parser not available, Chemistry::ReactionGrammar installed?")
      if $Chemistry::GGL::no_gml;
    $outstr = Chemistry::ReactionGrammar::File::GML->write_string($mol);
    $outstr =~ s/#.*//;
  }
  else {
    $outstr = $mol->sprintf("%".$Chemistry::GGL::mol_mode);
  }
  return $outstr;
}


=item canonical_smiles($mol, $noProtons)

takes either a Chemistry::Mol object or a SMILES string and returns the canonical SMILES for that molecule. Set $noProtons to a true value to remove protons ([H]) from the resulting SMILES.

=cut

sub canonical_smiles {
  my ($mol, $noProtons) = @_;
  my $ret;
  $noProtons||=0;
  eval {
    $ret = canonicalSmiles(_convert_mol($mol),$noProtons);
  };
  if ($@) {
    die("canonicalize: $@");
  }
  return $ret;

}

=item SV* toyChem(SV* substrates, SV* rules, int outmode)

called by apply_rules

=item SV* canonicalSmiles(SV*)

called by canonical_smiles

=item SV* checkRules(SV* ruleref)

called by check_rules

=back

=head1 VERSION

0.01

=head1 AUTHOR

Heinz Ekker E<lt>hekker@tbi.univie.ac.atE<gt>

=head1 COPYRIGHT

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

Since the GGL is statically linked into this module, it is also subject to
the licence agreement of the GGL itself.

=cut

1;

