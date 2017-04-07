#!perl

use Test::More;
use Chemistry::GGL;

my %tests = (
  # simple
  'CC(=O)O' => 'CC(O)=O',
  # aromatize
  'C1=CC=CC=C1' => 'c1ccccc1',
  # keep aromatization
  'c1ccccc1' => 'c1ccccc1',
  # complex ring systems
  'O=C(N)c1ccc[n+](c1)C2OC(C(O)C2O)COP([O-])(=O)OP(=O)(O)OCC5OC(n4cnc3c(ncnc34)N)C(O)C5O' => 'NC(=O)c1ccc[n+](c1)C2OC(COP([O-])(=O)OP(O)(=O)OCC3OC(C(O)C3O)n4cnc5c(N)ncnc54)C(O)C2O',
  # Complex labels
  'C[Mg]' => 'C[Mg]',
  # Invalid labels
  'C[Mx]' => -1,
);

plan tests => (scalar keys %tests);

for my $in (keys %tests) {
  eval {
    my $ret = Chemistry::GGL::canonical_smiles($in, 0);
    is($ret, $tests{$in});
  };
  if ($@) {
    if ($tests{$in} == -1) {
      pass($in);
    }
    else {
      fail($in);
    }
  }
}

