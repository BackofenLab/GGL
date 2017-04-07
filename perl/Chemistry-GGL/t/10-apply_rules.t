#!perl

use Chemistry::GGL;

use Test::More;
use File::Find;

my @tests;
find(sub { m#\.test$# && push @tests, $File::Find::name }, 
  qw{t/10-apply_rules});

plan tests => scalar @tests;

for my $f (@tests) {
  local $/=undef;
  open(IF, "<$f") or die("Open $f: $!");
  my $tmp = <IF>;
  close(IF);
  my @p = split(qr#^//.*$#m, $tmp);
  unless(@p == 4) {
    warn "Invalid input file: $f";
    next;
  }
  my @in;
  for my $m (split("\n", $p[1])) {
    next if $m=~m#^$#;
    push @in, $m;
  }
  my $rule = $p[2];
  my %want;
  for my $m (split("\n", $p[3])) {
    next if $m=~m#^$#;
    $want{$m}++;
  }
  my $ret = Chemistry::GGL::apply_rules(\@in, [$rule], 2);

  my $ok=1;
  for my $prod (@$ret) {
    $ok=0 unless exists $want{$prod};
    if (! exists $want{$prod}) {warn "\nCannot find in target products: $prod\n";}
    $want{$prod}--;
  }
  for my $prod (keys %want) {
    $ok=0 unless $want{$prod}==0;
    if ($want{$prod} > 0) {warn "Missing product: $prod\n";}
  }

  ok($ok, "$f: some products missing");
}

  
