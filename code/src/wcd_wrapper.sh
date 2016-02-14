#!/bin/sh

SEQFILE=$2
OUTFILE=${SEQFILE}X.clt
WCD=wcd
COMP2EXT=comp2ext.pl
$WCD --show_clusters $SEQFILE | $COMP2EXT >$OUTFILE

