#!/usr/bin/tclsh

# -------------------------------------------------------------------

proc DNA {} {

  set input AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC

  set a 0 ; set g 0
  set c 0 ; set t 0

  set length [string length $input]
  for {set i 0} {$i < $length} {incr i} {
    incr [string tolower [string index $input $i]]
  }

  puts "$a $c $g $t"

}

# -------------------------------------------------------------------

proc RNA {} {

  set input GATGGAACTTGACTACGTAAATT

  set length [string length $input]
  for {set i 0} {$i < $length} {incr i} {
    set s [string index $input $i]
    puts -nonewline [expr {$s == "T" ? "U" : $s}]
  }

  puts ""

}

# -------------------------------------------------------------------

proc REVC {} {

  set input AAAACCCGGT

  set rv [string reverse $input]

  set comp(A) T ; set comp(T) A
  set comp(G) C ; set comp(C) G

  set length [string length $rv]
  for {set i 0} {$i < $length} {incr i} {
    puts -nonewline $comp([string index $rv $i])
  }

  puts ""

}

# -------------------------------------------------------------------

proc FIB {} {

  set lst {1 1} ; set n 5 ;  set k 3

  for {set i 2} {$i < $n} {incr i} {
    set f1 [lindex $lst [expr {$i - 1}]]
    set f2 [lindex $lst [expr {$i - 2}]]
    lappend lst [expr {$f1 + $f2 * $k}]
  }

  puts $lst

}

# -------------------------------------------------------------------

proc GC {} {

  set lines [readLines "gc.fasta"]

  gcsAndLengths $lines gcs lengths
  set rslt [maxPercentageGC gcs lengths]
  puts "[lindex $rslt 0]\n[lindex $rslt 1]"

}

proc readLines {fname} {

  # Read data
  set f [open $fname r]
  set lines [split [read $f] \n]
  close $f

  # Remove empty line
  set le [expr {[llength $lines] - 1}]
  return [lreplace $lines $le $le]

}

proc gcsAndLengths {lines gcs lengths} {

  upvar $gcs gg ; upvar $lengths ll

  set cString ""
  foreach l $lines {

    if {">" == [string index $l 0]} {

      set cString [string range $l 1 end]
      set gg($cString) 0
      set ll($cString) 0

    } else {

      incr gg($cString) [lineGC $l]
      incr ll($cString) [string length $l]

    }

  }

}

proc lineGC {line} {

  set gc 0 ; set len [string length $line]

  for {set i 0} {$i < $len} {incr i} {
    set c [string index $line $i]
    if {"G" == $c || "C" == $c} { incr gc }
  }

  return $gc

}

proc maxPercentageGC {gcs lengths} {

  upvar $gcs gg ; upvar $lengths ll

  set maxV -1 ; set maxS ""
  foreach k [array names gg] {

    set g $gg($k) ; set l $ll($k)
    set p [expr {double($g) / $l * 100}]

    if {$p > $maxV} { set maxV $p ; set maxS $k }

  }

  return "$maxS $maxV"

}

# -------------------------------------------------------------------

proc HAMM {} {

  set str1 GAGCCTACTAACGGGAT
  set str2 CATCGTAATGACGGCCT

  set le [string length $str1] ; set distance 0
  for {set i 0} {$i < $le} {incr i} {
    if {[string index $str1 $i] !=\
        [string index $str2 $i]} {
      incr distance
    }
  }

  puts $distance

}

# -------------------------------------------------------------------

proc IPRB {} {

  set types 3
  set k 2 ; set m 2 ; set n 2
  set orgs "$k $m $n"
  set tot [expr {double($k + $m + $n)}]
  set tot1 [expr {$tot - 1}]

  # Pairs and their probabilities
  set pairs [calPairs $types]

  set pro 0
  for {set i 0} {$i < $types} {incr i} {

    # First organism
    set fProb [expr {[lindex $orgs $i] / $tot}]

    for {set j 0} {$j < $types} {incr j} {

      # Second organism
      set sProb [expr {([lindex $orgs $j] - ($i == $j ? 1 : 0)) / $tot1}]

      # Probability of dominant allele
      set pro [expr {$pro + $fProb * $sProb * [getProb "$i $j" $pairs]}]

    }
  }

  # Sum of probabilities
  puts $pro

}

proc calPairs {types} {

  set pairs {}
  for {set i 0} {$i < $types} {incr i} {
    for {set j $i} {$j < $types} {incr j} {
      lappend pairs "{$i $j} [probDom $i $j]"
    }
  }

  return $pairs

}

# Calculate probability for pair.
proc probDom {i1 i2} {

  set alleles {YY Yy yy}

  set l 2 ; set dom 0 ; set d "Y"
  for {set i 0} {$i < $l} {incr i} {
    for {set j 0} {$j < $l} {incr j} {
      if {[string index [lindex $alleles $i1] $i] == $d ||
          [string index [lindex $alleles $i2] $j] == $d} {
        incr dom
      }
    }
  }

  return [expr {$dom / 4.0}]

}

# Find probability for pair
proc getProb {pair pairs} {

  set pair [lsort -integer $pair]
  set l [llength $pairs]

  for {set i 0} {$i < $l} {incr i} {
    set cpair [lindex $pairs $i]
    if {$pair == [lindex $cpair 0]} {
      return [lindex $cpair 1]
    }
  }

}

# -------------------------------------------------------------------

proc PROT {} {

  set input AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

  set f [open "rna_codon.txt" r]
  set codons [regexp -all -inline {\S+} [read $f]]
  close $f

  set le [string length $input]
  for {set i 0} {$i < $le} {incr i 3} {
    set codon [string range $input $i [expr {$i + 2}]]
    set aa [findAminoAcid $codon $codons]
    if {$aa != "Stop"} { puts -nonewline $aa }
  }

  puts ""

}

proc findAminoAcid {codon table} {
  set le [llength $table]
  for {set i 0} {$i < $le} {incr i 2} {
    if {$codon == [lindex $table $i]} {
      return [lindex $table [expr {$i + 1}]]
    }
  }
}

# -------------------------------------------------------------------

proc SUBS {} {

  set input "GATATATGCATATACTT
ATAT"

  set str [lindex $input 0]
  set mot [lindex $input 1]

  set strLen [string length $str]
  set motLen [string length $mot]

  for {set i 0} {$i <= $strLen - $motLen} {incr i} {
    set sub [string range $str $i [expr {$i + $motLen - 1}]]
    if {$sub == $mot} { puts -nonewline "[expr {$i + 1}] " }
  }

  puts ""

}

# -------------------------------------------------------------------

proc CONS {} {

  set lines [readLines "cons.fasta"]

  # Collect strings
  set dnaStrings {} ; set currentDna ""
  foreach line $lines {
    if {">" == [string index $line 0]} {
      if {"" != $currentDna} {
        lappend dnaStrings $currentDna
        set currentDna ""
      }
    } else {
      set currentDna $currentDna$line
    }
  }

  lappend dnaStrings $currentDna ; unset currentDna

  # Display consensus string and profile matrix
  conStringMatrix [calFreqs $dnaStrings]

}

proc conStringMatrix {frequencies} {

  set s {A C G T}

  # Display consensus string
  conString $frequencies $s

  # Display profile matrix
  profileMatrix $frequencies $s

}

proc conString {frequencies s} {

  foreach f $frequencies {

    set maxI 0
    for {set i 1} {$i < 4} {incr i} {
      if {[lindex $f $i] > [lindex $f $maxI]} {
        set maxI $i
      }
    }

    puts -nonewline [lindex $s $maxI]

  }

  puts ""
}

proc profileMatrix {frequencies s} {
  for {set i 0} {$i < 4} {incr i} {
    puts -nonewline "[lindex $s $i]: "
    foreach f $frequencies { puts -nonewline "[lindex $f $i] " }
    puts ""
  }
}

proc calFreqs {dnaStrings} {

  set frequencies {}
  set dnaLen [string length [lindex $dnaStrings 0]]
  for {set i 0} {$i < $dnaLen} {incr i} {

    set fqs(A) 0 ; set fqs(C) 0
    set fqs(G) 0 ; set fqs(T) 0

    foreach dna $dnaStrings {
      incr fqs([string index $dna $i])
    }

    lappend frequencies "$fqs(A) $fqs(C) $fqs(G) $fqs(T)"

  }

  return $frequencies

}

# -------------------------------------------------------------------

proc FIBD {} {

  set lst {1 1} ; set n 10 ; set m 3

  for {set i 2} {$i < $n} {incr i} {

    set ind1 [expr {$m > $i ? 0 : $i - $m}]
    set ind2 [expr {$ind1 + 1}]

    set fc [lindex $lst $ind2]
    set ff [lindex $lst $ind1]

    lappend lst [expr {$fc + $ff}]

  }

  puts $lst

}

# -------------------------------------------------------------------

proc GRPH {} {

  set lines [readLines "grph.fasta"]

  set cKey "" ; set cStr ""
  foreach l $lines {

    if {">" == [string index $l 0]} {

      if {"" != $cStr} { set arr($cKey) $cStr ; set cStr "" }
      set cKey [string range $l 1 end]

    } else {
      set cStr $cStr$l
    }
  }

  set arr($cKey) $cStr ; unset cStr l

  # O3
  set k 3
  set keys [array names arr]

  foreach i $keys {
    foreach j $keys {
      if {$i != $j &&\
          [string range $arr($i) [expr {[string length $arr($i)] - $k}] end] ==\
          [string range $arr($j) 0 [expr {$k - 1}]]} {
        puts "$i $j"
      }
    }
  }

}

# -------------------------------------------------------------------

proc IEV {} {

  set input {1 0 0 1 0 1}
  set le [llength $input]
  set off 2 ; set prbs [calPairs 3]

  set s 0
  for {set i 0} {$i < $le} {incr i} {
    set s [expr {$s + [lindex $input $i] *\
                      [lindex $prbs $i 1]}]
  }

  puts [expr {$s * $off}]

}

# -------------------------------------------------------------------

proc LCSM {} {
}

# -------------------------------------------------------------------

#DNA
#RNA
#REVC
#FIB
#GC
#HAMM
#IPRB
#PROT
#SUBS
#CONS
#FIBD
#GRPH
#IEV
LCSM

# -------------------------------------------------------------------

