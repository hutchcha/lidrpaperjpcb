# This is a set of scripts for finding the lipid order parameter from
# coordinate data.  The lipid order parameter is defined as
#     S_CD = < 3/2 cos^2 w - 1/2 >
# where w is the angle made by the bond between an aliphatic carbon
# and a hydrogen. 

# Author: Justin Gullingsrud
# Email:  justinrocks@gmail.com


# Find order param for the c3 (Sn1 chain) tail.
#proc orderparam-c3 { result { resname $x all } } {
  #upvar result arr
set res [atomselect top "resname POPC POPE SAPI PSM POPS"]
set resn [lsort -u [$res get resname]]
set n [molinfo top get numframes]
set memcen [lindex [measure center $res] 2]

set bot [[atomselect top "name P and z < $memcen"] get resid]
set fourth [expr {$n / 4}]
set beg 0
set end $fourth
set n 1
puts $n
while {$n < 5 } {
  puts $n
  foreach x $resn {
  set top [lsort -u [[atomselect top "resname $x and name P and z > $memcen"] get resid]]
  puts $x
  set outfile [open "rheborderparams/sn1orderparam${x}${n}top.dat" w]
    for { set i 2 } { $i <= 20 } { incr i } {
      puts $i
      set cp [atomselect top "(resname $x and resid $top) and name C3$i"]
      if { [$cp num] == 0 } {
        puts "skipping $i"
        continue
      }
      set hx [atomselect top "(resname $x and resid $top) and name H${i}X"]
      set hy [atomselect top "(resname $x and resid $top) and name H${i}Y"]
      set hz [atomselect top "(resname $x and resid $top) and name H${i}Z"]

      set sum 0.0
      set nh 0
      set nres [$cp num]
      for { set frame $beg } { $frame < $end } { incr frame } {
        $cp frame $frame
        $hx frame $frame
        $hy frame $frame
        $hz frame $frame
        set cpx [$cp get x]
        set cpy [$cp get y]
        set cpz [$cp get z]
        set hxx [vecsub $cpx [$hx get x]]
        set hxy [vecsub $cpy [$hx get y]]
        set hxz [vecsub $cpz [$hx get z]]
        foreach dx $hxx dy $hxy dz $hxz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
        if { [$hy num] != 0 } {
          set hyx [vecsub $cpx [$hy get x]]
          set hyy [vecsub $cpy [$hy get y]]
          set hyz [vecsub $cpz [$hy get z]]
          foreach dx $hyx dy $hyy dz $hyz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }
        if { [$hz num] != 0 } {
          set hzx [vecsub $cpx [$hz get x]]
          set hzy [vecsub $cpy [$hz get y]]
          set hzz [vecsub $cpz [$hz get z]]
          foreach dx $hzx dy $hzy dz $hzz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

      }
      set arr($i) [expr {-1.5*$sum/$nh + 0.5}]
    puts $outfile "${i} [expr {-1.5*$sum/$nh + 0.5}]"
    }
    close $outfile
    
    

    }
  
  set beg [expr {$beg + $fourth}]
  set end [ expr {$end + $fourth}]
  set n [expr {$n + 1}]
}
foreach x $resn {
  exec awk -f columnscript.awk rheborderparams/sn1orderparam${x}1top.dat rheborderparams/sn1orderparam${x}2top.dat rheborderparams/sn1orderparam${x}3top.dat rheborderparams/sn1orderparam${x}4top.dat > rheborderparams/sn1orderparam${x}top.dat
  
  for {set i 1} {$i <= 4} {incr i} {
    exec rm rheborderparams/sn1orderparam${x}${i}top.dat
  }
}
set beg 0
set end $fourth
set n 1

while {$n < 5 } {
  puts $n
  foreach x $resn {
  set bot [lsort -u [[atomselect top "resname $x and name P and z < $memcen"] get resid]]
  puts "$x $n"
  set outfile [open "rheborderparams/sn1orderparam${x}${n}bot.dat" w]
    for { set i 2 } { $i <= 20 } { incr i } {
      puts $i
      set cp [atomselect top "(resname $x and resid $bot) and name C3$i"]
      if { [$cp num] == 0 } {
        puts "skipping $i"
        continue
      }
      set hx [atomselect top "(resname $x and resid $bot) and name H${i}X"]
      set hy [atomselect top "(resname $x and resid $bot) and name H${i}Y"]
      set hz [atomselect top "(resname $x and resid $bot) and name H${i}Z"]

      set sum 0.0
      set nh 0
      set nres [$cp num]
      for { set frame $beg } { $frame < $end } { incr frame } {
        $cp frame $frame
        $hx frame $frame
        $hy frame $frame
        $hz frame $frame
        set cpx [$cp get x]
        set cpy [$cp get y]
        set cpz [$cp get z]
        set hxx [vecsub $cpx [$hx get x]]
        set hxy [vecsub $cpy [$hx get y]]
        set hxz [vecsub $cpz [$hx get z]]
        foreach dx $hxx dy $hxy dz $hxz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
        if { [$hy num] != 0 } {
          set hyx [vecsub $cpx [$hy get x]]
          set hyy [vecsub $cpy [$hy get y]]
          set hyz [vecsub $cpz [$hy get z]]
          foreach dx $hyx dy $hyy dz $hyz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }
        if { [$hz num] != 0 } {
          set hzx [vecsub $cpx [$hz get x]]
          set hzy [vecsub $cpy [$hz get y]]
          set hzz [vecsub $cpz [$hz get z]]
          foreach dx $hzx dy $hzy dz $hzz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

      }
      set arr($i) [expr {-1.5*$sum/$nh + 0.5}]
    puts $outfile "${i} [expr {-1.5*$sum/$nh + 0.5}]"
    }
    close $outfile
    
    

    }
  
  set beg [expr {$beg + $fourth}]
  set end [ expr {$end + $fourth}]
  set n [expr {$n + 1}]
}
foreach x $resn {
  exec awk -f columnscript.awk rheborderparams/sn1orderparam${x}1bot.dat rheborderparams/sn1orderparam${x}2bot.dat rheborderparams/sn1orderparam${x}3bot.dat rheborderparams/sn1orderparam${x}4bot.dat > rheborderparams/sn1orderparam${x}bot.dat
  
  for {set i 1} {$i <= 4} {incr i} {
    exec rm rheborderparams/sn1orderparam${x}${i}bot.dat
  }
}

set beg 0
set end $fourth
set n 1
puts "sn2 ${n}"
while {$n < 5} {
  foreach x $resn {
    set top [lsort -u [[atomselect top "chain M and resname $x and name P and z > $memcen"] get resid]]
    set outfile [open rheborderparams/sn2orderparam${x}${n}top.dat w]
    if {$x == "SAPI"} {
      set c 20
    } else {
      set c 18
    }
    for { set i 2 } { $i <= $c } { incr i } {
      puts "$i $x"
      set cp [atomselect top "(resname $x and resid $top)and name C2$i"]
      if { [$cp num] == 0 } {
        puts "skipping $i"
      }
      set hx [atomselect top "(resname $x and resid $top) and name H${i}R"]
      set hy [atomselect top "(resname $x and resid $top) and name H${i}S"]
      set hz [atomselect top "(resname $x and resid $top) and name H${i}T"]
      set h9 [atomselect top "(resname $x and resid $top) and name H${i}1"]
      set sum 0.0
      set nh 0
      set nres [$cp num]
      for { set frame $beg } { $frame < $end } { incr frame } {
        $cp frame $frame
        $hx frame $frame
        $hy frame $frame
        $hz frame $frame
        $h9 frame $frame
        set cpx [$cp get x]
        set cpy [$cp get y]
        set cpz [$cp get z]

        if { [$hx num] != 0 } {
          set hxx [vecsub $cpx [$hx get x]]
          set hxy [vecsub $cpy [$hx get y]]
          set hxz [vecsub $cpz [$hx get z]]
          foreach dx $hxx dy $hxy dz $hxz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

        if { [$hy num] != 0 } {
          set hyx [vecsub $cpx [$hy get x]]
          set hyy [vecsub $cpy [$hy get y]]
          set hyz [vecsub $cpz [$hy get z]]
          foreach dx $hyx dy $hyy dz $hyz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

        if { [$hz num] != 0 } {
          set hzx [vecsub $cpx [$hz get x]]
          set hzy [vecsub $cpy [$hz get y]]
          set hzz [vecsub $cpz [$hz get z]]
          foreach dx $hzx dy $hzy dz $hzz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

        if { [$h9 num] != 0 } {
          set h9x [vecsub $cpx [$h9 get x]]
          set h9y [vecsub $cpy [$h9 get y]]
          set h9z [vecsub $cpz [$h9 get z]]
          foreach dx $h9x dy $h9y dz $h9z {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

      }
      set arr($i) [expr {-1.5*$sum/$nh + 0.5}]
      
      puts $outfile "${i} [expr {-1.5*$sum/$nh + 0.5}]"
    }
    close $outfile
  }
  set beg [expr {$beg + $fourth}]
  set end [ expr {$end + $fourth}]
  set n [expr {$n + 1}]
}

foreach x $resn {
  exec awk -f columnscript.awk rheborderparams/sn2orderparam${x}1top.dat rheborderparams/sn2orderparam${x}2top.dat rheborderparams/sn2orderparam${x}3top.dat rheborderparams/sn2orderparam${x}4top.dat > rheborderparams/sn2orderparam${x}top.dat
  
  for {set i 1} {$i <= 4} {incr i} {
    exec rm rheborderparams/sn2orderparam${x}${i}top.dat
  }
}
set beg 0
set end $fourth
set n 1
puts "sn2 ${n}"
while {$n < 5} {
  foreach x $resn {
    set bot [lsort -u [[atomselect top "chain M and resname $x and name P and z < $memcen"] get resid]]
    set outfile [open rheborderparams/sn2orderparam${x}${n}bot.dat w]
    if {$x == "SAPI"} {
      set c 20
    } else {
      set c 18
    }
    for { set i 2 } { $i <= $c } { incr i } {
      puts "$i $x"
      set cp [atomselect top "(resname $x and resid $bot)and name C2$i"]
      if { [$cp num] == 0 } {
        puts "skipping $i"
      }
      set hx [atomselect top "(resname $x and resid $bot) and name H${i}R"]
      set hy [atomselect top "(resname $x and resid $bot) and name H${i}S"]
      set hz [atomselect top "(resname $x and resid $bot) and name H${i}T"]
      set h9 [atomselect top "(resname $x and resid $bot) and name H${i}1"]
      set sum 0.0
      set nh 0
      set nres [$cp num]
      for { set frame $beg } { $frame < $end } { incr frame } {
        $cp frame $frame
        $hx frame $frame
        $hy frame $frame
        $hz frame $frame
        $h9 frame $frame
        set cpx [$cp get x]
        set cpy [$cp get y]
        set cpz [$cp get z]

        if { [$hx num] != 0 } {
          set hxx [vecsub $cpx [$hx get x]]
          set hxy [vecsub $cpy [$hx get y]]
          set hxz [vecsub $cpz [$hx get z]]
          foreach dx $hxx dy $hxy dz $hxz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

        if { [$hy num] != 0 } {
          set hyx [vecsub $cpx [$hy get x]]
          set hyy [vecsub $cpy [$hy get y]]
          set hyz [vecsub $cpz [$hy get z]]
          foreach dx $hyx dy $hyy dz $hyz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

        if { [$hz num] != 0 } {
          set hzx [vecsub $cpx [$hz get x]]
          set hzy [vecsub $cpy [$hz get y]]
          set hzz [vecsub $cpz [$hz get z]]
          foreach dx $hzx dy $hzy dz $hzz {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

        if { [$h9 num] != 0 } {
          set h9x [vecsub $cpx [$h9 get x]]
          set h9y [vecsub $cpy [$h9 get y]]
          set h9z [vecsub $cpz [$h9 get z]]
          foreach dx $h9x dy $h9y dz $h9z {
            set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            set sum [expr {$sum + $dz*$dz/$norm2}]
          }
          incr nh $nres
        }

      }
      set arr($i) [expr {-1.5*$sum/$nh + 0.5}]
      
      puts $outfile "${i} [expr {-1.5*$sum/$nh + 0.5}]"
    }
    close $outfile
  }
    
    
    # unset hxx [vecsub $cpx [$hx get x]]
    # unset hxy [vecsub $cpy [$hx get y]]
    # unset hxz [vecsub $cpz [$hx get z]]
    # unset cpx [$cp get x]
    # unset cpy [$cp get y]
    # unset cpz [$cp get z]
    # unset hx 
    # unset hy 
    # unset hz 
    # unset h9 
  set beg [expr {$beg + $fourth}]
  set end [ expr {$end + $fourth}]
  set n [expr {$n + 1}]
}

foreach x $resn {
  exec awk -f columnscript.awk rheborderparams/sn2orderparam${x}1bot.dat rheborderparams/sn2orderparam${x}2bot.dat rheborderparams/sn2orderparam${x}3bot.dat rheborderparams/sn2orderparam${x}4bot.dat > rheborderparams/sn2orderparam${x}bot.dat
  
  for {set i 1} {$i <= 4} {incr i} {
    exec rm rheborderparams/sn2orderparam${x}${i}bot.dat
  }
}


exec python3.10 /home/chutchins@uthouston.edu/projects/rheb/anchoranalysis/plotordpar.py
