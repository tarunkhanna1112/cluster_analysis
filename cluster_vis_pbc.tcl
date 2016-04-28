set inp [open "input" "r"]
set dinp [read $inp]
close $inp

set n_res [lindex $dinp 1] 

set start_frame [lindex $dinp 5 0]
set end_frame [lindex $dinp 5 1]
set total_frame [expr { $end_frame - $start_frame }]
set step [lindex $dinp 5 2]

set R [lindex $dinp 8]

set first [lindex $dinp 9 0]
set last [lindex $dinp 9 1]

set jj 1

set emp [open "empty_file" "w"]
puts $emp "0 0 0"
close $emp

for {set j $start_frame } {$j <= $end_frame} {incr j $step} {
	puts "			**** STEP :: $j ****"
	set f [open "mg_coord.$j" "r"]
	set data [read $f]
	close $f

	set box_x [lindex $data 1]
	set box_y [lindex $data 2]
	set box_z [lindex $data 3]

	# GETTING THE NEAREST NEIGHBOURS TO EACH RESIDUE

	for {set i 0} {$i < $n_res} {incr i} {
		set m 0
		set k 4
		set k1 [expr { ($i * 4 ) + 4 }]
		set x1 [lindex $data $k1]
		set y1 [lindex $data [expr { $k1 + 1 }]]
		set z1 [lindex $data [expr { $k1 + 2 }]]
	 
		while { $k < [llength $data]} {
			if { $k != $k1 } {

				# ALL NEAREST PERIODIC IMAGES

				for {set px -$box_x} {$px <= $box_x} {set px [expr { $px + $box_x }]} {
					for {set py -$box_y} {$py <= $box_y} {set py [expr { $py + $box_y }]} {
						for {set pz -$box_z} {$pz <= $box_z} {set pz [expr { $pz + $box_z }]} {

							set x2 [expr { [lindex $data $k] + $px }]
							set y2 [expr { [lindex $data [expr { $k + 1 }]] + $py }]
							set z2 [expr { [lindex $data [expr { $k + 2 }]] + $pz }]
		
							set delx [expr { $x1 - $x2 }]
							set delx [expr { $delx*$delx }]

							set dely [expr { $y1 - $y2 }]
							set dely [expr { $dely*$dely }]

							set delz [expr { $z1 - $z2 }]
							set delz [expr { $delz*$delz }]

							set del [expr { sqrt($delx + $dely + $delz) }]

							if { $del < [expr { 2 * $R }] } {
								set n($i,$m) $k
								incr m
							}
						}
					}
				}
			}
		incr k 4
		}
		set n($i,$m) -99
	}

	# CALCULATING THE GEOMETRIC CENTRE OF THE SPHERE

	set g [open "plot.gp" "w"]
	set label [expr { ($j * 50.0) / $total_frame }]
	set label [format "%.3f" $label]
	puts $g "set label \"TIME $label nsec\" at [expr { $box_x + 5.0 }],[expr { $box_y + 5.0 }],[expr { $box_z + 5.0 }] font \"LiberationSans-Regular,24\""
	puts $g "set parametric"
	puts $g "set isosamples 10.0"
	puts $g "set hidden3d"
	#puts $g "set view 0,0"
	puts $g "set xlabel 'X'"
	puts $g "set ylabel 'Y'"
	puts $g "set zlabel 'Z'"
	puts $g "set terminal jpeg size 1920,1024"
	puts $g "r = 5.0"
	puts $g "set urange \[0:2*pi\]"
	puts $g "set vrange \[-pi/2:pi/2\]"

	for {set i 0} {$i < $n_res} {incr i} {
		set m 0
		#puts "			**** PLOTING NEIGHBOURS OF RESIDUE $i ****"
		set p [expr { ($i * 4) + 4 }]
		set comx [lindex $data $p]
		set comy [lindex $data [expr { $p + 1 }]]
		set comz [lindex $data [expr { $p + 2 }]]
		while { $n($i,$m) != -99 } {
			set comx [expr { $comx + [lindex $data $n($i,$m)] }]
			set comy [expr { $comy + [lindex $data [expr { $n($i,$m) + 1 }]] }]
			set comz [expr { $comz + [lindex $data [expr { $n($i,$m) + 2 }]] }]
			incr m
		}
		incr m
		if { $m > 0 } {
			set cmx($i) [expr { $comx / $m }]
			set cmy($i) [expr { $comy / $m }]
			set cmz($i) [expr { $comz / $m }]

			puts $g "fx$i\(v,u) = r*cos(v)*cos(u) + $cmx($i)"
			puts $g "fy$i\(v,u) = r*cos(v)*sin(u) + $cmy($i)"
			puts $g "fz$i\(v)   = r*sin(v) + $cmz($i)"
		}
	}

	set m 0
	#if { $n($first,$m) != -99 } {
		puts $g "set output \"$jj.jpeg\""
		incr jj
	#}

	# LABELLING THE PLOT

	for {set i $first} {$i < $last} {incr i} {
		set m 0
		set p [expr { ($i * 4) + 4 }]
		set comx [lindex $data $p]
		set comy [lindex $data [expr { $p + 1 }]]
		set comz [lindex $data [expr { $p + 2 }]]
		#puts $g "set label \"[expr {$i + 1 }]\,\" at 0,0,$comz font \"LiberationSans-Regular,12\""
		set lab "[expr {$i + 1}]"
		while { $n($i,$m) != -99 } {
			set neig [lindex $data [expr { $n($i,$m) + 3 }]]
			set comx [lindex $data $n($i,$m)]
			set comy [lindex $data [expr { $n($i,$m) + 1 }]] 
			set comz [lindex $data [expr { $n($i,$m) + 2 }]] 
			#puts $g "set label \"$neig\,\" at 0,0,$comz font \"LiberationSans-Regular,12\""
			set lab [linsert $lab end  $neig]
			incr m
		}
		puts $g "set label \"\[$lab\]\" at $comx,$comy,$comz font \"LiberationSans-Regular,16\""
	}
	for {set i $first} {$i < $last} {incr i} {
		#puts "		** NEIGHBORS OF RESIDUE [expr { $i + 1 }] ARE: ** " 
		set m 0
		if { $n($i,$m) != -99 } {
			if { $i == $first && $first != [expr { $last - 1 }]} {
				puts $g "splot fx$i\(v,u),fy$i\(v,u),fz$i\(v) notitle,\\"
			} elseif { $i == $first && $first == [expr { $last - 1 }] } {
					puts $g "splot fx$i\(v,u),fy$i\(v,u),fz$i\(v) notitle,'mg_coord.$j' pt 7 ps 2 notitle"
			} elseif { $i < [expr { $last - 1 }] } {
					puts $g "fx$i\(v,u),fy$i\(v,u),fz$i\(v) notitle,\\"
			} elseif { $first != [expr { $last - 1 }] } {
					puts $g "fx$i\(v,u),fy$i\(v,u),fz$i\(v) notitle,'mg_coord.$j' pt 7 ps 2 notitle"
			}	
		} elseif { $i == $first && $i != [expr { $last - 1 }]} {
				puts $g "splot 'empty_file' u 1:2:3 w p notitle,\\"
		} elseif { $i == $first && $i == [expr { $last - 1 }] } {
				puts $g "splot 'empty_file' u 1:2:3 w p"
		} elseif { $i == [expr { $last - 1 }] && $i != $first} {
				#puts $g "'empty_file' u 1:2:3 w p"
				puts $g "'mg_coord.$j' pt 7 ps 2 notitle"
		}
		while { $n($i,$m) != -99 } {
			set neig [lindex $data [expr { $n($i,$m) + 3 }]]
			puts "$neig"
			incr m
		}
	}
	#puts $g "pause 2"
	close $g
	exec gnuplot plot.gp	
}

set i $jj

exec mkdir -p images_all

for {set i 1} {$i < $jj} {incr i} {
	exec mv $i.jpeg ./images_all/
}

file delete empty_file
file delete input

for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {
	puts "$fr"
	file delete mg_coord.$fr
}


# ffmpeg -framerate 3 -i "%d.jpeg" out3.mp4
	

	
