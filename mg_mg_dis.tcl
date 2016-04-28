# calculating the out of plane distance of mg in the final structure

# matrix order

proc mg {} {

	# GETTING THE COORDINATES FROM THE TRAJECTORY FILE TO GET THE GOOD ENSEMBLE OF VALUES

	set f [open "frame.mdcrd" "r"]
	set data [read $f]
	close $f

	set num_cla [input2]
	set num_per_res [input3]

	# MG COORDINATES

	set atom_pos [input10]
	set elem ""

	set start [expr { 3 + $atom_pos }]
	set end [expr { $start + (3*$num_cla*$num_per_res) }]
	set count $start
	set res 1
	for {set i $start} {$i < $end} {incr i [expr { 3*$num_per_res }]} {
		set coordx [lindex $data $i]
		set coordy [lindex $data [expr { $i + 1 }]]
		set coordz [lindex $data [expr { $i + 2 }]]
		set elem2 [list $coordx $coordy $coordz $res]
		set elem [linsert $elem end $elem2]
		incr res
	}
	return $elem
}

proc com {} {

		# GETTING THE COORDINATES FROM THE TRAJECTORY FILE TO GET THE GOOD ENSEMBLE OF VALUES

	set f [open "frame.mdcrd" "r"]
	set data [read $f]
	close $f

	set g [open "[input4]" "r"]
	set data1 [read $g]
	close $g

	set k1 0

	# READING THE MASS FROM THE PRMTOP FILE

	while { [lindex $data1 $k1] != "MASS" } {
		incr k1
	}
	incr k1 2
	set k2 $k1

	set num_cla [input2]
	set num_per_res [input3]
	set num_frames 100

	set nres_atoms [expr { $num_cla * $num_per_res }]

	# CENTRE OF MASS COORDINATES

	set atom_pos [input8]
	set elem ""

	set start [expr { 3 + $atom_pos }]
	set end [expr { $start + (3*$num_cla*$num_per_res) }]
	set k1 $k2
	set comx 0.0
	set comy 0.0
	set comz 0.0
	set total_mass 0.0
	set count $start
	set res 1
	for {set i $start} {$i < $end} {incr i 3} {
		set mass [lindex $data1 $k1]
		set total_mass [expr {$total_mass + $mass}]
		incr k1
		set coordx [lindex $data $i]
		set coordy [lindex $data [expr { $i + 1 }]]
		set coordz [lindex $data [expr { $i + 2 }]]
		set comx [expr { $comx + ($mass*$coordx) }]
		set comy [expr { $comy + ($mass*$coordy) }]
		set comz [expr { $comz + ($mass*$coordz) }]
		if { [expr { $count % (3*$num_per_res) }] == 0 } {
			set comx [expr { $comx / $total_mass }]
			set comy [expr { $comy / $total_mass }]
			set comz [expr { $comz / $total_mass }]
			set elem2 [list $comx $comy $comz $res]
			set elem [linsert $elem end $elem2]
			set comx 0.0
			set comy 0.0
			set comz 0.0
			set total_mass 0.0
			set count $start
			incr res
			set k1 $k2
		} else { 
			incr count 3
		}		
	}
	return $elem
}
	
proc mg_pair {} {
	
	# IN THIS PROCEDURE WE WILL CALCULATE THE NEAREST NEIGBOUR OF EACH MG AND CALCULATE THE VORNOI POLYHEDRA OF EACH ATOM IN THAT RESIDUE WRT TO CENTRAL MG OF THE RESIDUE UNDER CONSIDERATION

	package require math::linearalgebra 

	set n_cla [input2]

	set num [open "[input4]" "r"]	
	set natom [read $num]
	close $num

	set k 0
	while { [lindex $natom $k] != "%FORMAT(10I8)" } {
		incr k
	}
	incr k
	set n_atoms [lindex $natom $k]
	set num_at_pres [input3]

	for {set i 0} {$i <= $n_cla} {incr i} {
		set d($i) 0.0
		set l($i) ""
		set n_dl($i) ""
		for {set j 0} {$j <= $n_cla} {incr j} {
			set dpair($i,$j) 0.0
		}
	}
	
	set m [open "res_res_distance" "w"]

	set mp [open "res_distance_all" "w"]

	set start_frame [input5]
	set end_frame [input6]
	set step [input7]
	set n_frames [expr { (($end_frame - $start_frame)/$step) + 1 }]

	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {

		set mgc [open "mg_coord.$fr" "w"]

		# EXECUTING CPPTRAJ

		set inp [open "ai.dat" "w"]
		puts $inp "trajin [cpptraj_traj] $fr $fr"
		puts $inp "trajout frame.mdcrd mdcrd nobox"
		puts $inp "trajout frame.rst rst"
		puts $inp "trajout frame.pdb pdb"
		puts $inp "go"
		close $inp

		set inp [open "ai.dat" "r"]

		exec cpptraj -p [input4] -i ai.dat 

		close $inp

		set h [open "frame.rst" "r"]
		set data2 [read $h]
		close $h

		set end 0

		while { $end < [llength $data2] } {	
			incr end
		}

		set box_x "[lindex $data2 [expr { $end - 3 }]]"
		set box_y "[lindex $data2 [expr { $end - 2 }]]"
		set box_z "[lindex $data2 [expr { $end - 1 }]]"

		puts $mgc " # $box_x $box_y $box_z"

		puts "				**** FRAME :: $fr"
		set f [open "$fr" "w"]
		if { [input9] == 0 } {
			set coord [com]
		} else { 
			set coord [mg]
		}
		puts $f "{$coord}"
		close $f

		# CALCULATING THE DISTANCE OF EACH RESIDUE FROM OTHER RESIDUES

		set f [open "$fr" "r"]
		set data [read $f]
		close $f
		for {set i 1} {$i <= $n_cla} {incr i} {
			set res1x [lindex $data 0 [expr { $i - 1 }] 0]
			set res1y [lindex $data 0 [expr { $i - 1 }] 1]
			set res1z [lindex $data 0 [expr { $i - 1 }] 2]
		
			puts $mgc "$res1x $res1y $res1z $i"

			set dis_temp 1000.0

			for {set j 1} {$j <= $n_cla} {incr j} {
				if {$j != $i } {
					set res2x [lindex $data 0 [expr { $j - 1 }] 0]
					set res2y [lindex $data 0 [expr { $j - 1 }] 1]
					set res2z [lindex $data 0 [expr { $j - 1 }] 2]

					# CALCULATING THE DISTANCE BETWEEN THE TWO RESIDUES

					set difx [expr { $res1x-$res2x }]
					set difx [expr { $difx * $difx }]
					set dify [expr { $res1y-$res2y }]
					set dify [expr { $dify * $dify }]
					set difz [expr { $res1z-$res2z }]
					set difz [expr { $difz * $difz }]

					set dis [expr { sqrt($difx + $dify + $difz) }]

					set dpair($i,$j) [expr { $dpair($i,$j) + $dis }]

					if { $dis < $dis_temp } {
						set dis_temp $dis
						set pair $j
					}			
				}
			}
			set d($i) [expr { $d($i) + $dis_temp }]	
			if { [expr { $fr % 10 }] == 0 } {
				set l($i) [linsert $l($i) end $pair]
				set n_dl($i) [linsert $n_dl($i) end $dis_temp]
			}
		}
	close $mgc
	} 
	# AVERAGING

	for {set i 1} {$i <= $n_cla} {incr i} {
		set value [expr { $d($i) / $n_frames }]
		puts $m " $i $value { $l($i) $n_dl($i) }"
		for {set j 1} {$j <= $n_cla} {incr j} {
			if { $j != $i } {
				set value [expr {$dpair($i,$j) / $n_frames }]
				puts $mp " $i $j $value "
			}
		}
	}
	close $m
}

proc delete_files {} {
	
	# THIS PROCEDURE WILL DELETE ALL THE FILES REQUIRED DURING CALCULATIONS ONLY

	for {set i 0} {$i < [input2]} {incr i} {
		file delete $i
		file delete $i.avg
	}
	
	set start_frame [input5]
	set end_frame [input6]
	set step [input7]
	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {
		file delete $fr
	}

	file delete dummy
	file delete dummy_ester
	file delete frame.pdb 
	file delete frame.rst 
	file delete frame.mdcrd
	file delete ai.dat

}

proc input1 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 0]
}

proc input2 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 1]
}

proc input3 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 2]
}

proc input4 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 3]
}

proc input5 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 5 0]
}

proc input6 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 5 1]
}

proc input7 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 5 2]
} 

proc cpptraj_traj {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 6]
}

proc input8 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 4]
}

proc input9 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 7 0]
}

proc input10 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 7 1]
}

mg_pair
delete_files





















