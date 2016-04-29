puts ""
puts ""
puts "				************************************************************************"
puts "				CLUSTER ANALYSIS BASED ON THE CLOSEST DISTANCE BETWEEN ATOMS OR RESIDUES"
puts "				************************************************************************"
puts ""
puts "						WRITTEN BY: TARUN KHANNA"
puts "						IMPERIAL COLLEGE LONDON,U.K."
puts ""
puts ""

puts "		NOTE: FOR RUNNING THIS CODE YOU NEED THE CPPTRAJ FUNCTIONALITY OF AMBERTOOLS"
puts ""
puts ""

puts "ENTER THE AMBER PRMTOP"
set opt1 [gets stdin]
exec ls $opt1

puts ""

puts "ENTER THE AMBER TRAJECTORY"
set opt2 [gets stdin]
exec ls $opt2

puts ""

puts "ENTER THE STARTING FRAME NUMBER"
set opt21 [gets stdin]

puts ""

puts "ENTER THE LAST FRAME NUMBER"
set opt22 [gets stdin]

puts ""

puts "ENTER THE FRAME STEP SIZE"
set opt23 [gets stdin]

puts ""


puts "ENTER THE RESIDUE NAME FOR WHICH YOU WANT TO DO THE CLUSTER ANALYSIS"
set opt3 [gets stdin]

puts ""

puts "ENTER '0' TO DO THE CLUSTER ANALYSIS BASED ON CENTER OF MASS OR '1' TO BE BASED ON CERTAIN ATOM TYPE"
set opt4 [gets stdin]

puts ""

if { $opt4 == 1 } {
	puts "ENTER THE ATOM TYPE (PDB NAME OF THE ATOM)"
	set opt5 [gets stdin]
	
	puts ""
}

set f [open "input" "w"]

puts $f "INPUT"

set g [open "$opt1" "r"]
set data [read $g]
close $g

if { $opt4 == 1 } {
	set k 0
	while { [lindex $data $k] != "ATOM_NAME" } {
		incr k
	}
	incr k 2
	
	set atompos 0
	while { [lindex $data $k] != "$opt5" } {
		if { $k > [llength $data] } {
			puts "		## ERROR :: ATOM TYPE NOT FOUND IN THE PRMTOP"
			exit
		}
		incr k
		incr atompos
	}
}
set k 0

while { [lindex $data $k] != "RESIDUE_LABEL" } {
	incr k
}
incr k 2

set fres $k

while { [lindex $data $k] != "$opt3" } {
	if { $k > [llength $data] } {
		puts "		## ERROR :: RESIDUE NAME NOT FOUND IN THE PRMTOP"
		exit
	}
	incr k
}
set yres $k
set respos [expr { $yres - $fres }]

set nres 0
while { [lindex $data $k] == "$opt3" } {
	incr nres
	incr k
}
puts $f "$nres"

set k 0

while { [lindex $data $k] != "RESIDUE_POINTER" } {
	incr k
}
incr k 2 

set k [expr { $k + $respos }]

set natoms [expr { [lindex $data [expr { $k + 1 }]] - [lindex $data $k] }]

puts $f "$natoms"

puts $f "$opt1"

puts $f "$respos"

puts $f "{$opt21 $opt22 $opt23}"

puts $f "$opt2"

if { $opt4 == 1 } {
	puts $f "{ 1 $atompos }"
} else {
	puts $f "0"
}

puts "ENTER THE RADIUS OF THE SPHERE WHICH DEFINES A CLUSTER OF MOLEULES"
set opt6 [gets stdin]
puts ""

puts $f "$opt6"

puts "ENTER THE NUMBER OF RESIDUES YOU WANT TO SEE ON THE SAME GRAPH"
puts "NOTE: WITH ALL ($nres) THE GRAPH MIGHT LOOK MESSY"
set opt7 [gets stdin]
puts ""

if { $opt7 == $nres } {
	puts $f "{0 $nres}"
} else { 
	puts "ENTER THE STARTING RESIDUE"
	set opt8 [gets stdin]
	puts ""
	puts "ENTER THE END RESIDUE"
	set opt9 [gets stdin]
	puts ""

	puts "NOTE:: ONLY FOR RESIDUES $opt8 TO $opt9 THE CLUSTER ANALYSIS WILL BE DONE"
	puts $f "{[expr { $opt8 - 1 }] $opt9}"
}



close $f

exec tclsh mg_mg_dis.tcl | tee log1
exec tclsh cluster_vis_pbc.tcl | tee log2

puts ""
puts ""
puts "				### DONE"
puts ""
puts "				### FILE : 'res_res_distance' CONTAINS THE NEAREST NEIGHBOUR OF EACH RESIDUE"
puts "				### FILE : 'res_distance_all' CONTAINS THE DISTANCE OF EACH RESIDUE FROM OTHER RESIDUE IN THE MASK"
puts "				### FOLDER 'images_all' CONTAINS THE IMAGES OF THE CLUSTER"
puts "				### USE ffmpeg TO CONVERT THESE IMAGES INTO A VIDEO"
puts ""
puts "				### (SYNTAX: ffmpeg -framerate 3 -i \"%d.jpeg\" out3.mp4"
puts ""
puts ""


