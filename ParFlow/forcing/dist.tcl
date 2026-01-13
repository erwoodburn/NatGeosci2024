lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*


pfset Process.Topology.P        8
pfset Process.Topology.Q        8
pfset Process.Topology.R        1

pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

# Define the number of grid blocks in domain
pfset ComputationalGrid.NX               150
pfset ComputationalGrid.NY               170
pfset ComputationalGrid.NZ               1

# Define size of domain grid -- length units equals hydraulic conductivity -- m/h
pfset ComputationalGrid.DX              100.0
pfset ComputationalGrid.DY              100.0
pfset ComputationalGrid.DZ               2.0




set name "NLDAS"
for {set t 0} {$t <= 8759} {incr t 1} {

	#set var [list "APCP" "DLWR" "DSWR" "Press" "SPFH" "UGRD" "VGRD"]
	set var [list "Temp"]

	foreach v $var {
	set filein [format "%s.%s.%06d.pfb" $name $v $t]

	pfdist $filein
	}

puts $t

#set filein [format "%s.%05d.pfb" $name $t]
#puts $filein
#set cond [pfload $filein]
#pfsave $cond -silo [format "%s.%06d.silo" $name $t]

#pfdelete $cond
}
