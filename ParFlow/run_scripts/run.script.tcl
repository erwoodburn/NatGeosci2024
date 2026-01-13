### 100m East River Model 
### This version compiled after parameter determination in Parameter Sensitivity Study, February 2018 ###
### Permeability values represent 'effective K' (sin(theta)*K) to incorporate slope ###
### mannings n values are chosen by landcover, stream cells have their own values however, then empirical*50 to best match behavior ###

set tcl_precision 17
set runname [lindex $argv 0]

#---------------------------------------------------------
# Import the ParFlow TCL package
#---------------------------------------------------------

lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*
pfset FileVersion 4

#---------------------------------------------------------
# Process Topology -- Sets the parallel process topology - P*Q*R must equal nodes*ntasks per node
#---------------------------------------------------------
# DX 
pfset Process.Topology.P        8
# DY
pfset Process.Topology.Q        8
# DZ - usually leave as one processor, don't want to split up dz
pfset Process.Topology.R        1

set  nproc [expr [pfget Process.Topology.P]*[pfget Process.Topology.Q]*[pfget Process.Topology.R]]

#-----------------------------------------------------------------------------
### Initialize variables for restart
###-----------------------------------------------------------------------------
set log_file [open "restart.log" a 0600]
set rst_file clm_restart.tcl

if { [file exists $rst_file] == 1} {
        source clm_restart.tcl
#        if { $istep == 8760 } {
#                puts "Last run completed full year, starting new run from beginning"
#                set systemTime [clock seconds]

 #               puts $log_file "istep = pfstep = 8760, starting a new run from the begining"
 #               puts $log_file  "The time is: [clock format $systemTime -format %H:%M:%S]"
 #               puts $log_file "The date is: [clock format $systemTime -format %D]"
 #               close $log_file
        #Set the ParFlow counter to 0 and clm to 1
 #               set istep 0
 #               set clmstep [expr int($istep+1)]

        #copy in the intial pressure file and set the initial conditions
#                set ip "$runname.out.press.08760.pfb"

        #Copy inputs into the run directory
 #               file copy -force "../pf_inputs/ER_slope_x.pfb" .
 #               file copy -force "../pf_inputs/ER_slope_y.pfb" .
 #               file copy -force "../pf_inputs/ER_subsurface_indicator.pfb" .
 #               file copy -force "../pf_inputs/ER_mannings_100m.pfb" .
 #               file copy -force "../clm_inputs/drv_vegm.dat" .
 #               file copy -force "../clm_inputs/drv_vegp.dat" .
 #               file copy -force "../clm_inputs/drv_clmin.startup.dat" .
 #               file copy -force "../clm_inputs/drv_clmin.restart.dat" .

        #Copy the clm driver file with the restart flags turned off
#                file copy -force "./drv_clmin.startup.dat" ./drv_clmin.dat
#        } else {

        puts $istep
        set clmstep [expr int($istep+1)]
	#ERW 8-27-18 modified ip:
	set ip_time [expr int($istep-1)]

        #Record the restart
        set systemTime [clock seconds]
        puts $log_file  "---------------------------------------------------------"
        puts $log_file  "The time is: [clock format $systemTime -format %H:%M:%S]"
        puts $log_file "The date is: [clock format $systemTime -format %D]"
        puts $log_file "Restarting istep = pfstep = $istep"
        close $log_file

        #Set the intial pressure file
#       set ip [format $runname.out.press.%05d.pfb $istep]
        set ip [format $runname.out.press.%05d.pfb $ip_time]


        #Copy the clm driver with restart flags
        file copy -force ./drv_clmin.restart.dat ./drv_clmin.dat

#ERW changed on 2-23-18

        #Copy the restart files to the correct time stamp
#        for {set i 0 } { $i < $nproc } { incr i 1 } {
#                set fname_rst [format "clm.rst.%05d.$i" $istep]
#                exec cp clm.rst.00000.$i $fname_rst
#        }

        #Copy the kinsol log and the clm log so you have a copy saved
        set fname_kin2 [format $runname.out.%05d.kinsol.log $istep]
        exec cp $runname.out.kinsol.log $fname_kin2

        set fname_clm2 [format CLM.out.clm.%05d.log $istep]
        exec cp CLM.out.clm.log $fname_clm2
#        }

#If there isn't a restart file then start from scratch
} else {
        puts $log_file "CLM restart tcl not found starting a new run from the begining"
        set systemTime [clock seconds]
        puts $log_file  "The time is: [clock format $systemTime -format %H:%M:%S]"
        puts $log_file "The date is: [clock format $systemTime -format %D]"
        puts $log_file "Starting New Run"

        #Copy inputs into the run directory
        file copy -force "pf_inputs/ER_slope_x.pfb" .
        file copy -force "pf_inputs/ER_slope_y.pfb" .
        file copy -force "pf_inputs/ER_subsurface_indicator.pfb" .
        file copy -force "pf_inputs/ER_mannings_100m.pfb" .
        
	#10-4-22 ERW changed to NEON
	#file copy -force "clm_inputs/drv_vegm.dat" .
	file copy -force "clm_inputs/drv_vegm_NEON_ERW_10-4-22.dat" .
        file copy -force "clm_inputs/drv_vegp.dat" .
        file copy -force "clm_inputs/drv_clmin.startup.dat" .
        file copy -force "clm_inputs/drv_clmin.restart.dat" .

        #Set the ParFlow counter to 0 and clm to 1
        set istep 0
        set clmstep [expr int($istep+1)]

        set ip "BaseR8.out.press.08760.pfb"

        #Copy the clm driver file with the restart flags turned off
        file copy -force "./drv_clmin.startup.dat" ./drv_clmin.dat
}

set mannings ER_mannings_100m.pfb
set slopex ER_slope_x.pfb
set slopey ER_slope_y.pfb
set sub ER_subsurface_indicator.pfb

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------

# Locate the origin
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

# Define the number of grid blocks in domain
pfset ComputationalGrid.NX               150 
pfset ComputationalGrid.NY               170
pfset ComputationalGrid.NZ               5

# Define size of domain grid -- length units equals hydraulic conductivity -- m/h
pfset ComputationalGrid.DX		100.0
pfset ComputationalGrid.DY		100.0
pfset ComputationalGrid.DZ		 2.0     

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------

pfset GeomInput.Names                 	"domain_input indi_input"
pfset GeomInput.domain_input.GeomName  		domain
pfset GeomInput.domain_input.InputType  	Box


#---------------------------------------------------------
# Domain Geometry 
#---------------------------------------------------------
pfset Geom.domain.Lower.X                        0.0
pfset Geom.domain.Lower.Y                        0.0
pfset Geom.domain.Lower.Z                        0.0
 
pfset Geom.domain.Upper.X                        15000.0
pfset Geom.domain.Upper.Y                        17000.0
pfset Geom.domain.Upper.Z                        10.0
pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"

#---------------------------------------------------------
# variable dz assignments
#---------------------------------------------------------

pfset Solver.Nonlinear.VariableDz	True
pfset dzScale.GeomNames				domain
pfset dzScale.Type            		nzList
pfset dzScale.nzListNumber       	5

pfset Cell.0.dzScale.Value 10.50
pfset Cell.1.dzScale.Value 4.00
pfset Cell.2.dzScale.Value 0.300
pfset Cell.3.dzScale.Value 0.150
pfset Cell.4.dzScale.Value 0.050

#---------------------------------------------------------
# Indicator Geometry Input
#---------------------------------------------------------
pfset GeomInput.indi_input.InputType    IndicatorField
pfset GeomInput.indi_input.GeomNames    "s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15"
pfset Geom.indi_input.FileName        	$sub

pfset GeomInput.s1.Value        1
pfset GeomInput.s2.Value        2
pfset GeomInput.s3.Value        3
pfset GeomInput.s4.Value        4
pfset GeomInput.s5.Value        5
pfset GeomInput.s6.Value        6
pfset GeomInput.s7.Value        7
pfset GeomInput.s8.Value        8
pfset GeomInput.s9.Value        9
pfset GeomInput.s10.Value       10
pfset GeomInput.s11.Value       11
pfset GeomInput.s12.Value       12
pfset GeomInput.s13.Value       13
pfset GeomInput.g1.Value        21
pfset GeomInput.g2.Value        22
pfset GeomInput.g3.Value        23
pfset GeomInput.g4.Value        24
pfset GeomInput.g5.Value        25
pfset GeomInput.g6.Value        26
pfset GeomInput.g7.Value        27
pfset GeomInput.g8.Value        28
pfset GeomInput.g9.Value		29
pfset GeomInput.g10.Value		30
pfset GeomInput.g11.Value		31
pfset GeomInput.g12.Value		32
pfset GeomInput.g13.Value		33
pfset GeomInput.g14.Value		34
pfset GeomInput.g15.Value		35

#-----------------------------------------------------------------------------
# Perm -- Ksat is what is specified here -- m/h units
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                  "domain s3 s6 s9 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15"

pfset Geom.domain.Perm.Type             Constant
pfset Geom.domain.Perm.Value            0.00069988

pfset Geom.s3.Perm.Type                 Constant
pfset Geom.s3.Perm.Value                0.02

pfset Geom.s6.Perm.Type                 Constant
pfset Geom.s6.Perm.Value                0.015

pfset Geom.s9.Perm.Type                 Constant
pfset Geom.s9.Perm.Value                0.01

pfset Geom.g1.Perm.Type                 Constant
pfset Geom.g1.Perm.Value                0.00122772

pfset Geom.g2.Perm.Type                 Constant
pfset Geom.g2.Perm.Value                0.00736632

pfset Geom.g3.Perm.Type                 Constant
pfset Geom.g3.Perm.Value                0.00245544

pfset Geom.g4.Perm.Type                 Constant
pfset Geom.g4.Perm.Value                0.0030693

pfset Geom.g5.Perm.Type                 Constant
pfset Geom.g5.Perm.Value                0.0061386

pfset Geom.g6.Perm.Type                 Constant
pfset Geom.g6.Perm.Value                0.00491088

pfset Geom.g7.Perm.Type                 Constant
pfset Geom.g7.Perm.Value                0.00552474

pfset Geom.g8.Perm.Type                 Constant
pfset Geom.g8.Perm.Value                0.0122772

pfset Geom.g9.Perm.Type                 Constant
pfset Geom.g9.Perm.Value                0.00368316

pfset Geom.g10.Perm.Type                Constant
pfset Geom.g10.Perm.Value               0.0122772

pfset Geom.g11.Perm.Type                Constant
pfset Geom.g11.Perm.Value               0.0122772

pfset Geom.g12.Perm.Type                Constant
pfset Geom.g12.Perm.Value               0.0122772

pfset Geom.g13.Perm.Type                Constant
pfset Geom.g13.Perm.Value               0.0122772

pfset Geom.g14.Perm.Type                Constant
pfset Geom.g14.Perm.Value               0.0122772

pfset Geom.g15.Perm.Type                Constant
pfset Geom.g15.Perm.Value               0.0122772

#-----------------------------------------------------------------------------
# Permeability Tensors
#-----------------------------------------------------------------------------

pfset Perm.TensorType			TensorByGeom
pfset Geom.Perm.TensorByGeom.Names	"domain"
pfset Geom.domain.Perm.TensorValX  	1.0d0
pfset Geom.domain.Perm.TensorValY  	1.0d0
pfset Geom.domain.Perm.TensorValZ  	1.0d0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type			Constant
pfset SpecificStorage.GeomNames			"domain"
pfset Geom.domain.SpecificStorage.Value 	1.0e-5


#-----------------------------------------------------------------------------
# Phases  -- Setting to 1.0 allows for the calculation of K instead of k
#-----------------------------------------------------------------------------

pfset Phase.Names					"water"
pfset Phase.water.Density.Type	    Constant
pfset Phase.water.Density.Value	    1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names		""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------
pfset Geom.Retardation.GeomNames	""

#-----------------------------------------------------------------------------
# Gravity  -- Setting to 1.0 allows for the calculation of K instead of k
#-----------------------------------------------------------------------------
pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info -- units set by permeability -- m/h
#-----------------------------------------------------------------------------

# Base unit of time for all time values entered. All time should be expressed as multiples of this value.  Perm = m/hr, base unit = hr
pfset TimingInfo.BaseUnit        1.0
pfset TimingInfo.StartCount      $istep 
pfset TimingInfo.StartTime       $istep
pfset TimingInfo.StopTime        8759

# Interval in which output files will be written
pfset TimingInfo.DumpInterval    1.0

#Setup time step type and value
pfset TimeStep.Type              Constant
pfset TimeStep.Value             1.0

#-----------------------------------------------------------------------------
# Porosity -- Assign for soil units, geologic units will receive domain porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames		"domain s3 s6 s9 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15"
#pfset Geom.Porosity.GeomNames		"domain"

pfset Geom.domain.Porosity.Type		Constant
pfset Geom.domain.Porosity.Value	0.2

pfset Geom.s3.Porosity.Type			Constant
pfset Geom.s3.Porosity.Value		0.387

pfset Geom.s6.Porosity.Type			Constant
pfset Geom.s6.Porosity.Value		0.399

pfset Geom.s9.Porosity.Type			Constant
pfset Geom.s9.Porosity.Value		0.442

pfset Geom.g1.Porosity.Type         Constant
pfset Geom.g1.Porosity.Value        0.2

pfset Geom.g2.Porosity.Type         Constant
pfset Geom.g2.Porosity.Value        0.3

pfset Geom.g3.Porosity.Type         Constant
pfset Geom.g3.Porosity.Value        0.15

pfset Geom.g4.Porosity.Type         Constant
pfset Geom.g4.Porosity.Value        0.3

pfset Geom.g5.Porosity.Type         Constant
pfset Geom.g5.Porosity.Value        0.35

pfset Geom.g6.Porosity.Type         Constant
pfset Geom.g6.Porosity.Value        0.3

pfset Geom.g7.Porosity.Type         Constant
pfset Geom.g7.Porosity.Value        0.3

pfset Geom.g8.Porosity.Type         Constant
pfset Geom.g8.Porosity.Value        0.4

pfset Geom.g9.Porosity.Type         Constant
pfset Geom.g9.Porosity.Value        0.2

pfset Geom.g10.Porosity.Type        Constant
pfset Geom.g10.Porosity.Value       0.35

pfset Geom.g11.Porosity.Type        Constant
pfset Geom.g11.Porosity.Value       0.35

pfset Geom.g12.Porosity.Type        Constant
pfset Geom.g12.Porosity.Value       0.35

pfset Geom.g13.Porosity.Type        Constant
pfset Geom.g13.Porosity.Value       0.35

pfset Geom.g14.Porosity.Type        Constant
pfset Geom.g14.Porosity.Value       0.35

pfset Geom.g15.Porosity.Type        Constant
pfset Geom.g15.Porosity.Value       0.35

#-----------------------------------------------------------------------------
# Domain -- Report/Upload the domain problem that is defined above
#-----------------------------------------------------------------------------
pfset Domain.GeomName			"domain"


#-----------------------------------------------------------------------------
# Mobility -- Mobility between phases -- Only one phase in this problem
#-----------------------------------------------------------------------------

pfset Phase.water.Mobility.Type		Constant
pfset Phase.water.Mobility.Value	1.0

#-----------------------------------------------------------------------------
# Wells
#----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------

pfset Cycle.Names 				constant
pfset Cycle.constant.Names			"alltime"
pfset Cycle.constant.alltime.Length    		  1
pfset Cycle.constant.Repeat           		  -1
 
#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------

pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"  
pfset Patch.z-lower.BCPressure.alltime.Value	       0.0 

pfset Patch.x-upper.BCPressure.Type		      FluxConst
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

pfset Patch.z-upper.BCPressure.Type		      OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle		      "constant"  
pfset Patch.z-upper.BCPressure.alltime.Value		0.0

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type 			"PFBFile"
pfset TopoSlopesX.GeomNames 	"domain"
pfset TopoSlopesX.FileName		$slopex

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type			"PFBFile"
pfset TopoSlopesY.GeomNames		"domain"
pfset TopoSlopesY.FileName  	$slopey

#---------------------------------------------------------
# Mannings Coefficient
#---------------------------------------------------------
pfset Mannings.Type			"PFBFile"
pfset Mannings.GeomNames	"domain"
pfset Mannings.FileName		$mannings

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type            VanGenuchten
pfset Phase.RelPerm.GeomNames       "domain s3 s6 s9 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15"

pfset Geom.domain.RelPerm.Alpha     3.5
pfset Geom.domain.RelPerm.N         2.0

pfset Geom.s3.RelPerm.Alpha			2.692
pfset Geom.s3.RelPerm.N				2.445
	
pfset Geom.s6.RelPerm.Alpha			1.122
pfset Geom.s6.RelPerm.N				2.479
	
pfset Geom.s9.RelPerm.Alpha			1.585
pfset Geom.s9.RelPerm.N				2.413

pfset Geom.g1.RelPerm.Alpha         1.585
pfset Geom.g1.RelPerm.N             2.413

pfset Geom.g2.RelPerm.Alpha         1.585
pfset Geom.g2.RelPerm.N             2.413

pfset Geom.g3.RelPerm.Alpha         1.585
pfset Geom.g3.RelPerm.N             2.413

pfset Geom.g4.RelPerm.Alpha         1.585
pfset Geom.g4.RelPerm.N             2.413

pfset Geom.g5.RelPerm.Alpha         1.585
pfset Geom.g5.RelPerm.N             2.413

pfset Geom.g6.RelPerm.Alpha         1.585
pfset Geom.g6.RelPerm.N             2.413

pfset Geom.g7.RelPerm.Alpha         1.585
pfset Geom.g7.RelPerm.N             2.413

pfset Geom.g8.RelPerm.Alpha         1.585
pfset Geom.g8.RelPerm.N             2.413

pfset Geom.g9.RelPerm.Alpha         1.585
pfset Geom.g9.RelPerm.N             2.413

pfset Geom.g10.RelPerm.Alpha        1.585
pfset Geom.g10.RelPerm.N            2.413

pfset Geom.g11.RelPerm.Alpha        1.585
pfset Geom.g11.RelPerm.N            2.413

pfset Geom.g12.RelPerm.Alpha        1.585
pfset Geom.g12.RelPerm.N            2.413

pfset Geom.g13.RelPerm.Alpha        1.585
pfset Geom.g13.RelPerm.N            2.413

pfset Geom.g14.RelPerm.Alpha        1.585
pfset Geom.g14.RelPerm.N            2.413

pfset Geom.g15.RelPerm.Alpha        1.585
pfset Geom.g15.RelPerm.N			2.413

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type               VanGenuchten
pfset Phase.Saturation.GeomNames          "domain s3 s6 s9 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15"

pfset Geom.domain.Saturation.Alpha		3.5
pfset Geom.domain.Saturation.N 			2.0
pfset Geom.domain.Saturation.SRes		0.2
pfset Geom.domain.Saturation.SSat		1.0

pfset Geom.s3.Saturation.Alpha       	2.692
pfset Geom.s3.Saturation.N				2.445
pfset Geom.s3.Saturation.SRes			0.1
pfset Geom.s3.Saturation.SSat			1.0

pfset Geom.s6.Saturation.Alpha			1.122
pfset Geom.s6.Saturation.N				2.479
pfset Geom.s6.Saturation.SRes			0.15
pfset Geom.s6.Saturation.SSat			1.0

pfset Geom.s9.Saturation.Alpha			1.585
pfset Geom.s9.Saturation.N				2.413
pfset Geom.s9.Saturation.SRes			0.18
pfset Geom.s9.Saturation.SSat			1.0

pfset Geom.g1.Saturation.Alpha          1.585
pfset Geom.g1.Saturation.N              2.413
pfset Geom.g1.Saturation.SRes           0.18
pfset Geom.g1.Saturation.SSat           1.0

pfset Geom.g2.Saturation.Alpha          1.585
pfset Geom.g2.Saturation.N              2.413
pfset Geom.g2.Saturation.SRes           0.18
pfset Geom.g2.Saturation.SSat           1.0

pfset Geom.g3.Saturation.Alpha          1.585
pfset Geom.g3.Saturation.N              2.413
pfset Geom.g3.Saturation.SRes           0.18
pfset Geom.g3.Saturation.SSat           1.0

pfset Geom.g4.Saturation.Alpha          1.585
pfset Geom.g4.Saturation.N              2.413
pfset Geom.g4.Saturation.SRes           0.18
pfset Geom.g4.Saturation.SSat           1.0

pfset Geom.g5.Saturation.Alpha          1.585
pfset Geom.g5.Saturation.N              2.413
pfset Geom.g5.Saturation.SRes           0.18
pfset Geom.g5.Saturation.SSat           1.0

pfset Geom.g6.Saturation.Alpha          1.585
pfset Geom.g6.Saturation.N              2.413
pfset Geom.g6.Saturation.SRes           0.18
pfset Geom.g6.Saturation.SSat           1.0

pfset Geom.g7.Saturation.Alpha          1.585
pfset Geom.g7.Saturation.N              2.413
pfset Geom.g7.Saturation.SRes           0.18
pfset Geom.g7.Saturation.SSat           1.0

pfset Geom.g8.Saturation.Alpha          1.585
pfset Geom.g8.Saturation.N              2.413
pfset Geom.g8.Saturation.SRes           0.18
pfset Geom.g8.Saturation.SSat           1.0

pfset Geom.g9.Saturation.Alpha          1.585
pfset Geom.g9.Saturation.N              2.413
pfset Geom.g9.Saturation.SRes           0.18
pfset Geom.g9.Saturation.SSat           1.0

pfset Geom.g10.Saturation.Alpha         1.585
pfset Geom.g10.Saturation.N             2.413
pfset Geom.g10.Saturation.SRes          0.18
pfset Geom.g10.Saturation.SSat          1.0

pfset Geom.g11.Saturation.Alpha         1.585
pfset Geom.g11.Saturation.N             2.413
pfset Geom.g11.Saturation.SRes          0.18
pfset Geom.g11.Saturation.SSat          1.0

pfset Geom.g12.Saturation.Alpha         1.585
pfset Geom.g12.Saturation.N             2.413
pfset Geom.g12.Saturation.SRes          0.18
pfset Geom.g12.Saturation.SSat          1.0

pfset Geom.g13.Saturation.Alpha         1.585
pfset Geom.g13.Saturation.N             2.413
pfset Geom.g13.Saturation.SRes          0.18
pfset Geom.g13.Saturation.SSat          1.0

pfset Geom.g14.Saturation.Alpha         1.585
pfset Geom.g14.Saturation.N             2.413
pfset Geom.g14.Saturation.SRes          0.18
pfset Geom.g14.Saturation.SSat          1.0

pfset Geom.g15.Saturation.Alpha         1.585
pfset Geom.g15.Saturation.N             2.413
pfset Geom.g15.Saturation.SRes          0.18
pfset Geom.g15.Saturation.SSat			1.0

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type				Constant
pfset PhaseSources.water.GeomNames			domain
pfset PhaseSources.water.Geom.domain.Value	0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
pfset KnownSolution		NoKnownSolution

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------

pfset Solver										Richards
pfset Solver.MaxIter								1000000000
			
pfset Solver.TerrainFollowingGrid					True
			
pfset Solver.Nonlinear.MaxIter						1500
pfset Solver.Nonlinear.ResidualTol					0.00001
pfset Solver.Nonlinear.EtaValue						0.001
			
pfset Solver.PrintSubsurf							False
pfset Solver.Drop									1E-15
pfset Solver.AbsTol									1E-9
pfset Solver.MaxConvergenceFailures					7
			
pfset Solver.Nonlinear.UseJacobian					True 
pfset Solver.Nonlinear.StepTol						1e-19
pfset Solver.Nonlinear.AbsTol						1e-9
			
pfset Solver.Linear.MaxIter							15000
pfset Solver.Nonlinear.Globalization				LineSearch
pfset Solver.Linear.KrylovDimension					150
pfset Solver.Linear.MaxRestarts						6
			
pfset Solver.Linear.Preconditioner					PFMG
pfset Solver.Linear.Preconditioner.PCMatrixType		FullJacobian

#-----------------------------------------------------------------------------
#CLM: Output setup
#-----------------------------------------------------------------------------
  pfset Solver.LSM                          CLM
  pfset Solver.CLM.CLMFileDir               "clm_output/"
  pfset Solver.CLM.Print1dOut               False
  pfset Solver.BinaryOutDir                 False
  pfset Solver.CLM.CLMDumpInterval          1
  pfset Solver.CLM.EvapBeta                 Linear
  pfset Solver.CLM.VegWaterStress           Saturation
  pfset Solver.CLM.ResSat                   0.2
  pfset Solver.CLM.WiltingPoint             0.2
  pfset Solver.CLM.FieldCapacity            1.00

#------Met forcing and timestep setup
  pfset Solver.CLM.MetForcing				2D
  #pfset Solver.CLM.MetFileName              "ER"
  #pfset Solver.CLM.MetFilePath              ../Forcing_WY06_LF

  pfset Solver.CLM.MetFileName              "NLDAS"
  pfset Solver.CLM.MetFilePath              ../forcing/wy2015/2.5/

#  pfset Solver.CLM.MetFileNT                1
  pfset Solver.CLM.IstepStart               $clmstep
  pfset Solver.CLM.ReuseCount               1
  #pfset Solver.CLM.MetForcing              1D
  #pfset Solver.CLM.MetFileName             "1D_forcing_castnet.txt"
  #pfset Solver.CLM.MetFilePath             ./

#------Irrigation setup
  pfset Solver.CLM.IrrigationType			none
  pfset Solver.CLM.ForceVegetation			False
 
#------Writing output:
  pfset Solver.WriteSiloEvapTrans			False
  pfset Solver.WriteSiloOverlandBCFlux		False
  pfset Solver.PrintSubsurfData				True
  pfset Solver.PrintPressure				True
  pfset Solver.PrintSaturation				True
  pfset Solver.WriteCLMBinary				False
  pfset Solver.WriteSiloCLM					False
  pfset Solver.PrintCLM						True
  pfset Solver.CLM.RootZoneNZ				4
  pfset Solver.CLM.SoiLayer					4
  pfset Solver.CLM.WriteLogs				False
  pfset Solver.CLM.WriteLastRST				False
  pfset Solver.CLM.DailyRST				True
  pfset Solver.CLM.SingleFile				True

#------Write velocities and EvapTrans for use with SLIM
  pfset Solver.PrintEvapTrans				True
  pfset Solver.PrintVelocities              True
#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

# Restart from previous pressure file
pfset ICPressure.Type					PFBFile
pfset ICPressure.GeomNames				domain
pfset Geom.domain.ICPressure.FileName	$ip
pfdist 									$ip
#--------------------------------------------------------
# Distribute Inputs
#--------------------------------------------------------

pfset ComputationalGrid.NZ		1	 
pfdist 							$slopex
pfdist 							$slopey
pfdist 							$mannings
# Be sure to reset computational grid to actual dz values (already done below for use of indicator file)	
pfset ComputationalGrid.NZ      5
pfdist 							$sub

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
#pfrun 		$runname
pfwritedb	$runname

