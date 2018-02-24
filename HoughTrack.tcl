#!/bin/sh
# \
exec ./intiwish "$0" "$@"
###############################################

#Update TCL librairy:
auto_mkindex ./LIBTCL *.tcl

lappend auto_path ./LIBTCL/

# Load package tkimg, if present
catch {package require Img}

if {$argc != 1} {
 puts "Usage: HoughTrack.tcl <Test_directory>"
 exit
}

set arg(rep_test) [lindex $argv 0]
set sequence [lsort [glob $arg(rep_test)/*]]
set first_image [lindex $sequence 0]

# Two image types: prototype and test

image create photo testimage -file $first_image -palette 256/256/256
# 5 images in the interface, 3 with type "proto", 2 with type "test"
Intinew CMD 5

image create photo testsaliency -width [image width testimage] \
    -height [image height testimage] -palette 256/256/256
 
CMD change_bord 1

label .prototestparam 

label .protoparam

label .prototype -relief ridge -bg red
label .proto0
canvas .protocnv0 -width 50 -height 50 -bg black 
label .proto0txt -text "Prototype"
pack .protocnv0 .proto0txt -side top -in .proto0 
label .proto1
canvas .protocnv1 -width 50 -height 50 -bg black 
label .proto1txt -text "Index"
pack .protocnv1 .proto1txt -side top -in .proto1
label .proto2
canvas .protocnv2 -width 50 -height 50 -bg black 
label .proto2txt -text "Weight"
pack .protocnv2 .proto2txt -side top -in .proto2
pack .proto0 .proto1 .proto2 -side left -in .prototype

label .param1
##################################################################
label .differentielles -relief ridge -borderwidth 3 -bg blue
scale .sigma -from 0.5 -to 10 -resolution 0.5 -label "Sigma Differentials"\
                             -variable sigma -orient horizontal
scale .nbetiq -from 10 -to 100 -resolution 1 -label "Number of labels"\
                             -variable nbetiq -orient horizontal
pack .sigma .nbetiq -in .differentielles -side top -expand true -fill both
label .choix_ordre -relief ridge -borderwidth 3 -bg blue
radiobutton .ord0 -variable ord_diff -value 0 -text "Order 0"
radiobutton .ord1 -variable ord_diff -value 1 -text "Order 1"
radiobutton .ord2 -variable ord_diff -value 2 -text "Order 2"
pack .ord0 .ord1 .ord2 -side left -in .choix_ordre
##################################################################
button .initRtable -text "Init R-Table" -bg lightgreen -command "Calcul_RTable"
##################################################################
label .gamma -relief ridge -borderwidth 5 -bg green
label .gammalab -text "Gamma correction"
scale .gammascale -from 0.1 -to 5 -variable gamma \
   -orient horizontal -command "change_gamma_map" \
    -resolution 0.1
pack .gammalab .gammascale -in .gamma -side top -expand true -fill x
##################################################################
label .nbest -relief ridge -borderwidth 5 -bg blue
label .nbestlab -text "Number of objects"
scale .nbestscale -from 0 -to 8 -variable nbest \
    -orient horizontal -resolution 1
label .maj -relief ridge -borderwidth 3 -bg blue
label .maj_lab -text "Update Prototype"
label .maj_choix
radiobutton .maj_non -variable maj_proto -value 0 -text "No"
radiobutton .maj_oui -variable maj_proto -value 1 -text "Yes"
pack .maj_non .maj_oui -side left -in .maj_choix -expand true -fill x
pack .maj_lab .maj_choix -in .maj -side top -expand true -fill x
pack .nbestlab .nbestscale .maj -in .nbest -side top -expand true -fill x
##################################################################
button .tracebest -text "TRACKING" -bg lightblue -command "Go_Track"

pack .differentielles .choix_ordre .initRtable .gamma .nbest .tracebest -in .param1 -side top -expand true -fill both

pack .prototype .param1 -side top -in .protoparam -expand true -fill both

label .testparam
label .tests 

Imgnew testimage 3
Imgnew testsaliency 4

canvas .testimgcanvas -width [image width testimage] -height [image height testimage] -bg black
.testimgcanvas create image [expr [image width testimage]/2] [expr [image height testimage]/2] -image testimage
bind .testimgcanvas <Button-1> "debutproto %x %y"
bind .testimgcanvas <B1-Motion> "dessineproto %x %y"
bind .testimgcanvas <ButtonRelease-1> "finproto %x %y"
canvas .testsaliencycanvas -width [image width testsaliency] -height [image height testsaliency] -bg black
.testsaliencycanvas create image [expr [image width testsaliency]/2] [expr [image height testsaliency]/2] -image testsaliency
pack .testimgcanvas .testsaliencycanvas -side top -in .tests

label .param2

pack .tests .param2 -side top -in .testparam -expand true -fill both

set sigma 2.0
set nbetiq 10
set gamma 1
set nbest 0
set ord_diff 1
set maj_proto 0

set nbetiq_courant $nbetiq 
set sigma_courant $sigma

set init_proto 0

pack .protoparam .testparam -in .prototestparam -side left -expand true -fill both

update

label .menubar

menubutton .gene -text "File" -menu .gene.menu
set m [menu .gene.menu]
$m add command -label "Quit" -command "exit"

menubutton .opt -text "Options" -menu .opt.menu
set m [menu .opt.menu]
$m add cascade -label "Line Colour" -menu $m.sub1
set m1 [menu $m.sub1 -tearoff 0]
$m1 add radio -label "black" -variable line_color -value 0 -foreground black
$m1 add radio -label "red" -variable line_color -value 1 -foreground red
$m1 add radio -label "green" -variable line_color -value 2 -foreground green
$m1 add radio -label "yellow" -variable line_color -value 3 -foreground yellow
$m1 add radio -label "blue" -variable line_color -value 4 -foreground blue
$m1 add radio -label "magenta" -variable line_color -value 5 -foreground magenta
$m1 add radio -label "cyan" -variable line_color -value 6 -foreground cyan
$m1 add radio -label "white" -variable line_color -value 7 -foreground white
$m1 invoke 5
$m add command -label "About..." -command "pub"


pack .gene .opt -in .menubar -side left

pack .menubar .prototestparam -side top
 
