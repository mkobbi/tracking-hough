#!/bin/sh
# \
exec ./intiwish "$0" "$@"
###############################################

#Update TCL librairy:
auto_mkindex ./LIBTCL *.tcl

lappend auto_path ./LIBTCL/

# Load package tkimg, if present
catch {package require Img}

if {$argc != 2} {
    puts "Usage: HoughTrack_noGUI.tcl <Img_Proto> <Test_directory>"
    exit
}

set arg(image_proto) [lindex $argv 0]
set arg(rep_test) [lindex $argv 1]

# 5 images, 3 with type "proto", 2 with type "test"
Intinew CMD 5 

# Create the 3 images with type "proto"
image create photo proto -file $arg(image_proto) -palette 256/256/256
image create photo etiquette -width [image width proto] \
    -height [image height proto] -palette 256/256/256
image create photo poids -width [image width proto] \
    -height [image height proto] -palette 256/256/256
# Attribute their numbers
Imgnew proto 0
Imgnew etiquette 1
Imgnew poids 2

# Global parameters
set sigma 2.0
set nbetiq 10
set nbest 5

CMD change_bord 1
# Allocate and create the R-Table
CMD Allocate_RTable $nbetiq
CMD Calcul_RTable_Proto_ord1 0 1 2 $sigma $nbetiq
# Allocate the "best position" array
CMD Allocate_BestPos $nbest

set num_image 0
# Loop on test directory images (alphanumeric order)
foreach fich [lsort [glob $arg(rep_test)/*]] {
  set nom [file tail $fich]
  puts $nom
  set nom2 [file rootname $nom]
  # Create the 2 images with type "test" and attribute their numbers
  image create photo testimage -file $fich -palette 256/256/256
  Imgnew testimage 3
  if {$num_image == 0} {
    image create photo testsaliency -width [image width testimage] \
      -height [image height testimage] -palette 256/256/256
    Imgnew testsaliency 4
    # Allocate the GHT 
    CMD Allocate_Hough_2d [image width testimage] [image height testimage]
    # Default position of the proto: image centre 
    set position_proto_x [expr [image width testimage]/2]
    set position_proto_y [expr [image height testimage]/2]
    }
  CMD Reset_Hough2d [image width testimage] [image height testimage]
  CMD Hough_General_ord1 3 $sigma $nbetiq
  CMD MaJ_Hough2d 4
  # Uncomment to save the transform sequence
    #set outimg TH_$nom2.jpg
    #testsaliency write $outimg -format {JPEG -quality 85}
  set pos [CMD Find_Best_THG [image width testimage] [image height testimage] \
                                 $nbest [image width proto] [image height proto]]
  set new_pos [CMD Update_Proto 0 1 2 3 $position_proto_x $position_proto_y $nbetiq $sigma] 
  set position_proto_x [lindex $new_pos 0]
  set position_proto_y [lindex $new_pos 1]
  puts "New Prototype Position: x = $position_proto_x ; y = $position_proto_y"
  # Uncomment to save the proto sequence
     #set outimg PROTO_$nom2.jpg
     #proto write $outimg -format {JPEG -quality 85}
  # Uncomment to save the result sequence
    for {set i 0} {$i < $nbest} {incr i} {
      set x [lindex $pos [expr 2*$i]]
      set y [lindex $pos [expr 2*$i +1]]
      puts "Detection n°$i : x = $x ; y = $y"
      CMD Dessine_rectangle 3 $x $y \
      [image width proto] [image height proto] [expr $i + 1]	
      }
    #set outimg OUT_$nom2.jpg
    #testimage write $outimg -format {JPEG -quality 85}
}

