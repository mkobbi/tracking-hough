proc debutproto {x y} {
global x1_proto 
global y1_proto
global box_anchor
global box_last

set x1_proto $x 
set y1_proto $y
catch {.testimgcanvas delete $box_last}
set box_anchor [list $x $y]
catch {unset box_last}
}

proc dessineproto {x y} {
global box_anchor
global box_last
global box_tag

catch {.testimgcanvas delete $box_last}
set box_last [eval {.testimgcanvas create rect} $box_anchor \
			{$x $y -tag box_tag -outline red}]
}

proc finproto {x y} {
global x1_proto 
global y1_proto
global x2_proto 
global y2_proto
global box_last
global position_proto_x
global position_proto_y

set x2_proto $x 
set y2_proto $y

if {[expr ($x2_proto - $x1_proto)%2] == 0} {incr x2_proto -1}
if {[expr ($y2_proto - $y1_proto)%2] == 0} {incr y2_proto -1}
# Les dimensions du proto sont toujours impaires
image create photo proto -palette 256/256/256
proto copy testimage -from $x1_proto $y1_proto $x2_proto $y2_proto
image create photo protoreco -width [image width proto] \
    -height [image height proto] -palette 256/256/256
image create photo protopoids -width [image width proto] \
    -height [image height proto] -palette 256/256/256
for {set i 0} {$i < 3} {incr i} {
    .protocnv$i configure -width [image width proto] \
	-height [image height proto] -bg black }
Imgnew proto 0
Imgnew protoreco 1
Imgnew protopoids 2
.protocnv0 create image [expr [image width proto]/2] \
                        [expr [image height proto]/2] -image proto
.protocnv1 create image [expr [image width protoreco]/2] \
                        [expr [image height protoreco]/2] -image protoreco
.protocnv2 create image [expr [image width protoreco]/2] \
                        [expr [image height protoreco]/2] -image protopoids
set position_proto_x [expr ($x1_proto + $x2_proto)/2]
set position_proto_y [expr ($y1_proto + $y2_proto)/2]
puts "Initial Prototype Position: x = $position_proto_x ; y = $position_proto_y"
}

proc change_gamma_map {value} {

testsaliency configure -gamma $value
}

proc Calcul_RTable {} {
global sigma
global nbetiq
global sigma_courant
global nbetiq_courant
global ord_diff
global init_proto

CMD Allocate_RTable $nbetiq
CMD Calcul_RTable_Proto_ord$ord_diff 0 1 2 $sigma $nbetiq

set nbetiq_courant $nbetiq
set sigma_courant $sigma
set init_proto 1
}

proc Go_Track {} {
   global ord_diff
   global arg
   global nbetiq
   global sigma
   global nbest
   global maj_proto
   global rectangle
   global sequence
   global box_last
   global init_proto
   global position_proto_x
   global position_proto_y

    #Nettoyage première image
    catch {.testimgcanvas delete $box_last}
    #Message d'erreur si proto non initialisé
    if {$init_proto == 0} {
        set error_msg "******Error: Prototype does not exist******"
   	Put_message $error_msg
        return
	} 
    CMD Allocate_Hough_2d [image width testimage] [image height testimage]
    Imgnew testsaliency 4

    CMD Allocate_BestPos $nbest
set i 0
foreach fich $sequence {
   set nom [file tail $fich]
    #Creation test image and initialise GHT
    image create photo testimage -file $fich -palette 256/256/256
    Imgnew testimage 3
    CMD Reset_Hough2d [image width testimage] [image height testimage]
    CMD Hough_General_ord$ord_diff 3 $sigma $nbetiq
    CMD MaJ_Hough2d 4
    update
	# Uncomment to save the transform sequence
         #set outimg [format TH_%03d.jpg $i]
         #testsaliency write $outimg -format {JPEG -quality 85}
      set pos [CMD Find_Best_THG [image width testimage] [image height testimage] \
                                 $nbest [image width proto] [image height proto]]
      nettoie_rectangles
      for {set k 0} {$k < $nbest} {incr k} {
        set x [lindex $pos [expr 2*$k]]
        set y [lindex $pos [expr 2*$k +1]]
        puts "Detection n°$k: x = $x ; y = $y"
        #CMD Dessine_rectangle 3 $x $y [image width proto] [image height proto] $k
        trace_rectangle $k $x $y [image width proto] [image height proto]
         }
    #Uncommenter to check the effect of non maxima suppression
    #CMD MaJ_Hough2d 4
    update
    # Uncomment to save the result sequence
      #set outimg [format OUT_%03d.jpg $i]
      #image create photo mon_image -format window -data .testimgcanvas
      #mon_image write $outimg -format {JPEG -quality 100}
    if {($nbest > 0)&&($maj_proto==1)} {
      set new_pos [CMD Update_Proto 0 1 2 3 $position_proto_x $position_proto_y $nbetiq $sigma] 
      set position_proto_x [lindex $new_pos 0]
      set position_proto_y [lindex $new_pos 1]
      puts "Nouvelle Position Prototype : x = $position_proto_x ; y = $position_proto_y"
       }
    # Uncomment to save the proto sequence
     #set outimg [format PROTO_%03d.jpg $i]
     #proto write $outimg -format {JPEG -quality 85}
    incr i
    }
}

proc trace_rectangle {num_colour x y width_proto height_proto} {
global rectangle 

  set colour_index {"red" "green" "yellow" "blue" "magenta" "cyan" "white" "black"}
  set colour [lindex $colour_index $num_colour]
  set x1 [expr $x - $width_proto/2]
  set y1 [expr $y - $height_proto/2]
  set x2 [expr $x + $width_proto/2]
  set y2 [expr $y + $height_proto/2]
  set rect_debut [list $x1 $y1]
  set rect_fin [list $x2 $y2]
  set rectangle($num_colour) [eval {.testimgcanvas create rect} $rect_debut \
	       {$x2 $y2 -tag rect_tag -outline $colour}]
} 

proc nettoie_rectangles {} {
global rectangle
global nbest

  for {set k 0} {$k < $nbest} {incr k} {
    catch {.testimgcanvas delete $rectangle($k)}
  }
}


