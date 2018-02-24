#!/bin/sh
# \
exec tclsh "$0" "$@"
###############################################

if {$argc != 1} {
 puts "Usage : ComposeImg.tcl <Directory>"
 exit
}

set directory [lindex $argv 0]
set i 0

set list_proto [lsort [glob $directory/PROTO_*.jpg]]
set list_result [lsort [glob $directory/OUT_*.jpg]]
set list_thg [lsort [glob $directory/TH_*.jpg]]
foreach fich_result $list_result {
    set num [format "%03d" $i]
    set fich_proto [lindex $list_proto $i]
    set fich_thg [lindex $list_thg $i]
    exec composite $fich_proto $fich_thg TMP.jpg 
    exec convert $fich_result TMP.jpg -append Trame_$num.jpg
    puts $num
    incr i
}
# Création vidéo : à adapter
#exec ffmpeg -framerate 25 -i Trame_%03d.jpg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
