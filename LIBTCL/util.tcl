proc change_fmt {} {
global entree 
global f

regsub gif|pgm|ppm|tif|jpg|bmp|png $entree(result) $entree(format) entree(result)
}

proc Select_fichier {action} {
global entree
global f

set w .filebox
toplevel $w
wm title $w "Selection de fichiers"
wm iconname $w "filebox"

label $w.msg -wraplength 4i -justify left -text "Nom du fichier ?"
pack $w.msg -side top
 
frame $w.buttons
pack $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.cancel -text "Annuler" -command {set entree(ok) 0}
button $w.buttons.ok -text "OK" -command {set entree(ok) 1}
pack $w.buttons.cancel $w.buttons.ok -side left -expand 1

   set f [frame $w.$action]
if {$action == "sauvegarder"} {
    set entree(result) "sans_titre.png"
} else { 
    set entree(result) ""
}
    label $f.lab -text "$action : " -anchor e
    entry $f.ent -width 20 -textvar entree(result)
    bind $f.ent <Return> {set entree(ok) 1}
    button $f.but -text "Parcourir ..." -command "fileDialog $w $f.ent $action"
    pack $f.lab -side left
    pack $f.ent -side left -expand yes -fill x
    pack $f.but -side left
    pack $f -fill x -padx 1c -pady 3

    checkbutton $w.strict -text "Look Motif" \
	-variable tk_strictMotif -onvalue 1 -offvalue 0
    pack $w.strict -anchor c

if {$action == "sauvegarder"} {
    set g [label $w.form_choice]
   label $g.texte -text "Format de sauvegarde"
   radiobutton $g.gif -variable entree(format) -value "gif" -text "gif" -command "change_fmt"
   radiobutton $g.pgm -variable entree(format) -value "pgm" -text "pgm" -command "change_fmt"
   radiobutton $g.ppm -variable entree(format) -value "ppm" -text "ppm" -command "change_fmt"
   radiobutton $g.png -variable entree(format) -value "png" -text "png" -command "change_fmt"
   radiobutton $g.tiff -variable entree(format) -value "tif" -text "tiff" -command "change_fmt"
   radiobutton $g.jpg -variable entree(format) -value "jpg" -text "jpg" -command "change_fmt"
   radiobutton $g.bmp -variable entree(format) -value "bmp" -text "bmp" -command "change_fmt"
   set entree(format) "png"
   pack $g.texte -in $w.form_choice -side top
   pack $g.gif $g.pgm $g.ppm $g.png $g.tiff $g.jpg $g.bmp -in $w.form_choice -side left
   pack $g -side top
}
    
focus $f.ent
grab $w
tkwait variable entree(ok)
grab release $w
destroy $w

}

proc fileDialog {w ent operation} {
global entree

    switch $operation {
"sauvegarder" {
     set types {
	 {"Image"   {.gif .pgm .ppm .png .tif .jpg .bmp}     }
	 {"Autres fichiers"		*}
     }
 set entree(result) [tk_getSaveFile -filetypes $types -parent $w \
 -initialfile "sans_titre" -defaultextension .$entree(format) -title "sauver"]
 }
 "ouvrir" {
     set types {
	 {"Images gif"		{.gif}		}
	 {"Images pgm"		{.pgm}	        }
	 {"Images ppm"		{.ppm}		}
         {"Images png"		{.png}		}
	 {"Images jpeg"		{.jpg .jpeg}	        }
	 {"Images tiff"		{.tif .tiff}		}
	 {"Images bmp"		{.bmp}	        }
	 {"Autres fichiers"		*}
     }
set entree(result) [tk_getOpenFile -filetypes $types -parent $w \
			  -title $operation]
    }
     "enregistrer_macro" {
	 set types {
	     {"Fichiers macro" {.mac}}
	     {"Autre fichiers" *}
	 }
  set entree(result) [tk_getSaveFile -filetypes $types -parent $w \
 -initialfile "sans_titre" -defaultextension .mac -title "Enregistrer macro"]
 }
      "charger_macro" {
       set types {
	     {"Fichiers macro" {.mac}}
	     {"Autre fichiers" *}
	 }
set entree(result) [tk_getOpenFile -filetypes $types -parent $w \
			  -title $operation]
    }
}
}

proc new {} {
    message "Commande indisponible pour l'instant ;\n Veuillez relancer l'interface manuellement"
}
    

proc set_titre {image numero} {
global version

    set Imginfo [list Image [image width toto($numero)] x [image height toto($numero)]]
    set Imgnom [file tail $image]
    wm title .inti($numero) [concat Inti $version : $Imginfo - $Imgnom]
}

proc save_canvas {} {
global entree

Select_fichier sauvegarder
  if {$entree(ok)} {
    image create photo monimage -format window -data .imagentrada
    switch $entree(format) {
      "ppm" {monimage write $entree(result) -format PPM/PGM}
      "gif" {monimage write $entree(result) -format GIF}
      "tif" {monimage write $entree(result) -format TIFF}
      "bmp" {monimage write $entree(result) -format BMP}
      "png" {monimage write $entree(result) -format PNG}
      "jpg" {monimage write $entree(result) -format JPEG}
    }
  }
}


proc save {num} {
global entree

Select_fichier sauvegarder
    if {$entree(ok)} {
	if {$num == 0} {
	    set img entree
	} else {
	    set img sortie
	}
  switch $entree(format) {
      "ppm" {$img write $entree(result) -format PPM/PGM}
      "pgm" {CMD save $num $entree(result)}
      "gif" {$img write $entree(result) -format GIF}
      "tif" {$img write $entree(result) -format TIFF}
      "bmp" {$img write $entree(result) -format BMP}
      "png" {$img write $entree(result) -format PNG}
      "jpg" {$img write $entree(result) -format JPEG}
}
}
}

proc save2 {num} {
global entree
global radius

Select_fichier sauvegarder
    if {$entree(ok)} {
	if {$num == 0} {
	    set img entree
	} else {
	    set img transfo($radius)
	}
  switch $entree(format) {
      "ppm" {$img write $entree(result) -format PPM/PGM}
      "pgm" {$img write $entree(result) -format PPM/PGM}
      "gif" {$img write $entree(result) -format GIF}
      "tif" {$img write $entree(result) -format TIFF}
      "bmp" {$img write $entree(result) -format BMP}
      "png" {$img write $entree(result) -format PNG}
      "jpg" {$img write $entree(result) -format JPEG}
}
}
}

proc savealltsf {} {
global entree
global arg
    
    Select_fichier sauvegarder
    if {$entree(ok)} {
	switch $entree(format) {
	    "ppm" {set format PPM/PGM; set suffix ppm}
	    "pgm" {set format PPM/PGM; set suffix pgm}
	    "gif" {set format GIF; set suffix gif}
	    "tif" {set format TIFF; set suffix tif}
	    "bmp" {set format BMP; set suffix bmp}
            "png" {set format PNG; set suffix png}
	    "jpg" {set format JPEG; set suffix jpg}
	}
	for {set i 1} {$i <= $arg(resol_radius)} {incr i} {
	    set img transfo($i)
	    set prefix [file rootname $entree(result)]
	    set file [format $prefix%02d.$suffix $i]
	    $img write $file -format $format
	}
    }
}

proc savegeneral {num} {
global entree

    Select_fichier sauvegarder
    if {$entree(ok)} {
	switch $num {
	    0 {set img etiquette}
	    1 {set img poids}
	    2 {set img testimage}
	    3 {set img testhough}
	}
	switch $entree(format) {
	    "ppm" {$img write $entree(result) -format PPM/PGM}
	    "pgm" {CMD save $num $entree(result)}
	    "gif" {$img write $entree(result) -format GIF}
	    "tif" {$img write $entree(result) -format TIFF}
	    "bmp" {$img write $entree(result) -format BMP}
            "png" {$img write $entree(result) -format PNG}
	    "jpg" {$img write $entree(result) -format JPEG}
	}
    }
}


proc Clean_Image {} {
global arg

image create photo entree -file $arg(image) -palette 256/256/256 
Imgnew entree 0
update
}

proc Clean_Transfo {} {
global arg

    CMD Reset_Hough $arg(resol_radius) [image width entree] [image height entree]
    for {set i 1} {$i < $arg(resol_radius)} {incr i} {
	CMD MaJ_Hough $i
    }
update
}

proc Clean_Image2 {} {
global arg

image create photo testimage -file $arg(image_test) -palette 256/256/256
Imgnew testimage 3
update
}

proc Clean_Transfo2 {} {

CMD Reset_Hough2d [image width testimage] [image height testimage]
CMD MaJ_Hough2d 4   
update
}

proc Put_message {text} {

toplevel .msg
wm title .msg "Message"

label .msg.corps -wraplength 4i -justify left -text $text
pack .msg.corps -side top
 
button .msg.button -text "OK" -command {set message_lu 1}
pack .msg.button -side left -expand 1
 
grab .msg
tkwait variable message_lu
grab release .msg
destroy .msg

}

proc saisie_variables {num lscale lradio} {
global parametres

set p [toplevel .param($num)]
wm title $p "Saisie de paramètres"
#set haut [expr 100 * ([llength $lscale]+[llength $lradio])]
#wm geometry $p =200x$haut
foreach var $lscale {
	set varname [lindex $var 0]
        set varmin [lindex $var 1]
	set varmax [lindex $var 2]
        if {[llength $var]>3} {
	    set resol [lindex $var 3]
	} else  {
	    set resol 1
	}
	scale $p.sc_$varname -from $varmin -to $varmax\
 	                  -label $varname \
                          -variable parametres($varname) \
	                  -resolution $resol \
                          -orient horizontal
        pack $p.sc_$varname -side top -expand 1
        }
foreach var $lradio {
    set varname [lindex $var 0]
    label $p.lab_$varname -text $varname
    frame $p.fr_$varname
    for {set i 1} {$i < [llength $var]} {incr i} {
	radiobutton $p.but_$varname$i -variable parametres($varname) \
	    -value [expr $i -1] -text [lindex $var $i]
        pack $p.but_$varname$i -in $p.fr_$varname -side left -expand 1
    }
    pack $p.lab_$varname  $p.fr_$varname -side top -expand 1
}
label $p.rampe
button $p.vazy -text "OK" -command {set parametres(fini) glop}
button $p.vazypas -text "Annuler" -command {set parametres(fini) pasglop}
pack $p.vazy $p.vazypas -in $p.rampe -side left -expand 1
pack $p.rampe -side top -expand 1
grab $p
tkwait variable parametres(fini)
grab release $p
destroy $p

}


proc pub {} {

toplevel .pub
wm title .pub "Info Inti2.0"

image create photo logo -file ./.intilogo.gif -palette 256/256/256

label .pub.img -image logo
pack .pub.img -side top
 
button .pub.button -text "Fermer" -command {set pub_vue 1}
pack .pub.button -side left -expand 1
 
grab .pub
tkwait variable pub_vue
grab release .pub
destroy .pub

}

proc debutsegment {x y} {

global x1_profil 
global y1_profil
global segment_anchor
global segment_last

if {[winfo exists .plotprof]} {
    destroy .plotprof
}
set x1_profil $x 
set y1_profil $y
set xmax [expr [image width sortie] - 1]
set ymax [expr [image height sortie] - 1]
if {$x1_profil > $xmax} {set x1_profil $xmax}
if {$y1_profil > $ymax} {set y1_profil $ymax}
if {$x1_profil < 0} {set x1_profil 0}
if {$y1_profil < 0} {set y1_profil 0}
catch {.imagensaida delete $segment_last}
set segment_anchor [list $x $y]
catch {unset segment_last}
}

proc dessinesegment {x y} {
global segment_anchor
global segment_last
global segment_tag

catch {.imagensaida delete $segment_last}
set segment_last [eval {.imagensaida create line} $segment_anchor {$x $y -tag segment_tag -fill yellow -arrow last}]
}

proc finsegment {x y} {

global x1_profil 
global y1_profil
global x2_profil 
global y2_profil

set x2_profil $x
set y2_profil $y
set xmax [expr [image width sortie] - 1]
set ymax [expr [image height sortie] - 1]
if {$x2_profil > [expr $xmax-1]} {set x2_profil [expr $xmax-1]}
if {$y2_profil > [expr $ymax-1]} {set y2_profil [expr $ymax-1]}
if {$x2_profil < 1} {set x2_profil 1}
if {$y2_profil < 1} {set y2_profil 1}
#puts "($x2_profil,$y2_profil)"
if {($x2_profil != $x1_profil) || \
	($y2_profil != $y1_profil)} {
    draw_profile
}
}

proc draw_profile {} {

global profil_color
global x1_profil 
global y1_profil
global x2_profil 
global y2_profil

set profil_color 0
set winpr [toplevel .plotprof]
wm title $winpr "Profil ($x1_profil,$y1_profil) -> ($x2_profil,$y2_profil)"
# Calcul de la largeur de fenetre
if {[expr abs($x2_profil-$x1_profil)]>[expr abs($y2_profil-$y1_profil)]} {
    set largeur [expr abs($x2_profil-$x1_profil)]
} else {
    set largeur [expr abs($y2_profil-$y1_profil)]
}
puts "Profil : $largeur points"
canvas $winpr.cadre -bg lightgrey
canvas $winpr.body -width $largeur -height 256
pack $winpr.body -in $winpr.cadre -padx 15 -pady 15
pack $winpr.cadre -side top
$winpr.body create line 0 256 $largeur 256 -fill lightgrey -dash 3
$winpr.body create line 0 206 $largeur 206 -fill lightgrey -dash 3
$winpr.body create line 0 156 $largeur 156 -fill lightgrey -dash 3
$winpr.body create line 0 106 $largeur 106 -fill lightgrey -dash 3
$winpr.body create line 0 55 $largeur 56 -fill lightgrey -dash 3
$winpr.body create line 0 6 $largeur 6 -fill lightgrey -dash 3
$winpr.cadre create line 14 286 14 0 -fill black -arrow last
$winpr.cadre create line 0 273 [expr $largeur + 30] 273  \
    -fill black -arrow last
button $winpr.but1 -text "Fermer" -command "destroy $winpr"
button $winpr.but2 -text "Update" -command "update_profile"
pack $winpr.but1 $winpr.but2 -side top

update_profile
}

proc update_profile {} {

global profil_color
global x1_profil 
global y1_profil
global x2_profil 
global y2_profil

set winpr .plotprof
switch $profil_color {
    0 {set plot_color "black"}
    1 {set plot_color "magenta"}
    2 {set plot_color "green"}
    3 {set plot_color "orange"}
    4 {set plot_color "red"}
    5 {set plot_color "blue"}
    6 {set plot_color "purple"}
    7 {set plot_color "cyan"}
    8 {set plot_color "brown"}
    9 {set plot_color "darkgrey"}
}
incr profil_color
if {$profil_color>9} {set profil_color 0}
# calcul du quadrant
if {$x1_profil < $x2_profil} {set sens_x 1} else {set sens_x -1}
if {$y1_profil < $y2_profil} {set sens_y 1} else {set sens_y -1}
# Calcul de l'octant
set absy [expr abs($y2_profil-$y1_profil)]
set absx [expr abs($x2_profil-$x1_profil)]
# Initialisations diverses...
set index 1
set x_profil $x1_profil  
set y_profil $y1_profil
set last_val [lindex [sortie get $x1_profil $y1_profil] 0]
set indic_x [expr $sens_x*($x2_profil-$x_profil)]
set indic_y [expr $sens_y*($y2_profil-$y_profil)]
if {$absx > $absy} {
    set pente [expr double($y2_profil-$y1_profil)/ \
		   double($x2_profil-$x1_profil)]
    #puts $pente
    while {($indic_x>=0)&&($indic_y>=0)} {
	set denom [expr double($x_profil+$sens_x-$x1_profil)]
	set match_horiz [expr abs((double($y_profil-$y1_profil)/$denom)-$pente)]
	set match_both [expr abs((double($y_profil+$sens_y-$y1_profil)/$denom)-$pente)]
	incr x_profil $sens_x
	if {$match_both < $match_horiz} {
	    incr y_profil $sens_y
	}
	set val [lindex [sortie get $x_profil $y_profil] 0]
	set indic_x [expr $sens_x*($x2_profil-$x_profil)]
	set indic_y [expr $sens_y*($y2_profil-$y_profil)]
	$winpr.body create line $index [expr 256 -$last_val] \
	    [expr $index + 1] [expr 256 - $val] -fill $plot_color
	incr index
	set last_val $val
	#puts "($x_profil,$y_profil):$val"
    }
} else {
    set antipente [expr double($x2_profil-$x1_profil)/ \
		       double($y2_profil-$y1_profil)]
    #puts $antipente
    while {($indic_x>=0)&&($indic_y>=0)} {
	set denom [expr double($y_profil+$sens_y-$y1_profil)]
	set match_vert [expr abs((double($x_profil-$x1_profil)/$denom)-$antipente)]
	set match_both [expr abs((double($x_profil+$sens_x-$x1_profil)/$denom)-$antipente)]
	incr y_profil $sens_y
	if {$match_both < $match_vert} {
	    incr x_profil $sens_x
	}
	set val [lindex [sortie get $x_profil $y_profil] 0]
	set indic_x [expr $sens_x*($x2_profil-$x_profil)]
	set indic_y [expr $sens_y*($y2_profil-$y_profil)]
	$winpr.body create line $index [expr 256 -$last_val] \
	    [expr $index + 1] [expr 256 - $val] -fill $plot_color
	incr index
	set last_val $val
	#puts "($x_profil,$y_profil):$val"
    }
}
}
