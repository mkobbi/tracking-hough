Inti 2.0 - Version 2016
*************************************
Instructions de configuration, installation et exécution
*********************************************************

(1a) Configuration pour LINUX
------------------------------
Packages requis : 
* gcc, g++, libstdc++, make, cmake
* tcl, tcl-dev, tk, tk-dev
Package recommandé
* libtk-img
(1b) Configuration pour Windows/Cygwin
---------------------------------------
Installer Cygwin avec les packages suivants :
* gcc-g++, make, cmake
* tcl-tk, tcl-tk-devel
* xorg-server, xinit
(2) Installation et compilation
--------------------------------
- Copier et désarchiver le logiciel sur votre répertoire
de travail (ou votre répertoire Cygwin)
- lancer la commande ./build.sh
(3) Exécution du logiciel
--------------------------
Le logiciel se lance en invoquant l'un des script tcl... 

Attention, sous Cygwin, il faut avoir préalablement démarré 
le serveur X avec la commande :
$ startxwin &
puis avoir positionné le DISPLAY à la valeur fournie, par ex :
$ export DISPLAY=:0.0
